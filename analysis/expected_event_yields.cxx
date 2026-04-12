#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TLeaf.h>
#include <TObject.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <initializer_list>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace AnaExpectedYields {

const int max_particles = 500;
const double units = 1.0e38;
const double default_argon_a = 40.0;
const double mass_lambda = 1.115683;
const double mass_sigma0 = 1.192642;
const double mass_xi0 = 1.31486;
const double mass_xim = 1.32171;
const double mass_omega = 1.67245;
const double mass_proton = 0.9382720813;
const double mass_pion_charged = 0.13957039;
const double mass_pion0 = 0.1349768;
const double mass_kaon_charged = 0.493677;
const double lambda_to_p_piminus_br = 0.639;
const int direct_decay_steps = 48;
const int cascade_decay_steps = 24;

struct StatusRow {
    int plan_line = 0;
    TString generator;
    TString version;
    TString variation;
    TString knob;
    TString beam_mode;
    TString beam_species;
    TString interaction;
    TString fsi_state;
    TString sample;
    TString input_file;
    TString primary_status;
    TString nuclear_exit_status;
    int order = 0;
};

struct YieldCategory {
    TString name;
    TString definition;
};

struct YieldAcc {
    Long64_t raw = 0;
    double expected = 0.0;
    double expected_sum_squares = 0.0;
    double xsec = 0.0;
    double xsec_sum_squares = 0.0;
};

struct YieldRow {
    StatusRow status;
    TString category;
    TString definition;
    YieldAcc acc;
};

struct FluxEnv {
    bool initialized = false;
    bool active = false;
    double emin = 0.0;
    double emax = 0.0;
    double proposal_emin = 0.0;
    double proposal_emax = 10.0;
    double mean_bin_flux = 1.0;
    TString source;
    std::vector<TH1*> hists;
};

struct TreeInputs {
    TLeaf* enu = nullptr;
    TLeaf* weight = nullptr;
    TLeaf* scale = nullptr;
    TLeaf* target_a = nullptr;
    TLeaf* analysis_weight = nullptr;
    TLeaf* xsec_weight = nullptr;
    TLeaf* n_fsp = nullptr;
    TLeaf* pdg = nullptr;
    TLeaf* px = nullptr;
    TLeaf* py = nullptr;
    TLeaf* pz = nullptr;
    TLeaf* energy = nullptr;
};

struct FourVector {
    double px = 0.0;
    double py = 0.0;
    double pz = 0.0;
    double e = 0.0;
};

struct EventSummary {
    bool has_lambda = false;
    bool has_sigma0 = false;
    bool has_xi0 = false;
    bool has_xim = false;
    bool has_omegam = false;
    bool has_visible_proton = false;
    bool has_visible_piminus = false;
    double detector_branching = 0.0;
    double detector_visible = 0.0;
    double visible_from_lambda = 0.0;
    double visible_from_sigma0 = 0.0;
    double visible_from_xi0 = 0.0;
    double visible_from_xim = 0.0;
    double visible_from_omegam = 0.0;
};

std::vector<std::string> split_csv_line(const std::string& line)
{
    std::vector<std::string> fields;
    std::string field;
    bool in_quote = false;
    for (size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (c == '"') {
            if (in_quote && i + 1 < line.size() && line[i + 1] == '"') {
                field += '"';
                ++i;
            } else {
                in_quote = !in_quote;
            }
        } else if (c == ',' && !in_quote) {
            fields.push_back(field);
            field.clear();
        } else {
            field += c;
        }
    }
    fields.push_back(field);
    return fields;
}

TString csv(TString text)
{
    text.ReplaceAll("\"", "\"\"");
    TString out = "\"";
    out += text;
    out += "\"";
    return out;
}

TString lower(TString text)
{
    text.ToLower();
    return text;
}

TLeaf* leaf(TTree* tree, std::initializer_list<const char*> names)
{
    for (const char* name : names) {
        if (TLeaf* found = tree->GetLeaf(name)) return found;
    }
    return nullptr;
}

double value(TLeaf* leaf_in, int i = 0, double fallback = 0.0)
{
    return leaf_in ? leaf_in->GetValue(i) : fallback;
}

int int_value(TLeaf* leaf_in, int i = 0, int fallback = 0)
{
    return leaf_in ? int(std::lround(leaf_in->GetValue(i))) : fallback;
}

double momentum(const TreeInputs& in, int i)
{
    const double x = value(in.px, i);
    const double y = value(in.py, i);
    const double z = value(in.pz, i);
    return std::sqrt(x * x + y * y + z * z);
}

double momentum(const FourVector& p4)
{
    return std::sqrt(p4.px * p4.px + p4.py * p4.py + p4.pz * p4.pz);
}

double mass_gev(int pdg)
{
    if (pdg == 3122) return mass_lambda;
    if (pdg == 3212) return mass_sigma0;
    if (pdg == 3322) return mass_xi0;
    if (pdg == 3312) return mass_xim;
    if (pdg == 3334) return mass_omega;
    if (pdg == 2212) return mass_proton;
    if (std::abs(pdg) == 211) return mass_pion_charged;
    if (pdg == 111) return mass_pion0;
    if (std::abs(pdg) == 321) return mass_kaon_charged;
    return 0.0;
}

FourVector particle_p4(const TreeInputs& in, int i, double mass)
{
    FourVector p4;
    p4.px = value(in.px, i);
    p4.py = value(in.py, i);
    p4.pz = value(in.pz, i);
    const double p = momentum(p4);
    const double nominal_energy = std::sqrt(std::max(0.0, p * p + mass * mass));
    const double tree_energy = value(in.energy, i, nominal_energy);
    p4.e = tree_energy >= 0.999 * nominal_energy ? tree_energy : nominal_energy;
    return p4;
}

int final_count(const TreeInputs& in)
{
    return std::max(0, std::min(int_value(in.n_fsp), max_particles));
}

double clamp_probability(double value_in)
{
    if (!std::isfinite(value_in)) return 0.0;
    return std::max(0.0, std::min(1.0, value_in));
}

void combine_or_probability(double& event_probability, double contribution)
{
    contribution = clamp_probability(contribution);
    event_probability = 1.0 - (1.0 - clamp_probability(event_probability)) * (1.0 - contribution);
}

double two_body_momentum(double parent_mass, double child_mass, double companion_mass)
{
    const double m2 = parent_mass * parent_mass;
    const double plus = (child_mass + companion_mass) * (child_mass + companion_mass);
    const double minus = (child_mass - companion_mass) * (child_mass - companion_mass);
    const double lambda = (m2 - plus) * (m2 - minus);
    if (parent_mass <= 0.0 || lambda <= 0.0) return 0.0;
    return std::sqrt(lambda) / (2.0 * parent_mass);
}

double child_momentum_after_isotropic_decay(double parent_p, double parent_e,
                                            double parent_mass, double child_mass,
                                            double companion_mass, double costheta)
{
    const double pstar = two_body_momentum(parent_mass, child_mass, companion_mass);
    const double estar = std::sqrt(std::max(0.0, pstar * pstar + child_mass * child_mass));
    const double beta = parent_e > 0.0 ? std::min(0.999999, parent_p / parent_e) : 0.0;
    const double gamma = parent_mass > 0.0 ? std::max(1.0, parent_e / parent_mass) : 1.0;
    const double sintheta = std::sqrt(std::max(0.0, 1.0 - costheta * costheta));
    const double p_parallel = gamma * (pstar * costheta + beta * estar);
    const double p_perp = pstar * sintheta;
    return std::sqrt(std::max(0.0, p_parallel * p_parallel + p_perp * p_perp));
}

double child_energy_after_isotropic_decay(double parent_p, double parent_e,
                                          double parent_mass, double child_mass,
                                          double companion_mass, double costheta)
{
    const double pstar = two_body_momentum(parent_mass, child_mass, companion_mass);
    const double estar = std::sqrt(std::max(0.0, pstar * pstar + child_mass * child_mass));
    const double beta = parent_e > 0.0 ? std::min(0.999999, parent_p / parent_e) : 0.0;
    const double gamma = parent_mass > 0.0 ? std::max(1.0, parent_e / parent_mass) : 1.0;
    return gamma * (estar + beta * pstar * costheta);
}

double lambda_daughter_threshold_acceptance(double lambda_p, double lambda_e,
                                            double proton_threshold, double piminus_threshold,
                                            int steps = direct_decay_steps)
{
    if (lambda_e <= 0.0 || steps <= 0) return 0.0;
    const double pstar = two_body_momentum(mass_lambda, mass_proton, mass_pion_charged);
    const double e_proton_star = std::sqrt(pstar * pstar + mass_proton * mass_proton);
    const double e_piminus_star = std::sqrt(pstar * pstar + mass_pion_charged * mass_pion_charged);
    const double beta = std::min(0.999999, std::max(0.0, lambda_p / lambda_e));
    const double gamma = std::max(1.0, lambda_e / mass_lambda);

    int pass = 0;
    for (int i = 0; i < steps; ++i) {
        const double costheta = -1.0 + (2.0 * (i + 0.5) / double(steps));
        const double sintheta = std::sqrt(std::max(0.0, 1.0 - costheta * costheta));
        const double proton_parallel = gamma * (pstar * costheta + beta * e_proton_star);
        const double piminus_parallel = gamma * (-pstar * costheta + beta * e_piminus_star);
        const double p_perp = pstar * sintheta;
        const double proton_p = std::sqrt(proton_parallel * proton_parallel + p_perp * p_perp);
        const double piminus_p = std::sqrt(piminus_parallel * piminus_parallel + p_perp * p_perp);
        if (proton_p > proton_threshold && piminus_p > piminus_threshold) ++pass;
    }
    return double(pass) / double(steps);
}

double lambda_daughter_threshold_acceptance(const FourVector& lambda_p4,
                                            double proton_threshold,
                                            double piminus_threshold)
{
    return lambda_daughter_threshold_acceptance(momentum(lambda_p4), lambda_p4.e,
                                               proton_threshold, piminus_threshold);
}

double two_body_feeddown_visible_acceptance(const FourVector& parent_p4,
                                           double parent_mass,
                                           double lambda_companion_mass,
                                           double proton_threshold,
                                           double piminus_threshold,
                                           int steps = direct_decay_steps)
{
    const double parent_p = momentum(parent_p4);
    const double parent_e = parent_p4.e > parent_p ? parent_p4.e :
        std::sqrt(std::max(0.0, parent_p * parent_p + parent_mass * parent_mass));
    double total = 0.0;
    for (int i = 0; i < steps; ++i) {
        const double costheta = -1.0 + (2.0 * (i + 0.5) / double(steps));
        const double lambda_p = child_momentum_after_isotropic_decay(
            parent_p, parent_e, parent_mass, mass_lambda, lambda_companion_mass, costheta);
        const double lambda_e = child_energy_after_isotropic_decay(
            parent_p, parent_e, parent_mass, mass_lambda, lambda_companion_mass, costheta);
        total += lambda_daughter_threshold_acceptance(lambda_p, lambda_e, proton_threshold,
                                                     piminus_threshold, direct_decay_steps);
    }
    return steps > 0 ? total / double(steps) : 0.0;
}

double cascade_feeddown_visible_acceptance(const FourVector& parent_p4,
                                          double parent_mass,
                                          double intermediate_mass,
                                          double first_companion_mass,
                                          double second_companion_mass,
                                          double proton_threshold,
                                          double piminus_threshold)
{
    const double parent_p = momentum(parent_p4);
    const double parent_e = parent_p4.e > parent_p ? parent_p4.e :
        std::sqrt(std::max(0.0, parent_p * parent_p + parent_mass * parent_mass));
    double total = 0.0;
    for (int i = 0; i < cascade_decay_steps; ++i) {
        const double cos_parent = -1.0 + (2.0 * (i + 0.5) / double(cascade_decay_steps));
        const double intermediate_p = child_momentum_after_isotropic_decay(
            parent_p, parent_e, parent_mass, intermediate_mass, first_companion_mass, cos_parent);
        const double intermediate_e = child_energy_after_isotropic_decay(
            parent_p, parent_e, parent_mass, intermediate_mass, first_companion_mass, cos_parent);
        for (int j = 0; j < cascade_decay_steps; ++j) {
            const double cos_intermediate = -1.0 + (2.0 * (j + 0.5) / double(cascade_decay_steps));
            const double lambda_p = child_momentum_after_isotropic_decay(
                intermediate_p, intermediate_e, intermediate_mass, mass_lambda,
                second_companion_mass, cos_intermediate);
            const double lambda_e = child_energy_after_isotropic_decay(
                intermediate_p, intermediate_e, intermediate_mass, mass_lambda,
                second_companion_mass, cos_intermediate);
            total += lambda_daughter_threshold_acceptance(lambda_p, lambda_e, proton_threshold,
                                                         piminus_threshold, cascade_decay_steps);
        }
    }
    return total / double(cascade_decay_steps * cascade_decay_steps);
}

double visible_feeddown_factor(int parent_pdg, const FourVector& parent_p4,
                               double proton_threshold, double piminus_threshold)
{
    if (parent_pdg == 3122) {
        return lambda_to_p_piminus_br *
            lambda_daughter_threshold_acceptance(parent_p4, proton_threshold, piminus_threshold);
    }
    if (parent_pdg == 3212) {
        return lambda_to_p_piminus_br *
            two_body_feeddown_visible_acceptance(parent_p4, mass_sigma0, 0.0,
                                                proton_threshold, piminus_threshold);
    }
    if (parent_pdg == 3322) {
        return lambda_to_p_piminus_br *
            two_body_feeddown_visible_acceptance(parent_p4, mass_xi0, mass_pion0,
                                                proton_threshold, piminus_threshold);
    }
    if (parent_pdg == 3312) {
        return lambda_to_p_piminus_br *
            two_body_feeddown_visible_acceptance(parent_p4, mass_xim, mass_pion_charged,
                                                proton_threshold, piminus_threshold);
    }
    if (parent_pdg == 3334) {
        const double direct_lambda_k = 0.678 * two_body_feeddown_visible_acceptance(
            parent_p4, mass_omega, mass_kaon_charged, proton_threshold, piminus_threshold);
        const double via_xi0 = 0.236 * cascade_feeddown_visible_acceptance(
            parent_p4, mass_omega, mass_xi0, mass_pion_charged, mass_pion0,
            proton_threshold, piminus_threshold);
        const double via_xim = 0.086 * cascade_feeddown_visible_acceptance(
            parent_p4, mass_omega, mass_xim, mass_pion0, mass_pion_charged,
            proton_threshold, piminus_threshold);
        return lambda_to_p_piminus_br * (direct_lambda_k + via_xi0 + via_xim);
    }
    return 0.0;
}

EventSummary summarize_event(const TreeInputs& in, double proton_threshold, double piminus_threshold)
{
    EventSummary summary;
    const int n = final_count(in);
    for (int i = 0; i < n; ++i) {
        const int code = int_value(in.pdg, i);
        summary.has_lambda = summary.has_lambda || code == 3122;
        summary.has_sigma0 = summary.has_sigma0 || code == 3212;
        summary.has_xi0 = summary.has_xi0 || code == 3322;
        summary.has_xim = summary.has_xim || code == 3312;
        summary.has_omegam = summary.has_omegam || code == 3334;
        if (in.px && in.py && in.pz) {
            const double p = momentum(in, i);
            summary.has_visible_proton = summary.has_visible_proton || (code == 2212 && p > proton_threshold);
            summary.has_visible_piminus = summary.has_visible_piminus || (code == -211 && p > piminus_threshold);
        }

        const double parent_mass = mass_gev(code);
        if (parent_mass <= 0.0) continue;
        const FourVector parent = particle_p4(in, i, parent_mass);
        const double branching_factor = lambda_to_p_piminus_br;
        const double visible_factor = visible_feeddown_factor(code, parent, proton_threshold, piminus_threshold);
        if (!std::isfinite(visible_factor)) continue;

        if (code == 3122) {
            combine_or_probability(summary.visible_from_lambda, visible_factor);
        } else if (code == 3212) {
            combine_or_probability(summary.visible_from_sigma0, visible_factor);
        } else if (code == 3322) {
            combine_or_probability(summary.visible_from_xi0, visible_factor);
        } else if (code == 3312) {
            combine_or_probability(summary.visible_from_xim, visible_factor);
        } else if (code == 3334) {
            combine_or_probability(summary.visible_from_omegam, visible_factor);
        } else {
            continue;
        }
        combine_or_probability(summary.detector_branching, branching_factor);
        combine_or_probability(summary.detector_visible, visible_factor);
    }
    return summary;
}

std::vector<YieldCategory> categories()
{
    return {
        {"ana_final_hyperon_inclusive",
         "generator final-state contains Lambda, Sigma0, Xi0, Xi-, or Omega- with no momentum requirement"},
        {"final_Lambda",
         "generator final-state contains Lambda, PDG 3122, with no momentum requirement"},
        {"final_Sigma0",
         "generator final-state contains Sigma0, PDG 3212, with no momentum requirement"},
        {"final_Xi0",
         "generator final-state contains Xi0, PDG 3322, with no momentum requirement"},
        {"final_Xim",
         "generator final-state contains Xi-, PDG 3312, with no momentum requirement"},
        {"final_Omegam",
         "generator final-state contains Omega-, PDG 3334, with no momentum requirement"},
        {"final_Lambda_visible_p_piminus_proxy",
         "final-state Lambda plus proxy visible proton and pi- above the configured thresholds"},
        {"detector_Lambda_to_p_piminus_branching",
         "fast detector envelope: final-state Lambda/Sigma0/Xi0/Xi-/Omega- reach a Lambda -> p pi- decay, branching ratios only"},
        {"detector_visible_Lambda_to_p_piminus",
         "fast detector envelope: branching ratios plus two-body Lambda decay probability for visible proton and pi- thresholds"},
        {"detector_visible_from_final_Lambda",
         "fast detector envelope contribution seeded by final-state Lambda"},
        {"detector_visible_from_final_Sigma0",
         "fast detector envelope contribution seeded by final-state Sigma0 -> Lambda gamma"},
        {"detector_visible_from_final_Xi0",
         "fast detector envelope contribution seeded by final-state Xi0 -> Lambda pi0"},
        {"detector_visible_from_final_Xim",
         "fast detector envelope contribution seeded by final-state Xi- -> Lambda pi-"},
        {"detector_visible_from_final_Omegam",
         "fast detector envelope contribution seeded by final-state Omega- -> Lambda K or Xi pi feed-down"}
    };
}

double category_factor(const EventSummary& summary, const TString& category)
{
    if (category == "ana_final_hyperon_inclusive") {
        return summary.has_lambda || summary.has_sigma0 || summary.has_xi0 ||
            summary.has_xim || summary.has_omegam ? 1.0 : 0.0;
    }
    if (category == "final_Lambda") return summary.has_lambda ? 1.0 : 0.0;
    if (category == "final_Sigma0") return summary.has_sigma0 ? 1.0 : 0.0;
    if (category == "final_Xi0") return summary.has_xi0 ? 1.0 : 0.0;
    if (category == "final_Xim") return summary.has_xim ? 1.0 : 0.0;
    if (category == "final_Omegam") return summary.has_omegam ? 1.0 : 0.0;
    if (category == "final_Lambda_visible_p_piminus_proxy") {
        return summary.has_lambda && summary.has_visible_proton && summary.has_visible_piminus ? 1.0 : 0.0;
    }
    if (category == "detector_Lambda_to_p_piminus_branching") return summary.detector_branching;
    if (category == "detector_visible_Lambda_to_p_piminus") return summary.detector_visible;
    if (category == "detector_visible_from_final_Lambda") return summary.visible_from_lambda;
    if (category == "detector_visible_from_final_Sigma0") return summary.visible_from_sigma0;
    if (category == "detector_visible_from_final_Xi0") return summary.visible_from_xi0;
    if (category == "detector_visible_from_final_Xim") return summary.visible_from_xim;
    if (category == "detector_visible_from_final_Omegam") return summary.visible_from_omegam;
    return 0.0;
}

void add(YieldAcc& acc, double expected_weight, double xsec_weight)
{
    ++acc.raw;
    acc.expected += expected_weight;
    acc.expected_sum_squares += expected_weight * expected_weight;
    acc.xsec += xsec_weight;
    acc.xsec_sum_squares += xsec_weight * xsec_weight;
}

TString normalize_beam_mode(TString beam_mode)
{
    beam_mode.ToLower();
    if (beam_mode.Contains("rhc")) return "rhc";
    if (beam_mode.Contains("fhc")) return "fhc";
    if (beam_mode == "" || beam_mode == "combined" || beam_mode == "both" || beam_mode == "all") return "combined";
    return "";
}

TString normalize_beam_species(TString beam_species)
{
    beam_species.ToLower();
    if (beam_species == "" || beam_species == "combined" || beam_species == "both" ||
        beam_species == "all") return "";
    if (beam_species == "numu" || beam_species == "nu_mu" || beam_species == "14") return "numu";
    if (beam_species == "numubar" || beam_species == "anti_numu" || beam_species == "-14") return "numubar";
    return beam_species;
}

TString infer_beam_mode(TString text)
{
    text.ToLower();
    if (text.Contains("rhc")) return "rhc";
    if (text.Contains("fhc")) return "fhc";
    return "combined";
}

TString infer_beam_species(TString text)
{
    text.ToLower();
    if (text.Contains("numubar") || text.Contains("anti_numu")) return "numubar";
    if (text.Contains("numu")) return "numu";
    return "";
}

TString flux_hist_path(const TString& beam_mode, const TString& beam_species)
{
    TString path = beam_mode;
    path += "/";
    path += beam_species;
    path += "/Detsmear/";
    path += (beam_species == "numubar" ? "numubar_CV_AV_TPC_5MeV_bin" : "numu_CV_AV_TPC_5MeV_bin");
    return path;
}

TString hist_leaf_name(TString path)
{
    const Ssiz_t slash = path.Last('/');
    return slash < 0 ? path : path(slash + 1, path.Length() - slash - 1);
}

std::vector<TString> flux_hist_paths(TString beam_mode, TString beam_species)
{
    beam_mode = normalize_beam_mode(beam_mode);
    beam_species = normalize_beam_species(beam_species);
    std::vector<TString> paths;
    auto add = [&](const TString& mode, const TString& species) {
        paths.push_back(flux_hist_path(mode, species));
    };
    if (beam_mode == "combined") {
        if (beam_species == "") {
            add("fhc", "numu");
            add("fhc", "numubar");
            add("rhc", "numu");
            add("rhc", "numubar");
        } else {
            add("fhc", beam_species);
            add("rhc", beam_species);
        }
    } else if (beam_mode == "fhc" || beam_mode == "rhc") {
        add(beam_mode, beam_species == "" ? TString(beam_mode == "rhc" ? "numubar" : "numu") : beam_species);
    }
    return paths;
}

TH1* find_hist(TDirectory* dir, const TString& name)
{
    if (!dir) return nullptr;
    if (TH1* direct = dynamic_cast<TH1*>(dir->Get(name))) return direct;
    TIter next(dir->GetListOfKeys());
    while (TKey* key = static_cast<TKey*>(next())) {
        TObject* obj = key->ReadObj();
        if (!obj) continue;
        if (obj->InheritsFrom(TH1::Class()) && name == obj->GetName()) return static_cast<TH1*>(obj);
        if (obj->InheritsFrom(TDirectory::Class())) {
            if (TH1* found = find_hist(static_cast<TDirectory*>(obj), name)) return found;
        }
    }
    return nullptr;
}

TH1* find_flux_hist(TFile* file, const TString& path)
{
    TH1* hist = dynamic_cast<TH1*>(file->Get(path));
    return hist ? hist : find_hist(file, hist_leaf_name(path));
}

bool extend_flux(TH1* hist, double floor_fraction, FluxEnv& env)
{
    if (!hist || hist->GetMaximum() <= 0.0) return false;
    const double floor = std::max(0.0, floor_fraction) * hist->GetMaximum();
    bool found = false;
    double emin = 0.0;
    double emax = 0.0;
    for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
        if (hist->GetBinContent(bin) <= floor) continue;
        const double low = hist->GetXaxis()->GetBinLowEdge(bin);
        const double high = hist->GetXaxis()->GetBinUpEdge(bin);
        emin = found ? std::min(emin, low) : low;
        emax = found ? std::max(emax, high) : high;
        found = true;
    }
    if (!found) return false;
    TH1* clone = static_cast<TH1*>(hist->Clone(TString::Format("%s_expected_yield_flux_%d", hist->GetName(), int(env.hists.size()))));
    clone->SetDirectory(nullptr);
    env.hists.push_back(clone);
    env.emin = env.active ? std::min(env.emin, emin) : emin;
    env.emax = env.active ? std::max(env.emax, emax) : emax;
    env.active = true;
    return true;
}

double mean_flux_in_proposal(const FluxEnv& env)
{
    double total = 0.0;
    int bins = 0;
    for (TH1* hist : env.hists) {
        for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
            const double center = hist->GetXaxis()->GetBinCenter(bin);
            if (center < env.proposal_emin || center >= env.proposal_emax) continue;
            total += std::max(0.0, hist->GetBinContent(bin));
            ++bins;
        }
    }
    return bins > 0 && !env.hists.empty() ? total / (double(bins) / double(env.hists.size())) : 0.0;
}

FluxEnv read_flux(TString path, TString beam_mode, TString beam_species, double floor_fraction,
                  double proposal_emin, double proposal_emax)
{
    FluxEnv env;
    env.initialized = true;
    env.proposal_emin = proposal_emin;
    env.proposal_emax = proposal_emax;

    if (path == "") {
        env.active = false;
        env.source = "no_flux_reweight";
        return env;
    }

    beam_mode = normalize_beam_mode(beam_mode);
    if (beam_mode == "") {
        std::cerr << "Beam mode must be FHC, RHC, or combined." << std::endl;
        return env;
    }

    TFile* file = TFile::Open(path, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Could not open NuMI flux file: " << path.Data() << std::endl;
        if (file) file->Close();
        return env;
    }

    std::vector<TString> used;
    const std::vector<TString> names = flux_hist_paths(beam_mode, beam_species);
    for (const TString& name : names) {
        if (extend_flux(find_flux_hist(file, name), floor_fraction, env)) used.push_back(name);
        else std::cerr << "Could not use NuMI flux histogram: " << name.Data() << std::endl;
    }
    file->Close();
    delete file;

    if (!env.active || env.emax <= env.emin) {
        env.active = false;
        std::cerr << "No valid non-zero NuMI flux support was found in " << path.Data() << std::endl;
        return env;
    }
    env.emin = std::max(env.emin, env.proposal_emin);
    env.emax = std::min(env.emax, env.proposal_emax);
    env.mean_bin_flux = mean_flux_in_proposal(env);
    if (env.emax <= env.emin || env.mean_bin_flux <= 0.0) {
        env.active = false;
        std::cerr << "No valid NuMI flux overlap with flat proposal range "
                  << env.proposal_emin << " to " << env.proposal_emax << " GeV." << std::endl;
        return env;
    }

    env.source = path;
    env.source += " [";
    for (int i = 0; i < int(used.size()); ++i) {
        if (i) env.source += ", ";
        env.source += used[i];
    }
    env.source += "]";
    return env;
}

double flux_reweight(const FluxEnv& env, double enu)
{
    if (!env.initialized || !env.active) return 1.0;
    if (enu < env.emin || enu >= env.emax || env.mean_bin_flux <= 0.0) return 0.0;
    double flux = 0.0;
    for (TH1* hist : env.hists) {
        const int bin = hist->FindBin(enu);
        if (bin < 1 || bin > hist->GetNbinsX()) continue;
        flux += std::max(0.0, hist->GetBinContent(bin));
    }
    return flux / env.mean_bin_flux;
}

void clear_flux_cache(std::map<std::string, FluxEnv>& cache)
{
    for (auto& item : cache) {
        for (TH1* hist : item.second.hists) delete hist;
        item.second.hists.clear();
    }
    cache.clear();
}

FluxEnv& flux_for_row(const StatusRow& row, const TString& flux_file, double floor_fraction,
                      double proposal_emin, double proposal_emax,
                      std::map<std::string, FluxEnv>& cache)
{
    TString beam_mode = normalize_beam_mode(row.beam_mode);
    TString beam_species = normalize_beam_species(row.beam_species);
    if (beam_mode == "") beam_mode = infer_beam_mode(row.input_file + " " + row.sample);
    if (beam_species == "") beam_species = infer_beam_species(row.input_file + " " + row.sample);
    const std::string key = std::string(flux_file.Data()) + "|" + beam_mode.Data() + "|" + beam_species.Data();
    FluxEnv& env = cache[key];
    if (!env.initialized) env = read_flux(flux_file, beam_mode, beam_species, floor_fraction, proposal_emin, proposal_emax);
    return env;
}

TreeInputs bind_inputs(TTree* tree)
{
    TreeInputs in;
    in.enu = leaf(tree, {"Enu_true", "enu_true"});
    in.weight = leaf(tree, {"Weight", "weight"});
    in.scale = leaf(tree, {"fScaleFactor", "scale_factor"});
    in.target_a = leaf(tree, {"tgta", "target_a"});
    in.analysis_weight = leaf(tree, {"analysis_weight"});
    in.xsec_weight = leaf(tree, {"xsec_weight_1e38_per_Ar", "xsec_weight_1e38_cm2_per_Ar",
                                 "xsec_weight_1e38_per_target", "xsec_weight"});
    in.n_fsp = leaf(tree, {"nfsp", "n_fsp"});
    in.pdg = leaf(tree, {"pdg"});
    in.px = leaf(tree, {"px"});
    in.py = leaf(tree, {"py"});
    in.pz = leaf(tree, {"pz"});
    in.energy = leaf(tree, {"E", "energy"});
    return in;
}

double base_analysis_weight(const TreeInputs& in)
{
    if (in.analysis_weight) return value(in.analysis_weight, 0, 1.0);
    return value(in.weight, 0, 1.0) * value(in.scale, 0, 1.0);
}

double xsec_weight_1e38_per_Ar(const TreeInputs& in, double analysis_weight)
{
    if (in.xsec_weight) return value(in.xsec_weight, 0, 0.0);
    const double target_a = std::max(1.0, value(in.target_a, 0, default_argon_a));
    return analysis_weight * target_a * units;
}

bool scan_file(const StatusRow& row, const TString& flux_file, double flux_floor_fraction,
               double proposal_emin, double proposal_emax, double exposure_scale,
               double proton_threshold, double piminus_threshold,
               std::map<std::string, FluxEnv>& flux_cache, std::vector<YieldRow>& out)
{
    TFile* file = TFile::Open(row.input_file, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Could not open yield input: " << row.input_file.Data() << std::endl;
        if (file) file->Close();
        return false;
    }

    TTree* tree = nullptr;
    file->GetObject("FlatTree_VARS", tree);
    if (!tree) {
        std::cerr << "Could not find FlatTree_VARS in " << row.input_file.Data() << std::endl;
        file->Close();
        delete file;
        return false;
    }

    TreeInputs in = bind_inputs(tree);
    if (!in.enu && flux_file != "") {
        std::cerr << "Missing Enu_true/enu_true in " << row.input_file.Data()
                  << "; cannot apply NuMI flux reweighting." << std::endl;
        file->Close();
        delete file;
        return false;
    }
    if (!in.n_fsp || !in.pdg) {
        std::cerr << "Missing exact FlatTree final-state nfsp/pdg arrays in "
                  << row.input_file.Data() << "; cannot apply ana final-state hyperon rule."
                  << std::endl;
        file->Close();
        delete file;
        return false;
    }

    std::vector<YieldCategory> cats = categories();
    std::vector<YieldRow> local;
    for (const YieldCategory& cat : cats) {
        YieldRow yield;
        yield.status = row;
        yield.category = cat.name;
        yield.definition = cat.definition;
        local.push_back(yield);
    }

    FluxEnv& flux = flux_for_row(row, flux_file, flux_floor_fraction, proposal_emin, proposal_emax, flux_cache);
    if (flux_file != "" && !flux.active) {
        std::cerr << "Skipping " << row.input_file.Data() << " because NuMI flux reweighting is inactive." << std::endl;
        file->Close();
        delete file;
        return false;
    }

    const Long64_t entries = tree->GetEntries();
    for (Long64_t entry = 0; entry < entries; ++entry) {
        tree->GetEntry(entry);
        const double flux_weight = in.enu ? flux_reweight(flux, value(in.enu)) : 1.0;
        if (!std::isfinite(flux_weight) || flux_weight <= 0.0) continue;

        const double base_weight = base_analysis_weight(in);
        const double expected_weight = base_weight * flux_weight * exposure_scale;
        const double xsec_weight = xsec_weight_1e38_per_Ar(in, base_weight) * flux_weight;
        if (!std::isfinite(expected_weight) || !std::isfinite(xsec_weight)) continue;

        const EventSummary summary = summarize_event(in, proton_threshold, piminus_threshold);
        for (YieldRow& yield : local) {
            const double factor = category_factor(summary, yield.category);
            if (factor > 0.0) add(yield.acc, expected_weight * factor, xsec_weight * factor);
        }
    }

    out.insert(out.end(), local.begin(), local.end());
    file->Close();
    delete file;
    return true;
}

std::vector<StatusRow> read_status_csv(const TString& status_csv)
{
    std::ifstream input(status_csv.Data());
    std::vector<StatusRow> rows;
    if (!input) {
        std::cerr << "Could not open generator matrix status CSV: " << status_csv.Data() << std::endl;
        return rows;
    }

    std::string line;
    bool first = true;
    int order = 0;
    while (std::getline(input, line)) {
        if (first) {
            first = false;
            continue;
        }
        const std::vector<std::string> f = split_csv_line(line);
        if (f.size() < 13) continue;
        StatusRow row;
        row.plan_line = std::atoi(f[0].c_str());
        row.generator = f[1].c_str();
        row.version = f[2].c_str();
        row.variation = f[3].c_str();
        row.knob = f[4].c_str();
        row.beam_mode = f[5].c_str();
        row.beam_species = f[6].c_str();
        row.interaction = f[7].c_str();
        row.fsi_state = f[8].c_str();
        row.sample = f[9].c_str();
        row.input_file = f[10].c_str();
        row.primary_status = f[11].c_str();
        row.nuclear_exit_status = f[12].c_str();
        row.order = order++;
        if (!gSystem->AccessPathName(row.input_file.Data())) rows.push_back(row);
    }
    return rows;
}

int central_score(const YieldRow& row)
{
    const TString labels = lower(row.status.variation + " " + row.status.knob + " " +
                                 row.status.version + " " + row.status.fsi_state);
    int score = 1000 + row.status.order;
    if (labels.Contains("nominal") || labels.Contains("central") ||
        labels.Contains("ar23_20i_00_000")) {
        score = 0;
    } else if (labels.Contains("all_strange")) {
        score = 100;
    } else if (labels.Contains("dis_only")) {
        score = 200;
    } else if (labels.Contains("hyp_all")) {
        score = 300;
    }
    if (labels.Contains("fsi_off")) score += 25;
    return score + row.status.order;
}

const YieldRow& central_row(const std::vector<YieldRow>& rows)
{
    int best = 0;
    int best_score = central_score(rows[0]);
    for (int i = 1; i < int(rows.size()); ++i) {
        const int score = central_score(rows[i]);
        if (score < best_score) {
            best = i;
            best_score = score;
        }
    }
    return rows[best];
}

TString group_key(const YieldRow& row, bool include_fsi_state)
{
    TString key = row.status.generator + "|" + row.status.beam_mode + "|" +
        row.status.beam_species + "|" + row.status.interaction + "|";
    if (include_fsi_state) key += row.status.fsi_state + "|";
    key += row.category;
    return key;
}

TString generator_group_key(const YieldRow& row, bool include_fsi_state)
{
    TString key = row.status.beam_mode + "|" + row.status.beam_species + "|" +
        row.status.interaction + "|";
    if (include_fsi_state) key += row.status.fsi_state + "|";
    key += row.category;
    return key;
}

TString variation_label(const YieldRow& row)
{
    TString label = row.status.generator + " ";
    label += row.status.version;
    label += " ";
    label += row.status.variation;
    if (row.status.knob != "" && !row.status.variation.Contains(row.status.knob)) {
        label += " ";
        label += row.status.knob;
    }
    if (row.status.fsi_state != "") {
        label += " ";
        label += row.status.fsi_state;
    }
    return label;
}

void write_variation_rows(const TString& path, const std::vector<YieldRow>& rows,
                          double exposure_scale, double proton_threshold,
                          double piminus_threshold)
{
    std::ofstream output(path.Data());
    if (!output) {
        std::cerr << "Could not open yield variation CSV: " << path.Data() << std::endl;
        return;
    }

    output << "sample,generator,version,variation,knob,beam_mode,beam_species,interaction,"
           << "fsi_state,category,raw_events,expected_events,expected_stat_uncertainty,"
           << "xsec_weighted_1e38_cm2_per_Ar,xsec_weighted_stat_uncertainty_1e38_cm2_per_Ar,"
           << "exposure_scale,lambda_to_p_piminus_br,proton_threshold_gev,"
           << "piminus_threshold_gev,definition\n";
    output << std::setprecision(12);
    for (const YieldRow& row : rows) {
        output << csv(row.status.sample).Data() << ','
               << csv(row.status.generator).Data() << ','
               << csv(row.status.version).Data() << ','
               << csv(row.status.variation).Data() << ','
               << csv(row.status.knob).Data() << ','
               << csv(row.status.beam_mode).Data() << ','
               << csv(row.status.beam_species).Data() << ','
               << csv(row.status.interaction).Data() << ','
               << csv(row.status.fsi_state).Data() << ','
               << csv(row.category).Data() << ','
               << row.acc.raw << ','
               << row.acc.expected << ','
               << std::sqrt(row.acc.expected_sum_squares) << ','
               << row.acc.xsec << ','
               << std::sqrt(row.acc.xsec_sum_squares) << ','
               << exposure_scale << ','
               << lambda_to_p_piminus_br << ','
               << proton_threshold << ','
               << piminus_threshold << ','
               << csv(row.definition).Data() << '\n';
    }
}

void write_uncertainty_rows(const TString& path, const std::vector<YieldRow>& rows,
                            bool include_fsi_state, double exposure_scale,
                            double proton_threshold, double piminus_threshold)
{
    std::map<std::string, std::vector<YieldRow>> grouped;
    for (const YieldRow& row : rows) grouped[group_key(row, include_fsi_state).Data()].push_back(row);

    std::ofstream output(path.Data());
    if (!output) {
        std::cerr << "Could not open yield uncertainty CSV: " << path.Data() << std::endl;
        return;
    }

    output << "generator,beam_mode,beam_species,interaction";
    if (include_fsi_state) output << ",fsi_state";
    output << ",category,central_sample,central_version,central_variation,central_knob,"
           << "central_fsi_state,n_variations,central_raw_events,expected_events,"
           << "expected_stat_uncertainty,expected_systematic_uncertainty,"
           << "expected_total_uncertainty,xsec_weighted_1e38_cm2_per_Ar,"
           << "xsec_weighted_stat_uncertainty_1e38_cm2_per_Ar,"
           << "xsec_weighted_systematic_uncertainty_1e38_cm2_per_Ar,"
           << "xsec_weighted_total_uncertainty_1e38_cm2_per_Ar,"
           << "min_expected_events,max_expected_events,"
           << "min_xsec_weighted_1e38_cm2_per_Ar,max_xsec_weighted_1e38_cm2_per_Ar,"
           << "exposure_scale,lambda_to_p_piminus_br,proton_threshold_gev,"
           << "piminus_threshold_gev,definition\n";
    output << std::setprecision(12);

    for (const auto& item : grouped) {
        const std::vector<YieldRow>& group = item.second;
        if (group.empty()) continue;
        const YieldRow& central = central_row(group);
        double expected_syst = 0.0;
        double xsec_syst = 0.0;
        double min_expected = central.acc.expected;
        double max_expected = central.acc.expected;
        double min_xsec = central.acc.xsec;
        double max_xsec = central.acc.xsec;

        for (const YieldRow& row : group) {
            expected_syst = std::max(expected_syst, std::abs(row.acc.expected - central.acc.expected));
            xsec_syst = std::max(xsec_syst, std::abs(row.acc.xsec - central.acc.xsec));
            min_expected = std::min(min_expected, row.acc.expected);
            max_expected = std::max(max_expected, row.acc.expected);
            min_xsec = std::min(min_xsec, row.acc.xsec);
            max_xsec = std::max(max_xsec, row.acc.xsec);
        }

        const double expected_stat = std::sqrt(central.acc.expected_sum_squares);
        const double xsec_stat = std::sqrt(central.acc.xsec_sum_squares);
        const double expected_total = std::sqrt(expected_stat * expected_stat + expected_syst * expected_syst);
        const double xsec_total = std::sqrt(xsec_stat * xsec_stat + xsec_syst * xsec_syst);

        output << csv(central.status.generator).Data() << ','
               << csv(central.status.beam_mode).Data() << ','
               << csv(central.status.beam_species).Data() << ','
               << csv(central.status.interaction).Data() << ',';
        if (include_fsi_state) output << csv(central.status.fsi_state).Data() << ',';
        output << csv(central.category).Data() << ','
               << csv(central.status.sample).Data() << ','
               << csv(central.status.version).Data() << ','
               << csv(central.status.variation).Data() << ','
               << csv(central.status.knob).Data() << ','
               << csv(central.status.fsi_state).Data() << ','
               << group.size() << ','
               << central.acc.raw << ','
               << central.acc.expected << ','
               << expected_stat << ','
               << expected_syst << ','
               << expected_total << ','
               << central.acc.xsec << ','
               << xsec_stat << ','
               << xsec_syst << ','
               << xsec_total << ','
               << min_expected << ','
               << max_expected << ','
               << min_xsec << ','
               << max_xsec << ','
               << exposure_scale << ','
               << lambda_to_p_piminus_br << ','
               << proton_threshold << ','
               << piminus_threshold << ','
               << csv(central.definition).Data() << '\n';
    }
}

void write_generator_envelope_rows(const TString& path, const std::vector<YieldRow>& rows,
                                   bool include_fsi_state, double exposure_scale,
                                   double proton_threshold, double piminus_threshold)
{
    std::map<std::string, std::vector<YieldRow>> grouped;
    for (const YieldRow& row : rows) {
        grouped[generator_group_key(row, include_fsi_state).Data()].push_back(row);
    }

    std::ofstream output(path.Data());
    if (!output) {
        std::cerr << "Could not open generator-envelope yield CSV: " << path.Data() << std::endl;
        return;
    }

    output << "beam_mode,beam_species,interaction";
    if (include_fsi_state) output << ",fsi_state";
    output << ",category,reference_sample,reference_generator,reference_version,"
           << "reference_variation,reference_knob,reference_fsi_state,n_members,"
           << "reference_raw_events,expected_events,expected_stat_uncertainty,"
           << "expected_generator_envelope_uncertainty,expected_total_uncertainty,"
           << "xsec_weighted_1e38_cm2_per_Ar,"
           << "xsec_weighted_stat_uncertainty_1e38_cm2_per_Ar,"
           << "xsec_weighted_generator_envelope_uncertainty_1e38_cm2_per_Ar,"
           << "xsec_weighted_total_uncertainty_1e38_cm2_per_Ar,"
           << "min_expected_events,max_expected_events,min_expected_label,"
           << "max_expected_label,min_xsec_weighted_1e38_cm2_per_Ar,"
           << "max_xsec_weighted_1e38_cm2_per_Ar,min_xsec_label,max_xsec_label,"
           << "exposure_scale,lambda_to_p_piminus_br,proton_threshold_gev,"
           << "piminus_threshold_gev,definition\n";
    output << std::setprecision(12);

    for (const auto& item : grouped) {
        const std::vector<YieldRow>& group = item.second;
        if (group.empty()) continue;
        const YieldRow& central = central_row(group);
        double expected_syst = 0.0;
        double xsec_syst = 0.0;
        double min_expected = central.acc.expected;
        double max_expected = central.acc.expected;
        double min_xsec = central.acc.xsec;
        double max_xsec = central.acc.xsec;
        const YieldRow* min_expected_row = &central;
        const YieldRow* max_expected_row = &central;
        const YieldRow* min_xsec_row = &central;
        const YieldRow* max_xsec_row = &central;

        for (const YieldRow& row : group) {
            expected_syst = std::max(expected_syst, std::abs(row.acc.expected - central.acc.expected));
            xsec_syst = std::max(xsec_syst, std::abs(row.acc.xsec - central.acc.xsec));
            if (row.acc.expected < min_expected) {
                min_expected = row.acc.expected;
                min_expected_row = &row;
            }
            if (row.acc.expected > max_expected) {
                max_expected = row.acc.expected;
                max_expected_row = &row;
            }
            if (row.acc.xsec < min_xsec) {
                min_xsec = row.acc.xsec;
                min_xsec_row = &row;
            }
            if (row.acc.xsec > max_xsec) {
                max_xsec = row.acc.xsec;
                max_xsec_row = &row;
            }
        }

        const double expected_stat = std::sqrt(central.acc.expected_sum_squares);
        const double xsec_stat = std::sqrt(central.acc.xsec_sum_squares);
        const double expected_total = std::sqrt(expected_stat * expected_stat + expected_syst * expected_syst);
        const double xsec_total = std::sqrt(xsec_stat * xsec_stat + xsec_syst * xsec_syst);

        output << csv(central.status.beam_mode).Data() << ','
               << csv(central.status.beam_species).Data() << ','
               << csv(central.status.interaction).Data() << ',';
        if (include_fsi_state) output << csv(central.status.fsi_state).Data() << ',';
        output << csv(central.category).Data() << ','
               << csv(central.status.sample).Data() << ','
               << csv(central.status.generator).Data() << ','
               << csv(central.status.version).Data() << ','
               << csv(central.status.variation).Data() << ','
               << csv(central.status.knob).Data() << ','
               << csv(central.status.fsi_state).Data() << ','
               << group.size() << ','
               << central.acc.raw << ','
               << central.acc.expected << ','
               << expected_stat << ','
               << expected_syst << ','
               << expected_total << ','
               << central.acc.xsec << ','
               << xsec_stat << ','
               << xsec_syst << ','
               << xsec_total << ','
               << min_expected << ','
               << max_expected << ','
               << csv(variation_label(*min_expected_row)).Data() << ','
               << csv(variation_label(*max_expected_row)).Data() << ','
               << min_xsec << ','
               << max_xsec << ','
               << csv(variation_label(*min_xsec_row)).Data() << ','
               << csv(variation_label(*max_xsec_row)).Data() << ','
               << exposure_scale << ','
               << lambda_to_p_piminus_br << ','
               << proton_threshold << ','
               << piminus_threshold << ','
               << csv(central.definition).Data() << '\n';
    }
}

void write_summary_md(const TString& path, const TString& status_csv, const TString& flux_file,
                      double proposal_emin, double proposal_emax, double exposure_scale,
                      double proton_threshold, double piminus_threshold,
                      int samples_processed, int rows_written)
{
    std::ofstream out(path.Data());
    if (!out) return;
    out << "Expected event yields\n"
        << "=====================\n\n"
        << "Input status CSV: `" << status_csv.Data() << "`\n\n"
        << "NuMI flux file: `" << (flux_file == "" ? TString("none") : flux_file).Data() << "`\n\n"
        << "Flat proposal range: " << proposal_emin << " to " << proposal_emax << " GeV\n\n"
        << "Exposure scale applied to `expected_events`: " << exposure_scale << "\n\n"
        << "Fast detector-envelope Lambda -> p pi- branching ratio: "
        << lambda_to_p_piminus_br << "\n\n"
        << "Visible daughter thresholds: proton p > " << proton_threshold
        << " GeV, pi- p > " << piminus_threshold << " GeV\n\n"
        << "Samples processed: " << samples_processed << "\n\n"
        << "Per-variation rows written: " << rows_written << "\n\n"
        << "Definitions:\n\n"
        << "* `expected_events` is `analysis_weight * normalized_NuMI_flux_weight * exposure_scale`.\n"
        << "  With the default exposure scale of 1, this is a unit-normalized analysis yield;\n"
        << "  pass an exposure/POT/target scale if absolute detector event counts are needed.\n"
        << "* `xsec_weighted_1e38_cm2_per_Ar` uses an existing xsec-weight branch when present;\n"
        << "  otherwise it uses `analysis_weight * A * 1e38`, with A taken from `tgta/target_a`\n"
        << "  and defaulting to argon A=40.\n"
        << "* Statistical uncertainties are from `sqrt(sum weights^2)`.\n"
        << "* Systematic uncertainties are variation-envelope uncertainties around the chosen\n"
        << "  central row. The fixed-FSI CSV groups by generator, beam, interaction, FSI state,\n"
        << "  and category; the all-variations CSV also envelopes FSI-state differences.\n"
        << "* The generator-envelope CSVs drop generator from the grouping and envelope over\n"
        << "  generator, version, configuration/variation, knob, and optionally FSI state.\n"
        << "* Detector-envelope rows are fast truth-level corrections. They keep the generator\n"
        << "  sample definition unchanged, then weight eligible final-state hyperons by\n"
        << "  Lambda reachability, the Lambda -> p pi- branching ratio, and the probability\n"
        << "  for the two visible daughter momenta to pass threshold. This does not replace\n"
        << "  a full geometry/containment envelope or a reduced Stage C/D ntuple.\n"
        << "* The inclusive detector-envelope row uses an event-level OR probability when\n"
        << "  multiple eligible final-state hyperons are present. Parent-seeded detector\n"
        << "  rows are non-exclusive composition rows and should not be summed as an event\n"
        << "  yield without accounting for overlap.\n"
        << "* The main selection category is exactly the `ana` final-state rule: Lambda,\n"
        << "  Sigma0, Xi0, Xi-, or Omega- in the generator final-state PDG array, with no\n"
        << "  momentum threshold.\n";
}

}  // namespace AnaExpectedYields

void expected_event_yields(
    TString status_csv = "analysis/output/matrix/generator_matrix_status.csv",
    TString output_dir = "analysis/output/expected_event_yields",
    TString flux_file = "analysis/flux/microboone_numi_flux_5mev.root",
    double flux_floor_fraction = 0.0,
    double proposal_emin_gev = 0.0,
    double proposal_emax_gev = 10.0,
    double exposure_scale = 1.0,
    double proton_threshold_gev = 0.30,
    double piminus_threshold_gev = 0.07
)
{
    using namespace AnaExpectedYields;

    const std::vector<StatusRow> rows = read_status_csv(status_csv);
    if (rows.empty()) {
        std::cerr << "No present matrix inputs found in " << status_csv.Data() << std::endl;
        return;
    }

    gSystem->mkdir(output_dir.Data(), true);
    std::vector<YieldRow> yields;
    std::map<std::string, FluxEnv> flux_cache;
    int samples_processed = 0;
    for (const StatusRow& row : rows) {
        if (scan_file(row, flux_file, flux_floor_fraction, proposal_emin_gev,
                      proposal_emax_gev, exposure_scale, proton_threshold_gev,
                      piminus_threshold_gev, flux_cache, yields)) {
            ++samples_processed;
        }
    }

    if (yields.empty()) {
        std::cerr << "No yield rows were produced." << std::endl;
        clear_flux_cache(flux_cache);
        return;
    }

    TString variation_csv = output_dir;
    if (!variation_csv.EndsWith("/")) variation_csv += "/";
    variation_csv += "expected_event_yields_by_variation.csv";
    TString fixed_fsi_csv = output_dir;
    if (!fixed_fsi_csv.EndsWith("/")) fixed_fsi_csv += "/";
    fixed_fsi_csv += "expected_event_yields_uncertainties.csv";
    TString all_variation_csv = output_dir;
    if (!all_variation_csv.EndsWith("/")) all_variation_csv += "/";
    all_variation_csv += "expected_event_yields_uncertainties_all_variations.csv";
    TString generator_envelope_fixed_fsi_csv = output_dir;
    if (!generator_envelope_fixed_fsi_csv.EndsWith("/")) generator_envelope_fixed_fsi_csv += "/";
    generator_envelope_fixed_fsi_csv += "expected_event_yields_generator_envelope_fixed_fsi.csv";
    TString generator_envelope_csv = output_dir;
    if (!generator_envelope_csv.EndsWith("/")) generator_envelope_csv += "/";
    generator_envelope_csv += "expected_event_yields_generator_envelope.csv";
    TString summary_md = output_dir;
    if (!summary_md.EndsWith("/")) summary_md += "/";
    summary_md += "expected_event_yields_summary.md";

    write_variation_rows(variation_csv, yields, exposure_scale,
                         proton_threshold_gev, piminus_threshold_gev);
    write_uncertainty_rows(fixed_fsi_csv, yields, true, exposure_scale,
                           proton_threshold_gev, piminus_threshold_gev);
    write_uncertainty_rows(all_variation_csv, yields, false, exposure_scale,
                           proton_threshold_gev, piminus_threshold_gev);
    write_generator_envelope_rows(generator_envelope_fixed_fsi_csv, yields, true, exposure_scale,
                                  proton_threshold_gev, piminus_threshold_gev);
    write_generator_envelope_rows(generator_envelope_csv, yields, false, exposure_scale,
                                  proton_threshold_gev, piminus_threshold_gev);
    write_summary_md(summary_md, status_csv, flux_file, proposal_emin_gev, proposal_emax_gev,
                     exposure_scale, proton_threshold_gev, piminus_threshold_gev,
                     samples_processed, int(yields.size()));

    std::cout << "Expected-yield samples processed: " << samples_processed
              << "\nPer-variation yields: " << variation_csv.Data()
              << "\nFixed-FSI uncertainty summary: " << fixed_fsi_csv.Data()
              << "\nAll-variation uncertainty summary: " << all_variation_csv.Data()
              << "\nFixed-FSI generator envelope: " << generator_envelope_fixed_fsi_csv.Data()
              << "\nAll-variation generator envelope: " << generator_envelope_csv.Data()
              << "\nNotes: " << summary_md.Data() << std::endl;

    clear_flux_cache(flux_cache);
}
