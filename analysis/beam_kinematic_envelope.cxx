#include <TAxis.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TKey.h>
#include <TLeaf.h>
#include <TLegend.h>
#include <TObject.h>
#include <TPad.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <initializer_list>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace AnaBeamEnvelope {

const int max_particles = 500;
const double units = 1.0e38;
const double default_argon_a = 40.0;
const double muon_mass_gev = 0.1056583745;
const double electron_mass_gev = 0.00051099895;
const double nucleon_mass_gev = 0.9395654133;
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
const int direct_decay_steps = 36;
const int cascade_decay_steps = 18;

struct StatusRow {
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
};

struct VariableSpec {
    TString id;
    TString axis_title;
    int bins = 0;
    double min = 0.0;
    double max = 0.0;
    bool valid = false;
};

struct FluxEnv {
    bool initialized = false;
    bool active = false;
    double emin = 0.0;
    double emax = 0.0;
    double proposal_emin = 0.0;
    double proposal_emax = 10.0;
    double mean_bin_flux = 1.0;
    std::vector<TH1*> hists;
};

struct TreeInputs {
    TLeaf* enu = nullptr;
    TLeaf* q2 = nullptr;
    TLeaf* w = nullptr;
    TLeaf* elep = nullptr;
    TLeaf* coslep = nullptr;
    TLeaf* pdglep = nullptr;
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
    bool has_hyperon = false;
    double detector_branching = 0.0;
    double detector_visible = 0.0;
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

std::vector<TString> split_list(TString text)
{
    text.ReplaceAll(",", " ");
    std::vector<TString> out;
    int start = 0;
    while (start < text.Length()) {
        while (start < text.Length() && std::isspace(static_cast<unsigned char>(text[start]))) ++start;
        int end = start;
        while (end < text.Length() && !std::isspace(static_cast<unsigned char>(text[end]))) ++end;
        if (end > start) out.push_back(text(start, end - start));
        start = end + 1;
    }
    return out;
}

TString lower(TString text)
{
    text.ToLower();
    return text;
}

TString csv(TString text)
{
    text.ReplaceAll("\"", "\"\"");
    TString out = "\"";
    out += text;
    out += "\"";
    return out;
}

TString safe_name(TString text)
{
    text.ReplaceAll("/", "_");
    text.ReplaceAll(" ", "_");
    text.ReplaceAll(":", "_");
    text.ReplaceAll(";", "_");
    text.ReplaceAll(",", "_");
    text.ReplaceAll("|", "_");
    text.ReplaceAll("(", "_");
    text.ReplaceAll(")", "_");
    text.ReplaceAll("#", "");
    text.ReplaceAll("{", "");
    text.ReplaceAll("}", "");
    text.ReplaceAll("^", "");
    text.ReplaceAll("+", "plus");
    text.ReplaceAll("-", "minus");
    return text;
}

TString join_path(TString dir, TString name)
{
    if (!dir.EndsWith("/")) dir += "/";
    return dir + name;
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

bool is_ana_hyperon(int pdg)
{
    return pdg == 3122 || pdg == 3212 || pdg == 3322 || pdg == 3312 || pdg == 3334;
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

void combine_or_probability(double& event_probability, double contribution)
{
    if (!std::isfinite(contribution)) contribution = 0.0;
    contribution = std::max(0.0, std::min(1.0, contribution));
    event_probability = 1.0 - (1.0 - event_probability) * (1.0 - contribution);
    event_probability = std::max(0.0, std::min(1.0, event_probability));
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
            lambda_daughter_threshold_acceptance(momentum(parent_p4), parent_p4.e,
                                                proton_threshold, piminus_threshold);
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
    for (int i = 0; i < final_count(in); ++i) {
        const int code = int_value(in.pdg, i);
        if (!is_ana_hyperon(code)) continue;
        summary.has_hyperon = true;
        const double parent_mass = mass_gev(code);
        if (parent_mass <= 0.0) continue;
        const FourVector parent = particle_p4(in, i, parent_mass);
        combine_or_probability(summary.detector_branching, lambda_to_p_piminus_br);
        combine_or_probability(summary.detector_visible,
                               visible_feeddown_factor(code, parent, proton_threshold, piminus_threshold));
    }
    return summary;
}

double selection_factor(const TreeInputs& in, TString selection,
                        double proton_threshold, double piminus_threshold)
{
    selection = lower(selection);
    if (selection == "" || selection == "detector_visible") selection = "detector_visible_lambda";
    const EventSummary summary = summarize_event(in, proton_threshold, piminus_threshold);
    if (selection == "all") return 1.0;
    if (selection == "final_hyperon" || selection == "ana_final_hyperon") return summary.has_hyperon ? 1.0 : 0.0;
    if (selection == "detector_branching") return summary.detector_branching;
    if (selection == "detector_visible_lambda" || selection == "visible_lambda") return summary.detector_visible;
    return summary.detector_visible;
}

double lepton_mass_gev(int pdg)
{
    const int apdg = std::abs(pdg);
    if (apdg == 13) return muon_mass_gev;
    if (apdg == 11) return electron_mass_gev;
    return 0.0;
}

double compute_q2(double enu, double elep, double coslep, int pdglep)
{
    const double mass = lepton_mass_gev(pdglep);
    const double p2 = std::max(0.0, elep * elep - mass * mass);
    const double plep = std::sqrt(p2);
    return 2.0 * enu * (elep - plep * coslep) - mass * mass;
}

double compute_w(double enu, double elep, double q2)
{
    const double w2 = nucleon_mass_gev * nucleon_mass_gev +
        2.0 * nucleon_mass_gev * (enu - elep) - q2;
    return w2 > 0.0 ? std::sqrt(w2) : -1.0;
}

double leading_hyperon_p(const TreeInputs& in)
{
    double best = -1.0;
    for (int i = 0; i < final_count(in); ++i) {
        if (!is_ana_hyperon(int_value(in.pdg, i))) continue;
        const FourVector p4 = particle_p4(in, i, mass_gev(int_value(in.pdg, i)));
        best = std::max(best, momentum(p4));
    }
    return best;
}

VariableSpec variable_for_id(TString id)
{
    id = lower(id);
    if (id == "enu") return {"enu", "E_{#nu} [GeV]", 80, 0.0, 10.0, true};
    if (id == "q2") return {"q2", "Q^{2} [GeV^{2}]", 80, 0.0, 5.0, true};
    if (id == "w") return {"w", "W [GeV]", 80, 0.0, 5.0, true};
    if (id == "hyperon_p") return {"hyperon_p", "leading selected hyperon p [GeV]", 80, 0.0, 4.0, true};
    return {"", "", 0, 0.0, 0.0, false};
}

std::vector<VariableSpec> variables_from_arg(TString arg)
{
    if (arg == "" || lower(arg) == "all") arg = "enu q2 w hyperon_p";
    std::vector<VariableSpec> out;
    for (const TString& id : split_list(arg)) {
        const VariableSpec spec = variable_for_id(id);
        if (spec.valid) out.push_back(spec);
        else std::cerr << "Skipping unknown beam-envelope variable: " << id.Data() << std::endl;
    }
    return out;
}

bool value_for_variable(const TreeInputs& in, const VariableSpec& spec, double& out)
{
    if (spec.id == "enu") {
        if (!in.enu) return false;
        out = value(in.enu);
        return true;
    }
    if (spec.id == "q2") {
        if (in.q2) {
            out = value(in.q2);
            return true;
        }
        if (!in.enu || !in.elep || !in.coslep || !in.pdglep) return false;
        out = compute_q2(value(in.enu), value(in.elep), value(in.coslep), int_value(in.pdglep));
        return true;
    }
    if (spec.id == "w") {
        if (in.w) {
            out = value(in.w);
            return true;
        }
        double q2 = 0.0;
        if (!value_for_variable(in, variable_for_id("q2"), q2) || !in.enu || !in.elep) return false;
        out = compute_w(value(in.enu), value(in.elep), q2);
        return true;
    }
    if (spec.id == "hyperon_p") {
        out = leading_hyperon_p(in);
        return out >= 0.0;
    }
    return false;
}

TString normalize_beam_mode(TString beam_mode)
{
    beam_mode.ToLower();
    if (beam_mode.Contains("rhc")) return "rhc";
    if (beam_mode.Contains("fhc")) return "fhc";
    return "";
}

TString normalize_beam_species(TString beam_species)
{
    beam_species.ToLower();
    if (beam_species == "numu" || beam_species == "nu_mu" || beam_species == "14") return "numu";
    if (beam_species == "numubar" || beam_species == "anti_numu" || beam_species == "-14") return "numubar";
    return beam_species;
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

bool add_flux_hist(TH1* hist, FluxEnv& env, double floor_fraction)
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
    TH1* clone = static_cast<TH1*>(hist->Clone(TString::Format("%s_beam_envelope_flux_%d", hist->GetName(), int(env.hists.size()))));
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

    TFile* file = TFile::Open(path, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Could not open NuMI flux file: " << path.Data() << std::endl;
        if (file) file->Close();
        return env;
    }

    TH1* hist = dynamic_cast<TH1*>(file->Get(flux_hist_path(beam_mode, beam_species)));
    if (!hist) hist = find_hist(file, hist_leaf_name(flux_hist_path(beam_mode, beam_species)));
    add_flux_hist(hist, env, floor_fraction);
    file->Close();
    delete file;

    env.emin = std::max(env.emin, env.proposal_emin);
    env.emax = std::min(env.emax, env.proposal_emax);
    env.mean_bin_flux = mean_flux_in_proposal(env);
    env.active = !env.hists.empty() && env.mean_bin_flux > 0.0 && env.emax > env.emin;
    if (!env.active) {
        std::cerr << "No valid NuMI flux envelope for " << beam_mode.Data()
                  << " " << beam_species.Data() << std::endl;
    }
    return env;
}

FluxEnv& flux_for(const StatusRow& row, const TString& flux_file, double floor_fraction,
                  double proposal_emin, double proposal_emax,
                  std::map<std::string, FluxEnv>& cache)
{
    const TString beam_mode = normalize_beam_mode(row.beam_mode);
    const TString beam_species = normalize_beam_species(row.beam_species);
    const std::string key = std::string(flux_file.Data()) + "|" + beam_mode.Data() + "|" + beam_species.Data();
    FluxEnv& env = cache[key];
    if (!env.initialized) env = read_flux(flux_file, beam_mode, beam_species, floor_fraction, proposal_emin, proposal_emax);
    return env;
}

double flux_weight(const FluxEnv& env, double enu)
{
    if (!env.active || env.mean_bin_flux <= 0.0) return 0.0;
    if (enu < env.emin || enu >= env.emax) return 0.0;
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

TreeInputs bind_inputs(TTree* tree)
{
    TreeInputs in;
    in.enu = leaf(tree, {"Enu_true", "enu_true"});
    in.q2 = leaf(tree, {"Q2", "Q2_true", "q2", "q2_true"});
    in.w = leaf(tree, {"W", "W_true", "w", "w_true"});
    in.elep = leaf(tree, {"ELep", "e_lep", "E_lep"});
    in.coslep = leaf(tree, {"CosLep", "cos_lep", "cth_lep"});
    in.pdglep = leaf(tree, {"PDGLep", "pdg_lep"});
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

double event_weight(const TreeInputs& in)
{
    if (in.xsec_weight) return value(in.xsec_weight, 0, 1.0);
    const double analysis_weight = in.analysis_weight ? value(in.analysis_weight, 0, 1.0) :
        value(in.weight, 0, 1.0) * value(in.scale, 0, 1.0);
    const double target_a = std::max(1.0, value(in.target_a, 0, default_argon_a));
    return analysis_weight * target_a * units;
}

std::vector<StatusRow> read_status(const TString& status_csv)
{
    std::ifstream input(status_csv.Data());
    std::vector<StatusRow> rows;
    if (!input) {
        std::cerr << "Could not open generator matrix status CSV: " << status_csv.Data() << std::endl;
        return rows;
    }
    std::string line;
    bool first = true;
    while (std::getline(input, line)) {
        if (first) {
            first = false;
            continue;
        }
        const std::vector<std::string> f = split_csv_line(line);
        if (f.size() < 13) continue;
        StatusRow row;
        row.generator = f[1].c_str();
        row.version = f[2].c_str();
        row.variation = f[3].c_str();
        row.knob = f[4].c_str();
        row.beam_mode = normalize_beam_mode(f[5].c_str());
        row.beam_species = normalize_beam_species(f[6].c_str());
        row.interaction = f[7].c_str();
        row.fsi_state = f[8].c_str();
        row.sample = f[9].c_str();
        row.input_file = f[10].c_str();
        if ((row.beam_mode == "fhc" || row.beam_mode == "rhc") &&
            !gSystem->AccessPathName(row.input_file.Data())) rows.push_back(row);
    }
    return rows;
}

TString group_key(const StatusRow& row)
{
    return row.generator + "|" + row.version + "|" + row.variation + "|" + row.knob + "|" +
        row.beam_species + "|" + row.interaction + "|" + row.fsi_state;
}

TString group_title(TString key)
{
    key.ReplaceAll("|", " ");
    return key;
}

bool fill_hist(const StatusRow& row, TH1D* hist, const VariableSpec& variable,
               const TString& selection, const TString& flux_file, double flux_floor_fraction,
               double proposal_emin, double proposal_emax, double proton_threshold,
               double piminus_threshold, std::map<std::string, FluxEnv>& flux_cache)
{
    TFile* file = TFile::Open(row.input_file, "READ");
    if (!file || file->IsZombie()) {
        if (file) file->Close();
        return false;
    }
    TTree* tree = nullptr;
    file->GetObject("FlatTree_VARS", tree);
    if (!tree) {
        file->Close();
        delete file;
        return false;
    }
    const TreeInputs in = bind_inputs(tree);
    if (!in.enu || !in.n_fsp || !in.pdg || !in.px || !in.py || !in.pz) {
        std::cerr << "Missing required kinematic beam-envelope branches in "
                  << row.input_file.Data() << std::endl;
        file->Close();
        delete file;
        return false;
    }
    FluxEnv& flux = flux_for(row, flux_file, flux_floor_fraction, proposal_emin, proposal_emax, flux_cache);
    if (!flux.active) {
        file->Close();
        delete file;
        return false;
    }

    const Long64_t entries = tree->GetEntries();
    for (Long64_t entry = 0; entry < entries; ++entry) {
        tree->GetEntry(entry);
        const double factor = selection_factor(in, selection, proton_threshold, piminus_threshold);
        if (factor <= 0.0) continue;
        double x = 0.0;
        if (!value_for_variable(in, variable, x)) continue;
        if (!std::isfinite(x) || x < variable.min || x >= variable.max) continue;
        const double flux_w = flux_weight(flux, value(in.enu));
        const double weight = event_weight(in) * flux_w * factor;
        if (!std::isfinite(weight) || weight <= 0.0) continue;
        hist->Fill(x, weight);
    }
    file->Close();
    delete file;
    return true;
}

void write_bin_csv(const TString& path, const TString& key, const VariableSpec& variable,
                   const TString& selection, TH1D* fhc, TH1D* rhc)
{
    std::ofstream out(path.Data());
    if (!out) return;
    out << "group,selection,variable,bin,low,high,fhc_yield,fhc_stat,"
        << "rhc_yield,rhc_stat,envelope_min,envelope_max,rhc_over_fhc\n";
    out << std::setprecision(12);
    for (int bin = 1; bin <= fhc->GetNbinsX(); ++bin) {
        const double f = fhc->GetBinContent(bin);
        const double r = rhc->GetBinContent(bin);
        out << csv(key).Data() << ',' << csv(selection).Data() << ','
            << csv(variable.id).Data() << ',' << bin << ','
            << fhc->GetXaxis()->GetBinLowEdge(bin) << ','
            << fhc->GetXaxis()->GetBinUpEdge(bin) << ','
            << f << ',' << fhc->GetBinError(bin) << ','
            << r << ',' << rhc->GetBinError(bin) << ','
            << std::min(f, r) << ',' << std::max(f, r) << ','
            << (f > 0.0 ? r / f : 0.0) << '\n';
    }
}

TGraphAsymmErrors* make_band(TH1D* a, TH1D* b)
{
    TGraphAsymmErrors* band = new TGraphAsymmErrors(a->GetNbinsX());
    band->SetFillColorAlpha(kGray + 1, 0.30);
    band->SetLineColor(kGray + 1);
    for (int bin = 1; bin <= a->GetNbinsX(); ++bin) {
        const double low = std::min(a->GetBinContent(bin), b->GetBinContent(bin));
        const double high = std::max(a->GetBinContent(bin), b->GetBinContent(bin));
        const double center = 0.5 * (low + high);
        const double x = a->GetXaxis()->GetBinCenter(bin);
        const double ex = 0.5 * a->GetXaxis()->GetBinWidth(bin);
        band->SetPoint(bin - 1, x, center);
        band->SetPointError(bin - 1, ex, ex, center - low, high - center);
    }
    return band;
}

void style_hist(TH1D* hist, int color)
{
    hist->SetStats(false);
    hist->SetLineColor(color);
    hist->SetMarkerColor(color);
    hist->SetLineWidth(3);
    hist->SetMarkerStyle(20);
}

void draw_plot(const TString& path_base, const TString& key, const VariableSpec& variable,
               const TString& selection, TH1D* fhc, TH1D* rhc)
{
    TString dir = gSystem->DirName(path_base.Data());
    gSystem->mkdir(dir.Data(), true);
    TCanvas canvas("beam_envelope_canvas", "beam_envelope_canvas", 980, 720);
    canvas.SetLeftMargin(0.14);
    canvas.SetBottomMargin(0.14);
    canvas.SetTopMargin(0.08);
    canvas.SetRightMargin(0.05);

    style_hist(fhc, kBlue + 1);
    style_hist(rhc, kRed + 1);
    TH1D* frame = static_cast<TH1D*>(fhc->Clone("beam_envelope_frame"));
    frame->Reset();
    frame->SetDirectory(nullptr);
    frame->SetStats(false);
    TString title = group_title(key);
    title += ";";
    title += variable.axis_title;
    title += ";weighted yield [10^{-38} cm^{2}/Ar]";
    frame->SetTitle(title);
    const double ymax = std::max(fhc->GetMaximum(), rhc->GetMaximum());
    frame->GetYaxis()->SetRangeUser(0.0, ymax > 0.0 ? 1.35 * ymax : 1.0);
    frame->Draw("axis");

    TGraphAsymmErrors* band = make_band(fhc, rhc);
    band->Draw("2 same");
    fhc->Draw("hist e same");
    rhc->Draw("hist e same");
    frame->Draw("axis same");

    TLegend legend(0.58, 0.72, 0.92, 0.90);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetHeader(selection.Data());
    legend.AddEntry(fhc, "FHC reweight", "l");
    legend.AddEntry(rhc, "RHC reweight", "l");
    legend.AddEntry(band, "FHC/RHC envelope", "f");
    legend.Draw();

    canvas.SaveAs(path_base + ".pdf");
    canvas.SaveAs(path_base + ".png");
    delete band;
    delete frame;
}

}  // namespace AnaBeamEnvelope

void beam_kinematic_envelope(
    TString status_csv = "analysis/output/matrix/generator_matrix_status.csv",
    TString output_dir = "analysis/output/beam_kinematic_envelope",
    TString variables = "enu q2 w hyperon_p",
    TString selection = "detector_visible_lambda",
    TString flux_file = "example/numi/flux/microboone_numi_flux_5mev.root",
    double flux_floor_fraction = 0.0,
    double proposal_emin_gev = 0.0,
    double proposal_emax_gev = 10.0,
    double proton_threshold_gev = 0.30,
    double piminus_threshold_gev = 0.07
)
{
    using namespace AnaBeamEnvelope;
    gStyle->SetOptStat(0);

    const std::vector<StatusRow> rows = read_status(status_csv);
    const std::vector<VariableSpec> specs = variables_from_arg(variables);
    if (rows.empty() || specs.empty()) {
        std::cerr << "No beam-envelope rows or variables found." << std::endl;
        return;
    }
    gSystem->mkdir(output_dir.Data(), true);

    std::map<std::string, std::vector<StatusRow>> groups;
    for (const StatusRow& row : rows) groups[group_key(row).Data()].push_back(row);

    std::map<std::string, FluxEnv> flux_cache;
    std::ofstream summary(join_path(output_dir, "beam_kinematic_envelope_summary.csv").Data());
    summary << "group,selection,variable,fhc_integral,rhc_integral,envelope_min_integral,"
            << "envelope_max_integral,rhc_over_fhc,plot_base,csv_file\n";
    summary << std::setprecision(12);

    int outputs = 0;
    for (const auto& item : groups) {
        const TString key = item.first.c_str();
        for (const VariableSpec& variable : specs) {
            TString fhc_name = "fhc_";
            fhc_name += key;
            fhc_name += "_";
            fhc_name += variable.id;
            TString rhc_name = "rhc_";
            rhc_name += key;
            rhc_name += "_";
            rhc_name += variable.id;
            TString hist_title = ";";
            hist_title += variable.axis_title;
            hist_title += ";weighted yield";
            TH1D* fhc = new TH1D(safe_name(fhc_name),
                                 hist_title,
                                 variable.bins, variable.min, variable.max);
            TH1D* rhc = new TH1D(safe_name(rhc_name),
                                 hist_title,
                                 variable.bins, variable.min, variable.max);
            fhc->SetDirectory(nullptr);
            rhc->SetDirectory(nullptr);
            fhc->Sumw2();
            rhc->Sumw2();

            for (const StatusRow& row : item.second) {
                TH1D* target = row.beam_mode == "fhc" ? fhc : rhc;
                fill_hist(row, target, variable, selection, flux_file, flux_floor_fraction,
                          proposal_emin_gev, proposal_emax_gev, proton_threshold_gev,
                          piminus_threshold_gev, flux_cache);
            }

            if (fhc->Integral() <= 0.0 && rhc->Integral() <= 0.0) {
                delete fhc;
                delete rhc;
                continue;
            }

            TString base_name = "beam_envelope_";
            base_name += key;
            base_name += "_";
            base_name += variable.id;
            base_name += "_";
            base_name += selection;
            const TString base = join_path(output_dir, safe_name(base_name));
            const TString csv_path = base + ".csv";
            write_bin_csv(csv_path, key, variable, selection, fhc, rhc);
            draw_plot(base, key, variable, selection, fhc, rhc);

            const double fhc_int = fhc->Integral();
            const double rhc_int = rhc->Integral();
            summary << csv(key).Data() << ',' << csv(selection).Data() << ','
                    << csv(variable.id).Data() << ','
                    << fhc_int << ',' << rhc_int << ','
                    << std::min(fhc_int, rhc_int) << ',' << std::max(fhc_int, rhc_int) << ','
                    << (fhc_int > 0.0 ? rhc_int / fhc_int : 0.0) << ','
                    << csv(base).Data() << ',' << csv(csv_path).Data() << '\n';
            ++outputs;
            delete fhc;
            delete rhc;
        }
    }

    clear_flux_cache(flux_cache);
    std::cout << "Wrote " << outputs << " beam kinematic envelope outputs under "
              << output_dir.Data() << std::endl;
}
