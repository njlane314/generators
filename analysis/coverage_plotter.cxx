#include <TAxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2D.h>
#include <TLeaf.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace AnaCoverage {

const int max_particles = 500;
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
const double units = 1.0e38;
const double containment_fractions[] = {0.70, 0.90, 0.99};
const int containment_line_styles[] = {1, 7, 3};
const int direct_decay_steps = 36;
const int cascade_decay_steps = 18;

struct VariableSpec {
    TString id;
    TString axis_title;
    int bins = 0;
    double min = 0.0;
    double max = 0.0;
};

struct PairSpec {
    TString id;
    TString title;
    VariableSpec x;
    VariableSpec y;
    bool valid = false;
};

struct InputSpec {
    TString path;
    TString label;
    TString generator;
    TString version;
    TString variation;
    TString knob;
    TString beam_mode;
    TString beam_species;
    TString interaction;
    TString fsi_state;
    TString sample;
};

struct FourVector {
    double px = 0.0;
    double py = 0.0;
    double pz = 0.0;
    double e = 0.0;
};

struct EventSummary {
    bool has_lambda = false;
    bool has_hyperon = false;
    double detector_branching = 0.0;
    double detector_visible = 0.0;
};

struct TreeInputs {
    TTree* tree = nullptr;
    TString tree_name;
    TLeaf* enu = nullptr;
    TLeaf* q2 = nullptr;
    TLeaf* w = nullptr;
    TLeaf* elep = nullptr;
    TLeaf* coslep = nullptr;
    TLeaf* pdglep = nullptr;
    TLeaf* weight = nullptr;
    TLeaf* scale = nullptr;
    TLeaf* target_a = nullptr;
    TLeaf* xsec_weight = nullptr;
    TLeaf* analysis_weight = nullptr;
    TLeaf* n_fsp = nullptr;
    TLeaf* pdg = nullptr;
    TLeaf* px = nullptr;
    TLeaf* py = nullptr;
    TLeaf* pz = nullptr;
    TLeaf* energy = nullptr;
    TLeaf* leading_lambda_p = nullptr;
    TLeaf* leading_lambda_costheta = nullptr;
    TLeaf* finalstate_lambda = nullptr;
    TLeaf* finalstate_lambda_count = nullptr;
    TLeaf* finalstate_strange_baryon_count = nullptr;
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

std::vector<TString> split_list(TString list)
{
    list.ReplaceAll(",", " ");
    TObjArray* tokens = list.Tokenize(" ");
    std::vector<TString> values;
    for (int i = 0; i < tokens->GetEntriesFast(); ++i) {
        TString value = static_cast<TObjString*>(tokens->At(i))->GetString();
        value = value.Strip(TString::kBoth);
        if (value != "") values.push_back(value);
    }
    delete tokens;
    return values;
}

std::vector<TString> split_comma_list(TString list)
{
    TObjArray* tokens = list.Tokenize(",");
    std::vector<TString> values;
    for (int i = 0; i < tokens->GetEntriesFast(); ++i) {
        TString value = static_cast<TObjString*>(tokens->At(i))->GetString();
        value = value.Strip(TString::kBoth);
        if (value != "") values.push_back(value);
    }
    delete tokens;
    return values;
}

std::vector<TString> split_root_file_list(TString list)
{
    std::vector<TString> values = split_comma_list(list);
    if (values.size() <= 1) values = split_list(list);
    return values;
}

bool has_token(const std::vector<TString>& values, const TString& target)
{
    return std::find(values.begin(), values.end(), target) != values.end();
}

TString lower(TString text)
{
    text.ToLower();
    return text;
}

TString join_path(TString dir, TString name)
{
    if (!dir.EndsWith("/")) dir += "/";
    return dir + name;
}

TString safe_name(TString name)
{
    name.ReplaceAll("/", "_");
    name.ReplaceAll(" ", "_");
    name.ReplaceAll(":", "_");
    name.ReplaceAll(";", "_");
    name.ReplaceAll(",", "_");
    name.ReplaceAll("|", "_");
    name.ReplaceAll("(", "_");
    name.ReplaceAll(")", "_");
    name.ReplaceAll("#", "");
    name.ReplaceAll("{", "");
    name.ReplaceAll("}", "");
    name.ReplaceAll("^", "");
    name.ReplaceAll("=", "");
    name.ReplaceAll("+", "plus");
    name.ReplaceAll("-", "minus");
    return name;
}

TString base_name(TString path)
{
    const Ssiz_t slash = path.Last('/');
    if (slash >= 0) path.Remove(0, slash + 1);
    path.ReplaceAll(".strange_taxonomy.root", "");
    path.ReplaceAll(".flat.root", "");
    path.ReplaceAll(".root", "");
    return path;
}

TLeaf* leaf(TTree* tree, std::initializer_list<const char*> names)
{
    for (const char* name : names) {
        if (TLeaf* found = tree->GetLeaf(name)) return found;
    }
    return nullptr;
}

double leaf_value(TLeaf* leaf, int index = 0, double fallback = 0.0)
{
    return leaf ? leaf->GetValue(index) : fallback;
}

int int_leaf_value(TLeaf* leaf, int index = 0, int fallback = 0)
{
    return leaf ? int(std::lround(leaf->GetValue(index))) : fallback;
}

double momentum(const FourVector& p4)
{
    return std::sqrt(p4.px * p4.px + p4.py * p4.py + p4.pz * p4.pz);
}

bool is_ana_final_hyperon(int pdg)
{
    return pdg == 3122 || pdg == 3212 || pdg == 3322 || pdg == 3312 || pdg == 3334;
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

bool is_kaon(int pdg)
{
    const int apdg = std::abs(pdg);
    return apdg == 321 || apdg == 311 || pdg == 310 || pdg == 130;
}

FourVector particle_p4(const TreeInputs& in, int i, double mass)
{
    FourVector p4;
    p4.px = leaf_value(in.px, i);
    p4.py = leaf_value(in.py, i);
    p4.pz = leaf_value(in.pz, i);
    const double p = momentum(p4);
    const double nominal_energy = std::sqrt(std::max(0.0, p * p + mass * mass));
    const double tree_energy = leaf_value(in.energy, i, nominal_energy);
    p4.e = tree_energy >= 0.999 * nominal_energy ? tree_energy : nominal_energy;
    return p4;
}

void combine_or_probability(double& event_probability, double contribution)
{
    if (!std::isfinite(contribution)) contribution = 0.0;
    contribution = std::max(0.0, std::min(1.0, contribution));
    event_probability = 1.0 - (1.0 - std::max(0.0, std::min(1.0, event_probability))) *
        (1.0 - contribution);
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
                                            double proton_threshold,
                                            double piminus_threshold,
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

VariableSpec variable_for_id(const TString& id)
{
    if (id == "enu") return {"enu", "E_{#nu} [GeV]", 80, 0.0, 10.0};
    if (id == "q2") return {"q2", "Q^{2} [GeV^{2}]", 80, 0.0, 5.0};
    if (id == "w") return {"w", "W [GeV]", 80, 0.0, 5.0};
    if (id == "lambda_p") return {"lambda_p", "p_{#Lambda} [GeV]", 80, 0.0, 4.0};
    if (id == "lambda_costheta") return {"lambda_costheta", "cos#theta_{#Lambda}", 60, -1.0, 1.0};
    return {"", "", 0, 0.0, 0.0};
}

PairSpec pair_for_id(const TString& id)
{
    if (id == "enu_q2") return {"enu_q2", "E_{#nu} vs Q^{2}", variable_for_id("enu"), variable_for_id("q2"), true};
    if (id == "enu_w") return {"enu_w", "E_{#nu} vs W", variable_for_id("enu"), variable_for_id("w"), true};
    if (id == "q2_w") return {"q2_w", "Q^{2} vs W", variable_for_id("q2"), variable_for_id("w"), true};
    if (id == "enu_lambda_p") return {"enu_lambda_p", "E_{#nu} vs p_{#Lambda}", variable_for_id("enu"), variable_for_id("lambda_p"), true};
    if (id == "w_lambda_p") return {"w_lambda_p", "W vs p_{#Lambda}", variable_for_id("w"), variable_for_id("lambda_p"), true};
    if (id == "lambda_p_costheta") {
        return {"lambda_p_costheta", "p_{#Lambda} vs cos#theta_{#Lambda}",
                variable_for_id("lambda_p"), variable_for_id("lambda_costheta"), true};
    }
    if (id == "q2_lambda_p") return {"q2_lambda_p", "Q^{2} vs p_{#Lambda}", variable_for_id("q2"), variable_for_id("lambda_p"), true};
    return {"", "", variable_for_id(""), variable_for_id(""), false};
}

std::vector<PairSpec> pair_specs_from_arg(TString pair_ids_arg)
{
    const TString default_pairs = "enu_q2 enu_w q2_w enu_lambda_p w_lambda_p lambda_p_costheta";
    std::vector<TString> ids = split_list(pair_ids_arg == "" ? default_pairs : pair_ids_arg);
    if (ids.empty() || has_token(ids, "all")) ids = split_list(default_pairs);

    std::vector<PairSpec> pairs;
    for (const TString& id : ids) {
        const PairSpec pair = pair_for_id(id);
        if (!pair.valid) {
            std::cerr << "Skipping unknown coverage pair: " << id.Data() << std::endl;
            continue;
        }
        pairs.push_back(pair);
    }
    return pairs;
}

int coverage_color(size_t index)
{
    const int colors[] = {
        TColor::GetColor("#D55E00"),
        TColor::GetColor("#0072B2"),
        TColor::GetColor("#009E73"),
        TColor::GetColor("#CC79A7"),
        TColor::GetColor("#E69F00"),
        TColor::GetColor("#56B4E9"),
        TColor::GetColor("#000000"),
        TColor::GetColor("#999999"),
    };
    return colors[index % (sizeof(colors) / sizeof(colors[0]))];
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
    if (beam_species == "" || beam_species == "combined" || beam_species == "both" || beam_species == "all") return "";
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

double containment_threshold(TH2D* hist, double fraction)
{
    std::vector<double> contents;
    contents.reserve(hist->GetNbinsX() * hist->GetNbinsY());
    double total = 0.0;
    for (int x = 1; x <= hist->GetNbinsX(); ++x) {
        for (int y = 1; y <= hist->GetNbinsY(); ++y) {
            const double value = hist->GetBinContent(x, y);
            if (value <= 0.0) continue;
            contents.push_back(value);
            total += value;
        }
    }
    if (total <= 0.0 || contents.empty()) return 0.0;
    std::sort(contents.begin(), contents.end(), std::greater<double>());
    const double target = fraction * total;
    double cumulative = 0.0;
    for (double value : contents) {
        cumulative += value;
        if (cumulative >= target) return value;
    }
    return contents.back();
}

double positive_bin_sum(TH2D* hist)
{
    double total = 0.0;
    for (int x = 1; x <= hist->GetNbinsX(); ++x) {
        for (int y = 1; y <= hist->GetNbinsY(); ++y) {
            total += std::max(0.0, hist->GetBinContent(x, y));
        }
    }
    return total;
}

TString selection_label(TString selection, double proton_threshold, double piminus_threshold)
{
    selection.ToLower();
    if (selection == "all") return "all events";
    if (selection == "final_lambda" || selection == "lambda") return "final-state #Lambda";
    if (selection == "detector_branching") return "#Lambda reachability #times B(#Lambda #rightarrow p#pi^{-})";
    if (selection == "detector_visible_lambda" || selection == "visible_lambda" ||
        selection == "detector_visible") {
        return Form("detector visible #Lambda #rightarrow p#pi^{-}: p_{p}>%.2f, p_{#pi^{-}}>%.2f GeV",
                    proton_threshold, piminus_threshold);
    }
    return "final-state #Lambda/#Sigma^{0}/#Xi^{0}/#Xi^{-}/#Omega^{-}";
}

TString visibility_stem(TString label, TString selection,
                        double proton_threshold, double piminus_threshold)
{
    TString stem = label == "" ? selection : label;
    selection = lower(selection);
    if (selection == "detector_visible_lambda" || selection == "visible_lambda" ||
        selection == "detector_visible") {
        stem += Form("_pp%.2f_pim%.2f", proton_threshold, piminus_threshold);
    }
    return safe_name(stem);
}

class CoveragePlotter {
public:
    explicit CoveragePlotter(TString output_dir = "analysis/output/coverage")
        : output_dir_(output_dir)
    {
        gROOT->SetBatch(true);
        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);
        gStyle->SetLabelFont(42, "xyz");
        gStyle->SetTitleFont(42, "xyz");
        gStyle->SetLegendFont(42);
    }

    void SetPairIds(TString pair_ids)
    {
        pairs_ = pair_specs_from_arg(pair_ids);
    }

    void SetSelection(TString selection)
    {
        selection_ = lower(selection);
        if (selection_ == "") selection_ = "final_hyperon";
    }

    void SetDetectorVisibility(double proton_threshold, double piminus_threshold)
    {
        proton_threshold_ = proton_threshold;
        piminus_threshold_ = piminus_threshold;
    }

    void DrawExplicitRootFiles(TString input_source, TString labels_arg, TString output_stem)
    {
        std::vector<InputSpec> specs;
        const std::vector<TString> paths = split_root_file_list(input_source);
        const std::vector<TString> labels = split_comma_list(labels_arg);
        for (size_t i = 0; i < paths.size(); ++i) {
            InputSpec spec;
            spec.path = paths[i];
            spec.label = i < labels.size() ? labels[i] : base_name(paths[i]);
            spec.sample = base_name(paths[i]);
            spec.beam_mode = infer_beam_mode(paths[i] + " " + spec.label);
            spec.beam_species = infer_beam_species(paths[i] + " " + spec.label);
            specs.push_back(spec);
        }
        draw_specs_for_all_pairs(specs, join_path(output_dir_, safe_name(output_stem)), output_stem);
    }

    void DrawStatusCsv(TString status_csv, TString groupings_arg)
    {
        const std::vector<InputSpec> rows = read_status_csv(status_csv);
        if (rows.empty()) {
            std::cerr << "No coverage inputs found in " << status_csv.Data() << std::endl;
            return;
        }
        std::vector<TString> groupings = split_list(groupings_arg);
        if (groupings.empty()) groupings = split_list("variation beam generator");
        if (has_token(groupings, "all")) groupings = split_list("variation beam generator");

        for (const TString& grouping : groupings) {
            std::map<std::string, std::vector<InputSpec>> grouped;
            for (const InputSpec& row : rows) grouped[group_key(row, grouping).Data()].push_back(row);
            const TString group_dir = join_path(output_dir_, grouping);
            gSystem->mkdir(group_dir, true);
            for (const auto& item : grouped) {
                if (item.second.empty()) continue;
                const TString key = item.first.c_str();
                const TString title = group_title(key, grouping);
                const TString output_base = join_path(group_dir, safe_name("coverage_" + grouping + "_" + key));
                std::vector<InputSpec> specs = item.second;
                for (InputSpec& spec : specs) spec.label = row_label(spec, grouping);
                draw_specs_for_all_pairs(specs, output_base, title);
            }
        }
    }

private:
    TString output_dir_;
    TString selection_ = "final_hyperon";
    double proton_threshold_ = 0.30;
    double piminus_threshold_ = 0.07;
    std::vector<PairSpec> pairs_ = pair_specs_from_arg("");

    TString selection_stem() const
    {
        return visibility_stem("", selection_, proton_threshold_, piminus_threshold_);
    }

    std::vector<InputSpec> read_status_csv(const TString& status_csv)
    {
        std::ifstream input(status_csv.Data());
        std::vector<InputSpec> rows;
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
            InputSpec row;
            row.generator = f[1].c_str();
            row.version = f[2].c_str();
            row.variation = f[3].c_str();
            row.knob = f[4].c_str();
            row.beam_mode = f[5].c_str();
            row.beam_species = f[6].c_str();
            row.interaction = f[7].c_str();
            row.fsi_state = f[8].c_str();
            row.sample = f[9].c_str();
            row.path = f[10].c_str();
            row.label = row_label(row, "generator");
            if (!gSystem->AccessPathName(row.path.Data())) rows.push_back(row);
        }
        return rows;
    }

    TString row_label(const InputSpec& row, const TString& grouping)
    {
        if (grouping == "variation") return row.variation;
        if (grouping == "beam") return row.beam_mode + " " + row.beam_species;
        if (grouping == "generator") return row.generator + " " + row.variation;
        return row.sample == "" ? base_name(row.path) : row.sample;
    }

    TString group_key(const InputSpec& row, const TString& grouping)
    {
        if (grouping == "variation") {
            return row.generator + "|" + row.beam_mode + "|" + row.beam_species + "|" + row.interaction + "|" + row.fsi_state;
        }
        if (grouping == "beam") {
            return row.generator + "|" + row.variation + "|" + row.interaction + "|" + row.fsi_state;
        }
        if (grouping == "generator") {
            return row.beam_mode + "|" + row.beam_species + "|" + row.interaction + "|" + row.fsi_state;
        }
        return row.generator;
    }

    TString group_title(TString key, const TString& grouping)
    {
        key.ReplaceAll("|", " ");
        if (grouping == "variation") return "Variation coverage: " + key;
        if (grouping == "beam") return "Beam coverage: " + key;
        if (grouping == "generator") return "Generator coverage: " + key;
        return "Coverage: " + key;
    }

    TTree* get_tree(TFile* file, TString& tree_name)
    {
        TTree* tree = nullptr;
        file->GetObject("FlatTree_VARS", tree);
        if (tree) {
            tree_name = "FlatTree_VARS";
            return tree;
        }
        file->GetObject("strange_taxonomy", tree);
        if (tree) {
            tree_name = "strange_taxonomy";
            return tree;
        }
        tree_name = "";
        return nullptr;
    }

    TreeInputs bind_tree(TTree* tree, const TString& tree_name)
    {
        TreeInputs in;
        in.tree = tree;
        in.tree_name = tree_name;
        in.enu = leaf(tree, {"Enu_true", "enu_true"});
        in.q2 = leaf(tree, {"Q2", "Q2_true", "q2", "q2_true"});
        in.w = leaf(tree, {"W", "W_true", "w", "w_true"});
        in.elep = leaf(tree, {"ELep", "e_lep", "E_lep"});
        in.coslep = leaf(tree, {"CosLep", "cos_lep", "cth_lep"});
        in.pdglep = leaf(tree, {"PDGLep", "pdg_lep"});
        in.weight = leaf(tree, {"Weight", "weight"});
        in.scale = leaf(tree, {"fScaleFactor", "scale_factor", "analysis_scale"});
        in.target_a = leaf(tree, {"tgta", "target_a"});
        in.xsec_weight = leaf(tree, {"xsec_weight_1e38_per_Ar"});
        in.analysis_weight = leaf(tree, {"analysis_weight"});
        in.n_fsp = leaf(tree, {"nfsp", "n_fsp"});
        in.pdg = leaf(tree, {"pdg"});
        in.px = leaf(tree, {"px"});
        in.py = leaf(tree, {"py"});
        in.pz = leaf(tree, {"pz"});
        in.energy = leaf(tree, {"E", "energy"});
        in.leading_lambda_p = leaf(tree, {"leading_Lambda_p"});
        in.leading_lambda_costheta = leaf(tree, {"leading_Lambda_costheta"});
        in.finalstate_lambda = leaf(tree, {"finalstate_Lambda"});
        in.finalstate_lambda_count = leaf(tree, {"finalstate_Lambda_count"});
        in.finalstate_strange_baryon_count = leaf(tree, {"finalstate_strange_baryon_count"});
        return in;
    }

    int final_count(const TreeInputs& in) const
    {
        return std::max(0, std::min(int_leaf_value(in.n_fsp), max_particles));
    }

    bool has_final_state_pdg(const TreeInputs& in) const
    {
        return in.n_fsp && in.pdg;
    }

    bool event_has_final_lambda(const TreeInputs& in) const
    {
        if (has_final_state_pdg(in)) {
            for (int i = 0; i < final_count(in); ++i) {
                if (int_leaf_value(in.pdg, i) == 3122) return true;
            }
            return false;
        }
        if (in.finalstate_lambda) return int_leaf_value(in.finalstate_lambda) != 0;
        if (in.finalstate_lambda_count) return int_leaf_value(in.finalstate_lambda_count) > 0;
        return true;
    }

    bool event_has_ana_final_hyperon(const TreeInputs& in) const
    {
        if (has_final_state_pdg(in)) {
            for (int i = 0; i < final_count(in); ++i) {
                if (is_ana_final_hyperon(int_leaf_value(in.pdg, i))) return true;
            }
            return false;
        }
        if (in.finalstate_strange_baryon_count) return int_leaf_value(in.finalstate_strange_baryon_count) > 0;
        if (in.finalstate_lambda_count) return int_leaf_value(in.finalstate_lambda_count) > 0;
        if (in.finalstate_lambda) return int_leaf_value(in.finalstate_lambda) != 0;
        return true;
    }

    EventSummary summarize_event(const TreeInputs& in) const
    {
        EventSummary summary;
        if (has_final_state_pdg(in)) {
            for (int i = 0; i < final_count(in); ++i) {
                const int code = int_leaf_value(in.pdg, i);
                if (!is_ana_final_hyperon(code)) continue;
                summary.has_hyperon = true;
                summary.has_lambda = summary.has_lambda || code == 3122;
                combine_or_probability(summary.detector_branching, lambda_to_p_piminus_br);
                if (!in.px || !in.py || !in.pz) continue;
                const double parent_mass = mass_gev(code);
                if (parent_mass <= 0.0) continue;
                combine_or_probability(summary.detector_visible,
                                       visible_feeddown_factor(code,
                                                               particle_p4(in, i, parent_mass),
                                                               proton_threshold_,
                                                               piminus_threshold_));
            }
            return summary;
        }

        summary.has_lambda = event_has_final_lambda(in);
        summary.has_hyperon = event_has_ana_final_hyperon(in);
        if (summary.has_hyperon) {
            summary.detector_branching = lambda_to_p_piminus_br;
        }
        double lambda_p = -1.0;
        double lambda_costheta = -2.0;
        if (summary.has_lambda && leading_lambda(in, lambda_p, lambda_costheta)) {
            const double lambda_e = std::sqrt(lambda_p * lambda_p + mass_lambda * mass_lambda);
            summary.detector_visible = lambda_to_p_piminus_br *
                lambda_daughter_threshold_acceptance(lambda_p, lambda_e,
                                                     proton_threshold_,
                                                     piminus_threshold_);
        }
        return summary;
    }

    double selection_factor(const TreeInputs& in) const
    {
        if (selection_ == "all") return 1.0;
        const EventSummary summary = summarize_event(in);
        if (selection_ == "final_lambda" || selection_ == "lambda") return summary.has_lambda ? 1.0 : 0.0;
        if (selection_ == "final_hyperon" || selection_ == "ana_final_hyperon") return summary.has_hyperon ? 1.0 : 0.0;
        if (selection_ == "detector_branching") return summary.detector_branching;
        if (selection_ == "detector_visible_lambda" || selection_ == "visible_lambda" ||
            selection_ == "detector_visible") {
            return summary.detector_visible;
        }
        return summary.has_hyperon ? 1.0 : 0.0;
    }

    bool leading_lambda(const TreeInputs& in, double& p, double& costheta) const
    {
        if (in.leading_lambda_p) {
            p = leaf_value(in.leading_lambda_p, 0, -1.0);
            costheta = in.leading_lambda_costheta ? leaf_value(in.leading_lambda_costheta, 0, -2.0) : -2.0;
            return p >= 0.0;
        }
        if (!in.n_fsp || !in.pdg || !in.px || !in.py || !in.pz) return false;
        p = -1.0;
        costheta = -2.0;
        for (int i = 0; i < final_count(in); ++i) {
            if (int_leaf_value(in.pdg, i) != 3122) continue;
            const double px = leaf_value(in.px, i);
            const double py = leaf_value(in.py, i);
            const double pz = leaf_value(in.pz, i);
            const double candidate = std::sqrt(px * px + py * py + pz * pz);
            if (candidate > p) {
                p = candidate;
                costheta = candidate > 0.0 ? pz / candidate : -2.0;
            }
        }
        return p >= 0.0;
    }

    bool value_for_variable(const TreeInputs& in, const VariableSpec& spec, double& value) const
    {
        if (spec.id == "enu") {
            if (!in.enu) return false;
            value = leaf_value(in.enu);
            return true;
        }
        if (spec.id == "q2") {
            if (in.q2) {
                value = leaf_value(in.q2);
                return true;
            }
            if (!in.enu || !in.elep || !in.coslep || !in.pdglep) return false;
            value = compute_q2(leaf_value(in.enu), leaf_value(in.elep), leaf_value(in.coslep), int_leaf_value(in.pdglep));
            return true;
        }
        if (spec.id == "w") {
            if (in.w) {
                value = leaf_value(in.w);
                return true;
            }
            double q2 = 0.0;
            if (!value_for_variable(in, variable_for_id("q2"), q2) || !in.enu || !in.elep) return false;
            value = compute_w(leaf_value(in.enu), leaf_value(in.elep), q2);
            return true;
        }
        if (spec.id == "lambda_p" || spec.id == "lambda_costheta") {
            double p = -1.0;
            double costheta = -2.0;
            if (!leading_lambda(in, p, costheta)) return false;
            value = spec.id == "lambda_p" ? p : costheta;
            return true;
        }
        return false;
    }

    double event_weight(const TreeInputs& in) const
    {
        if (in.xsec_weight) return leaf_value(in.xsec_weight, 0, 1.0);
        if (in.analysis_weight) return leaf_value(in.analysis_weight, 0, 1.0) * units * 40.0;
        const double weight = leaf_value(in.weight, 0, 1.0);
        const double scale = leaf_value(in.scale, 0, 1.0);
        const double target_a = std::max(1.0, leaf_value(in.target_a, 0, 40.0));
        return weight * scale * target_a * units;
    }

    TH2D* make_coverage_hist(const InputSpec& spec, const PairSpec& pair)
    {
        TFile* file = TFile::Open(spec.path, "READ");
        if (!file || file->IsZombie()) {
            std::cerr << "Could not open coverage input: " << spec.path.Data() << std::endl;
            if (file) file->Close();
            return nullptr;
        }

        TString tree_name;
        TTree* tree = get_tree(file, tree_name);
        if (!tree) {
            std::cerr << "Missing FlatTree_VARS or strange_taxonomy tree in " << spec.path.Data() << std::endl;
            file->Close();
            delete file;
            return nullptr;
        }

        const TreeInputs in = bind_tree(tree, tree_name);
        const TString hist_title = TString(";") + pair.x.axis_title + ";" + pair.y.axis_title + ";weighted yield";
        TH2D* hist = new TH2D(safe_name(spec.label + "_" + pair.id),
                              hist_title,
                              pair.x.bins, pair.x.min, pair.x.max,
                              pair.y.bins, pair.y.min, pair.y.max);
        hist->SetDirectory(nullptr);
        hist->SetStats(false);
        hist->Sumw2();

        const Long64_t entries = tree->GetEntries();
        for (Long64_t entry = 0; entry < entries; ++entry) {
            tree->GetEntry(entry);
            const double detector_factor = selection_factor(in);
            if (detector_factor <= 0.0) continue;

            double x = 0.0;
            double y = 0.0;
            if (!value_for_variable(in, pair.x, x) || !value_for_variable(in, pair.y, y)) continue;
            if (!std::isfinite(x) || !std::isfinite(y)) continue;
            if (x < pair.x.min || x >= pair.x.max || y < pair.y.min || y >= pair.y.max) continue;

            const double weight = event_weight(in) * detector_factor;
            if (!std::isfinite(weight) || weight <= 0.0) continue;
            hist->Fill(x, y, weight);
        }

        file->Close();
        delete file;

        if (positive_bin_sum(hist) <= 0.0) {
            delete hist;
            return nullptr;
        }
        return hist;
    }

    void draw_specs_for_all_pairs(const std::vector<InputSpec>& specs, const TString& output_base, const TString& title)
    {
        if (specs.empty()) return;
        for (const PairSpec& pair : pairs_) {
            draw_coverage_contours(specs, pair, output_base + "_" + pair.id + "_" + selection_stem(), title);
        }
    }

    void draw_coverage_contours(const std::vector<InputSpec>& specs,
                                const PairSpec& pair,
                                const TString& output_base,
                                const TString& title)
    {
        std::vector<TH2D*> hists;
        std::vector<TString> labels;
        for (const InputSpec& spec : specs) {
            TH2D* hist = make_coverage_hist(spec, pair);
            if (!hist) {
                std::cerr << "No coverage entries for " << spec.label.Data()
                          << " pair " << pair.id.Data() << " from " << spec.path.Data() << std::endl;
                continue;
            }
            hists.push_back(hist);
            labels.push_back(spec.label);
        }
        if (hists.empty()) {
            std::cerr << "No coverage histograms for " << output_base.Data() << std::endl;
            return;
        }

        TCanvas canvas("coverage_canvas", "coverage_canvas", 940, 720);
        canvas.SetLeftMargin(0.13);
        canvas.SetRightMargin(0.04);
        canvas.SetTopMargin(0.08);
        canvas.SetBottomMargin(0.13);

        TH2D frame("coverage_frame",
                   TString(";") + pair.x.axis_title + ";" + pair.y.axis_title,
                   pair.x.bins, pair.x.min, pair.x.max,
                   pair.y.bins, pair.y.min, pair.y.max);
        frame.SetDirectory(nullptr);
        frame.SetStats(false);
        frame.GetXaxis()->SetTitleSize(0.046);
        frame.GetYaxis()->SetTitleSize(0.046);
        frame.GetXaxis()->SetLabelSize(0.038);
        frame.GetYaxis()->SetLabelSize(0.038);
        frame.GetYaxis()->SetTitleOffset(1.15);
        frame.Draw("axis");

        std::vector<TH2D*> contours;
        for (size_t sample_index = 0; sample_index < hists.size(); ++sample_index) {
            const int color = coverage_color(sample_index);
            for (size_t level = 0; level < sizeof(containment_fractions) / sizeof(containment_fractions[0]); ++level) {
                const double threshold = containment_threshold(hists[sample_index], containment_fractions[level]);
                if (threshold <= 0.0) continue;
                TH2D* contour = static_cast<TH2D*>(hists[sample_index]->Clone(
                    safe_name(labels[sample_index] + "_" + pair.id + Form("_contour_%zu", level))));
                contour->SetDirectory(nullptr);
                contour->SetStats(false);
                contour->SetLineColor(color);
                contour->SetLineStyle(containment_line_styles[level]);
                contour->SetLineWidth(2);
                double levels[1] = {threshold};
                contour->SetContour(1, levels);
                contour->Draw("cont3 same");
                contours.push_back(contour);
            }
        }
        frame.Draw("axis same");

        TLegend containment_legend(0.15, 0.74, 0.37, 0.90);
        containment_legend.SetBorderSize(0);
        containment_legend.SetFillStyle(0);
        containment_legend.SetTextSize(0.034);
        std::vector<TLine*> containment_lines;
        for (size_t level = 0; level < sizeof(containment_fractions) / sizeof(containment_fractions[0]); ++level) {
            TLine* line = new TLine();
            line->SetLineColor(TColor::GetColor("#333333"));
            line->SetLineStyle(containment_line_styles[level]);
            line->SetLineWidth(2);
            containment_lines.push_back(line);
            containment_legend.AddEntry(line, Form("%.0f%% of yield", 100.0 * containment_fractions[level]), "l");
        }
        containment_legend.Draw();

        TLegend sample_legend(labels.size() > 3 ? 0.48 : 0.60, 0.70, 0.93, 0.90);
        sample_legend.SetBorderSize(0);
        sample_legend.SetFillStyle(0);
        sample_legend.SetTextSize(0.034);
        sample_legend.SetNColumns(labels.size() > 3 ? 2 : 1);
        std::vector<TLine*> sample_lines;
        for (size_t sample_index = 0; sample_index < labels.size(); ++sample_index) {
            TLine* line = new TLine();
            line->SetLineColor(coverage_color(sample_index));
            line->SetLineStyle(1);
            line->SetLineWidth(3);
            sample_lines.push_back(line);
            sample_legend.AddEntry(line, labels[sample_index], "l");
        }
        sample_legend.Draw();

        TLatex label;
        label.SetNDC(true);
        label.SetTextFont(42);
        label.SetTextSize(0.038);
        label.DrawLatex(0.13, 0.94,
                        title + ", " + selection_label(selection_,
                                                       proton_threshold_,
                                                       piminus_threshold_));

        TString dir = output_base;
        const Ssiz_t slash = dir.Last('/');
        if (slash >= 0) {
            dir.Remove(slash);
            gSystem->mkdir(dir, true);
        }
        canvas.SaveAs(output_base + ".pdf");
        canvas.SaveAs(output_base + ".png");

        for (TH2D* hist : hists) delete hist;
        for (TH2D* contour : contours) delete contour;
        for (TLine* line : containment_lines) delete line;
        for (TLine* line : sample_lines) delete line;
    }
};

}  // namespace AnaCoverage

void plot_kinematic_coverage(
    TString input_source = "analysis/output/matrix/generator_matrix_status.csv",
    TString output_dir = "analysis/output/coverage",
    TString labels_or_groupings = "variation beam generator",
    TString output_stem = "coverage",
    TString pair_ids = "enu_q2 enu_w q2_w enu_lambda_p w_lambda_p lambda_p_costheta",
    TString selection = "final_hyperon",
    double proton_threshold_gev = 0.30,
    double piminus_threshold_gev = 0.07
)
{
    AnaCoverage::CoveragePlotter plotter(output_dir);
    plotter.SetPairIds(pair_ids);
    plotter.SetSelection(selection);
    plotter.SetDetectorVisibility(proton_threshold_gev, piminus_threshold_gev);
    if (input_source.EndsWith(".csv")) {
        plotter.DrawStatusCsv(input_source, labels_or_groupings);
    } else {
        plotter.DrawExplicitRootFiles(input_source, labels_or_groupings, output_stem);
    }
}

void coverage_plotter(
    TString input_source = "analysis/output/matrix/generator_matrix_status.csv",
    TString output_dir = "analysis/output/coverage",
    TString labels_or_groupings = "variation beam generator",
    TString output_stem = "coverage",
    TString pair_ids = "enu_q2 enu_w q2_w enu_lambda_p w_lambda_p lambda_p_costheta",
    TString selection = "final_hyperon",
    double proton_threshold_gev = 0.30,
    double piminus_threshold_gev = 0.07
)
{
    plot_kinematic_coverage(input_source,
                            output_dir,
                            labels_or_groupings,
                            output_stem,
                            pair_ids,
                            selection,
                            proton_threshold_gev,
                            piminus_threshold_gev);
}

void plot_enu_q2_coverage(
    TString input_source = "analysis/output/matrix/generator_matrix_status.csv",
    TString output_dir = "analysis/output/coverage/enu_q2",
    TString labels_or_groupings = "variation beam generator",
    TString output_stem = "enu_q2_coverage",
    TString selection = "final_hyperon",
    double proton_threshold_gev = 0.30,
    double piminus_threshold_gev = 0.07
)
{
    plot_kinematic_coverage(input_source,
                            output_dir,
                            labels_or_groupings,
                            output_stem,
                            "enu_q2",
                            selection,
                            proton_threshold_gev,
                            piminus_threshold_gev);
}
