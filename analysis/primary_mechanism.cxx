#include <TAxis.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TKey.h>
#include <TLeaf.h>
#include <TNamed.h>
#include <TObject.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <initializer_list>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

namespace {

const int max_particles = 100;
const double units = 1.0e38;
const double default_argon_a = 40.0;
const double muon_mass_gev = 0.1056583745;

struct acc {
    Long64_t raw = 0;
    double weight = 0.0;
    double weight_sum_squares = 0.0;
    double xsec = 0.0;
    double xsec_sum_squares = 0.0;
};

struct stage {
    TLeaf* n = nullptr;
    TLeaf* pdg = nullptr;
    TLeaf* px = nullptr;
    TLeaf* py = nullptr;
    TLeaf* pz = nullptr;
    TLeaf* energy = nullptr;
};

struct topology {
    TLeaf* flag = nullptr;
    TLeaf* fiducial = nullptr;
    TLeaf* lambda_active = nullptr;
    TLeaf* proton = nullptr;
    TLeaf* pion = nullptr;
    TLeaf* detached = nullptr;
    TString source = "flat-tree final Lambda + visible p pi- proxy";
};

struct flux_envelope {
    bool active = false;
    double emin = 0.0;
    double emax = 0.0;
    TString source;
};

struct primary_counts {
    int lambda = 0;
    int sigma_0 = 0;
    int charged_sigma = 0;
    int kaon = 0;
    int strange_baryon = 0;
    int xi_or_omega = 0;
    int pion = 0;
};

std::vector<TString> origin_labels()
{
    return {"O0_direct_Lambda", "O1_Sigma0_feedthrough", "O2_charged_Sigma_exchange",
            "O3_neutral_Sigma_exchange", "O4_associated_KY",
            "O5_inelastic_hyperon_or_Ypi", "O6_other_strange", "O7_no_true_Lambda"};
}

std::vector<TString> origin_defs()
{
    return {"primary Lambda at the interaction vertex",
            "primary Sigma0 at the interaction vertex",
            "primary Sigma+ or Sigma- at the interaction vertex",
            "reserved until nuclear-exit ancestry is available",
            "primary kaon plus strange baryon at the interaction vertex",
            "primary Xi/Omega or other hyperon/Ypi-like strange ancestry",
            "selected Lambda topology with other or incomplete strange ancestry",
            "selected topology with no true-Lambda evidence"};
}

std::vector<TString> mode_labels()
{
    return {"QE", "MEC", "RES", "DIS", "COH", "diffractive", "other", "missing_mode"};
}

std::vector<TString> mode_defs()
{
    return {"abs(Mode) == 1", "abs(Mode) == 2", "abs(Mode) in {10,11,12,13,17,22,23}",
            "abs(Mode) in {21,26}", "abs(Mode) == 16", "abs(Mode) == 15",
            "Mode present but outside the standard grouping", "no Mode/mode branch"};
}

std::vector<TString> category_labels()
{
    return {"S0_inclusive", "S1_CC_muon", "S2_visible_kaon", "S3_low_extra",
            "S4_visible_EM", "S5_no_visible_muon"};
}

std::vector<TString> category_defs()
{
    return {"inclusive selected Lambda topology", "S0 plus charged-current visible muon tag",
            "S0 plus visible final-state kaon", "S0 plus extra visible kinetic energy below 0.15 GeV",
            "S0 plus visible final-state gamma/EM tag", "S0 without a visible muon tag"};
}

TLeaf* leaf(TTree* tree, std::initializer_list<TString> names)
{
    for (const TString& name : names) {
        if (TLeaf* found = tree->GetLeaf(name.Data())) return found;
    }
    return nullptr;
}

double value(TLeaf* leaf, int i = 0, double fallback = 0.0)
{
    return leaf ? leaf->GetValue(i) : fallback;
}

int int_value(TLeaf* leaf, int i = 0, int fallback = 0)
{
    return leaf ? int(std::lround(leaf->GetValue(i))) : fallback;
}

int count(const stage& s)
{
    return std::max(0, std::min(int_value(s.n), max_particles));
}

void add(acc& total, double weight, double xsec)
{
    ++total.raw;
    total.weight += weight;
    total.weight_sum_squares += weight * weight;
    total.xsec += xsec;
    total.xsec_sum_squares += xsec * xsec;
}

TString csv(TString text)
{
    text.ReplaceAll("\"", "\"\"");
    TString out = "\"";
    out += text;
    out += "\"";
    return out;
}

TString clean(TString text)
{
    text.ReplaceAll(".flat.root", "");
    text.ReplaceAll(".root", "");
    text.ReplaceAll("/", "_");
    text.ReplaceAll(" ", "_");
    text.ReplaceAll(":", "_");
    text.ReplaceAll(";", "_");
    return text;
}

TString infer_generator(TString sample)
{
    sample.ToLower();
    if (sample.Contains("nuwro")) return "NuWro";
    if (sample.Contains("genie")) return "GENIE";
    if (sample.Contains("gibuu")) return "GiBUU";
    return "unspecified";
}

TString infer_knob(TString sample)
{
    sample.ReplaceAll(".flat.root", "");
    sample.ReplaceAll(".root", "");
    std::vector<TString> known = {"AR23_20i_00_000", "G18_10a_02_11a", "G18_10b_02_11a",
                                  "G18_10a_02_11b", "G18_10b_02_11b", "G18_10c_02_11b",
                                  "G18_10d_02_11b", "hyp_lambda_only", "hyp_sigma0_only",
                                  "hyp_sigmam_only", "hyp_no_effmass", "all_strange",
                                  "dis_only", "hyp_all", "fsi_on", "fsi_off"};
    TString knob = "";
    for (const TString& name : known) {
        if (!sample.Contains(name)) continue;
        if (knob != "") knob += "+";
        knob += name;
    }
    return knob == "" ? "nominal" : knob;
}

void label(TAxis* axis, const std::vector<TString>& names)
{
    for (int i = 0; i < int(names.size()); ++i) axis->SetBinLabel(i + 1, names[i]);
}

bool is_kaon(int pdg)
{
    const int a = std::abs(pdg);
    return a == 321 || a == 311 || pdg == 310 || pdg == 130;
}

bool is_sigma(int pdg)
{
    const int a = std::abs(pdg);
    return a == 3112 || a == 3212 || a == 3222;
}

bool is_xi_or_omega(int pdg)
{
    const int a = std::abs(pdg);
    return a == 3312 || a == 3322 || a == 3334;
}

bool is_strange_baryon(int pdg)
{
    return std::abs(pdg) == 3122 || is_sigma(pdg) || is_xi_or_omega(pdg);
}

int species_bin(int pdg)
{
    if (pdg == 3122) return 0;
    if (pdg == 3112) return 1;
    if (pdg == 3212) return 2;
    if (pdg == 3222) return 3;
    if (std::abs(pdg) == 3312 || std::abs(pdg) == 3322) return 4;
    if (std::abs(pdg) == 3334) return 5;
    if (is_kaon(pdg)) return 6;
    if (is_strange_baryon(pdg)) return 7;
    return -1;
}

double mass_gev(int pdg)
{
    const int a = std::abs(pdg);
    if (a == 2212) return 0.9382720813;
    if (a == 2112) return 0.9395654133;
    if (a == 211) return 0.13957039;
    if (a == 321) return 0.493677;
    if (a == 311 || a == 310 || a == 130) return 0.497611;
    if (a == 3122) return 1.115683;
    if (a == 13) return muon_mass_gev;
    if (a == 11) return 0.00051099895;
    return 0.0;
}

double momentum(const stage& s, int i)
{
    const double x = value(s.px, i);
    const double y = value(s.py, i);
    const double z = value(s.pz, i);
    return std::sqrt(x * x + y * y + z * z);
}

primary_counts read_primary(const stage& s)
{
    primary_counts out;
    for (int i = 0; i < count(s); ++i) {
        const int code = int_value(s.pdg, i);
        out.lambda += code == 3122;
        out.sigma_0 += code == 3212;
        out.charged_sigma += code == 3112 || code == 3222;
        out.kaon += is_kaon(code);
        out.strange_baryon += is_strange_baryon(code);
        out.xi_or_omega += is_xi_or_omega(code);
        out.pion += std::abs(code) == 211 || code == 111;
    }
    return out;
}

int origin_bin(const primary_counts& p, bool true_lambda)
{
    if (p.kaon && p.strange_baryon) return 4;
    if (p.sigma_0) return 1;
    if (p.charged_sigma) return 2;
    if (p.lambda) return 0;
    if (p.xi_or_omega || (p.strange_baryon && p.pion) || p.strange_baryon) return 5;
    if (p.kaon || true_lambda) return 6;
    return 7;
}

int mode_bin(int mode, bool has_mode)
{
    if (!has_mode) return 7;
    const int m = std::abs(mode);
    if (m == 1) return 0;
    if (m == 2) return 1;
    if (m == 10 || m == 11 || m == 12 || m == 13 || m == 17 || m == 22 || m == 23) return 2;
    if (m == 21 || m == 26) return 3;
    if (m == 16) return 4;
    if (m == 15) return 5;
    return 6;
}

TH1* find_hist(TDirectory* dir, const TString& name)
{
    if (!dir) return nullptr;
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

bool extend_flux(TH1* hist, double floor_fraction, flux_envelope& env)
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
    env.emin = env.active ? std::min(env.emin, emin) : emin;
    env.emax = env.active ? std::max(env.emax, emax) : emax;
    env.active = true;
    return true;
}

flux_envelope read_flux(TString path, double floor_fraction)
{
    flux_envelope env;
    TFile* file = TFile::Open(path, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Could not open NuMI flux file: " << path.Data() << std::endl;
        return env;
    }
    std::vector<TString> used;
    const std::vector<TString> names = {"numu_CV_AV_TPC_5MeV_bin", "numubar_CV_AV_TPC_5MeV_bin"};
    for (const TString& name : names) {
        if (extend_flux(find_hist(file, name), floor_fraction, env)) used.push_back(name);
        else std::cerr << "Could not use NuMI flux histogram: " << name.Data() << std::endl;
    }
    file->Close();
    if (!env.active || env.emax <= env.emin) {
        env.active = false;
        std::cerr << "No valid non-zero NuMI flux support was found in " << path.Data() << std::endl;
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

bool in_flux(const flux_envelope& env, double enu)
{
    return env.active && enu >= env.emin && enu <= env.emax;
}

topology find_topology(TTree* tree, TString wp)
{
    wp.ToLower();
    topology out;
    out.flag = leaf(tree, {TString("topology_") + wp + "_flag",
                           TString("lambda_topology_") + wp + "_flag",
                           wp + "_topology_flag",
                           wp + "_flag"});
    if (out.flag) {
        out.source = out.flag->GetName();
        return out;
    }
    out.fiducial = leaf(tree, {"fiducial_vertex_flag", "fiducial_flag", "is_fiducial"});
    out.lambda_active = leaf(tree, {"lambda_decays_in_active_flag", "lambda_decay_in_active_flag",
                                    "lambda_active_flag"});
    out.proton = leaf(tree, {"proton_visible_flag", "p_visible_flag"});
    out.pion = leaf(tree, {"pion_visible_flag", "piminus_visible_flag", "pi_visible_flag"});
    out.detached = leaf(tree, {"detached_vertex_flag", "lambda_detached_flag"});
    if (out.fiducial && out.lambda_active && out.proton && out.pion) {
        out.source = "fiducial/lambda-active/daughter-visible truth flags";
        if (out.detached) out.source += " plus detached_vertex_flag";
    }
    return out;
}

bool has_topology_flags(const topology& t)
{
    return t.flag || (t.fiducial && t.lambda_active && t.proton && t.pion);
}

bool pass_topology_flags(const topology& t)
{
    if (t.flag) return value(t.flag) > 0.5;
    bool pass = value(t.fiducial) > 0.5 && value(t.lambda_active) > 0.5 &&
                value(t.proton) > 0.5 && value(t.pion) > 0.5;
    return pass && (!t.detached || value(t.detached) > 0.5);
}

void scale_density(TH1D* hist)
{
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        const double width = hist->GetBinWidth(i);
        hist->SetBinContent(i, hist->GetBinContent(i) / width);
        hist->SetBinError(i, hist->GetBinError(i) / width);
    }
}

void write_row(std::ofstream& out, const TString& sample, const TString& gen, const TString& knob,
               const TString& section, const TString& name, const acc& row, double denom,
               const TString& definition)
{
    const double frac = denom > 0.0 ? row.weight / denom : 0.0;
    out << csv(sample).Data() << ',' << csv(gen).Data() << ',' << csv(knob).Data() << ','
        << csv(section).Data() << ',' << csv(name).Data() << ',' << row.raw << ','
        << row.weight << ',' << std::sqrt(row.weight_sum_squares) << ','
        << row.xsec << ',' << std::sqrt(row.xsec_sum_squares) << ','
        << frac << ',' << csv(definition).Data() << '\n';
}

void write_summary(const TString& path, const TString& sample, const TString& gen,
                   const TString& knob, const acc& selected,
                   const std::vector<acc>& origins, const std::vector<acc>& modes,
                   const std::vector<acc>& categories)
{
    std::ofstream out(path.Data());
    out << "sample,generator,knob,section,name,raw_events,weighted_yield,"
        << "weighted_stat_uncertainty,xsec_weighted_1e38_cm2_per_target,"
        << "xsec_weighted_stat_uncertainty_1e38_cm2_per_target,"
        << "fraction_of_selected,definition\n";
    out << std::setprecision(12);
    write_row(out, sample, gen, knob, "total", "selected_topology", selected, selected.weight,
              "inclusive fiducial Lambda topology inside the NuMI flux envelope");

    const auto o_labels = origin_labels();
    const auto o_defs = origin_defs();
    const auto m_labels = mode_labels();
    const auto m_defs = mode_defs();
    const auto c_labels = category_labels();
    const auto c_defs = category_defs();
    for (int i = 0; i < int(origins.size()); ++i) {
        write_row(out, sample, gen, knob, "primary_origin", o_labels[i], origins[i], selected.weight, o_defs[i]);
    }
    for (int i = 0; i < int(modes.size()); ++i) {
        write_row(out, sample, gen, knob, "interaction_mode", m_labels[i], modes[i], selected.weight, m_defs[i]);
    }
    for (int i = 0; i < int(categories.size()); ++i) {
        write_row(out, sample, gen, knob, "visible_category", c_labels[i], categories[i], selected.weight, c_defs[i]);
    }
}

void write_origin_mode(const TString& path, const TString& sample, const TString& gen,
                       const TString& knob, const std::vector<std::vector<acc>>& matrix,
                       double denom)
{
    std::ofstream out(path.Data());
    out << "sample,generator,knob,primary_origin,interaction_mode,raw_events,"
        << "weighted_yield,weighted_stat_uncertainty,xsec_weighted_1e38_cm2_per_target,"
        << "xsec_weighted_stat_uncertainty_1e38_cm2_per_target,fraction_of_selected\n";
    out << std::setprecision(12);
    const auto origins = origin_labels();
    const auto modes = mode_labels();
    for (int o = 0; o < int(origins.size()); ++o) {
        for (int m = 0; m < int(modes.size()); ++m) {
            const acc& row = matrix[o][m];
            const double frac = denom > 0.0 ? row.weight / denom : 0.0;
            out << csv(sample).Data() << ',' << csv(gen).Data() << ',' << csv(knob).Data() << ','
                << csv(origins[o]).Data() << ',' << csv(modes[m]).Data() << ',' << row.raw << ','
                << row.weight << ',' << std::sqrt(row.weight_sum_squares) << ','
                << row.xsec << ',' << std::sqrt(row.xsec_sum_squares) << ','
                << frac << '\n';
        }
    }
}

void write_raw_modes(const TString& path, const TString& sample, const TString& gen,
                     const TString& knob, const std::map<int, acc>& rows, double denom)
{
    std::ofstream out(path.Data());
    out << "sample,generator,knob,raw_mode,raw_events,weighted_yield,"
        << "weighted_stat_uncertainty,xsec_weighted_1e38_cm2_per_target,"
        << "xsec_weighted_stat_uncertainty_1e38_cm2_per_target,"
        << "fraction_of_selected\n";
    out << std::setprecision(12);
    for (const auto& item : rows) {
        const double frac = denom > 0.0 ? item.second.weight / denom : 0.0;
        out << csv(sample).Data() << ',' << csv(gen).Data() << ',' << csv(knob).Data() << ','
            << item.first << ',' << item.second.raw << ',' << item.second.weight << ','
            << std::sqrt(item.second.weight_sum_squares) << ','
            << item.second.xsec << ',' << std::sqrt(item.second.xsec_sum_squares) << ','
            << frac << '\n';
    }
}

}  // namespace

void primary_mechanism(TString input_file,
                       TString output_dir = "analysis/output",
                       TString sample_label = "",
                       TString generator = "",
                       TString knob = "",
                       TString working_point = "nominal",
                       TString flux_file = "example/numi/flux/microboone_numi_flux_5mev.root",
                       double flux_floor_fraction = 0.0)
{
    if (sample_label == "") sample_label = gSystem->BaseName(input_file.Data());
    if (generator == "") generator = infer_generator(sample_label);
    if (knob == "") knob = infer_knob(sample_label);
    const TString sample = clean(sample_label);

    TFile* input = TFile::Open(input_file, "READ");
    if (!input || input->IsZombie()) {
        std::cerr << "Could not open input file: " << input_file.Data() << std::endl;
        return;
    }

    TTree* tree = nullptr;
    input->GetObject("FlatTree_VARS", tree);
    if (!tree) {
        std::cerr << "Could not find FlatTree_VARS in " << input_file.Data() << std::endl;
        input->Close();
        return;
    }

    TLeaf* enu = leaf(tree, {"Enu_true", "enu_true"});
    if (!enu) {
        std::cerr << "Missing Enu_true/enu_true branch; cannot apply the NuMI flux envelope." << std::endl;
        input->Close();
        return;
    }

    const flux_envelope flux = read_flux(flux_file, flux_floor_fraction);
    if (!flux.active) {
        std::cerr << "Refusing to make primary-mechanism summaries without a NuMI flux envelope." << std::endl;
        input->Close();
        return;
    }

    stage fsp{leaf(tree, {"nfsp", "n_fsp"}), leaf(tree, {"pdg"}), leaf(tree, {"px"}),
              leaf(tree, {"py"}), leaf(tree, {"pz"}), leaf(tree, {"E", "energy"})};
    stage primary{leaf(tree, {"nvertp", "n_vertp"}), leaf(tree, {"pdg_vert"})};
    if (!fsp.n || !fsp.pdg || !fsp.px || !fsp.py || !fsp.pz || !fsp.energy || !primary.n || !primary.pdg) {
        std::cerr << "Missing required FlatTree final-state or primary-vertex branches." << std::endl;
        input->Close();
        return;
    }

    TLeaf* mode_leaf = leaf(tree, {"Mode", "mode"});
    TLeaf* cc = leaf(tree, {"cc"});
    TLeaf* pdg_lep = leaf(tree, {"PDGLep", "pdg_lep"});
    TLeaf* e_lep = leaf(tree, {"ELep", "e_lep"});
    TLeaf* weight = leaf(tree, {"Weight", "weight"});
    TLeaf* scale = leaf(tree, {"fScaleFactor", "scale_factor"});
    TLeaf* target_a = leaf(tree, {"tgta", "target_a"});
    const topology topo = find_topology(tree, working_point);
    const bool has_topo_flags = has_topology_flags(topo);

    const auto origins = origin_labels();
    const auto modes = mode_labels();
    const auto categories = category_labels();
    const std::vector<TString> species = {"Lambda", "Sigma-", "Sigma0", "Sigma+", "Xi", "Omega", "K", "other_Y"};

    TH1D h_cutflow("topology_cutflow", ";selection;events [10^{-38} cm^{2}/target]", 10, 0.5, 10.5);
    label(h_cutflow.GetXaxis(), {"numi_flux", "final_Lambda", "visible_p_pi", "flat_proxy", "selected_topology",
                                  "S1_CC_muon", "S2_K", "S3_lowExtra", "S4_EM", "S5_noMuon"});
    TH1D h_origin("primary_origin", ";primary origin;events [10^{-38} cm^{2}/target]",
                  origins.size(), -0.5, origins.size() - 0.5);
    TH1D h_mode("interaction_mode", ";interaction mode;events [10^{-38} cm^{2}/target]",
                modes.size(), -0.5, modes.size() - 0.5);
    TH1D h_raw_mode("raw_mode", ";raw mode;events [10^{-38} cm^{2}/target]", 121, -60.5, 60.5);
    TH1D h_primary_species("primary_strange_species",
                           ";primary strange species;particles [10^{-38} cm^{2}/target]",
                           species.size(), -0.5, species.size() - 0.5);
    TH1D h_enu("enu_selected_topology", ";E_{#nu} [GeV];events / GeV [10^{-38} cm^{2}/target]",
               80, flux.emin, flux.emax);
    TH1D h_p_lambda("p_lambda_selected_topology",
                    ";p_{#Lambda} [GeV];events / GeV [10^{-38} cm^{2}/target]", 80, 0.0, 4.0);
    TH2D h_origin_mode("primary_origin_vs_interaction_mode", ";primary origin;interaction mode",
                       origins.size(), -0.5, origins.size() - 0.5, modes.size(), -0.5, modes.size() - 0.5);
    TH2D h_origin_category("primary_origin_vs_visible_category", ";primary origin;visible category",
                           origins.size(), -0.5, origins.size() - 0.5,
                           categories.size(), -0.5, categories.size() - 0.5);
    TH2D h_p_lambda_origin("p_lambda_vs_primary_origin", ";p_{#Lambda} [GeV];primary origin",
                           80, 0.0, 4.0, origins.size(), -0.5, origins.size() - 0.5);
    label(h_origin.GetXaxis(), origins);
    label(h_mode.GetXaxis(), modes);
    label(h_primary_species.GetXaxis(), species);
    label(h_origin_mode.GetXaxis(), origins);
    label(h_origin_mode.GetYaxis(), modes);
    label(h_origin_category.GetXaxis(), origins);
    label(h_origin_category.GetYaxis(), categories);
    label(h_p_lambda_origin.GetYaxis(), origins);

    acc selected;
    std::vector<acc> origin_count(origins.size()), mode_count(modes.size()), category_count(categories.size());
    std::vector<std::vector<acc>> origin_mode(origins.size(), std::vector<acc>(modes.size()));
    std::map<int, acc> raw_mode;
    Long64_t entries_in_flux = 0;

    const Long64_t entries = tree->GetEntries();
    for (Long64_t entry = 0; entry < entries; ++entry) {
        tree->GetEntry(entry);
        const double enu_gev = value(enu);
        if (!in_flux(flux, enu_gev)) continue;
        ++entries_in_flux;

        const double evt_w = value(weight, 0, 1.0) * value(scale, 0, 1.0);
        const double xsec_w = evt_w * std::max(1, int_value(target_a, 0, int(default_argon_a))) * units;
        h_cutflow.Fill(1, xsec_w);

        int n_lambda = 0, n_proton = 0, n_pi = 0, n_kaon = 0, n_gamma = 0;
        int best_lambda = -1, best_proton = -1, best_pi = -1;
        double p_lambda = -1.0, p_proton = -1.0, p_pi = -1.0;
        bool visible_muon = false;
        for (int i = 0; i < count(fsp); ++i) {
            const int code = int_value(fsp.pdg, i);
            const double p = momentum(fsp, i);
            if (code == 3122 && p > p_lambda) {
                ++n_lambda;
                best_lambda = i;
                p_lambda = p;
            }
            if (code == 2212 && p > 0.30 && p > p_proton) {
                ++n_proton;
                best_proton = i;
                p_proton = p;
            }
            if (code == -211 && p > 0.07 && p > p_pi) {
                ++n_pi;
                best_pi = i;
                p_pi = p;
            }
            if (is_kaon(code) && p > 0.10) ++n_kaon;
            if (code == 22 && value(fsp.energy, i) > 0.03) ++n_gamma;
            if (std::abs(code) == 13 && p > 0.10) visible_muon = true;
        }
        if (!visible_muon && int_value(cc) == 1 && std::abs(int_value(pdg_lep)) == 13 &&
            value(e_lep) > muon_mass_gev) {
            visible_muon = std::sqrt(value(e_lep) * value(e_lep) - muon_mass_gev * muon_mass_gev) > 0.10;
        }

        double extra_energy = 0.0;
        for (int i = 0; i < count(fsp); ++i) {
            const int code = int_value(fsp.pdg, i);
            const int a = std::abs(code);
            if (i == best_lambda || i == best_proton || i == best_pi || a == 12 || a == 14 || a == 16) continue;
            extra_energy += std::max(0.0, value(fsp.energy, i) - mass_gev(code));
        }

        const bool final_lambda = n_lambda > 0;
        const bool visible_p_pi = n_proton > 0 && n_pi > 0;
        const bool flat_proxy = final_lambda && visible_p_pi;
        const bool pass_topo = has_topo_flags ? pass_topology_flags(topo) : flat_proxy;
        if (final_lambda) h_cutflow.Fill(2, xsec_w);
        if (visible_p_pi) h_cutflow.Fill(3, xsec_w);
        if (flat_proxy) h_cutflow.Fill(4, xsec_w);
        if (!pass_topo) continue;

        const primary_counts primary_data = read_primary(primary);
        const bool true_lambda = final_lambda || has_topo_flags || primary_data.lambda ||
                                 primary_data.sigma_0 || primary_data.charged_sigma ||
                                 primary_data.strange_baryon;
        const int origin = origin_bin(primary_data, true_lambda);
        const int mode = mode_bin(int_value(mode_leaf), mode_leaf != nullptr);
        const bool cat[] = {true, int_value(cc) == 1 && visible_muon, n_kaon > 0,
                            extra_energy < 0.15, n_gamma > 0, !visible_muon};

        h_cutflow.Fill(5, xsec_w);
        for (int i = 1; i < 6; ++i) {
            if (cat[i]) h_cutflow.Fill(5 + i, xsec_w);
        }
        add(selected, evt_w, xsec_w);
        add(origin_count[origin], evt_w, xsec_w);
        add(mode_count[mode], evt_w, xsec_w);
        add(origin_mode[origin][mode], evt_w, xsec_w);
        if (mode_leaf) add(raw_mode[int_value(mode_leaf)], evt_w, xsec_w);
        for (int i = 0; i < 6; ++i) {
            if (cat[i]) add(category_count[i], evt_w, xsec_w);
        }

        h_origin.Fill(origin, xsec_w);
        h_mode.Fill(mode, xsec_w);
        h_origin_mode.Fill(origin, mode, xsec_w);
        if (mode_leaf) h_raw_mode.Fill(int_value(mode_leaf), xsec_w);
        for (int i = 0; i < 6; ++i) {
            if (cat[i]) h_origin_category.Fill(origin, i, xsec_w);
        }
        if (best_lambda >= 0) {
            h_enu.Fill(enu_gev, xsec_w);
            h_p_lambda.Fill(p_lambda, xsec_w);
            h_p_lambda_origin.Fill(p_lambda, origin, xsec_w);
        }
        for (int i = 0; i < count(primary); ++i) {
            const int bin = species_bin(int_value(primary.pdg, i));
            if (bin >= 0) h_primary_species.Fill(bin, xsec_w);
        }
    }

    gSystem->mkdir(output_dir.Data(), true);
    if (!output_dir.EndsWith("/")) output_dir += "/";
    const TString root_path = output_dir + "primary_mechanism_" + sample + ".root";
    const TString summary_path = output_dir + "primary_mechanism_" + sample + "_summary.csv";
    const TString matrix_path = output_dir + "primary_mechanism_" + sample + "_origin_mode.csv";
    const TString raw_mode_path = output_dir + "primary_mechanism_" + sample + "_raw_mode.csv";

    write_summary(summary_path, sample_label, generator, knob, selected, origin_count, mode_count, category_count);
    write_origin_mode(matrix_path, sample_label, generator, knob, origin_mode, selected.weight);
    write_raw_modes(raw_mode_path, sample_label, generator, knob, raw_mode, selected.weight);
    scale_density(&h_enu);
    scale_density(&h_p_lambda);

    TFile output(root_path, "RECREATE");
    TString note = "selection_source=";
    note += topo.source;
    note += "; working_point=";
    note += working_point;
    note += "; reduced_topology=";
    note += (has_topo_flags ? "true" : "false");
    note += "; generator=";
    note += generator;
    note += "; knob=";
    note += knob;
    note += "; numi_flux_source=";
    note += flux.source;
    note += "; numi_flux_emin_GeV=";
    note += TString::Format("%.8g", flux.emin);
    note += "; numi_flux_emax_GeV=";
    note += TString::Format("%.8g", flux.emax);
    note += "; if reduced_topology=false this is not a true fiducial selection";
    TNamed("selection_note", note).Write();
    h_cutflow.Write();
    h_origin.Write();
    h_mode.Write();
    h_raw_mode.Write();
    h_primary_species.Write();
    h_enu.Write();
    h_p_lambda.Write();
    h_origin_mode.Write();
    h_origin_category.Write();
    h_p_lambda_origin.Write();
    output.Close();

    std::cout << "\nSample: " << sample_label.Data()
              << "\nGenerator: " << generator.Data()
              << "\nKnob: " << knob.Data()
              << "\nInput: " << input_file.Data()
              << "\nEntries: " << entries
              << "\nEntries inside NuMI flux envelope: " << entries_in_flux
              << "\nNuMI flux envelope: " << flux.emin << " to " << flux.emax
              << " GeV from " << flux.source.Data()
              << "\nTopology selection source: " << topo.source.Data()
              << "\nOutputs:\n  " << root_path.Data()
              << "\n  " << summary_path.Data()
              << "\n  " << matrix_path.Data()
              << "\n  " << raw_mode_path.Data() << std::endl;
    input->Close();
}
