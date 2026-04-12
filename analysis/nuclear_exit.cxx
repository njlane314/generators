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

struct acc {
    Long64_t raw = 0;
    double weight = 0.0;
    double xsec = 0.0;
};

struct flux_envelope {
    bool active = false;
    double emin = 0.0;
    double emax = 0.0;
    double proposal_emin = 0.0;
    double proposal_emax = 10.0;
    double mean_bin_flux = 1.0;
    TString source;
    std::vector<TH1*> hists;
};

struct stage {
    TLeaf* n = nullptr;
    TLeaf* pdg = nullptr;
    TLeaf* px = nullptr;
    TLeaf* py = nullptr;
    TLeaf* pz = nullptr;
    TLeaf* energy = nullptr;
    TString source;
};

struct strange_counts {
    int lambda = 0;
    int sigma_minus = 0;
    int sigma_0 = 0;
    int sigma_plus = 0;
    int xi = 0;
    int omega = 0;
    int kaon = 0;
    int other_y = 0;
    int pion = 0;

    bool strange_baryon() const { return lambda || sigma_minus || sigma_0 || sigma_plus || xi || omega || other_y; }
    bool strange() const { return strange_baryon() || kaon; }
    bool charged_sigma() const { return sigma_minus || sigma_plus; }
};

std::vector<TString> species_labels()
{
    return {"Lambda", "Sigma-", "Sigma0", "Sigma+", "Xi", "Omega", "K", "other_Y"};
}

std::vector<TString> exit_labels()
{
    return {"Lambda", "Sigma-", "Sigma0", "Sigma+", "Xi", "Omega", "K", "other_Y", "no_exit_strange"};
}

std::vector<TString> origin_labels()
{
    return {"O0_direct_Lambda", "O1_Sigma0_feedthrough", "O2_charged_Sigma_exchange",
            "O3_neutral_Sigma_exchange", "O4_associated_KY",
            "O5_inelastic_hyperon_or_Ypi", "O6_other_strange", "O7_no_primary_strange"};
}

TLeaf* leaf(TTree* tree, std::initializer_list<TString> names, TString* found_name = nullptr)
{
    for (const TString& name : names) {
        if (TLeaf* found = tree->GetLeaf(name.Data())) {
            if (found_name) *found_name = name;
            return found;
        }
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
    total.xsec += xsec;
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

TString normalize_beam_polarity(TString beam_polarity)
{
    beam_polarity.ToLower();
    if (beam_polarity == "" || beam_polarity == "combined" ||
        beam_polarity == "both" || beam_polarity == "all") return "combined";
    if (beam_polarity == "fhc" || beam_polarity == "rhc") return beam_polarity;
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

TString flux_hist_path(const TString& beam_polarity, const TString& beam_species)
{
    TString path = beam_polarity;
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

std::vector<TString> flux_hist_paths(const TString& beam_polarity, TString beam_species)
{
    beam_species = normalize_beam_species(beam_species);
    std::vector<TString> paths;
    auto add = [&](const TString& polarity, const TString& species) {
        paths.push_back(flux_hist_path(polarity, species));
    };
    if (beam_polarity == "combined") {
        if (beam_species == "") {
            add("fhc", "numu");
            add("fhc", "numubar");
            add("rhc", "numu");
            add("rhc", "numubar");
        } else {
            add("fhc", beam_species);
            add("rhc", beam_species);
        }
    } else if (beam_polarity == "fhc" || beam_polarity == "rhc") {
        add(beam_polarity, beam_species == "" ? TString(beam_polarity == "rhc" ? "numubar" : "numu") : beam_species);
    }
    return paths;
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
    std::vector<TString> known = {
        "AR23_20i_00_000", "G18_10a_02_11a", "G18_10b_02_11a",
        "G18_10a_02_11b", "G18_10b_02_11b", "G18_10c_02_11b",
        "G18_10d_02_11b", "hyp_lambda_only", "hyp_sigma0_only",
        "hyp_sigmam_only", "hyp_no_effmass", "all_strange",
        "dis_only", "hyp_all", "fsi_on", "fsi_off"
    };
    TString knob = "";
    for (const TString& name : known) {
        if (!sample.Contains(name)) continue;
        if (knob != "") knob += "+";
        knob += name;
    }
    return knob == "" ? TString("nominal") : knob;
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

int species_bin(int pdg)
{
    if (pdg == 3122) return 0;
    if (pdg == 3112) return 1;
    if (pdg == 3212) return 2;
    if (pdg == 3222) return 3;
    if (std::abs(pdg) == 3312 || std::abs(pdg) == 3322) return 4;
    if (std::abs(pdg) == 3334) return 5;
    if (is_kaon(pdg)) return 6;
    return -1;
}

void add_strange(strange_counts& out, int pdg)
{
    if (pdg == 3122) ++out.lambda;
    else if (pdg == 3112) ++out.sigma_minus;
    else if (pdg == 3212) ++out.sigma_0;
    else if (pdg == 3222) ++out.sigma_plus;
    else if (std::abs(pdg) == 3312 || std::abs(pdg) == 3322) ++out.xi;
    else if (std::abs(pdg) == 3334) ++out.omega;
    else if (is_kaon(pdg)) ++out.kaon;
    if (std::abs(pdg) == 211 || pdg == 111) ++out.pion;
}

strange_counts read_strange(const stage& s)
{
    strange_counts out;
    for (int i = 0; i < count(s); ++i) add_strange(out, int_value(s.pdg, i));
    return out;
}

int origin_bin(const strange_counts& p)
{
    if (p.kaon && p.strange_baryon()) return 4;
    if (p.sigma_0) return 1;
    if (p.charged_sigma()) return 2;
    if (p.lambda) return 0;
    if (p.xi || p.omega || (p.strange_baryon() && p.pion) || p.strange_baryon()) return 5;
    if (p.kaon) return 6;
    return 7;
}

int exit_bin(const strange_counts& x)
{
    if (x.lambda) return 0;
    if (x.sigma_minus) return 1;
    if (x.sigma_0) return 2;
    if (x.sigma_plus) return 3;
    if (x.xi) return 4;
    if (x.omega) return 5;
    if (x.kaon) return 6;
    if (x.other_y || x.strange_baryon()) return 7;
    return 8;
}

bool has_momentum(const stage& s)
{
    return s.px && s.py && s.pz;
}

double momentum(const stage& s, int i)
{
    const double x = value(s.px, i);
    const double y = value(s.py, i);
    const double z = value(s.pz, i);
    return std::sqrt(x * x + y * y + z * z);
}

stage make_stage(
    TTree* tree,
    std::initializer_list<TString> n_names,
    std::initializer_list<TString> pdg_names,
    std::initializer_list<TString> px_names,
    std::initializer_list<TString> py_names,
    std::initializer_list<TString> pz_names,
    std::initializer_list<TString> energy_names
)
{
    TString n_name, pdg_name, px_name, py_name, pz_name, energy_name;
    stage s;
    s.n = leaf(tree, n_names, &n_name);
    s.pdg = leaf(tree, pdg_names, &pdg_name);
    s.px = leaf(tree, px_names, &px_name);
    s.py = leaf(tree, py_names, &py_name);
    s.pz = leaf(tree, pz_names, &pz_name);
    s.energy = leaf(tree, energy_names, &energy_name);
    if (s.n && s.pdg) {
        s.source = n_name;
        s.source += ", ";
        s.source += pdg_name;
        if (has_momentum(s)) {
            s.source += ", ";
            s.source += px_name;
            s.source += "/";
            s.source += py_name;
            s.source += "/";
            s.source += pz_name;
        }
        if (s.energy) {
            s.source += ", ";
            s.source += energy_name;
        }
    }
    return s;
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

TH1* find_flux_hist(TFile* file, const TString& path)
{
    TH1* hist = dynamic_cast<TH1*>(file->Get(path));
    return hist ? hist : find_hist(file, hist_leaf_name(path));
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
    TH1* clone = static_cast<TH1*>(hist->Clone(TString::Format("%s_flat_proposal_reweight_%d", hist->GetName(), int(env.hists.size()))));
    clone->SetDirectory(nullptr);
    env.hists.push_back(clone);
    env.emin = env.active ? std::min(env.emin, emin) : emin;
    env.emax = env.active ? std::max(env.emax, emax) : emax;
    env.active = true;
    return true;
}

double mean_flux_in_proposal(const flux_envelope& env)
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

flux_envelope read_flux(TString path, TString beam_polarity, TString beam_species, double floor_fraction,
                        double proposal_emin, double proposal_emax)
{
    flux_envelope env;
    env.proposal_emin = proposal_emin;
    env.proposal_emax = proposal_emax;
    if (path == "") {
        env.active = true;
        env.emin = proposal_emin;
        env.emax = proposal_emax;
        env.mean_bin_flux = 1.0;
        env.source = "no_flux_reweight";
        return env;
    }
    beam_polarity = normalize_beam_polarity(beam_polarity);
    if (beam_polarity == "") {
        std::cerr << "Beam polarity must be FHC, RHC, or combined." << std::endl;
        return env;
    }
    TFile* file = TFile::Open(path, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Could not open NuMI flux file: " << path.Data() << std::endl;
        return env;
    }
    std::vector<TString> used;
    const std::vector<TString> names = flux_hist_paths(beam_polarity, beam_species);
    for (const TString& name : names) {
        if (extend_flux(find_flux_hist(file, name), floor_fraction, env)) used.push_back(name);
        else std::cerr << "Could not use NuMI flux histogram: " << name.Data() << std::endl;
    }
    file->Close();
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
    env.source += beam_polarity;
    env.source += ": ";
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

double flux_reweight(const flux_envelope& env, double enu)
{
    if (!in_flux(env, enu) || env.mean_bin_flux <= 0.0) return 0.0;
    if (env.hists.empty()) return 1.0;
    double flux = 0.0;
    for (TH1* hist : env.hists) {
        const int bin = hist->FindBin(enu);
        if (bin < 1 || bin > hist->GetNbinsX()) continue;
        flux += std::max(0.0, hist->GetBinContent(bin));
    }
    return flux / env.mean_bin_flux;
}

void scale_density(TH1D* hist)
{
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        const double width = hist->GetBinWidth(i);
        hist->SetBinContent(i, hist->GetBinContent(i) / width);
        hist->SetBinError(i, hist->GetBinError(i) / width);
    }
}

void write_row(
    std::ofstream& out, const TString& sample, const TString& gen, const TString& knob,
    const TString& beam_polarity, const TString& section, const TString& name,
    const acc& row, double denom
)
{
    const double frac = denom > 0.0 ? row.weight / denom : 0.0;
    out << csv(sample).Data() << ',' << csv(gen).Data() << ',' << csv(knob).Data() << ','
        << csv(beam_polarity).Data() << ',' << csv(section).Data() << ',' << csv(name).Data()
        << ',' << row.raw << ',' << row.weight << ',' << row.xsec << ','
        << frac << ',' << csv(name).Data() << '\n';
}

void write_summary(
    const TString& path, const TString& sample, const TString& gen, const TString& knob,
    const TString& beam_polarity, const acc& flux_total,
    const std::vector<std::pair<TString, acc>>& rows
)
{
    std::ofstream out(path.Data());
    out << "sample,generator,knob,beam_polarity,section,name,raw_events,weighted_yield,"
        << "xsec_weighted_1e38_cm2_per_target,fraction_of_numi_flux,definition\n";
    out << std::setprecision(12);
    write_row(
        out, sample, gen, knob, beam_polarity,
        "total", "numi_flux_envelope", flux_total, flux_total.weight
    );
    for (const auto& row : rows) {
        write_row(
            out, sample, gen, knob, beam_polarity,
            "nuclear_exit", row.first, row.second, flux_total.weight
        );
    }
}

void write_migration(
    const TString& path, const TString& sample, const TString& gen, const TString& knob,
    const TString& beam_polarity, const std::vector<std::vector<acc>>& matrix,
    const std::vector<acc>& origin_total, double flux_weight
)
{
    std::ofstream out(path.Data());
    out << "sample,generator,knob,beam_polarity,primary_origin,exit_class,raw_events,weighted_yield,"
        << "xsec_weighted_1e38_cm2_per_target,fraction_of_origin,fraction_of_numi_flux\n";
    out << std::setprecision(12);
    const auto origins = origin_labels();
    const auto exits = exit_labels();
    for (int o = 0; o < int(origins.size()); ++o) {
        for (int e = 0; e < int(exits.size()); ++e) {
            const acc& row = matrix[o][e];
            const double f_origin = origin_total[o].weight > 0.0 ? row.weight / origin_total[o].weight : 0.0;
            const double f_flux = flux_weight > 0.0 ? row.weight / flux_weight : 0.0;
            out << csv(sample).Data() << ',' << csv(gen).Data() << ',' << csv(knob).Data() << ','
                << csv(beam_polarity).Data() << ',' << csv(origins[o]).Data() << ','
                << csv(exits[e]).Data() << ',' << row.raw << ',' << row.weight << ','
                << row.xsec << ',' << f_origin << ',' << f_flux << '\n';
        }
    }
}

void write_raw_exit(
    const TString& path, const TString& sample, const TString& gen, const TString& knob,
    const TString& beam_polarity, const std::map<int, acc>& rows
)
{
    std::ofstream out(path.Data());
    out << "sample,generator,knob,beam_polarity,exit_pdg,raw_particles,weighted_particles,"
        << "xsec_weighted_particles_1e38_cm2_per_target\n";
    out << std::setprecision(12);
    for (const auto& row : rows) {
        out << csv(sample).Data() << ',' << csv(gen).Data() << ',' << csv(knob).Data() << ','
            << csv(beam_polarity).Data() << ',' << row.first << ',' << row.second.raw << ','
            << row.second.weight << ',' << row.second.xsec << '\n';
    }
}

}  // namespace

void nuclear_exit(
    TString input_file,
    TString output_dir = "analysis/output",
    TString sample_label = "",
    TString generator = "",
    TString knob = "",
    TString flux_file = "analysis/flux/microboone_numi_flux_5mev.root",
    double flux_floor_fraction = 0.0,
    TString beam_polarity = "combined",
    double proposal_emin_gev = 0.0,
    double proposal_emax_gev = 10.0,
    TString beam_species = ""
)
{
    if (sample_label == "") sample_label = gSystem->BaseName(input_file.Data());
    if (generator == "") generator = infer_generator(sample_label);
    if (knob == "") knob = infer_knob(sample_label);
    beam_polarity = normalize_beam_polarity(beam_polarity);
    if (beam_polarity == "") {
        std::cerr << "Beam polarity must be FHC, RHC, or combined." << std::endl;
        return;
    }
    const TString sample = clean(sample_label);
    const TString beam = clean(beam_polarity);

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
        std::cerr << "Missing Enu_true/enu_true branch; cannot apply the analysis energy window." << std::endl;
        input->Close();
        return;
    }

    beam_species = normalize_beam_species(beam_species);
    const flux_envelope flux = read_flux(flux_file, beam_polarity, beam_species, flux_floor_fraction,
                                         proposal_emin_gev, proposal_emax_gev);
    if (!flux.active) {
        std::cerr << "Refusing to make nuclear-exit summaries without an active flux or energy window." << std::endl;
        input->Close();
        return;
    }

    stage primary = make_stage(
        tree, {"nvertp", "n_vertp"}, {"pdg_vert"}, {"px_vert"}, {"py_vert"},
        {"pz_vert"}, {"E_vert", "energy_vert"}
    );
    stage exit = make_stage(
        tree,
        {"n_exitp", "nexitp", "n_nuclear_exit", "nuclear_exit_n",
        "n_exit", "n_strange_exit", "n_exit_strange", "n_exit_particles"},
        {"pdg_exit", "pdg_exitp", "pdg_nuclear_exit", "nuclear_exit_pdg",
        "exit_pdg", "pdg_strange_exit", "exit_strange_pdg"},
        {"px_exit", "px_exitp", "px_nuclear_exit", "nuclear_exit_px",
        "exit_px", "px_strange_exit"},
        {"py_exit", "py_exitp", "py_nuclear_exit", "nuclear_exit_py",
        "exit_py", "py_strange_exit"},
        {"pz_exit", "pz_exitp", "pz_nuclear_exit", "nuclear_exit_pz",
        "exit_pz", "pz_strange_exit"},
        {"E_exit", "energy_exit", "E_exitp", "E_nuclear_exit",
        "nuclear_exit_E", "exit_E", "energy_strange_exit"}
    );
    if (!primary.n || !primary.pdg) {
        std::cerr << "Missing primary ancestry branches: expected nvertp/n_vertp and pdg_vert." << std::endl;
        input->Close();
        return;
    }
    if (!exit.n || !exit.pdg) {
        std::cerr << "Missing nuclear-exit branches. Expected n_exitp/pdg_exit or "
            << "nuclear_exit_n/nuclear_exit_pdg. Final-state pdg is not used as an exit proxy."
            << std::endl;
        input->Close();
        return;
    }

    TLeaf* weight = leaf(tree, {"Weight", "weight"});
    TLeaf* scale = leaf(tree, {"fScaleFactor", "scale_factor"});
    TLeaf* target_a = leaf(tree, {"tgta", "target_a"});
    const auto species = species_labels();
    const auto origins = origin_labels();
    const auto exits = exit_labels();

    TH1D h_cutflow("nuclear_exit_cutflow", ";selection;events [10^{-38} cm^{2}/target]", 10, 0.5, 10.5);
    label(h_cutflow.GetXaxis(), {
        "numi_flux", "primary_strange", "exit_strange", "exit_Lambda",
        "primary_Lambda", "primary_Lambda_to_exit_Lambda",
        "Sigma0_to_exit_Lambda", "charged_Sigma_to_exit_Lambda",
        "KY_to_exit_Lambda", "primary_strange_lost"
    });
    TH1D h_primary_species("primary_strange_species",
        ";primary strange species;particles [10^{-38} cm^{2}/target]",
        species.size(), -0.5, species.size() - 0.5);
    TH1D h_exit_species("exit_strange_species",
        ";nuclear-exit strange species;particles [10^{-38} cm^{2}/target]",
        species.size(), -0.5, species.size() - 0.5);
    TH1D h_origin("primary_origin", ";primary origin;events [10^{-38} cm^{2}/target]",
        origins.size(), -0.5, origins.size() - 0.5);
    TH1D h_exit("exit_class", ";nuclear-exit class;events [10^{-38} cm^{2}/target]",
        exits.size(), -0.5, exits.size() - 0.5);
    TH1D h_enu_lambda("enu_exit_lambda", ";E_{#nu} [GeV];events / GeV [10^{-38} cm^{2}/target]",
        80, flux.emin, flux.emax);
    TH1D h_p_lambda("p_exit_lambda", ";p_{#Lambda}^{exit} [GeV];particles / GeV [10^{-38} cm^{2}/target]",
        80, 0.0, 4.0);
    TH2D h_origin_exit("primary_origin_vs_exit_class", ";primary origin;nuclear-exit class",
        origins.size(), -0.5, origins.size() - 0.5, exits.size(), -0.5, exits.size() - 0.5);
    label(h_primary_species.GetXaxis(), species);
    label(h_exit_species.GetXaxis(), species);
    label(h_origin.GetXaxis(), origins);
    label(h_exit.GetXaxis(), exits);
    label(h_origin_exit.GetXaxis(), origins);
    label(h_origin_exit.GetYaxis(), exits);

    acc flux_total, primary_strange, exit_strange, exit_lambda, primary_lambda;
    acc direct_lambda, sigma0_lambda, charged_sigma_lambda, ky_lambda;
    acc strange_lost, no_primary_exit_strange, exit_kaon, exit_sigma;
    std::vector<acc> origin_total(origins.size()), exit_total(exits.size());
    std::vector<std::vector<acc>> migration(origins.size(), std::vector<acc>(exits.size()));
    std::map<int, acc> raw_exit;
    Long64_t entries_in_flux = 0;

    const Long64_t entries = tree->GetEntries();
    for (Long64_t entry = 0; entry < entries; ++entry) {
        tree->GetEntry(entry);
        const double enu_gev = value(enu);
        if (!in_flux(flux, enu_gev)) continue;
        ++entries_in_flux;

        const double evt_w = value(weight, 0, 1.0) * value(scale, 0, 1.0) * flux_reweight(flux, enu_gev);
        const double xsec_w = evt_w * std::max(1, int_value(target_a, 0, int(default_argon_a))) * units;
        const strange_counts p = read_strange(primary);
        const strange_counts x = read_strange(exit);
        const int origin = origin_bin(p);
        const int xbin = exit_bin(x);

        add(flux_total, evt_w, xsec_w);
        add(origin_total[origin], evt_w, xsec_w);
        add(exit_total[xbin], evt_w, xsec_w);
        add(migration[origin][xbin], evt_w, xsec_w);
        h_cutflow.Fill(1, xsec_w);
        h_origin.Fill(origin, xsec_w);
        h_exit.Fill(xbin, xsec_w);
        h_origin_exit.Fill(origin, xbin, xsec_w);

        if (p.strange()) {
            add(primary_strange, evt_w, xsec_w);
            h_cutflow.Fill(2, xsec_w);
        }
        if (x.strange()) {
            add(exit_strange, evt_w, xsec_w);
            h_cutflow.Fill(3, xsec_w);
        }
        if (x.lambda) {
            add(exit_lambda, evt_w, xsec_w);
            h_cutflow.Fill(4, xsec_w);
            h_enu_lambda.Fill(enu_gev, xsec_w);
        }
        if (p.lambda) {
            add(primary_lambda, evt_w, xsec_w);
            h_cutflow.Fill(5, xsec_w);
        }
        if (p.lambda && x.lambda) {
            add(direct_lambda, evt_w, xsec_w);
            h_cutflow.Fill(6, xsec_w);
        }
        if (p.sigma_0 && x.lambda) {
            add(sigma0_lambda, evt_w, xsec_w);
            h_cutflow.Fill(7, xsec_w);
        }
        if (p.charged_sigma() && x.lambda) {
            add(charged_sigma_lambda, evt_w, xsec_w);
            h_cutflow.Fill(8, xsec_w);
        }
        if (p.kaon && p.strange_baryon() && x.lambda) {
            add(ky_lambda, evt_w, xsec_w);
            h_cutflow.Fill(9, xsec_w);
        }
        if (p.strange() && !x.strange()) {
            add(strange_lost, evt_w, xsec_w);
            h_cutflow.Fill(10, xsec_w);
        }
        if (!p.strange() && x.strange()) add(no_primary_exit_strange, evt_w, xsec_w);
        if (x.kaon) add(exit_kaon, evt_w, xsec_w);
        if (x.sigma_minus || x.sigma_0 || x.sigma_plus) add(exit_sigma, evt_w, xsec_w);

        for (int i = 0; i < count(primary); ++i) {
            const int bin = species_bin(int_value(primary.pdg, i));
            if (bin >= 0) h_primary_species.Fill(bin, xsec_w);
        }
        for (int i = 0; i < count(exit); ++i) {
            const int code = int_value(exit.pdg, i);
            const int bin = species_bin(code);
            if (bin < 0) continue;
            h_exit_species.Fill(bin, xsec_w);
            add(raw_exit[code], evt_w, xsec_w);
            if (code == 3122 && has_momentum(exit)) h_p_lambda.Fill(momentum(exit, i), xsec_w);
        }
    }

    gSystem->mkdir(output_dir.Data(), true);
    if (!output_dir.EndsWith("/")) output_dir += "/";
    const TString root_path = output_dir + "nuclear_exit_" + sample + "_" + beam + ".root";
    const TString summary_path = output_dir + "nuclear_exit_" + sample + "_" + beam + "_summary.csv";
    const TString migration_path = output_dir + "nuclear_exit_" + sample + "_" + beam + "_migration.csv";
    const TString raw_path = output_dir + "nuclear_exit_" + sample + "_" + beam + "_raw_exit_pdg.csv";

    std::vector<std::pair<TString, acc>> rows = {
        {"primary_strange", primary_strange}, {"exit_strange", exit_strange},
        {"exit_Lambda", exit_lambda}, {"primary_Lambda", primary_lambda},
        {"primary_Lambda_to_exit_Lambda", direct_lambda}, {"Sigma0_to_exit_Lambda", sigma0_lambda},
        {"charged_Sigma_to_exit_Lambda", charged_sigma_lambda}, {"associated_KY_to_exit_Lambda", ky_lambda},
        {"primary_strange_lost_before_exit", strange_lost}, {"no_primary_strange_but_exit_strange", no_primary_exit_strange},
        {"exit_kaon", exit_kaon}, {"exit_Sigma", exit_sigma}};
    write_summary(summary_path, sample_label, generator, knob, beam_polarity, flux_total, rows);
    write_migration(migration_path, sample_label, generator, knob, beam_polarity,
        migration, origin_total, flux_total.weight);
    write_raw_exit(raw_path, sample_label, generator, knob, beam_polarity, raw_exit);
    scale_density(&h_enu_lambda);
    scale_density(&h_p_lambda);

    TFile output(root_path, "RECREATE");
    TString note = "primary_source=";
    note += primary.source;
    note += "; nuclear_exit_source=";
    note += exit.source;
    note += "; generator=";
    note += generator;
    note += "; knob=";
    note += knob;
    note += "; beam_polarity=";
    note += beam_polarity;
    note += "; beam_species=";
    note += beam_species == "" ? TString("default") : beam_species;
    note += "; numi_flux_source=";
    note += flux.source;
    note += "; numi_flux_emin_GeV=";
    note += TString::Format("%.8g", flux.emin);
    note += "; numi_flux_emax_GeV=";
    note += TString::Format("%.8g", flux.emax);
    note += "; flat_proposal_emin_GeV=";
    note += TString::Format("%.8g", flux.proposal_emin);
    note += "; flat_proposal_emax_GeV=";
    note += TString::Format("%.8g", flux.proposal_emax);
    note += "; numi_flux_reweight_mean_bin_flux=";
    note += TString::Format("%.8g", flux.mean_bin_flux);
    note += "; final_state_pdg_not_used_as_exit_proxy=true";
    TNamed(TString("nuclear_exit_metadata"), note).Write();
    h_cutflow.Write();
    h_primary_species.Write();
    h_exit_species.Write();
    h_origin.Write();
    h_exit.Write();
    h_enu_lambda.Write();
    h_p_lambda.Write();
    h_origin_exit.Write();
    output.Close();

    std::cout << "\nSample: " << sample_label.Data()
        << "\nGenerator: " << generator.Data()
        << "\nKnob: " << knob.Data()
        << "\nBeam polarity: " << beam_polarity.Data()
        << "\nInput: " << input_file.Data()
        << "\nEntries: " << entries
        << "\nEntries inside NuMI flux envelope: " << entries_in_flux
        << "\nNuMI flux envelope: " << flux.emin << " to " << flux.emax
        << " GeV from " << flux.source.Data()
        << "\nPrimary source: " << primary.source.Data()
        << "\nNuclear-exit source: " << exit.source.Data()
        << "\nOutputs:\n  " << root_path.Data()
        << "\n  " << summary_path.Data()
        << "\n  " << migration_path.Data()
        << "\n  " << raw_path.Data() << std::endl;
    input->Close();
}
