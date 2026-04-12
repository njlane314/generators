#include <TAxis.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TKey.h>
#include <TLeaf.h>
#include <TNamed.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace {

const int max_particles = 100;
const double units = 1.0e38;
const double default_argon_a = 40.0;
const double muon_mass_gev = 0.1056583745;

struct accumulator {
  Long64_t raw = 0;
  double weighted = 0.0;
  double xsec_weighted = 0.0;
};

struct primary_info {
  int lambda = 0;
  int sigma0 = 0;
  int charged_sigma = 0;
  int kaon = 0;
  int strange_baryon = 0;
  int xi_or_omega = 0;
  int pion = 0;
};

struct topology_leaf_set {
  TLeaf* topology = nullptr;
  TLeaf* fiducial_vertex = nullptr;
  TLeaf* lambda_active = nullptr;
  TLeaf* proton_visible = nullptr;
  TLeaf* pion_visible = nullptr;
  TLeaf* detached_vertex = nullptr;
  TString source = "flat-tree final Lambda + visible p pi- proxy";
};

struct flux_envelope {
  bool active = false;
  double emin = 0.0;
  double emax = 0.0;
  TString source;
};

std::vector<TString> origin_label_names()
{
  return {"O0_direct_Lambda",
          "O1_Sigma0_feedthrough",
          "O2_charged_Sigma_exchange",
          "O3_neutral_Sigma_exchange",
          "O4_associated_KY",
          "O5_inelastic_hyperon_or_Ypi",
          "O6_other_strange",
          "O7_no_true_Lambda"};
}

std::vector<TString> origin_definition_texts()
{
  return {"primary Lambda at the interaction vertex",
          "primary Sigma0 at the interaction vertex",
          "primary Sigma+ or Sigma- at the interaction vertex",
          "reserved for neutral-Sigma migration when nuclear-exit ancestry is available",
          "primary kaon plus strange baryon at the interaction vertex",
          "primary Xi/Omega or other inelastic hyperon/Ypi-like strange ancestry",
          "selected Lambda topology with other or incomplete strange ancestry",
          "selected topology with no true-Lambda evidence; normally unused for truth-only samples"};
}

std::vector<TString> mode_label_names()
{
  return {"QE",
          "MEC",
          "RES",
          "DIS",
          "COH",
          "diffractive",
          "other",
          "missing_Mode"};
}

std::vector<TString> mode_definition_texts()
{
  return {"abs(Mode) == 1",
          "abs(Mode) == 2",
          "abs(Mode) in {10,11,12,13,17,22,23}",
          "abs(Mode) in {21,26}",
          "abs(Mode) == 16",
          "abs(Mode) == 15",
          "Mode is present but not in the standard NUISANCE grouping used here",
          "no Mode/mode branch was found"};
}

std::vector<TString> visible_category_label_names()
{
  return {"S0_inclusive",
          "S1_CC_muon",
          "S2_visible_kaon",
          "S3_low_extra",
          "S4_visible_EM",
          "S5_no_visible_muon"};
}

std::vector<TString> visible_category_definition_texts()
{
  return {"inclusive selected Lambda topology",
          "S0 plus charged-current visible muon tag",
          "S0 plus visible final-state kaon",
          "S0 plus extra visible kinetic energy below 0.15 GeV",
          "S0 plus visible final-state gamma/EM tag",
          "S0 without a visible muon tag"};
}

bool bind_branch(TTree* tree,
                const std::vector<const char*>& names,
                void* address,
                bool required,
                const char* description,
                TString* bound_name = nullptr)
{
  for (const char* name : names) {
    if (!tree->GetBranch(name)) continue;
    tree->SetBranchAddress(name, address);
    if (bound_name) *bound_name = name;
    return true;
  }

  if (required) {
    std::cerr << "Missing required FlatTree_VARS branch for " << description << ": ";
    for (size_t i = 0; i < names.size(); ++i) {
      std::cerr << names[i] << (i + 1 == names.size() ? "" : " or ");
    }
    std::cerr << std::endl;
  }
  return false;
}

TLeaf* find_leaf(TTree* tree, const std::vector<TString>& names, TString* found_name = nullptr)
{
  for (const TString& name : names) {
    TLeaf* leaf = tree->GetLeaf(name);
    if (!leaf) continue;
    if (found_name) *found_name = name;
    return leaf;
  }
  return nullptr;
}

void add(accumulator& acc, double weighted, double xsec_weighted)
{
  ++acc.raw;
  acc.weighted += weighted;
  acc.xsec_weighted += xsec_weighted;
}

bool flag(TLeaf* leaf)
{
  return leaf && leaf->GetValue() > 0.5;
}

bool has_reduced_topology_flags(const topology_leaf_set& leaves)
{
  return leaves.topology ||
         (leaves.fiducial_vertex && leaves.lambda_active &&
          leaves.proton_visible && leaves.pion_visible);
}

bool pass_reduced_topology(const topology_leaf_set& leaves)
{
  if (leaves.topology) return flag(leaves.topology);

  bool pass = flag(leaves.fiducial_vertex) &&
              flag(leaves.lambda_active) &&
              flag(leaves.proton_visible) &&
              flag(leaves.pion_visible);
  if (leaves.detached_vertex) pass = pass && flag(leaves.detached_vertex);
  return pass;
}

bool in_flux_envelope(const flux_envelope& envelope, double enu)
{
  return envelope.active && enu >= envelope.emin && enu <= envelope.emax;
}

TH1* find_hist_by_name(TDirectory* directory, const TString& hist_name)
{
  if (!directory) return nullptr;

  TIter next_key(directory->GetListOfKeys());
  TKey* key = nullptr;
  while ((key = static_cast<TKey*>(next_key()))) {
    TObject* object = key->ReadObj();
    if (!object) continue;

    if (object->InheritsFrom(TH1::Class()) && hist_name == object->GetName()) {
      return static_cast<TH1*>(object);
    }

    if (object->InheritsFrom(TDirectory::Class())) {
      TH1* hist = find_hist_by_name(static_cast<TDirectory*>(object), hist_name);
      if (hist) return hist;
    }
  }

  return nullptr;
}

bool add_flux_hist_support(TH1* hist, double floor_fraction, flux_envelope& envelope)
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
    if (!found) {
      emin = low;
      emax = high;
      found = true;
    } else {
      emin = std::min(emin, low);
      emax = std::max(emax, high);
    }
  }

  if (!found) return false;

  if (!envelope.active) {
    envelope.emin = emin;
    envelope.emax = emax;
    envelope.active = true;
  } else {
    envelope.emin = std::min(envelope.emin, emin);
    envelope.emax = std::max(envelope.emax, emax);
  }
  return true;
}

flux_envelope read_numi_flux_envelope(TString flux_file, double floor_fraction)
{
  flux_envelope envelope;

  TFile* flux_input = TFile::Open(flux_file, "READ");
  if (!flux_input || flux_input->IsZombie()) {
    std::cerr << "Could not open NuMI flux file: " << flux_file.Data() << std::endl;
    return envelope;
  }

  const std::vector<TString> hist_names = {
      "numu_CV_AV_TPC_5MeV_bin",
      "numubar_CV_AV_TPC_5MeV_bin"};

  std::vector<TString> used_hists;
  for (const TString& hist_name : hist_names) {
    TH1* hist = find_hist_by_name(flux_input, hist_name);
    if (!hist) {
      std::cerr << "Could not find NuMI flux histogram: " << hist_name.Data() << std::endl;
      continue;
    }
    if (add_flux_hist_support(hist, floor_fraction, envelope)) {
      used_hists.push_back(hist_name);
    }
  }

  flux_input->Close();

  if (!envelope.active) {
    std::cerr << "No non-zero NuMI flux support was found in " << flux_file.Data() << std::endl;
    return envelope;
  }
  if (envelope.emax <= envelope.emin) {
    std::cerr << "Invalid NuMI flux envelope from " << flux_file.Data() << ": "
              << envelope.emin << " to " << envelope.emax << " GeV" << std::endl;
    envelope.active = false;
    return envelope;
  }

  envelope.source = flux_file;
  envelope.source += " [";
  for (int i = 0; i < int(used_hists.size()); ++i) {
    if (i) envelope.source += ", ";
    envelope.source += used_hists[i];
  }
  envelope.source += "]";
  return envelope;
}

topology_leaf_set find_topology_leaves(TTree* tree, TString working_point)
{
  working_point.ToLower();
  topology_leaf_set leaves;
  TString found;

  leaves.topology = find_leaf(
      tree,
      {TString("topology_") + working_point + "_flag",
       TString("lambda_topology_") + working_point + "_flag",
       working_point + "_topology_flag",
       working_point + "_flag"},
      &found);
  if (leaves.topology) {
    leaves.source = found;
    return leaves;
  }

  leaves.fiducial_vertex =
      find_leaf(tree, {"fiducial_vertex_flag", "fiducial_flag", "is_fiducial"});
  leaves.lambda_active =
      find_leaf(tree,
               {"lambda_decays_in_active_flag",
                "lambda_decay_in_active_flag",
                "lambda_active_flag"});
  leaves.proton_visible =
      find_leaf(tree, {"proton_visible_flag", "p_visible_flag"});
  leaves.pion_visible =
      find_leaf(tree, {"pion_visible_flag", "piminus_visible_flag", "pi_visible_flag"});
  leaves.detached_vertex =
      find_leaf(tree, {"detached_vertex_flag", "lambda_detached_flag"});

  if (has_reduced_topology_flags(leaves)) {
    leaves.source = "fiducial/lambda-active/daughter-visible truth flags";
    if (leaves.detached_vertex) leaves.source += " plus detached_vertex_flag";
  }

  return leaves;
}

bool is_kaon(int pdg)
{
  const int apdg = std::abs(pdg);
  return apdg == 321 || apdg == 311 || pdg == 310 || pdg == 130;
}

bool is_sigma(int pdg)
{
  const int apdg = std::abs(pdg);
  return apdg == 3112 || apdg == 3212 || apdg == 3222;
}

bool is_xi_or_omega(int pdg)
{
  const int apdg = std::abs(pdg);
  return apdg == 3312 || apdg == 3322 || apdg == 3334;
}

bool is_strange_baryon(int pdg)
{
  return std::abs(pdg) == 3122 || is_sigma(pdg) || is_xi_or_omega(pdg);
}

double momentum(double px, double py, double pz)
{
  return std::sqrt(px * px + py * py + pz * pz);
}

double mass_gev(int pdg)
{
  const int apdg = std::abs(pdg);
  if (apdg == 2212) return 0.9382720813;
  if (apdg == 2112) return 0.9395654133;
  if (apdg == 211) return 0.13957039;
  if (apdg == 321) return 0.493677;
  if (apdg == 311 || apdg == 310 || apdg == 130) return 0.497611;
  if (apdg == 3122) return 1.115683;
  if (apdg == 3112) return 1.197449;
  if (apdg == 3212) return 1.192642;
  if (apdg == 3222) return 1.18937;
  if (apdg == 3312) return 1.32171;
  if (apdg == 3322) return 1.31486;
  if (apdg == 3334) return 1.67245;
  if (apdg == 13) return muon_mass_gev;
  if (apdg == 11) return 0.00051099895;
  return 0.0;
}

int strange_species_bin(int pdg)
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

primary_info get_primary_info(int n_vertp, const Int_t* pdg_vert)
{
  primary_info info;
  for (int i = 0; i < std::min(n_vertp, max_particles); ++i) {
    const int code = pdg_vert[i];
    info.lambda += code == 3122;
    info.sigma0 += code == 3212;
    info.charged_sigma += code == 3112 || code == 3222;
    info.kaon += is_kaon(code);
    info.strange_baryon += is_strange_baryon(code);
    info.xi_or_omega += is_xi_or_omega(code);
    info.pion += std::abs(code) == 211 || code == 111;
  }
  return info;
}

int origin_bin(const primary_info& primary, bool true_lambda_evidence)
{
  if (primary.kaon > 0 && primary.strange_baryon > 0) return 4;
  if (primary.sigma0 > 0) return 1;
  if (primary.charged_sigma > 0) return 2;
  if (primary.lambda > 0) return 0;
  if (primary.xi_or_omega > 0) return 5;
  if (primary.strange_baryon > 0 && primary.pion > 0) return 5;
  if (primary.strange_baryon > 0) return 5;
  if (primary.kaon > 0 || true_lambda_evidence) return 6;
  return 7;
}

int interaction_mode_bin(int mode, bool has_mode)
{
  if (!has_mode) return 7;

  const int abs_mode = std::abs(mode);
  if (abs_mode == 1) return 0;
  if (abs_mode == 2) return 1;
  if (abs_mode == 10 || abs_mode == 11 || abs_mode == 12 ||
      abs_mode == 13 || abs_mode == 17 || abs_mode == 22 ||
      abs_mode == 23) {
    return 2;
  }
  if (abs_mode == 21 || abs_mode == 26) return 3;
  if (abs_mode == 16) return 4;
  if (abs_mode == 15) return 5;

  return 6;
}

TString clean_label(TString text)
{
  text.ReplaceAll(".flat.root", "");
  text.ReplaceAll(".root", "");
  text.ReplaceAll("/", "_");
  text.ReplaceAll(" ", "_");
  text.ReplaceAll(":", "_");
  text.ReplaceAll(";", "_");
  return text;
}

TString csv_escape(TString value)
{
  value.ReplaceAll("\"", "\"\"");
  TString escaped = "\"";
  escaped += value;
  escaped += "\"";
  return escaped;
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

  const std::vector<TString> known_knobs = {
      "AR23_20i_00_000",
      "G18_10a_02_11a",
      "G18_10b_02_11a",
      "G18_10a_02_11b",
      "G18_10b_02_11b",
      "G18_10c_02_11b",
      "G18_10d_02_11b",
      "hyp_lambda_only",
      "hyp_sigma0_only",
      "hyp_sigmam_only",
      "hyp_no_effmass",
      "all_strange",
      "dis_only",
      "hyp_all",
      "fsi_on",
      "fsi_off",
  };

  TString knob = "";
  for (const TString& candidate : known_knobs) {
    if (sample.Contains(candidate)) {
      if (knob != "") knob += "+";
      knob += candidate;
    }
  }

  return knob == "" ? "nominal" : knob;
}

void label_axis(TAxis* axis, const std::vector<TString>& labels)
{
  for (int i = 0; i < int(labels.size()); ++i) {
    axis->SetBinLabel(i + 1, labels[i]);
  }
}

void scale_density(TH1D* hist)
{
  for (int i = 1; i <= hist->GetNbinsX(); ++i) {
    const double width = hist->GetBinWidth(i);
    hist->SetBinContent(i, hist->GetBinContent(i) / width);
    hist->SetBinError(i, hist->GetBinError(i) / width);
  }
}

void write_summary_row(std::ofstream& csv,
                     const TString& sample_label,
                     const TString& generator,
                     const TString& knob,
                     TString section,
                     TString name,
                     const accumulator& acc,
                     double selected_weighted,
                     TString definition)
{
  const double fraction =
      selected_weighted > 0.0 ? acc.weighted / selected_weighted : 0.0;
  csv << csv_escape(sample_label).Data() << ','
      << csv_escape(generator).Data() << ','
      << csv_escape(knob).Data() << ','
      << csv_escape(section).Data() << ','
      << csv_escape(name).Data() << ','
      << acc.raw << ','
      << acc.weighted << ','
      << acc.xsec_weighted << ','
      << fraction << ','
      << csv_escape(definition).Data() << '\n';
}

void write_summary(const TString& path,
                  const TString& sample_label,
                  const TString& generator,
                  const TString& knob,
                  const accumulator& selected_total,
                  const std::vector<accumulator>& origin_counts,
                  const std::vector<accumulator>& mode_counts,
                  const std::vector<accumulator>& category_counts)
{
  const std::vector<TString> origin_labels = origin_label_names();
  const std::vector<TString> origin_definitions = origin_definition_texts();
  const std::vector<TString> mode_labels = mode_label_names();
  const std::vector<TString> mode_definitions = mode_definition_texts();
  const std::vector<TString> category_labels = visible_category_label_names();
  const std::vector<TString> category_definitions = visible_category_definition_texts();

  std::ofstream csv(path.Data());
  csv << "sample,generator,knob,section,name,raw_events,weighted_yield,"
      << "xsec_weighted_1e38_cm2_per_target,fraction_of_selected,definition\n";
  csv << std::setprecision(12);

  write_summary_row(csv,
                  sample_label,
                  generator,
                  knob,
                  "total",
                  "selected_topology",
                  selected_total,
                  selected_total.weighted,
                  "inclusive fiducial Lambda topology selection inside the NuMI flux envelope");

  for (int i = 0; i < int(origin_counts.size()); ++i) {
    write_summary_row(csv,
                    sample_label,
                    generator,
                    knob,
                    "primary_origin",
                    origin_labels[i],
                    origin_counts[i],
                    selected_total.weighted,
                    origin_definitions[i]);
  }

  for (int i = 0; i < int(mode_counts.size()); ++i) {
    write_summary_row(csv,
                    sample_label,
                    generator,
                    knob,
                    "interaction_mode",
                    mode_labels[i],
                    mode_counts[i],
                    selected_total.weighted,
                    mode_definitions[i]);
  }

  for (int i = 0; i < int(category_counts.size()); ++i) {
    write_summary_row(csv,
                    sample_label,
                    generator,
                    knob,
                    "visible_category",
                    category_labels[i],
                    category_counts[i],
                    selected_total.weighted,
                    category_definitions[i]);
  }
}

void write_origin_mode_matrix(
    const TString& path,
    const TString& sample_label,
    const TString& generator,
    const TString& knob,
    const std::vector<std::vector<accumulator>>& origin_mode_counts,
    double selected_weighted)
{
  const std::vector<TString> origin_labels = origin_label_names();
  const std::vector<TString> mode_labels = mode_label_names();

  std::ofstream csv(path.Data());
  csv << "sample,generator,knob,primary_origin,interaction_mode,raw_events,weighted_yield,"
      << "xsec_weighted_1e38_cm2_per_target,fraction_of_selected\n";
  csv << std::setprecision(12);

  for (int origin = 0; origin < int(origin_labels.size()); ++origin) {
    for (int mode = 0; mode < int(mode_labels.size()); ++mode) {
      const accumulator& acc = origin_mode_counts[origin][mode];
      const double fraction =
          selected_weighted > 0.0 ? acc.weighted / selected_weighted : 0.0;
      csv << csv_escape(sample_label).Data() << ','
          << csv_escape(generator).Data() << ','
          << csv_escape(knob).Data() << ','
          << csv_escape(origin_labels[origin]).Data() << ','
          << csv_escape(mode_labels[mode]).Data() << ','
          << acc.raw << ','
          << acc.weighted << ','
          << acc.xsec_weighted << ','
          << fraction << '\n';
    }
  }
}

void write_raw_mode_codes(const TString& path,
                       const TString& sample_label,
                       const TString& generator,
                       const TString& knob,
                       const std::map<int, accumulator>& raw_mode_counts,
                       double selected_weighted)
{
  std::ofstream csv(path.Data());
  csv << "sample,generator,knob,raw_Mode,raw_events,weighted_yield,"
      << "xsec_weighted_1e38_cm2_per_target,fraction_of_selected\n";
  csv << std::setprecision(12);

  for (const auto& item : raw_mode_counts) {
    const double fraction =
        selected_weighted > 0.0 ? item.second.weighted / selected_weighted : 0.0;
    csv << csv_escape(sample_label).Data() << ','
        << csv_escape(generator).Data() << ','
        << csv_escape(knob).Data() << ','
        << item.first << ','
        << item.second.raw << ','
        << item.second.weighted << ','
        << item.second.xsec_weighted << ','
        << fraction << '\n';
  }
}

void print_block(const char* title,
                const std::vector<TString>& labels,
                const std::vector<accumulator>& counts,
                double selected_weighted)
{
  std::cout << '\n' << title << '\n';
  std::cout << std::left << std::setw(36) << "category"
            << std::right << std::setw(12) << "raw"
            << std::setw(18) << "weighted"
            << std::setw(14) << "fraction" << '\n';
  std::cout << std::string(80, '-') << '\n';
  std::cout << std::setprecision(8);
  for (int i = 0; i < int(labels.size()); ++i) {
    const double fraction =
        selected_weighted > 0.0 ? counts[i].weighted / selected_weighted : 0.0;
    std::cout << std::left << std::setw(36) << labels[i].Data()
              << std::right << std::setw(12) << counts[i].raw
              << std::setw(18) << counts[i].weighted
              << std::setw(14) << fraction << '\n';
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
  if (sample_label == "") {
    sample_label = gSystem->BaseName(input_file.Data());
  }
  if (generator == "") {
    generator = infer_generator(sample_label);
    if (generator == "unspecified") generator = infer_generator(input_file);
  }
  if (knob == "") {
    knob = infer_knob(sample_label);
    if (knob == "nominal") knob = infer_knob(input_file);
  }
  const TString clean_sample_label = clean_label(sample_label);

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

  Int_t mode = 0;
  Char_t cc = 0;
  Int_t pdg_nu = 0;
  Float_t enu_true = 0.0;
  Int_t target_a = int(default_argon_a);
  Int_t pdg_lep = 0;
  Float_t e_lep = 0.0;
  Int_t n_fsp = 0;
  Float_t px[max_particles] = {0.0};
  Float_t py[max_particles] = {0.0};
  Float_t pz[max_particles] = {0.0};
  Float_t energy[max_particles] = {0.0};
  Int_t pdg[max_particles] = {0};
  Int_t n_vertp = 0;
  Int_t pdg_vert[max_particles] = {0};
  Float_t weight = 1.0;
  Double_t scale_factor = 1.0;

  bool ok = true;
  const bool has_mode =
      bind_branch(tree, {"Mode", "mode"}, &mode, false, "interaction mode");
  const bool has_cc =
      bind_branch(tree, {"cc"}, &cc, false, "charged-current flag");
  bind_branch(tree, {"PDGnu", "pdg_nu"}, &pdg_nu, false, "incoming neutrino PDG");
  const bool has_enu =
      bind_branch(tree, {"Enu_true", "enu_true"}, &enu_true, false, "true neutrino energy");
  bind_branch(tree, {"tgta", "target_a"}, &target_a, false, "target A");
  const bool has_pdg_lep =
      bind_branch(tree, {"PDGLep", "pdg_lep"}, &pdg_lep, false, "outgoing lepton PDG");
  const bool has_e_lep =
      bind_branch(tree, {"ELep", "e_lep"}, &e_lep, false, "outgoing lepton energy");
  ok &= bind_branch(tree, {"nfsp", "n_fsp"}, &n_fsp, true, "final-state multiplicity");
  ok &= bind_branch(tree, {"px"}, px, true, "final-state px");
  ok &= bind_branch(tree, {"py"}, py, true, "final-state py");
  ok &= bind_branch(tree, {"pz"}, pz, true, "final-state pz");
  ok &= bind_branch(tree, {"E", "energy"}, energy, true, "final-state energy");
  ok &= bind_branch(tree, {"pdg"}, pdg, true, "final-state PDG");
  ok &= bind_branch(tree, {"nvertp", "n_vertp"}, &n_vertp, true, "primary-vertex multiplicity");
  ok &= bind_branch(tree, {"pdg_vert"}, pdg_vert, true, "primary-vertex PDG");
  const bool has_weight =
      bind_branch(tree, {"Weight", "weight"}, &weight, false, "event weight");
  const bool has_scale_factor =
      bind_branch(tree, {"fScaleFactor", "scale_factor"}, &scale_factor, false, "scale factor");

  if (!has_weight) {
    std::cout << "No Weight/weight branch found; using event weight = 1." << std::endl;
  }
  if (!has_scale_factor) {
    std::cout << "No fScaleFactor/scale_factor branch found; using scale factor = 1." << std::endl;
  }
  if (!has_mode) {
    std::cout << "No Mode/mode branch found; the interaction-mode summary will use missing_Mode." << std::endl;
  }
  if (!has_enu) {
    std::cerr << "Missing Enu_true/enu_true branch; cannot apply the NuMI flux envelope." << std::endl;
    input->Close();
    return;
  }

  if (!ok) {
    input->Close();
    return;
  }

  const flux_envelope numi_envelope = read_numi_flux_envelope(flux_file, flux_floor_fraction);
  if (!numi_envelope.active) {
    std::cerr << "Refusing to make primary-mechanism summaries without a NuMI flux envelope." << std::endl;
    input->Close();
    return;
  }

  const topology_leaf_set topology_leaves = find_topology_leaves(tree, working_point);
  const bool has_reduced_topology = has_reduced_topology_flags(topology_leaves);

  const std::vector<TString> origin_labels = origin_label_names();
  const std::vector<TString> mode_labels = mode_label_names();
  const std::vector<TString> category_labels = visible_category_label_names();

  TH1D h_cutflow("topology_cutflow",
                 ";selection;events [10^{-38} cm^{2}/target]",
                 10,
                 0.5,
                 10.5);
  label_axis(h_cutflow.GetXaxis(),
            {"numi_flux",
             "final_Lambda",
             "visible_p_pi",
             "flat_proxy",
             "selected_topology",
             "S1_CC_muon",
             "S2_K",
             "S3_lowExtra",
             "S4_EM",
             "S5_noMuon"});

  TH1D h_origin("primary_origin",
                ";primary origin;events [10^{-38} cm^{2}/target]",
                origin_labels.size(),
                -0.5,
                origin_labels.size() - 0.5);
  label_axis(h_origin.GetXaxis(), origin_labels);

  TH1D h_mode("interaction_mode",
              ";interaction mode;events [10^{-38} cm^{2}/target]",
              mode_labels.size(),
              -0.5,
              mode_labels.size() - 0.5);
  label_axis(h_mode.GetXaxis(), mode_labels);

  TH1D h_raw_mode("raw_Mode",
                  ";raw Mode;events [10^{-38} cm^{2}/target]",
                  121,
                  -60.5,
                  60.5);

  TH1D h_primary_species("primary_strange_species",
                         ";primary strange species;particles [10^{-38} cm^{2}/target]",
                         8,
                         -0.5,
                         7.5);
  label_axis(h_primary_species.GetXaxis(),
            {"Lambda", "Sigma-", "Sigma0", "Sigma+", "Xi", "Omega", "K", "other_Y"});

  TH1D h_enu("enu_selected_topology",
             ";E_{#nu} [GeV];events / GeV [10^{-38} cm^{2}/target]",
             80,
             numi_envelope.emin,
             numi_envelope.emax);

  TH1D h_p_lambda("p_lambda_selected_topology",
                  ";p_{#Lambda} [GeV];events / GeV [10^{-38} cm^{2}/target]",
                  80,
                  0.0,
                  4.0);

  TH2D h_origin_mode("primary_origin_vs_interaction_mode",
                     ";primary origin;interaction mode",
                     origin_labels.size(),
                     -0.5,
                     origin_labels.size() - 0.5,
                     mode_labels.size(),
                     -0.5,
                     mode_labels.size() - 0.5);
  label_axis(h_origin_mode.GetXaxis(), origin_labels);
  label_axis(h_origin_mode.GetYaxis(), mode_labels);

  TH2D h_origin_category("primary_origin_vs_visible_category",
                         ";primary origin;visible category",
                         origin_labels.size(),
                         -0.5,
                         origin_labels.size() - 0.5,
                         category_labels.size(),
                         -0.5,
                         category_labels.size() - 0.5);
  label_axis(h_origin_category.GetXaxis(), origin_labels);
  label_axis(h_origin_category.GetYaxis(), category_labels);

  TH2D h_p_lambda_origin("p_lambda_vs_primary_origin",
                         ";p_{#Lambda} [GeV];primary origin",
                         80,
                         0.0,
                         4.0,
                         origin_labels.size(),
                         -0.5,
                         origin_labels.size() - 0.5);
  label_axis(h_p_lambda_origin.GetYaxis(), origin_labels);

  accumulator selected_total;
  std::vector<accumulator> origin_counts(origin_labels.size());
  std::vector<accumulator> mode_counts(mode_labels.size());
  std::vector<accumulator> category_counts(category_labels.size());
  std::vector<std::vector<accumulator>> origin_mode_counts(
      origin_labels.size(), std::vector<accumulator>(mode_labels.size()));
  std::map<int, accumulator> raw_mode_counts;

  const Long64_t entries = tree->GetEntries();
  Long64_t entries_in_flux_envelope = 0;
  for (Long64_t entry = 0; entry < entries; ++entry) {
    tree->GetEntry(entry);
    if (!in_flux_envelope(numi_envelope, enu_true)) continue;
    ++entries_in_flux_envelope;

    const double event_weight = weight * scale_factor;
    const double target_scale = std::max(1, target_a) * units;
    const double xsec_weight = event_weight * target_scale;

    h_cutflow.Fill(1, xsec_weight);

    int n_final_lambda = 0;
    int n_visible_proton = 0;
    int n_visible_piminus = 0;
    int n_visible_kaon = 0;
    int n_visible_gamma = 0;
    int i_best_lambda = -1;
    int i_best_proton = -1;
    int i_best_piminus = -1;
    double p_best_lambda = -1.0;
    double p_best_proton = -1.0;
    double p_best_piminus = -1.0;
    bool visible_muon = false;

    for (int i = 0; i < std::min(n_fsp, max_particles); ++i) {
      const double p = momentum(px[i], py[i], pz[i]);
      const int apdg = std::abs(pdg[i]);

      if (pdg[i] == 3122) {
        ++n_final_lambda;
        if (p > p_best_lambda) {
          p_best_lambda = p;
          i_best_lambda = i;
        }
      }
      if (pdg[i] == 2212 && p > 0.30) {
        ++n_visible_proton;
        if (p > p_best_proton) {
          p_best_proton = p;
          i_best_proton = i;
        }
      }
      if (pdg[i] == -211 && p > 0.07) {
        ++n_visible_piminus;
        if (p > p_best_piminus) {
          p_best_piminus = p;
          i_best_piminus = i;
        }
      }
      if (is_kaon(pdg[i]) && p > 0.10) ++n_visible_kaon;
      if (pdg[i] == 22 && energy[i] > 0.03) ++n_visible_gamma;
      if (apdg == 13 && p > 0.10) visible_muon = true;
    }

    if (!visible_muon && has_cc && has_pdg_lep && has_e_lep &&
        cc == 1 && std::abs(pdg_lep) == 13 && e_lep > muon_mass_gev) {
      visible_muon =
          std::sqrt(e_lep * e_lep - muon_mass_gev * muon_mass_gev) > 0.10;
    }

    double extra_visible_energy = 0.0;
    for (int i = 0; i < std::min(n_fsp, max_particles); ++i) {
      if (i == i_best_proton || i == i_best_piminus || i == i_best_lambda) continue;
      const int apdg = std::abs(pdg[i]);
      if (apdg == 12 || apdg == 14 || apdg == 16) continue;
      extra_visible_energy +=
          std::max(0.0, static_cast<double>(energy[i]) - mass_gev(pdg[i]));
    }

    const bool final_lambda = n_final_lambda > 0;
    const bool visible_p_pi = n_visible_proton > 0 && n_visible_piminus > 0;
    const bool flat_proxy = final_lambda && visible_p_pi;
    const bool selected_topology =
        has_reduced_topology ? pass_reduced_topology(topology_leaves) : flat_proxy;

    if (final_lambda) h_cutflow.Fill(2, xsec_weight);
    if (visible_p_pi) h_cutflow.Fill(3, xsec_weight);
    if (flat_proxy) h_cutflow.Fill(4, xsec_weight);
    if (!selected_topology) continue;

    const bool s1_cc_muon = has_cc && cc == 1 && visible_muon;
    const bool s2_kaon = n_visible_kaon > 0;
    const bool s3_low_extra = extra_visible_energy < 0.15;
    const bool s4_em = n_visible_gamma > 0;
    const bool s5_no_muon = !visible_muon;

    const primary_info primary = get_primary_info(n_vertp, pdg_vert);
    const bool true_lambda_evidence =
        final_lambda || has_reduced_topology || primary.lambda > 0 ||
        primary.sigma0 > 0 || primary.charged_sigma > 0 ||
        primary.strange_baryon > 0;
    const int origin = origin_bin(primary, true_lambda_evidence);
    const int mode_bin = interaction_mode_bin(mode, has_mode);

    h_cutflow.Fill(5, xsec_weight);
    if (s1_cc_muon) h_cutflow.Fill(6, xsec_weight);
    if (s2_kaon) h_cutflow.Fill(7, xsec_weight);
    if (s3_low_extra) h_cutflow.Fill(8, xsec_weight);
    if (s4_em) h_cutflow.Fill(9, xsec_weight);
    if (s5_no_muon) h_cutflow.Fill(10, xsec_weight);

    add(selected_total, event_weight, xsec_weight);
    add(origin_counts[origin], event_weight, xsec_weight);
    add(mode_counts[mode_bin], event_weight, xsec_weight);
    add(category_counts[0], event_weight, xsec_weight);
    if (s1_cc_muon) add(category_counts[1], event_weight, xsec_weight);
    if (s2_kaon) add(category_counts[2], event_weight, xsec_weight);
    if (s3_low_extra) add(category_counts[3], event_weight, xsec_weight);
    if (s4_em) add(category_counts[4], event_weight, xsec_weight);
    if (s5_no_muon) add(category_counts[5], event_weight, xsec_weight);
    add(origin_mode_counts[origin][mode_bin], event_weight, xsec_weight);
    if (has_mode) add(raw_mode_counts[mode], event_weight, xsec_weight);

    h_origin.Fill(origin, xsec_weight);
    h_mode.Fill(mode_bin, xsec_weight);
    if (has_mode) h_raw_mode.Fill(mode, xsec_weight);
    h_origin_mode.Fill(origin, mode_bin, xsec_weight);
    h_origin_category.Fill(origin, 0, xsec_weight);
    if (s1_cc_muon) h_origin_category.Fill(origin, 1, xsec_weight);
    if (s2_kaon) h_origin_category.Fill(origin, 2, xsec_weight);
    if (s3_low_extra) h_origin_category.Fill(origin, 3, xsec_weight);
    if (s4_em) h_origin_category.Fill(origin, 4, xsec_weight);
    if (s5_no_muon) h_origin_category.Fill(origin, 5, xsec_weight);
    if (has_enu) h_enu.Fill(enu_true, xsec_weight);
    if (i_best_lambda >= 0) {
      h_p_lambda.Fill(p_best_lambda, xsec_weight);
      h_p_lambda_origin.Fill(p_best_lambda, origin, xsec_weight);
    }

    for (int i = 0; i < std::min(n_vertp, max_particles); ++i) {
      const int species = strange_species_bin(pdg_vert[i]);
      if (species >= 0) h_primary_species.Fill(species, xsec_weight);
    }
  }

  gSystem->mkdir(output_dir.Data(), true);
  if (!output_dir.EndsWith("/")) output_dir += "/";

  const TString root_path =
      output_dir + "primary_mechanism_" + clean_sample_label + ".root";
  const TString summary_path =
      output_dir + "primary_mechanism_" + clean_sample_label + "_summary.csv";
  const TString matrix_path =
      output_dir + "primary_mechanism_" + clean_sample_label + "_origin_mode.csv";
  const TString raw_mode_path =
      output_dir + "primary_mechanism_" + clean_sample_label + "_raw_mode.csv";

  write_summary(summary_path,
               sample_label,
               generator,
               knob,
               selected_total,
               origin_counts,
               mode_counts,
               category_counts);
  write_origin_mode_matrix(
      matrix_path, sample_label, generator, knob, origin_mode_counts, selected_total.weighted);
  write_raw_mode_codes(
      raw_mode_path, sample_label, generator, knob, raw_mode_counts, selected_total.weighted);

  scale_density(&h_enu);
  scale_density(&h_p_lambda);

  TFile output(root_path, "RECREATE");
  TString note = "selection_source=";
  note += topology_leaves.source;
  note += "; working_point=";
  note += working_point;
  note += "; generator=";
  note += generator;
  note += "; knob=";
  note += knob;
  note += "; reduced_topology=";
  note += (has_reduced_topology ? "true" : "false");
  note += "; numi_flux_source=";
  note += numi_envelope.source;
  note += "; numi_flux_emin_GeV=";
  note += TString::Format("%.8g", numi_envelope.emin);
  note += "; numi_flux_emax_GeV=";
  note += TString::Format("%.8g", numi_envelope.emax);
  note += "; if reduced_topology=false this is not a true fiducial selection";
  TNamed selection_note("selection_note", note);
  selection_note.Write();
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

  std::cout << '\n'
            << "Sample: " << sample_label.Data() << '\n'
            << "Generator: " << generator.Data() << '\n'
            << "Knob: " << knob.Data() << '\n'
            << "Input: " << input_file.Data() << '\n'
            << "Entries: " << entries << '\n'
            << "Entries inside NuMI flux envelope: " << entries_in_flux_envelope << '\n'
            << "NuMI flux envelope: " << numi_envelope.emin << " to "
            << numi_envelope.emax << " GeV from " << numi_envelope.source.Data() << '\n'
            << "Topology selection source: " << topology_leaves.source.Data() << '\n';
  if (!has_reduced_topology) {
    std::cout << "Detector note: no reduced fiducial topology flags were found, "
              << "so the macro used the flat-tree final Lambda + visible p pi- proxy."
              << '\n';
  }
  std::cout << "Outputs:\n"
            << "  " << root_path.Data() << '\n'
            << "  " << summary_path.Data() << '\n'
            << "  " << matrix_path.Data() << '\n'
            << "  " << raw_mode_path.Data() << '\n';

  print_block("Primary origin breakdown",
             origin_labels,
             origin_counts,
             selected_total.weighted);
  print_block("Interaction-mode breakdown",
             mode_labels,
             mode_counts,
             selected_total.weighted);
  print_block("Visible-category tags",
             category_labels,
             category_counts,
             selected_total.weighted);

  input->Close();
}
