#include <TFile.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

const int max_particles = 100;
const double argon_a = 40.0;
const double units = 1.0e38;
const double muon_mass_gev = 0.1056583745;

struct yield_row {
  TString name;
  TString definition;
  Long64_t raw = 0;
  double weighted = 0.0;
  double xsec_weighted = 0.0;
};

bool bind_required_branch(TTree* tree, const char* name, void* address)
{
  if (!tree->GetBranch(name)) {
    std::cerr << "Missing required FlatTree_VARS branch: " << name << std::endl;
    return false;
  }
  tree->SetBranchAddress(name, address);
  return true;
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

void fill(yield_row& row, bool pass, double weighted, double xsec_weighted)
{
  if (!pass) return;
  ++row.raw;
  row.weighted += weighted;
  row.xsec_weighted += xsec_weighted;
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

void write_rows(const TString& output_file,
               const TString& sample_label,
               const TString& generator,
               const TString& knob,
               const std::vector<yield_row>& rows)
{
  std::ofstream csv(output_file.Data());
  csv << "sample,generator,knob,category,raw_events,weighted_yield,"
      << "xsec_weighted_1e38_cm2_per_Ar,definition\n";
  csv << std::setprecision(12);
  for (const yield_row& row : rows) {
    const TString definition = csv_escape(row.definition);
    const TString escaped_sample = csv_escape(sample_label);
    const TString escaped_generator = csv_escape(generator);
    const TString escaped_knob = csv_escape(knob);
    csv << escaped_sample.Data() << ','
        << escaped_generator.Data() << ','
        << escaped_knob.Data() << ','
        << row.name.Data() << ','
        << row.raw << ','
        << row.weighted << ','
        << row.xsec_weighted << ','
        << definition.Data() << '\n';
  }
}

void print_rows(const TString& generator, const TString& knob, const std::vector<yield_row>& rows)
{
  std::cout << "Generator: " << generator.Data() << '\n'
            << "Knob: " << knob.Data() << '\n'
            << '\n';
  std::cout << std::left << std::setw(36) << "category"
            << std::right << std::setw(12) << "raw"
            << std::setw(18) << "weighted"
            << std::setw(24) << "xsec_1e38_per_Ar" << '\n';
  std::cout << std::string(90, '-') << '\n';
  std::cout << std::setprecision(8);
  for (const yield_row& row : rows) {
    std::cout << std::left << std::setw(36) << row.name.Data()
              << std::right << std::setw(12) << row.raw
              << std::setw(18) << row.weighted
              << std::setw(24) << row.xsec_weighted << '\n';
  }
}

}  // namespace

void topology_yields(TString input_file,
                     TString output_csv = "",
                     TString sample_label = "",
                     TString generator = "",
                     TString knob = "")
{
  if (output_csv == "") {
    output_csv = input_file;
    output_csv.ReplaceAll(".flat.root", ".topology_yields.csv");
    if (output_csv == input_file) {
      output_csv += ".topology_yields.csv";
    }
  }

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

  Char_t cc = 0;
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
  ok &= bind_required_branch(tree, "cc", &cc);
  ok &= bind_required_branch(tree, "PDGLep", &pdg_lep);
  ok &= bind_required_branch(tree, "ELep", &e_lep);
  ok &= bind_required_branch(tree, "nfsp", &n_fsp);
  ok &= bind_required_branch(tree, "px", px);
  ok &= bind_required_branch(tree, "py", py);
  ok &= bind_required_branch(tree, "pz", pz);
  ok &= bind_required_branch(tree, "E", energy);
  ok &= bind_required_branch(tree, "pdg", pdg);
  ok &= bind_required_branch(tree, "nvertp", &n_vertp);
  ok &= bind_required_branch(tree, "pdg_vert", pdg_vert);

  if (tree->GetBranch("Weight")) {
    tree->SetBranchAddress("Weight", &weight);
  } else {
    std::cout << "No Weight branch found; using weight = 1." << std::endl;
  }

  if (tree->GetBranch("fScaleFactor")) {
    tree->SetBranchAddress("fScaleFactor", &scale_factor);
  } else {
    std::cout << "No fScaleFactor branch found; using scale_factor = 1." << std::endl;
  }

  if (!ok) {
    input->Close();
    return;
  }

  std::vector<yield_row> rows = {
      {"t0_any_final_state_strangeness",
       "at least one final-state strange baryon or kaon"},
      {"t1_final_state_lambda",
       "at least one final-state Lambda"},
      {"t2_final_state_k_no_lambda",
       "final-state kaon and no final-state Lambda"},
      {"t3_final_state_sigma",
       "at least one final-state Sigma baryon"},
      {"t4_final_state_xi_or_omega",
       "at least one final-state Xi or Omega baryon"},
      {"s0_lambda_p_piminus_proxy",
       "final-state Lambda plus visible-proxy proton and pi-"},
      {"s1_cc_muon_tagged_lambda_proxy",
       "S0 plus charged-current visible-proxy muon"},
      {"s2_lambda_plus_visible_kaon_proxy",
       "S0 plus final-state visible-proxy kaon"},
      {"s3_lambda_low_extra_activity_proxy",
       "S0 plus extra visible kinetic energy below 0.15 GeV"},
      {"s4_lambda_visible_em_proxy",
       "S0 plus final-state gamma above 0.03 GeV"},
      {"s5_lambda_no_visible_muon_proxy",
       "S0 without a visible-proxy muon"},
      {"o0_primary_lambda_to_s0_proxy",
       "S0 event with primary vertex Lambda"},
      {"o1_primary_sigma0_to_s0_proxy",
       "S0 event with primary vertex Sigma0"},
      {"o2_primary_charged_sigma_to_s0_proxy",
       "S0 event with primary vertex Sigma+ or Sigma-"},
      {"o4_primary_ky_to_s0_proxy",
       "S0 event with primary vertex kaon and strange baryon"},
      {"o6_other_origin_to_s0_proxy",
       "S0 event not covered by the explicit primary labels"},
  };

  const Long64_t entries = tree->GetEntries();
  for (Long64_t entry = 0; entry < entries; ++entry) {
    tree->GetEntry(entry);

    const double event_weight = weight * scale_factor;
    const double xsec_weight = event_weight * units * argon_a;

    int n_final_lambda = 0;
    int n_final_kaon = 0;
    int n_final_sigma = 0;
    int n_final_xi_or_omega = 0;
    int n_final_strange_baryon = 0;
    int n_visible_proton = 0;
    int n_visible_piminus = 0;
    int n_visible_kaon = 0;
    int n_visible_gamma = 0;
    bool visible_muon = false;
    double extra_visible_energy = 0.0;

    int n_primary_lambda = 0;
    int n_primary_sigma0 = 0;
    int n_primary_charged_sigma = 0;
    int n_primary_kaon = 0;
    int n_primary_strange_baryon = 0;

    for (int i = 0; i < std::min(n_vertp, max_particles); ++i) {
      n_primary_lambda += pdg_vert[i] == 3122;
      n_primary_sigma0 += pdg_vert[i] == 3212;
      n_primary_charged_sigma += pdg_vert[i] == 3112 || pdg_vert[i] == 3222;
      n_primary_kaon += is_kaon(pdg_vert[i]);
      n_primary_strange_baryon += is_strange_baryon(pdg_vert[i]);
    }

    for (int i = 0; i < std::min(n_fsp, max_particles); ++i) {
      const double p = momentum(px[i], py[i], pz[i]);
      const int apdg = std::abs(pdg[i]);

      n_final_lambda += pdg[i] == 3122;
      n_final_kaon += is_kaon(pdg[i]);
      n_final_sigma += is_sigma(pdg[i]);
      n_final_xi_or_omega += is_xi_or_omega(pdg[i]);
      n_final_strange_baryon += is_strange_baryon(pdg[i]);
      n_visible_proton += pdg[i] == 2212 && p > 0.30;
      n_visible_piminus += pdg[i] == -211 && p > 0.07;
      n_visible_kaon += is_kaon(pdg[i]) && p > 0.10;
      n_visible_gamma += pdg[i] == 22 && energy[i] > 0.03;
      visible_muon = visible_muon || (apdg == 13 && p > 0.10);

      const bool is_invisible_neutrino = apdg == 12 || apdg == 14 || apdg == 16;
      const bool is_lambda_decay_proxy_daughter =
          (pdg[i] == 2212 && p > 0.30) || (pdg[i] == -211 && p > 0.07);
      if (!is_invisible_neutrino && !is_lambda_decay_proxy_daughter && pdg[i] != 3122) {
        extra_visible_energy += std::max(0.0, static_cast<double>(energy[i]) - mass_gev(pdg[i]));
      }
    }

    if (!visible_muon && std::abs(pdg_lep) == 13 && cc == 1 && e_lep > muon_mass_gev) {
      visible_muon = std::sqrt(e_lep * e_lep - muon_mass_gev * muon_mass_gev) > 0.10;
    }

    const bool any_final_strangeness = n_final_strange_baryon > 0 || n_final_kaon > 0;
    const bool final_lambda = n_final_lambda > 0;
    const bool s0_lambda_proxy = final_lambda && n_visible_proton > 0 && n_visible_piminus > 0;
    const bool s1_muon = s0_lambda_proxy && cc == 1 && visible_muon;
    const bool s2_kaon = s0_lambda_proxy && n_visible_kaon > 0;
    const bool s3_low_extra = s0_lambda_proxy && extra_visible_energy < 0.15;
    const bool s4_em = s0_lambda_proxy && n_visible_gamma > 0;
    const bool s5_no_muon = s0_lambda_proxy && !visible_muon;
    const bool primary_ky = n_primary_kaon > 0 && n_primary_strange_baryon > 0;
    const bool explicit_s0_origin =
        n_primary_lambda > 0 || n_primary_sigma0 > 0 || n_primary_charged_sigma > 0 || primary_ky;

    fill(rows[0], any_final_strangeness, event_weight, xsec_weight);
    fill(rows[1], final_lambda, event_weight, xsec_weight);
    fill(rows[2], n_final_kaon > 0 && !final_lambda, event_weight, xsec_weight);
    fill(rows[3], n_final_sigma > 0, event_weight, xsec_weight);
    fill(rows[4], n_final_xi_or_omega > 0, event_weight, xsec_weight);
    fill(rows[5], s0_lambda_proxy, event_weight, xsec_weight);
    fill(rows[6], s1_muon, event_weight, xsec_weight);
    fill(rows[7], s2_kaon, event_weight, xsec_weight);
    fill(rows[8], s3_low_extra, event_weight, xsec_weight);
    fill(rows[9], s4_em, event_weight, xsec_weight);
    fill(rows[10], s5_no_muon, event_weight, xsec_weight);
    fill(rows[11], s0_lambda_proxy && n_primary_lambda > 0, event_weight, xsec_weight);
    fill(rows[12], s0_lambda_proxy && n_primary_sigma0 > 0, event_weight, xsec_weight);
    fill(rows[13], s0_lambda_proxy && n_primary_charged_sigma > 0, event_weight, xsec_weight);
    fill(rows[14], s0_lambda_proxy && primary_ky, event_weight, xsec_weight);
    fill(rows[15], s0_lambda_proxy && !explicit_s0_origin, event_weight, xsec_weight);
  }

  TString output_dir = gSystem->DirName(output_csv.Data());
  if (output_dir != "." && output_dir != "") {
    gSystem->mkdir(output_dir.Data(), true);
  }

  write_rows(output_csv, sample_label, generator, knob, rows);

  std::cout << '\n'
            << "Sample: " << sample_label.Data() << '\n'
            << "Generator: " << generator.Data() << '\n'
            << "Knob: " << knob.Data() << '\n'
            << "Input: " << input_file.Data() << '\n'
            << "Output: " << output_csv.Data() << '\n'
            << "Entries: " << entries << '\n'
            << '\n'
            << "Weights: weighted_yield = weight * scale_factor; "
            << "xsec column also multiplies by 1e38 * A(Ar=40)." << '\n'
            << "Output: topology yield rows only; no histograms, canvases, ROOT output, or plots are produced." << '\n'
            << "Detector note: S0-S5 are truth-visible proxies from final-state PDGs. "
            << "No MicroBooNE containment, decay-vertex, or detached-vertex cuts are applied." << '\n'
            << '\n';
  print_rows(generator, knob, rows);

  input->Close();
}
