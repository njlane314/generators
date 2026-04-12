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

const int kMaxParticles = 100;
const double kArgonA = 40.0;
const double kUnits = 1.0e38;
const double kMuonMassGeV = 0.1056583745;

struct YieldRow {
  TString name;
  TString definition;
  Long64_t raw = 0;
  double weighted = 0.0;
  double xsec_weighted = 0.0;
};

bool BindRequiredBranch(TTree* tree, const char* name, void* address)
{
  if (!tree->GetBranch(name)) {
    std::cerr << "Missing required FlatTree_VARS branch: " << name << std::endl;
    return false;
  }
  tree->SetBranchAddress(name, address);
  return true;
}

bool IsKaon(int pdg)
{
  const int apdg = std::abs(pdg);
  return apdg == 321 || apdg == 311 || pdg == 310 || pdg == 130;
}

bool IsSigma(int pdg)
{
  const int apdg = std::abs(pdg);
  return apdg == 3112 || apdg == 3212 || apdg == 3222;
}

bool IsXiOrOmega(int pdg)
{
  const int apdg = std::abs(pdg);
  return apdg == 3312 || apdg == 3322 || apdg == 3334;
}

bool IsStrangeBaryon(int pdg)
{
  return std::abs(pdg) == 3122 || IsSigma(pdg) || IsXiOrOmega(pdg);
}

double Momentum(double px, double py, double pz)
{
  return std::sqrt(px * px + py * py + pz * pz);
}

double MassGeV(int pdg)
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
  if (apdg == 13) return kMuonMassGeV;
  if (apdg == 11) return 0.00051099895;
  return 0.0;
}

void Fill(YieldRow& row, bool pass, double weighted, double xsec_weighted)
{
  if (!pass) return;
  ++row.raw;
  row.weighted += weighted;
  row.xsec_weighted += xsec_weighted;
}

TString CsvEscape(TString value)
{
  value.ReplaceAll("\"", "\"\"");
  TString escaped = "\"";
  escaped += value;
  escaped += "\"";
  return escaped;
}

TString InferGenerator(TString sample)
{
  sample.ToLower();
  if (sample.Contains("nuwro")) return "NuWro";
  if (sample.Contains("genie")) return "GENIE";
  if (sample.Contains("gibuu")) return "GiBUU";
  return "unspecified";
}

TString InferKnob(TString sample)
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

void WriteRows(const TString& output_file,
               const TString& sample_label,
               const TString& generator,
               const TString& knob,
               const std::vector<YieldRow>& rows)
{
  std::ofstream csv(output_file.Data());
  csv << "sample,generator,knob,category,raw_events,weighted_yield,"
      << "xsec_weighted_1e38_cm2_per_Ar,definition\n";
  csv << std::setprecision(12);
  for (const YieldRow& row : rows) {
    const TString definition = CsvEscape(row.definition);
    const TString escaped_sample = CsvEscape(sample_label);
    const TString escaped_generator = CsvEscape(generator);
    const TString escaped_knob = CsvEscape(knob);
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

void PrintRows(const TString& generator, const TString& knob, const std::vector<YieldRow>& rows)
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
  for (const YieldRow& row : rows) {
    std::cout << std::left << std::setw(36) << row.name.Data()
              << std::right << std::setw(12) << row.raw
              << std::setw(18) << row.weighted
              << std::setw(24) << row.xsec_weighted << '\n';
  }
}

}  // namespace

void strangeness_topology_yields(TString input_file,
                                 TString output_csv = "",
                                 TString sample_label = "",
                                 TString generator = "",
                                 TString knob = "")
{
  if (output_csv == "") {
    output_csv = input_file;
    output_csv.ReplaceAll(".flat.root", ".strangeness_topology_yields.csv");
    if (output_csv == input_file) {
      output_csv += ".strangeness_topology_yields.csv";
    }
  }

  if (sample_label == "") {
    sample_label = gSystem->BaseName(input_file.Data());
  }
  if (generator == "") {
    generator = InferGenerator(sample_label);
  }
  if (knob == "") {
    knob = InferKnob(sample_label);
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
  Int_t PDGLep = 0;
  Float_t ELep = 0.0;
  Int_t nfsp = 0;
  Float_t px[kMaxParticles] = {0.0};
  Float_t py[kMaxParticles] = {0.0};
  Float_t pz[kMaxParticles] = {0.0};
  Float_t E[kMaxParticles] = {0.0};
  Int_t pdg[kMaxParticles] = {0};
  Int_t nvertp = 0;
  Int_t pdg_vert[kMaxParticles] = {0};
  Float_t Weight = 1.0;
  Double_t fScaleFactor = 1.0;

  bool ok = true;
  ok &= BindRequiredBranch(tree, "cc", &cc);
  ok &= BindRequiredBranch(tree, "PDGLep", &PDGLep);
  ok &= BindRequiredBranch(tree, "ELep", &ELep);
  ok &= BindRequiredBranch(tree, "nfsp", &nfsp);
  ok &= BindRequiredBranch(tree, "px", px);
  ok &= BindRequiredBranch(tree, "py", py);
  ok &= BindRequiredBranch(tree, "pz", pz);
  ok &= BindRequiredBranch(tree, "E", E);
  ok &= BindRequiredBranch(tree, "pdg", pdg);
  ok &= BindRequiredBranch(tree, "nvertp", &nvertp);
  ok &= BindRequiredBranch(tree, "pdg_vert", pdg_vert);

  if (tree->GetBranch("Weight")) {
    tree->SetBranchAddress("Weight", &Weight);
  } else {
    std::cout << "No Weight branch found; using Weight = 1." << std::endl;
  }

  if (tree->GetBranch("fScaleFactor")) {
    tree->SetBranchAddress("fScaleFactor", &fScaleFactor);
  } else {
    std::cout << "No fScaleFactor branch found; using fScaleFactor = 1." << std::endl;
  }

  if (!ok) {
    input->Close();
    return;
  }

  std::vector<YieldRow> rows = {
      {"T0_any_final_state_strangeness",
       "at least one final-state strange baryon or kaon"},
      {"T1_final_state_Lambda",
       "at least one final-state Lambda"},
      {"T2_final_state_K_no_Lambda",
       "final-state kaon and no final-state Lambda"},
      {"T3_final_state_Sigma",
       "at least one final-state Sigma baryon"},
      {"T4_final_state_Xi_or_Omega",
       "at least one final-state Xi or Omega baryon"},
      {"S0_Lambda_p_piminus_proxy",
       "final-state Lambda plus visible-proxy proton and pi-"},
      {"S1_CC_muon_tagged_Lambda_proxy",
       "S0 plus charged-current visible-proxy muon"},
      {"S2_Lambda_plus_visible_kaon_proxy",
       "S0 plus final-state visible-proxy kaon"},
      {"S3_Lambda_low_extra_activity_proxy",
       "S0 plus extra visible kinetic energy below 0.15 GeV"},
      {"S4_Lambda_visible_EM_proxy",
       "S0 plus final-state gamma above 0.03 GeV"},
      {"S5_Lambda_no_visible_muon_proxy",
       "S0 without a visible-proxy muon"},
      {"O0_primary_Lambda_to_S0_proxy",
       "S0 event with primary vertex Lambda"},
      {"O1_primary_Sigma0_to_S0_proxy",
       "S0 event with primary vertex Sigma0"},
      {"O2_primary_charged_Sigma_to_S0_proxy",
       "S0 event with primary vertex Sigma+ or Sigma-"},
      {"O4_primary_KY_to_S0_proxy",
       "S0 event with primary vertex kaon and strange baryon"},
      {"O6_other_origin_to_S0_proxy",
       "S0 event not covered by the explicit primary labels"},
  };

  const Long64_t entries = tree->GetEntries();
  for (Long64_t entry = 0; entry < entries; ++entry) {
    tree->GetEntry(entry);

    const double event_weight = Weight * fScaleFactor;
    const double xsec_weight = event_weight * kUnits * kArgonA;

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

    for (int i = 0; i < std::min(nvertp, kMaxParticles); ++i) {
      n_primary_lambda += pdg_vert[i] == 3122;
      n_primary_sigma0 += pdg_vert[i] == 3212;
      n_primary_charged_sigma += pdg_vert[i] == 3112 || pdg_vert[i] == 3222;
      n_primary_kaon += IsKaon(pdg_vert[i]);
      n_primary_strange_baryon += IsStrangeBaryon(pdg_vert[i]);
    }

    for (int i = 0; i < std::min(nfsp, kMaxParticles); ++i) {
      const double p = Momentum(px[i], py[i], pz[i]);
      const int apdg = std::abs(pdg[i]);

      n_final_lambda += pdg[i] == 3122;
      n_final_kaon += IsKaon(pdg[i]);
      n_final_sigma += IsSigma(pdg[i]);
      n_final_xi_or_omega += IsXiOrOmega(pdg[i]);
      n_final_strange_baryon += IsStrangeBaryon(pdg[i]);
      n_visible_proton += pdg[i] == 2212 && p > 0.30;
      n_visible_piminus += pdg[i] == -211 && p > 0.07;
      n_visible_kaon += IsKaon(pdg[i]) && p > 0.10;
      n_visible_gamma += pdg[i] == 22 && E[i] > 0.03;
      visible_muon = visible_muon || (apdg == 13 && p > 0.10);

      const bool is_invisible_neutrino = apdg == 12 || apdg == 14 || apdg == 16;
      const bool is_lambda_decay_proxy_daughter =
          (pdg[i] == 2212 && p > 0.30) || (pdg[i] == -211 && p > 0.07);
      if (!is_invisible_neutrino && !is_lambda_decay_proxy_daughter && pdg[i] != 3122) {
        extra_visible_energy += std::max(0.0, static_cast<double>(E[i]) - MassGeV(pdg[i]));
      }
    }

    if (!visible_muon && std::abs(PDGLep) == 13 && cc == 1 && ELep > kMuonMassGeV) {
      visible_muon = std::sqrt(ELep * ELep - kMuonMassGeV * kMuonMassGeV) > 0.10;
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

    Fill(rows[0], any_final_strangeness, event_weight, xsec_weight);
    Fill(rows[1], final_lambda, event_weight, xsec_weight);
    Fill(rows[2], n_final_kaon > 0 && !final_lambda, event_weight, xsec_weight);
    Fill(rows[3], n_final_sigma > 0, event_weight, xsec_weight);
    Fill(rows[4], n_final_xi_or_omega > 0, event_weight, xsec_weight);
    Fill(rows[5], s0_lambda_proxy, event_weight, xsec_weight);
    Fill(rows[6], s1_muon, event_weight, xsec_weight);
    Fill(rows[7], s2_kaon, event_weight, xsec_weight);
    Fill(rows[8], s3_low_extra, event_weight, xsec_weight);
    Fill(rows[9], s4_em, event_weight, xsec_weight);
    Fill(rows[10], s5_no_muon, event_weight, xsec_weight);
    Fill(rows[11], s0_lambda_proxy && n_primary_lambda > 0, event_weight, xsec_weight);
    Fill(rows[12], s0_lambda_proxy && n_primary_sigma0 > 0, event_weight, xsec_weight);
    Fill(rows[13], s0_lambda_proxy && n_primary_charged_sigma > 0, event_weight, xsec_weight);
    Fill(rows[14], s0_lambda_proxy && primary_ky, event_weight, xsec_weight);
    Fill(rows[15], s0_lambda_proxy && !explicit_s0_origin, event_weight, xsec_weight);
  }

  TString output_dir = gSystem->DirName(output_csv.Data());
  if (output_dir != "." && output_dir != "") {
    gSystem->mkdir(output_dir.Data(), true);
  }

  WriteRows(output_csv, sample_label, generator, knob, rows);

  std::cout << '\n'
            << "Sample: " << sample_label.Data() << '\n'
            << "Generator: " << generator.Data() << '\n'
            << "Knob: " << knob.Data() << '\n'
            << "Input: " << input_file.Data() << '\n'
            << "Output: " << output_csv.Data() << '\n'
            << "Entries: " << entries << '\n'
            << '\n'
            << "Weights: weighted_yield = Weight * fScaleFactor; "
            << "xsec column also multiplies by 1e38 * A(Ar=40)." << '\n'
            << "Output: topology yield rows only; no histograms, canvases, ROOT output, or plots are produced." << '\n'
            << "Detector note: S0-S5 are truth-visible proxies from final-state PDGs. "
            << "No MicroBooNE containment, decay-vertex, or detached-vertex cuts are applied." << '\n'
            << '\n';
  PrintRows(generator, knob, rows);

  input->Close();
}
