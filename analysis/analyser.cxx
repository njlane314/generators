#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVector3.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

class analyser {
public:
    analyser(TString input_file, TString tag, TString out_dir = "analysis/output")
        : input_path(input_file), sample_label(tag), output_dir(out_dir) {}

    void loop()
    {
        TFile *input_file = TFile::Open(input_path, "READ");
        if (!input_file || input_file->IsZombie()) return;

        TTree *tree = nullptr;
        input_file->GetObject("FlatTree_VARS", tree);
        if (!tree) {
            input_file->Close();
            return;
        }

        bind(tree);
        book();

        const Long64_t entries = tree->GetEntries();
        for (Long64_t entry = 0; entry < entries; ++entry) {
            tree->GetEntry(entry);
            fill_event();
        }

        write();
        input_file->Close();
    }

private:
    enum { max_particles = 100 };

    TString input_path, sample_label, output_dir;

    Int_t mode = 0, pdg_nu = 0, target = 0, target_a = 1, target_z = 0, pdg_lep = 0;
    Char_t cc = 0;
    Float_t enu_true = 0, e_lep = 0, cos_lep = 0, weight = 1;
    Double_t scale_factor = 1;
    Int_t n_fsp = 0, n_vertp = 0;
    Float_t px[max_particles] = {0}, py[max_particles] = {0}, pz[max_particles] = {0}, energy[max_particles] = {0};
    Float_t px_vert[max_particles] = {0}, py_vert[max_particles] = {0}, pz_vert[max_particles] = {0}, energy_vert[max_particles] = {0};
    Int_t pdg[max_particles] = {0}, pdg_vert[max_particles] = {0};

    TH1D *cutflow = nullptr, *origin = nullptr, *category = nullptr;
    TH1D *enu = nullptr, *enu_s0 = nullptr, *p_lam = nullptr, *cth_lam = nullptr, *opening_angle = nullptr;
    TH1D *primary = nullptr, *finals = nullptr, *extra = nullptr;
    TH2D *p_lam_cth = nullptr, *p_pip = nullptr, *origin_cat = nullptr;
    std::vector<TH1 *> hist_1d;
    std::vector<TH2 *> hist_2d;

    void bind(TTree *tree)
    {
        tree->SetBranchAddress("mode", &mode);
        tree->SetBranchAddress("cc", &cc);
        tree->SetBranchAddress("pdg_nu", &pdg_nu);
        tree->SetBranchAddress("enu_true", &enu_true);
        tree->SetBranchAddress("target", &target);
        tree->SetBranchAddress("target_a", &target_a);
        tree->SetBranchAddress("target_z", &target_z);
        tree->SetBranchAddress("pdg_lep", &pdg_lep);
        tree->SetBranchAddress("e_lep", &e_lep);
        tree->SetBranchAddress("cos_lep", &cos_lep);
        tree->SetBranchAddress("n_fsp", &n_fsp);
        tree->SetBranchAddress("px", px);
        tree->SetBranchAddress("py", py);
        tree->SetBranchAddress("pz", pz);
        tree->SetBranchAddress("energy", energy);
        tree->SetBranchAddress("pdg", pdg);
        tree->SetBranchAddress("n_vertp", &n_vertp);
        tree->SetBranchAddress("px_vert", px_vert);
        tree->SetBranchAddress("py_vert", py_vert);
        tree->SetBranchAddress("pz_vert", pz_vert);
        tree->SetBranchAddress("energy_vert", energy_vert);
        tree->SetBranchAddress("pdg_vert", pdg_vert);
        tree->SetBranchAddress("weight", &weight);
        tree->SetBranchAddress("scale_factor", &scale_factor);
    }

    TH1D *make_h1(const char *name, const char *title, int bins, double low, double high)
    {
        auto *hist = new TH1D(name, title, bins, low, high);
        hist->Sumw2();
        hist_1d.push_back(hist);
        return hist;
    }

    TH2D *make_h2(const char *name, const char *title, int x_bins, double x_low, double x_high, int y_bins, double y_low, double y_high)
    {
        auto *hist = new TH2D(name, title, x_bins, x_low, x_high, y_bins, y_low, y_high);
        hist->Sumw2();
        hist_2d.push_back(hist);
        return hist;
    }

    void book()
    {
        cutflow = make_h1("cutflow", ";selection;events", 8, 0.5, 8.5);
        label_bins(cutflow, {"all", "primary_strange", "final_lambda", "p_pi", "S0", "S1_mu", "S2_K", "S4_EM"});

        origin = make_h1("origin", ";origin;events", 8, -0.5, 7.5);
        label_bins(origin, {"O0_Lambda", "O1_Sigma0", "O2_Sigma_charged", "O3_Sigma_neutral", "O4_KY", "O5_Ypi", "O6_other", "O7_noLambda"});

        category = make_h1("category", ";category;events", 6, -0.5, 5.5);
        label_bins(category, {"S0", "S1_mu", "S2_K", "S3_lowExtra", "S4_EM", "S5_noMuon"});

        primary = make_h1("primary_strange", ";primary strange;particles", 8, -0.5, 7.5);
        finals = make_h1("final_strange", ";final strange;particles", 8, -0.5, 7.5);
        label_bins(primary, {"Lambda", "Sigma-", "Sigma0", "Sigma+", "Xi", "Omega", "K", "other"});
        label_bins(finals, {"Lambda", "Sigma-", "Sigma0", "Sigma+", "Xi", "Omega", "K", "other"});

        enu = make_h1("enu_all", ";E_{#nu} [GeV];events / GeV", 80, 0, 10);
        enu_s0 = make_h1("enu_S0", ";E_{#nu} [GeV];events / GeV", 80, 0, 10);
        p_lam = make_h1("p_lambda", ";p_{#Lambda} [GeV];particles / GeV", 80, 0, 4);
        cth_lam = make_h1("costheta_lambda", ";cos#theta_{#Lambda};particles", 50, -1, 1);
        opening_angle = make_h1("opening_angle_p_pi", ";opening angle p-#pi^{-} [rad];events / rad", 64, 0, 3.2);
        extra = make_h1("extra_visible_energy", ";extra visible energy [GeV];events / GeV", 80, 0, 2);
        p_lam_cth = make_h2("p_lambda_costheta", ";p_{#Lambda} [GeV];cos#theta_{#Lambda}", 60, 0, 4, 50, -1, 1);
        p_pip = make_h2("p_proton_p_piminus", ";p_{p} [GeV];p_{#pi^{-}} [GeV]", 60, 0, 2, 60, 0, 1.2);
        origin_cat = make_h2("origin_category", ";origin;category", 8, -0.5, 7.5, 6, -0.5, 5.5);
    }

    void fill_event()
    {
        const double event_weight = weight * scale_factor * std::max(1, target_a) * 1e38;
        int n_lam_0 = 0, n_sig_0 = 0, n_sig_c = 0, n_k_0 = 0, n_y_0 = 0;
        int n_lam = 0, n_p = 0, n_pi = 0, n_k = 0, n_gamma = 0, i_lam = -1, i_p = -1, i_pi = -1;
        double p_best_lam = -1, p_best_p = -1, p_best_pi = -1, e_extra = 0;
        bool mu = false;

        cutflow->Fill(1, event_weight);
        enu->Fill(enu_true, event_weight);

        for (int i = 0; i < std::min(n_vertp, max_particles); ++i) {
            const int species_id = species_bin(pdg_vert[i]);
            if (species_id < 0) continue;
            primary->Fill(species_id, event_weight);
            n_lam_0 += pdg_vert[i] == 3122;
            n_sig_0 += pdg_vert[i] == 3212;
            n_sig_c += pdg_vert[i] == 3112 || pdg_vert[i] == 3222;
            n_k_0 += is_k(pdg_vert[i]);
            n_y_0 += is_y(pdg_vert[i]);
        }

        for (int i = 0; i < std::min(n_fsp, max_particles); ++i) {
            const int species_id = species_bin(pdg[i]);
            const double p = momentum(i);
            if (species_id >= 0) finals->Fill(species_id, event_weight);
            if (pdg[i] == 3122 && p > p_best_lam) { n_lam++; i_lam = i; p_best_lam = p; }
            if (pdg[i] == 2212 && p > 0.30 && p > p_best_p) { n_p++; i_p = i; p_best_p = p; }
            if (pdg[i] == -211 && p > 0.07 && p > p_best_pi) { n_pi++; i_pi = i; p_best_pi = p; }
            if (is_k(pdg[i]) && p > 0.10) n_k++;
            if (pdg[i] == 22 && energy[i] > 0.03) n_gamma++;
            if (std::abs(pdg[i]) == 13 && p > 0.10) mu = true;
        }

        if (!mu && std::abs(pdg_lep) == 13 && cc == 1 && e_lep > 0.105658) {
            mu = std::sqrt(e_lep * e_lep - 0.105658 * 0.105658) > 0.10;
        }

        if (n_y_0 || n_k_0) cutflow->Fill(2, event_weight);
        if (n_lam) cutflow->Fill(3, event_weight);
        if (n_p && n_pi) cutflow->Fill(4, event_weight);

        if (i_lam >= 0) {
            TVector3 lam_vec(px[i_lam], py[i_lam], pz[i_lam]);
            p_lam->Fill(lam_vec.Mag(), event_weight);
            cth_lam->Fill(lam_vec.CosTheta(), event_weight);
            p_lam_cth->Fill(lam_vec.Mag(), lam_vec.CosTheta(), event_weight);
        }

        if (i_p >= 0 && i_pi >= 0) {
            TVector3 p_vec(px[i_p], py[i_p], pz[i_p]), pi_vec(px[i_pi], py[i_pi], pz[i_pi]);
            p_pip->Fill(p_vec.Mag(), pi_vec.Mag(), event_weight);
            opening_angle->Fill(p_vec.Angle(pi_vec), event_weight);
            e_extra = extra_energy(i_p, i_pi);
        }

        const bool s0 = n_p && n_pi && (n_lam || n_y_0);
        if (!s0) return;

        const int origin_id = origin_bin(n_lam_0, n_sig_0, n_sig_c, n_k_0, n_y_0, n_lam);
        cutflow->Fill(5, event_weight);
        origin->Fill(origin_id, event_weight);
        category->Fill(0, event_weight);
        origin_cat->Fill(origin_id, 0, event_weight);
        enu_s0->Fill(enu_true, event_weight);
        extra->Fill(e_extra, event_weight);

        if (mu) fill_cat(1, origin_id, event_weight);
        if (n_k) fill_cat(2, origin_id, event_weight);
        if (e_extra < 0.15) fill_cat(3, origin_id, event_weight);
        if (n_gamma) fill_cat(4, origin_id, event_weight);
        if (!mu) fill_cat(5, origin_id, event_weight);
    }

    void fill_cat(int cat, int origin_id, double event_weight)
    {
        category->Fill(cat, event_weight);
        origin_cat->Fill(origin_id, cat, event_weight);
        if (cat == 1) cutflow->Fill(6, event_weight);
        if (cat == 2) cutflow->Fill(7, event_weight);
        if (cat == 4) cutflow->Fill(8, event_weight);
    }

    double extra_energy(int proton_idx, int pion_idx)
    {
        double extra_e = 0;
        for (int i = 0; i < std::min(n_fsp, max_particles); ++i) {
            if (i == proton_idx || i == pion_idx) continue;
            if (std::abs(pdg[i]) == 12 || std::abs(pdg[i]) == 14 || std::abs(pdg[i]) == 16) continue;
            extra_e += std::max(0.0, double(energy[i]) - mass(pdg[i]));
        }
        return extra_e;
    }

    void write()
    {
        gSystem->mkdir(output_dir, true);
        TFile output_file(output_dir + "/lambda_analyser_" + clean_label(sample_label) + ".root", "RECREATE");
        for (auto *hist : hist_1d) {
            scale_density(hist);
            hist->Write();
        }
        for (auto *hist : hist_2d) hist->Write();
        output_file.Close();
    }

    double momentum(int i) { return std::sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i]); }
    static bool is_k(int code) { code = std::abs(code); return code == 321 || code == 311 || code == 310 || code == 130; }
    static bool is_y(int code) { code = std::abs(code); return code == 3122 || code == 3112 || code == 3212 || code == 3222 || code == 3312 || code == 3322 || code == 3334; }
    static int origin_bin(int lam, int sig_0, int sig_c, int kaon, int hyperon, int final_lam)
    {
        if (kaon && hyperon) return 4;
        if (sig_0) return 1;
        if (sig_c) return 2;
        if (lam) return 0;
        if (hyperon) return 5;
        if (kaon || final_lam) return 6;
        return 7;
    }
    static int species_bin(int code)
    {
        if (code == 3122) return 0;
        if (code == 3112) return 1;
        if (code == 3212) return 2;
        if (code == 3222) return 3;
        if (std::abs(code) == 3312 || std::abs(code) == 3322) return 4;
        if (std::abs(code) == 3334) return 5;
        if (is_k(code)) return 6;
        if (is_y(code)) return 7;
        return -1;
    }
    static double mass(int code)
    {
        code = std::abs(code);
        if (code == 2212) return 0.938272;
        if (code == 2112) return 0.939565;
        if (code == 211) return 0.139570;
        if (code == 321) return 0.493677;
        if (code == 311 || code == 310 || code == 130) return 0.497611;
        if (code == 3122) return 1.115683;
        if (code == 13) return 0.105658;
        if (code == 11) return 0.000511;
        return 0;
    }
    static TString clean_label(TString text) { text.ReplaceAll("/", "_"); text.ReplaceAll(" ", "_"); return text; }
    static void label_bins(TH1 *hist, const std::vector<const char *> &labels)
    {
        for (int i = 0; i < int(labels.size()); ++i) hist->GetXaxis()->SetBinLabel(i + 1, labels[i]);
    }
    static void scale_density(TH1 *hist)
    {
        TString name = hist->GetName();
        if (name.Contains("cutflow") || name.Contains("origin") || name.Contains("category") || name.Contains("strange")) return;
        for (int i = 1; i <= hist->GetNbinsX(); ++i) {
            hist->SetBinContent(i, hist->GetBinContent(i) / hist->GetBinWidth(i));
            hist->SetBinError(i, hist->GetBinError(i) / hist->GetBinWidth(i));
        }
    }
};
