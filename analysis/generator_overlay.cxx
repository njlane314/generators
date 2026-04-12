#include <TCanvas.h>
#include <TFile.h>
#include <TLegend.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TH1D.h>
#include <TH2D.h>

#include <algorithm>
#include <iostream>
#include <vector>

using namespace std;

bool is_count_hist(TString name)
{
    return name.Contains("cutflow") || name.Contains("origin") || name.Contains("category") || name.Contains("strange");
}

void generator_overlay()
{
    gStyle->SetOptStat(0);

    TString output_dir = "analysis/output/";
    vector<TString> tags = {"GENIE", "NuWro", "GiBUU"};
    vector<TString> labels = {"GENIE", "NuWro", "GiBUU"};
    vector<int> colors = {kBlue + 8, kRed + 1, kGreen + 2};
    vector<TFile *> files;
    vector<TString> active_labels;
    vector<int> active_colors;

    for (int i = 0; i < int(tags.size()); ++i) {
        TString path = output_dir + "lambda_analyser_" + tags[i] + ".root";
        if (gSystem->AccessPathName(path)) continue;
        TFile *file = TFile::Open(path);
        if (!file || file->IsZombie()) continue;
        files.push_back(file);
        active_labels.push_back(labels[i]);
        active_colors.push_back(colors[i]);
    }

    vector<TString> hist_1d_names = {
        "enu_all", "enu_S0", "p_lambda", "costheta_lambda", "opening_angle_p_pi",
        "extra_visible_energy", "primary_strange", "final_strange", "origin", "category", "cutflow"
    };

    for (auto name : hist_1d_names) {
        vector<TH1D *> hists;
        double y_max = 0;
        TCanvas canvas("c", name, 1000, 760);
        TLegend legend(0.18, 0.74, 0.84, 0.88);
        legend.SetBorderSize(0);
        legend.SetNColumns(3);

        for (int i = 0; i < int(files.size()); ++i) {
            TH1D *source = (TH1D *)files[i]->Get(name);
            if (!source) continue;
            TH1D *hist = (TH1D *)source->Clone(name + "_" + active_labels[i]);
            hist->SetDirectory(0);
            hist->SetLineColor(active_colors[i]);
            hist->SetLineWidth(4);
            if (!is_count_hist(name) && hist->Integral("width") > 0) hist->Scale(1.0 / hist->Integral("width"));
            hist->GetYaxis()->SetTitle(is_count_hist(name) ? "Weighted events" : "Normalized events");
            y_max = max(y_max, hist->GetMaximum());
            hists.push_back(hist);
            legend.AddEntry(hist, active_labels[i], "l");
        }

        if (!hists.size()) continue;
        canvas.SetLeftMargin(0.14);
        canvas.SetBottomMargin(0.14);
        for (int i = 0; i < int(hists.size()); ++i) {
            hists[i]->GetYaxis()->SetRangeUser(0, y_max > 0 ? 1.25 * y_max : 1);
            hists[i]->Draw(i ? "hist same" : "hist");
        }
        legend.Draw();
        canvas.SaveAs(output_dir + "overlay_" + name + ".pdf");
        for (auto *hist : hists) delete hist;
    }

    vector<TString> hist_2d_names = {"p_lambda_costheta", "p_proton_p_piminus", "origin_category"};
    for (auto name : hist_2d_names) {
        for (int i = 0; i < int(files.size()); ++i) {
            TH2D *hist = (TH2D *)files[i]->Get(name);
            if (!hist) continue;
            TCanvas canvas("c2", name + "_" + active_labels[i], 1000, 760);
            canvas.SetRightMargin(0.16);
            hist->Draw("colz text");
            canvas.SaveAs(output_dir + "map_" + active_labels[i] + "_" + name + ".pdf");
        }
    }

    for (auto *file : files) file->Close();
}
