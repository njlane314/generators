#include "beam_kinematic_envelope.cxx"

#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

namespace AnaGeneratorEnvelope {

using AnaBeamEnvelope::FluxEnv;
using AnaBeamEnvelope::StatusRow;
using AnaBeamEnvelope::VariableSpec;

TString group_key(const StatusRow& row, bool include_fsi_state)
{
    TString key = row.beam_mode + "|" + row.beam_species + "|" + row.interaction;
    if (include_fsi_state) {
        key += "|";
        key += row.fsi_state;
    }
    return key;
}

TString row_label(const StatusRow& row, bool include_fsi_state)
{
    TString label = row.generator + " ";
    label += row.version;
    label += " ";
    label += row.variation;
    if (row.knob != "" && !row.variation.Contains(row.knob)) {
        label += " ";
        label += row.knob;
    }
    if (!include_fsi_state && row.fsi_state != "") {
        label += " ";
        label += row.fsi_state;
    }
    return label;
}

TGraphAsymmErrors* make_envelope_band(const std::vector<TH1D*>& hists)
{
    if (hists.empty()) return nullptr;
    TH1D* ref = hists.front();
    TGraphAsymmErrors* band = new TGraphAsymmErrors(ref->GetNbinsX());
    band->SetFillColorAlpha(kGray + 1, 0.30);
    band->SetLineColor(kGray + 1);
    for (int bin = 1; bin <= ref->GetNbinsX(); ++bin) {
        double low = hists.front()->GetBinContent(bin);
        double high = low;
        for (TH1D* hist : hists) {
            low = std::min(low, hist->GetBinContent(bin));
            high = std::max(high, hist->GetBinContent(bin));
        }
        const double center = 0.5 * (low + high);
        const double x = ref->GetXaxis()->GetBinCenter(bin);
        const double ex = 0.5 * ref->GetXaxis()->GetBinWidth(bin);
        band->SetPoint(bin - 1, x, center);
        band->SetPointError(bin - 1, ex, ex, center - low, high - center);
    }
    return band;
}

void style_hist(TH1D* hist, int index)
{
    const int colors[] = {
        kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1, kOrange + 7,
        kCyan + 2, kViolet + 5, kGray + 2, kBlack, kSpring + 5
    };
    hist->SetStats(false);
    hist->SetLineColor(colors[index % (sizeof(colors) / sizeof(colors[0]))]);
    hist->SetMarkerColor(hist->GetLineColor());
    hist->SetLineWidth(2);
}

void write_bin_csv(const TString& path,
                   const TString& key,
                   const VariableSpec& variable,
                   const TString& selection,
                   const std::vector<TString>& labels,
                   const std::vector<TH1D*>& hists)
{
    std::ofstream out(path.Data());
    if (!out) return;
    out << "group,selection,variable,bin,low,high,envelope_min,envelope_max";
    for (const TString& label : labels) out << "," << AnaBeamEnvelope::csv(label).Data();
    out << "\n";
    out << std::setprecision(12);
    for (int bin = 1; bin <= hists.front()->GetNbinsX(); ++bin) {
        double low = hists.front()->GetBinContent(bin);
        double high = low;
        for (TH1D* hist : hists) {
            low = std::min(low, hist->GetBinContent(bin));
            high = std::max(high, hist->GetBinContent(bin));
        }
        out << AnaBeamEnvelope::csv(key).Data() << ','
            << AnaBeamEnvelope::csv(selection).Data() << ','
            << AnaBeamEnvelope::csv(variable.id).Data() << ','
            << bin << ','
            << hists.front()->GetXaxis()->GetBinLowEdge(bin) << ','
            << hists.front()->GetXaxis()->GetBinUpEdge(bin) << ','
            << low << ',' << high;
        for (TH1D* hist : hists) out << ',' << hist->GetBinContent(bin);
        out << '\n';
    }
}

void draw_plot(const TString& path_base,
               const TString& key,
               const VariableSpec& variable,
               const TString& selection,
               const std::vector<TString>& labels,
               const std::vector<TH1D*>& hists)
{
    TString dir = gSystem->DirName(path_base.Data());
    gSystem->mkdir(dir.Data(), true);
    TCanvas canvas("generator_envelope_canvas", "generator_envelope_canvas", 1050, 760);
    canvas.SetLeftMargin(0.14);
    canvas.SetBottomMargin(0.14);
    canvas.SetTopMargin(0.08);
    canvas.SetRightMargin(0.05);

    TH1D* frame = static_cast<TH1D*>(hists.front()->Clone("generator_envelope_frame"));
    frame->Reset();
    frame->SetDirectory(nullptr);
    frame->SetStats(false);
    TString title = key;
    title.ReplaceAll("|", " ");
    title += ";";
    title += variable.axis_title;
    title += ";weighted yield [10^{-38} cm^{2}/Ar]";
    frame->SetTitle(title);
    double ymax = 0.0;
    for (TH1D* hist : hists) ymax = std::max(ymax, hist->GetMaximum());
    frame->GetYaxis()->SetRangeUser(0.0, ymax > 0.0 ? 1.40 * ymax : 1.0);
    frame->Draw("axis");

    TGraphAsymmErrors* band = make_envelope_band(hists);
    if (band) band->Draw("2 same");
    for (int i = 0; i < int(hists.size()); ++i) {
        style_hist(hists[i], i);
        hists[i]->Draw("hist e same");
    }
    frame->Draw("axis same");

    TLegend legend(0.48, 0.64, 0.92, 0.90);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetHeader(selection.Data());
    if (band) legend.AddEntry(band, "generator/config envelope", "f");
    for (int i = 0; i < int(hists.size()); ++i) legend.AddEntry(hists[i], labels[i], "l");
    legend.Draw();

    canvas.SaveAs(path_base + ".pdf");
    canvas.SaveAs(path_base + ".png");
    delete band;
    delete frame;
}

}  // namespace AnaGeneratorEnvelope

void generator_kinematic_envelope(
    TString status_csv = "analysis/output/matrix/generator_matrix_status.csv",
    TString output_dir = "analysis/output/generator_kinematic_envelope",
    TString variables = "enu q2 w hyperon_p",
    TString selection = "detector_visible_lambda",
    TString flux_file = "analysis/flux/microboone_numi_flux_5mev.root",
    double flux_floor_fraction = 0.0,
    double proposal_emin_gev = 0.0,
    double proposal_emax_gev = 10.0,
    double proton_threshold_gev = 0.30,
    double piminus_threshold_gev = 0.07,
    bool include_fsi_state_in_group = false
)
{
    using namespace AnaGeneratorEnvelope;
    gStyle->SetOptStat(0);

    const std::vector<StatusRow> rows = AnaBeamEnvelope::read_status(status_csv);
    const std::vector<VariableSpec> specs = AnaBeamEnvelope::variables_from_arg(variables);
    if (rows.empty() || specs.empty()) {
        std::cerr << "No generator-envelope rows or variables found." << std::endl;
        return;
    }
    gSystem->mkdir(output_dir.Data(), true);

    std::map<std::string, std::vector<StatusRow>> groups;
    for (const StatusRow& row : rows) groups[group_key(row, include_fsi_state_in_group).Data()].push_back(row);

    std::map<std::string, FluxEnv> flux_cache;
    std::ofstream summary(AnaBeamEnvelope::join_path(output_dir, "generator_kinematic_envelope_summary.csv").Data());
    summary << "group,selection,variable,n_members,envelope_min_integral,envelope_max_integral,"
            << "plot_base,csv_file\n";
    summary << std::setprecision(12);

    int outputs = 0;
    for (const auto& item : groups) {
        const TString key = item.first.c_str();
        for (const VariableSpec& variable : specs) {
            std::vector<TH1D*> hists;
            std::vector<TString> labels;
            for (const StatusRow& row : item.second) {
                TString hist_name = "gen_";
                hist_name += row.sample;
                hist_name += "_";
                hist_name += variable.id;
                TString hist_title = ";";
                hist_title += variable.axis_title;
                hist_title += ";weighted yield";
                TH1D* hist = new TH1D(AnaBeamEnvelope::safe_name(hist_name),
                                      hist_title, variable.bins, variable.min, variable.max);
                hist->SetDirectory(nullptr);
                hist->Sumw2();
                AnaBeamEnvelope::fill_hist(row, hist, variable, selection, flux_file,
                                           flux_floor_fraction, proposal_emin_gev,
                                           proposal_emax_gev, proton_threshold_gev,
                                           piminus_threshold_gev, flux_cache);
                if (hist->Integral() > 0.0) {
                    hists.push_back(hist);
                    labels.push_back(row_label(row, include_fsi_state_in_group));
                } else {
                    delete hist;
                }
            }
            if (hists.size() < 2) {
                for (TH1D* hist : hists) delete hist;
                continue;
            }

            TString base_name = "generator_envelope_";
            base_name += key;
            base_name += "_";
            base_name += variable.id;
            base_name += "_";
            base_name += selection;
            const TString base = AnaBeamEnvelope::join_path(output_dir, AnaBeamEnvelope::safe_name(base_name));
            const TString csv_path = base + ".csv";
            write_bin_csv(csv_path, key, variable, selection, labels, hists);
            draw_plot(base, key, variable, selection, labels, hists);

            double min_integral = hists.front()->Integral();
            double max_integral = min_integral;
            for (TH1D* hist : hists) {
                min_integral = std::min(min_integral, hist->Integral());
                max_integral = std::max(max_integral, hist->Integral());
            }
            summary << AnaBeamEnvelope::csv(key).Data() << ','
                    << AnaBeamEnvelope::csv(selection).Data() << ','
                    << AnaBeamEnvelope::csv(variable.id).Data() << ','
                    << hists.size() << ','
                    << min_integral << ',' << max_integral << ','
                    << AnaBeamEnvelope::csv(base).Data() << ','
                    << AnaBeamEnvelope::csv(csv_path).Data() << '\n';
            ++outputs;
            for (TH1D* hist : hists) delete hist;
        }
    }

    AnaBeamEnvelope::clear_flux_cache(flux_cache);
    std::cout << "Wrote " << outputs << " generator kinematic envelope outputs under "
              << output_dir.Data() << std::endl;
}
