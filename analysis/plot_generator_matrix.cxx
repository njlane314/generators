#include <TAxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace {

struct status_row {
    int plan_line = 0;
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
    TString primary_status;
    TString nuclear_exit_status;
};

struct hist_spec {
    TString name;
    TString suffix;
    bool normalize = false;
};

std::vector<std::string> split_csv(const std::string& line)
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

TString beam_polarity_from_mode(TString beam_mode)
{
    beam_mode.ToLower();
    if (beam_mode.Contains("rhc")) return "rhc";
    if (beam_mode.Contains("fhc")) return "fhc";
    if (beam_mode == "combined" || beam_mode == "both" || beam_mode == "all") return "combined";
    return beam_mode;
}

TString stage_path(const status_row& row, const TString& input_dir, const TString& stage)
{
    const TString sample = clean(row.sample);
    const TString beam = beam_polarity_from_mode(row.beam_mode);
    if (stage == "primary") return join_path(input_dir, "primary_mechanism_" + sample + "_" + beam + ".root");
    if (stage == "exit") return join_path(input_dir, "nuclear_exit_" + sample + "_" + beam + ".root");
    return "";
}

std::vector<status_row> read_status(const TString& status_csv)
{
    std::ifstream input(status_csv.Data());
    std::vector<status_row> rows;
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
        const std::vector<std::string> f = split_csv(line);
        if (f.size() < 13) continue;
        status_row row;
        row.plan_line = std::atoi(f[0].c_str());
        row.generator = f[1].c_str();
        row.version = f[2].c_str();
        row.variation = f[3].c_str();
        row.knob = f[4].c_str();
        row.beam_mode = f[5].c_str();
        row.beam_species = f[6].c_str();
        row.interaction = f[7].c_str();
        row.fsi_state = f[8].c_str();
        row.sample = f[9].c_str();
        row.input_file = f[10].c_str();
        row.primary_status = f[11].c_str();
        row.nuclear_exit_status = f[12].c_str();
        rows.push_back(row);
    }
    return rows;
}

TH1* clone_hist1(const TString& path, const TString& hist_name, const TString& clone_name)
{
    TFile* file = TFile::Open(path, "READ");
    if (!file || file->IsZombie()) {
        if (file) file->Close();
        return nullptr;
    }
    TH1* source = dynamic_cast<TH1*>(file->Get(hist_name));
    if (!source) {
        file->Close();
        return nullptr;
    }
    TH1* clone = dynamic_cast<TH1*>(source->Clone(clone_name));
    if (clone) clone->SetDirectory(nullptr);
    file->Close();
    return clone;
}

TH2* clone_hist2(const TString& path, const TString& hist_name, const TString& clone_name)
{
    TFile* file = TFile::Open(path, "READ");
    if (!file || file->IsZombie()) {
        if (file) file->Close();
        return nullptr;
    }
    TH2* source = dynamic_cast<TH2*>(file->Get(hist_name));
    if (!source) {
        file->Close();
        return nullptr;
    }
    TH2* clone = dynamic_cast<TH2*>(source->Clone(clone_name));
    if (clone) clone->SetDirectory(nullptr);
    file->Close();
    return clone;
}

void style_hist(TH1* hist, int color, int marker)
{
    hist->SetStats(false);
    hist->SetLineColor(color);
    hist->SetMarkerColor(color);
    hist->SetMarkerStyle(marker);
    hist->SetLineWidth(3);
}

void normalize_shape(TH1* hist)
{
    const double integral = hist->Integral("width");
    if (integral > 0.0) hist->Scale(1.0 / integral);
}

std::vector<hist_spec> primary_hist_specs()
{
    return {
        {"topology_cutflow", "weighted", false},
        {"primary_origin", "weighted", false},
        {"interaction_mode", "weighted", false},
        {"primary_strange_species", "weighted", false},
        {"enu_selected_topology", "weighted", false},
        {"p_lambda_selected_topology", "weighted", false},
        {"enu_selected_topology", "shape", true},
        {"p_lambda_selected_topology", "shape", true},
    };
}

std::vector<hist_spec> exit_hist_specs()
{
    return {
        {"nuclear_exit_cutflow", "weighted", false},
        {"primary_strange_species", "weighted", false},
        {"exit_strange_species", "weighted", false},
        {"primary_origin", "weighted", false},
        {"exit_class", "weighted", false},
        {"enu_exit_lambda", "weighted", false},
        {"p_exit_lambda", "weighted", false},
        {"enu_exit_lambda", "shape", true},
        {"p_exit_lambda", "shape", true},
    };
}

std::vector<TString> primary_map_names()
{
    return {"primary_origin_vs_interaction_mode", "primary_origin_vs_visible_category", "p_lambda_vs_primary_origin"};
}

std::vector<TString> exit_map_names()
{
    return {"primary_origin_vs_exit_class"};
}

TString row_label(const status_row& row, const TString& group_mode)
{
    if (group_mode == "variation") return row.variation;
    if (group_mode == "beam") return row.beam_mode + " " + row.beam_species;
    if (group_mode == "generator") return row.generator + " " + row.variation;
    return row.sample;
}

TString group_key(const status_row& row, const TString& group_mode)
{
    if (group_mode == "variation") {
        return row.generator + "|" + row.beam_mode + "|" + row.beam_species + "|" + row.interaction + "|" + row.fsi_state;
    }
    if (group_mode == "beam") {
        return row.generator + "|" + row.variation + "|" + row.interaction + "|" + row.fsi_state;
    }
    if (group_mode == "generator") {
        return row.beam_mode + "|" + row.beam_species + "|" + row.interaction + "|" + row.fsi_state;
    }
    return row.generator;
}

TString group_title(const TString& key, const TString& group_mode)
{
    TString title = key;
    title.ReplaceAll("|", " ");
    if (group_mode == "variation") return "Variation loop: " + title;
    if (group_mode == "beam") return "Beam loop: " + title;
    if (group_mode == "generator") return "Generator loop: " + title;
    return title;
}

void save_canvas(TCanvas& canvas, const TString& output_base)
{
    canvas.SaveAs(output_base + ".pdf");
    canvas.SaveAs(output_base + ".png");
}

void draw_overlay_group(
    const std::vector<status_row>& rows,
    const TString& input_dir,
    const TString& output_dir,
    const TString& stage,
    const TString& group_mode,
    const TString& key,
    const hist_spec& spec
)
{
    const std::vector<int> colors = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1, kOrange + 7,
                                     kCyan + 2, kViolet + 5, kGray + 2, kBlack, kSpring + 5};
    std::vector<TH1*> hists;
    std::vector<TString> labels;
    double ymax = 0.0;

    for (const status_row& row : rows) {
        if (group_key(row, group_mode) != key) continue;
        if (stage == "primary" && row.primary_status != "ran") continue;
        if (stage == "exit" && row.nuclear_exit_status != "ran") continue;
        const TString path = stage_path(row, input_dir, stage);
        if (gSystem->AccessPathName(path.Data())) continue;

        TH1* hist = clone_hist1(path, spec.name, safe_name(stage + "_" + spec.name + "_" + row.sample));
        if (!hist) continue;
        if (spec.normalize) normalize_shape(hist);
        style_hist(hist, colors[hists.size() % colors.size()], 20 + (hists.size() % 10));
        hist->GetYaxis()->SetTitle(spec.normalize ? "Normalized events" : "Weighted events");
        ymax = std::max(ymax, hist->GetMaximum());
        hists.push_back(hist);
        labels.push_back(row_label(row, group_mode));
    }

    if (hists.empty()) return;

    gSystem->mkdir(output_dir.Data(), true);
    TCanvas canvas("c_overlay", spec.name, 1100, 820);
    canvas.SetLeftMargin(0.14);
    canvas.SetBottomMargin(0.16);
    TLegend legend(0.16, 0.74, 0.90, 0.90);
    legend.SetBorderSize(0);
    legend.SetNColumns(hists.size() > 4 ? 2 : 1);

    for (int i = 0; i < int(hists.size()); ++i) {
        hists[i]->SetTitle(group_title(key, group_mode));
        hists[i]->GetYaxis()->SetRangeUser(0.0, ymax > 0.0 ? 1.35 * ymax : 1.0);
        hists[i]->Draw(i == 0 ? "hist e" : "hist e same");
        legend.AddEntry(hists[i], labels[i], "l");
    }
    legend.Draw();

    TString output_base = join_path(output_dir, safe_name(stage + "_" + group_mode + "_" + key + "_" + spec.name + "_" + spec.suffix));
    save_canvas(canvas, output_base);
    for (TH1* hist : hists) delete hist;
}

void make_overlays(
    const std::vector<status_row>& rows,
    const TString& input_dir,
    const TString& output_dir,
    const TString& stage,
    const TString& group_mode,
    const std::vector<hist_spec>& specs
)
{
    std::set<TString> keys;
    for (const status_row& row : rows) {
        if (stage == "primary" && row.primary_status != "ran") continue;
        if (stage == "exit" && row.nuclear_exit_status != "ran") continue;
        keys.insert(group_key(row, group_mode));
    }
    for (const TString& key : keys) {
        for (const hist_spec& spec : specs) {
            draw_overlay_group(rows, input_dir, join_path(output_dir, stage + "/" + group_mode), stage, group_mode, key, spec);
        }
    }
}

void make_maps(
    const std::vector<status_row>& rows,
    const TString& input_dir,
    const TString& output_dir,
    const TString& stage,
    const std::vector<TString>& names
)
{
    gSystem->mkdir(output_dir.Data(), true);
    for (const status_row& row : rows) {
        if (stage == "primary" && row.primary_status != "ran") continue;
        if (stage == "exit" && row.nuclear_exit_status != "ran") continue;
        const TString path = stage_path(row, input_dir, stage);
        if (gSystem->AccessPathName(path.Data())) continue;

        for (const TString& name : names) {
            TH2* hist = clone_hist2(path, name, safe_name(stage + "_" + name + "_" + row.sample));
            if (!hist) continue;
            TCanvas canvas("c_map", name, 1100, 850);
            canvas.SetRightMargin(0.18);
            canvas.SetLeftMargin(0.16);
            canvas.SetBottomMargin(0.16);
            hist->SetStats(false);
            hist->SetTitle(row.generator + " " + row.variation + " " + row.beam_mode + " " + row.beam_species);
            hist->Draw(name.Contains("_vs_primary_origin") ? "colz" : "colz text");
            save_canvas(canvas, join_path(output_dir, safe_name(stage + "_map_" + row.sample + "_" + name)));
            delete hist;
        }
    }
}

}  // namespace

void plot_generator_matrix(
    TString status_csv = "analysis/output/matrix/generator_matrix_status.csv",
    TString input_dir = "analysis/output/matrix",
    TString output_dir = "analysis/output/matrix/plots"
)
{
    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);

    const std::vector<status_row> rows = read_status(status_csv);
    if (rows.empty()) {
        std::cerr << "No rows to plot from " << status_csv.Data() << std::endl;
        return;
    }

    make_overlays(rows, input_dir, output_dir, "primary", "variation", primary_hist_specs());
    make_overlays(rows, input_dir, output_dir, "primary", "beam", primary_hist_specs());
    make_overlays(rows, input_dir, output_dir, "primary", "generator", primary_hist_specs());
    make_maps(rows, input_dir, join_path(output_dir, "primary/maps"), "primary", primary_map_names());

    make_overlays(rows, input_dir, output_dir, "exit", "variation", exit_hist_specs());
    make_overlays(rows, input_dir, output_dir, "exit", "beam", exit_hist_specs());
    make_overlays(rows, input_dir, output_dir, "exit", "generator", exit_hist_specs());
    make_maps(rows, input_dir, join_path(output_dir, "exit/maps"), "exit", exit_map_names());

    std::cout << "Wrote generator matrix plots under " << output_dir.Data() << std::endl;
}
