#include <TFile.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct plan_row {
    int line = 0;
    TString generator;
    std::vector<TString> versions;
    std::vector<TString> knobs;
    std::vector<TString> beam_modes;
    std::vector<TString> beam_species;
    std::vector<TString> interactions;
    std::vector<TString> fsi_states;
    TString input_template;
    TString sample_template;
    TString variation_template;
    TString run_primary;
    TString run_nuclear_exit;
    TString notes;
};

struct expanded_row {
    int line = 0;
    TString generator;
    TString version;
    TString knob;
    TString beam_mode;
    TString beam_species;
    TString interaction;
    TString fsi_state;
    TString input_file;
    TString sample_label;
    TString variation;
    TString run_primary;
    TString run_nuclear_exit;
};

std::string trim_copy(std::string text)
{
    const std::string whitespace = " \t\r\n";
    const size_t first = text.find_first_not_of(whitespace);
    if (first == std::string::npos) return "";
    const size_t last = text.find_last_not_of(whitespace);
    return text.substr(first, last - first + 1);
}

TString trim_tstring(TString text)
{
    text = text.Strip(TString::kBoth);
    return text;
}

TString lower(TString text)
{
    text.ToLower();
    return text;
}

bool truthy(TString text)
{
    text = lower(trim_tstring(text));
    return text == "1" || text == "true" || text == "yes" || text == "y" || text == "run";
}

std::vector<std::string> split_tsv(const std::string& line)
{
    std::vector<std::string> fields;
    size_t start = 0;
    while (start <= line.size()) {
        const size_t tab = line.find('\t', start);
        if (tab == std::string::npos) {
            fields.push_back(line.substr(start));
            break;
        }
        fields.push_back(line.substr(start, tab - start));
        start = tab + 1;
    }
    return fields;
}

std::vector<TString> split_list(TString text)
{
    std::vector<TString> out;
    std::string raw = text.Data();
    size_t start = 0;
    while (start <= raw.size()) {
        const size_t comma = raw.find(',', start);
        const std::string item = comma == std::string::npos ? raw.substr(start) : raw.substr(start, comma - start);
        TString value = trim_copy(item).c_str();
        if (value != "") out.push_back(value);
        if (comma == std::string::npos) break;
        start = comma + 1;
    }
    if (out.empty()) out.push_back("");
    return out;
}

TString csv_escape(TString text)
{
    text.ReplaceAll("\"", "\"\"");
    TString out = "\"";
    out += text;
    out += "\"";
    return out;
}

TString root_quote(TString text)
{
    text.ReplaceAll("\\", "\\\\");
    text.ReplaceAll("\"", "\\\"");
    TString out = "\"";
    out += text;
    out += "\"";
    return out;
}

TString beam_polarity_from_mode(TString beam_mode)
{
    beam_mode.ToLower();
    if (beam_mode.Contains("rhc")) return "rhc";
    if (beam_mode.Contains("fhc")) return "fhc";
    if (beam_mode == "combined" || beam_mode == "both" || beam_mode == "all") return "combined";
    return beam_mode;
}

TString replace_tokens(
    TString text,
    const TString& generator,
    const TString& version,
    const TString& knob,
    const TString& beam_mode,
    const TString& beam_species,
    const TString& interaction,
    const TString& fsi_state
)
{
    text.ReplaceAll("{generator}", generator);
    text.ReplaceAll("{version}", version);
    text.ReplaceAll("{knob}", knob);
    text.ReplaceAll("{beam_mode}", beam_mode);
    text.ReplaceAll("{beam_species}", beam_species);
    text.ReplaceAll("{interaction}", interaction);
    text.ReplaceAll("{fsi_state}", fsi_state);
    return text;
}

bool has_any_leaf(TTree* tree, std::initializer_list<const char*> names)
{
    for (const char* name : names) {
        if (tree->GetLeaf(name)) return true;
    }
    return false;
}

bool has_exit_branches(const TString& input_file)
{
    TFile* file = TFile::Open(input_file, "READ");
    if (!file || file->IsZombie()) return false;
    TTree* tree = nullptr;
    file->GetObject("FlatTree_VARS", tree);
    if (!tree) {
        file->Close();
        return false;
    }
    const bool has_n = has_any_leaf(tree, {
        "n_exitp", "nexitp", "n_nuclear_exit", "nuclear_exit_n",
        "n_exit", "n_strange_exit", "n_exit_strange", "n_exit_particles"
    });
    const bool has_pdg = has_any_leaf(tree, {
        "pdg_exit", "pdg_exitp", "pdg_nuclear_exit", "nuclear_exit_pdg",
        "exit_pdg", "pdg_strange_exit", "exit_strange_pdg"
    });
    file->Close();
    return has_n && has_pdg;
}

std::vector<plan_row> read_plan(const TString& plan_file)
{
    std::ifstream input(plan_file.Data());
    std::vector<plan_row> rows;
    if (!input) {
        std::cerr << "Could not open generator loop plan: " << plan_file.Data() << std::endl;
        return rows;
    }

    std::map<std::string, int> header;
    std::string line;
    int line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        const std::string trimmed = trim_copy(line);
        if (trimmed.empty() || trimmed[0] == '#') continue;
        const std::vector<std::string> fields = split_tsv(line);
        if (header.empty()) {
            for (int i = 0; i < int(fields.size()); ++i) header[trim_copy(fields[i])] = i;
            continue;
        }

        auto field = [&](const char* name) -> TString {
            const auto found = header.find(name);
            if (found == header.end() || found->second >= int(fields.size())) return "";
            return trim_copy(fields[found->second]).c_str();
        };

        if (!truthy(field("enabled"))) continue;

        plan_row row;
        row.line = line_number;
        row.generator = field("generator");
        row.versions = split_list(field("versions"));
        row.knobs = split_list(field("knobs"));
        row.beam_modes = split_list(field("beam_modes"));
        row.beam_species = split_list(field("beam_species"));
        row.interactions = split_list(field("interactions"));
        row.fsi_states = split_list(field("fsi_states"));
        row.input_template = field("input_template");
        row.sample_template = field("sample_template");
        row.variation_template = field("variation_template");
        row.run_primary = field("run_primary");
        row.run_nuclear_exit = field("run_nuclear_exit");
        row.notes = field("notes");
        rows.push_back(row);
    }
    return rows;
}

std::vector<expanded_row> expand_plan(const std::vector<plan_row>& plan)
{
    std::vector<expanded_row> out;
    for (const plan_row& row : plan) {
        for (const TString& version : row.versions) {
            for (const TString& knob : row.knobs) {
                for (const TString& beam_mode : row.beam_modes) {
                    for (const TString& beam_species : row.beam_species) {
                        for (const TString& interaction : row.interactions) {
                            for (const TString& fsi_state : row.fsi_states) {
                                expanded_row item;
                                item.line = row.line;
                                item.generator = row.generator;
                                item.version = version;
                                item.knob = knob;
                                item.beam_mode = beam_mode;
                                item.beam_species = beam_species;
                                item.interaction = interaction;
                                item.fsi_state = fsi_state;
                                item.input_file = replace_tokens(row.input_template, row.generator, version, knob, beam_mode, beam_species, interaction, fsi_state);
                                item.sample_label = replace_tokens(row.sample_template, row.generator, version, knob, beam_mode, beam_species, interaction, fsi_state);
                                item.variation = replace_tokens(row.variation_template, row.generator, version, knob, beam_mode, beam_species, interaction, fsi_state);
                                if (item.variation == "") item.variation = knob;
                                item.run_primary = row.run_primary;
                                item.run_nuclear_exit = row.run_nuclear_exit;
                                out.push_back(item);
                            }
                        }
                    }
                }
            }
        }
    }
    return out;
}

void load_macro_once(const TString& path, bool& loaded)
{
    if (loaded) return;
    gROOT->ProcessLine(TString(".L ") + path + "+");
    loaded = true;
}

void run_primary_mechanism(
    const expanded_row& row,
    const TString& output_dir,
    const TString& working_point,
    const TString& flux_file,
    double flux_floor_fraction,
    double proposal_emin_gev,
    double proposal_emax_gev,
    bool& loaded
)
{
    load_macro_once("analysis/primary_mechanism.cxx", loaded);
    const TString beam_polarity = beam_polarity_from_mode(row.beam_mode);
    TString cmd = "primary_mechanism(";
    cmd += root_quote(row.input_file) + ",";
    cmd += root_quote(output_dir) + ",";
    cmd += root_quote(row.sample_label) + ",";
    cmd += root_quote(row.generator) + ",";
    cmd += root_quote(row.variation) + ",";
    cmd += root_quote(working_point) + ",";
    cmd += root_quote(flux_file) + ",";
    cmd += TString::Format("%.17g,", flux_floor_fraction);
    cmd += root_quote(beam_polarity) + ",";
    cmd += TString::Format("%.17g,%.17g,", proposal_emin_gev, proposal_emax_gev);
    cmd += root_quote(row.beam_species) + ")";
    gROOT->ProcessLine(cmd);
}

void run_nuclear_exit(
    const expanded_row& row,
    const TString& output_dir,
    const TString& flux_file,
    double flux_floor_fraction,
    double proposal_emin_gev,
    double proposal_emax_gev,
    bool& loaded
)
{
    load_macro_once("analysis/nuclear_exit.cxx", loaded);
    const TString beam_polarity = beam_polarity_from_mode(row.beam_mode);
    TString cmd = "nuclear_exit(";
    cmd += root_quote(row.input_file) + ",";
    cmd += root_quote(output_dir) + ",";
    cmd += root_quote(row.sample_label) + ",";
    cmd += root_quote(row.generator) + ",";
    cmd += root_quote(row.variation) + ",";
    cmd += root_quote(flux_file) + ",";
    cmd += TString::Format("%.17g,", flux_floor_fraction);
    cmd += root_quote(beam_polarity) + ",";
    cmd += TString::Format("%.17g,%.17g,", proposal_emin_gev, proposal_emax_gev);
    cmd += root_quote(row.beam_species) + ")";
    gROOT->ProcessLine(cmd);
}

}  // namespace

void run_generator_matrix(
    TString plan_file = "analysis/config/generator_loop_plan.tsv",
    TString output_dir = "analysis/output/matrix",
    TString working_point = "nominal",
    TString flux_file = "example/numi/flux/microboone_numi_flux_5mev.root",
    double flux_floor_fraction = 0.0,
    double proposal_emin_gev = 0.0,
    double proposal_emax_gev = 10.0
)
{
    const std::vector<plan_row> plan = read_plan(plan_file);
    const std::vector<expanded_row> rows = expand_plan(plan);
    if (rows.empty()) {
        std::cerr << "No enabled generator matrix rows in " << plan_file.Data() << std::endl;
        return;
    }

    gSystem->mkdir(output_dir.Data(), true);
    TString status_path = output_dir;
    if (!status_path.EndsWith("/")) status_path += "/";
    status_path += "generator_matrix_status.csv";
    std::ofstream status(status_path.Data());
    if (!status) {
        std::cerr << "Could not open status output: " << status_path.Data() << std::endl;
        return;
    }
    status << "plan_line,generator,version,variation,knob,beam_mode,beam_species,interaction,"
           << "fsi_state,sample,input_file,primary_status,nuclear_exit_status\n";

    bool primary_loaded = false;
    bool nuclear_loaded = false;
    int missing = 0;
    int primary_run = 0;
    int exit_run = 0;
    int exit_skipped = 0;

    for (const expanded_row& row : rows) {
        TString primary_status = "not_requested";
        TString exit_status = "not_requested";

        if (gSystem->AccessPathName(row.input_file.Data())) {
            primary_status = "missing_input";
            exit_status = "missing_input";
            ++missing;
        } else {
            if (truthy(row.run_primary)) {
                run_primary_mechanism(row, output_dir, working_point, flux_file, flux_floor_fraction,
                                      proposal_emin_gev, proposal_emax_gev, primary_loaded);
                primary_status = "ran";
                ++primary_run;
            }

            const TString run_exit = lower(trim_tstring(row.run_nuclear_exit));
            const bool exit_auto = run_exit == "" || run_exit == "auto";
            const bool exit_requested = exit_auto || truthy(run_exit);
            if (exit_requested) {
                if (has_exit_branches(row.input_file)) {
                    run_nuclear_exit(row, output_dir, flux_file, flux_floor_fraction,
                                     proposal_emin_gev, proposal_emax_gev, nuclear_loaded);
                    exit_status = "ran";
                    ++exit_run;
                } else {
                    exit_status = exit_auto ? "skipped_no_exit_branches" : "missing_exit_branches";
                    ++exit_skipped;
                }
            }
        }

        status << row.line << ','
               << csv_escape(row.generator).Data() << ','
               << csv_escape(row.version).Data() << ','
               << csv_escape(row.variation).Data() << ','
               << csv_escape(row.knob).Data() << ','
               << csv_escape(row.beam_mode).Data() << ','
               << csv_escape(row.beam_species).Data() << ','
               << csv_escape(row.interaction).Data() << ','
               << csv_escape(row.fsi_state).Data() << ','
               << csv_escape(row.sample_label).Data() << ','
               << csv_escape(row.input_file).Data() << ','
               << csv_escape(primary_status).Data() << ','
               << csv_escape(exit_status).Data() << '\n';
    }

    std::cout << "Generator matrix rows expanded: " << rows.size()
              << "\nMissing inputs: " << missing
              << "\nPrimary-mechanism runs: " << primary_run
              << "\nNuclear-exit runs: " << exit_run
              << "\nNuclear-exit skips: " << exit_skipped
              << "\nStatus: " << status_path.Data() << std::endl;
}
