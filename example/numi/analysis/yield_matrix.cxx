#include "topology_yields.cxx"

#include <TString.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

struct matrix_row {
  std::string sample;
  std::string generator;
  std::string knob;
  std::string category;
  double raw_events = 0.0;
  double weighted_yield = 0.0;
  double weighted_stat_uncertainty = 0.0;
  double xsec_weighted = 0.0;
  double xsec_weighted_stat_uncertainty = 0.0;
};

std::string trim_copy(std::string text)
{
  const std::string whitespace = " \t\r\n";
  const size_t first = text.find_first_not_of(whitespace);
  if (first == std::string::npos) return "";
  const size_t last = text.find_last_not_of(whitespace);
  return text.substr(first, last - first + 1);
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

double parse_double(const std::string& text)
{
  return std::strtod(text.c_str(), nullptr);
}

bool parse_yield_row(const std::string& line, matrix_row& row)
{
  const std::vector<std::string> fields = split_csv(line);
  if (fields.size() < 10) return false;

  row.sample = fields[0];
  row.generator = fields[1];
  row.knob = fields[2];
  row.category = fields[3];
  row.raw_events = parse_double(fields[4]);
  row.weighted_yield = parse_double(fields[5]);
  row.weighted_stat_uncertainty = parse_double(fields[6]);
  row.xsec_weighted = parse_double(fields[7]);
  row.xsec_weighted_stat_uncertainty = parse_double(fields[8]);
  return true;
}

bool append_csv(const TString& input_csv, std::ofstream& output, bool write_header)
{
  std::ifstream input(input_csv.Data());
  if (!input) {
    std::cerr << "Could not open temporary yield CSV: " << input_csv.Data() << std::endl;
    return false;
  }

  std::string line;
  bool first_line = true;
  while (std::getline(input, line)) {
    if (first_line) {
      first_line = false;
      if (!write_header) continue;
    }
    output << line << '\n';
  }
  return true;
}

bool collect_yield_rows(const TString& input_csv, std::vector<matrix_row>& rows)
{
  std::ifstream input(input_csv.Data());
  if (!input) {
    std::cerr << "Could not open temporary yield CSV: " << input_csv.Data() << std::endl;
    return false;
  }

  std::string line;
  bool first_line = true;
  while (std::getline(input, line)) {
    if (first_line) {
      first_line = false;
      continue;
    }

    matrix_row row;
    if (!parse_yield_row(line, row)) {
      std::cerr << "Could not parse yield CSV row: " << input_csv.Data() << std::endl;
      return false;
    }
    rows.push_back(row);
  }
  return true;
}

TString uncertainty_output_path(TString output_csv)
{
  if (output_csv.EndsWith(".csv")) {
    output_csv.Remove(output_csv.Length() - 4);
  }
  output_csv += "_uncertainties.csv";
  return output_csv;
}

void write_uncertainty_summary(const TString& output_csv, const std::vector<matrix_row>& rows)
{
  std::map<std::string, std::vector<matrix_row>> grouped_rows;
  for (const matrix_row& row : rows) {
    grouped_rows[row.generator + "\t" + row.category].push_back(row);
  }

  const TString uncertainty_csv = uncertainty_output_path(output_csv);
  std::ofstream output(uncertainty_csv.Data());
  if (!output) {
    std::cerr << "Could not open uncertainty CSV: " << uncertainty_csv.Data() << std::endl;
    return;
  }

  output << "generator,category,central_sample,central_knob,n_knobs,"
         << "weighted_yield,weighted_stat_uncertainty,weighted_systematic_uncertainty,"
         << "weighted_total_uncertainty,xsec_weighted_1e38_cm2_per_Ar,"
         << "xsec_weighted_stat_uncertainty_1e38_cm2_per_Ar,"
         << "xsec_weighted_systematic_uncertainty_1e38_cm2_per_Ar,"
         << "xsec_weighted_total_uncertainty_1e38_cm2_per_Ar,"
         << "minimum_xsec_weighted_1e38_cm2_per_Ar,"
         << "maximum_xsec_weighted_1e38_cm2_per_Ar\n";
  output << std::setprecision(12);

  for (const auto& item : grouped_rows) {
    const std::vector<matrix_row>& group = item.second;
    int central_index = 0;
    for (int i = 0; i < int(group.size()); ++i) {
      if (group[i].knob == "nominal") {
        central_index = i;
        break;
      }
    }

    const matrix_row& central = group[central_index];
    double weighted_systematic = 0.0;
    double xsec_systematic = 0.0;
    double min_xsec = central.xsec_weighted;
    double max_xsec = central.xsec_weighted;

    for (const matrix_row& row : group) {
      weighted_systematic =
          std::max(weighted_systematic, std::abs(row.weighted_yield - central.weighted_yield));
      xsec_systematic =
          std::max(xsec_systematic, std::abs(row.xsec_weighted - central.xsec_weighted));
      min_xsec = std::min(min_xsec, row.xsec_weighted);
      max_xsec = std::max(max_xsec, row.xsec_weighted);
    }

    const double weighted_total =
        std::sqrt(central.weighted_stat_uncertainty * central.weighted_stat_uncertainty +
                  weighted_systematic * weighted_systematic);
    const double xsec_total =
        std::sqrt(central.xsec_weighted_stat_uncertainty *
                      central.xsec_weighted_stat_uncertainty +
                  xsec_systematic * xsec_systematic);

    output << csv_escape(central.generator.c_str()).Data() << ','
           << csv_escape(central.category.c_str()).Data() << ','
           << csv_escape(central.sample.c_str()).Data() << ','
           << csv_escape(central.knob.c_str()).Data() << ','
           << group.size() << ','
           << central.weighted_yield << ','
           << central.weighted_stat_uncertainty << ','
           << weighted_systematic << ','
           << weighted_total << ','
           << central.xsec_weighted << ','
           << central.xsec_weighted_stat_uncertainty << ','
           << xsec_systematic << ','
           << xsec_total << ','
           << min_xsec << ','
           << max_xsec << '\n';
  }

  std::cout << "Wrote " << uncertainty_csv.Data() << std::endl;
}

}  // namespace

void yield_matrix(TString manifest_file, TString output_csv)
{
  std::ifstream manifest(manifest_file.Data());
  if (!manifest) {
    std::cerr << "Could not open manifest: " << manifest_file.Data() << std::endl;
    return;
  }

  TString output_dir = gSystem->DirName(output_csv.Data());
  if (output_dir != "." && output_dir != "") {
    gSystem->mkdir(output_dir.Data(), true);
  }

  std::ofstream output(output_csv.Data());
  if (!output) {
    std::cerr << "Could not open output CSV: " << output_csv.Data() << std::endl;
    return;
  }

  int line_number = 0;
  int sample_count = 0;
  std::vector<matrix_row> rows;
  std::string line;
  while (std::getline(manifest, line)) {
    ++line_number;

    const std::string trimmed_line = trim_copy(line);
    if (trimmed_line.empty() || trimmed_line[0] == '#') continue;

    std::vector<std::string> fields = split_tsv(line);
    if (fields.size() > 4) {
      std::cerr << "Too many columns in " << manifest_file.Data() << ":"
                << line_number << std::endl;
      return;
    }
    while (fields.size() < 4) fields.push_back("");

    const std::string input_text = trim_copy(fields[0]);
    const std::string generator_text = trim_copy(fields[1]);
    const std::string knob_text = trim_copy(fields[2]);
    const std::string sample_text = trim_copy(fields[3]);
    const TString input_file = input_text.c_str();
    const TString generator = generator_text.c_str();
    const TString knob = knob_text.c_str();
    const TString sample_label = sample_text.c_str();

    if (input_file == "") {
      std::cerr << "Missing input flat tree in " << manifest_file.Data() << ":"
                << line_number << std::endl;
      return;
    }
    if (gSystem->AccessPathName(input_file.Data())) {
      std::cerr << "Missing input flat tree in " << manifest_file.Data() << ":"
                << line_number << ": " << input_file.Data() << std::endl;
      return;
    }

    TString temp_csv = output_csv;
    temp_csv += TString::Format(".sample_%d.tmp.csv", line_number);
    topology_yields(input_file, temp_csv, sample_label, generator, knob);

    if (!collect_yield_rows(temp_csv, rows)) {
      gSystem->Unlink(temp_csv.Data());
      return;
    }
    if (!append_csv(temp_csv, output, sample_count == 0)) {
      gSystem->Unlink(temp_csv.Data());
      return;
    }
    gSystem->Unlink(temp_csv.Data());
    ++sample_count;
  }

  if (sample_count == 0) {
    std::cerr << "Manifest had no samples: " << manifest_file.Data() << std::endl;
    return;
  }

  std::cout << "Wrote " << output_csv.Data() << " for "
            << sample_count << " samples" << std::endl;
  write_uncertainty_summary(output_csv, rows);
}
