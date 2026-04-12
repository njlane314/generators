#include "topology_yields.cxx"

#include <TString.h>
#include <TSystem.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

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
}
