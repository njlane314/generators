#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TObject.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

namespace AnaSkim {
void fail(const TString& message)
{
    std::cerr << "ERROR: " << message << std::endl;
    gSystem->Exit(1);
}

std::vector<std::string> split_inputs(const char* input_paths)
{
    std::vector<std::string> paths;
    std::stringstream stream(input_paths ? input_paths : "");
    std::string path;
    while (std::getline(stream, path, ',')) {
        if (!path.empty()) paths.push_back(path);
    }
    return paths;
}
}

void skim_final_state_hyperon(const char* input_path,
                              const char* output_path,
                              const char* tree_name = "FlatTree_VARS",
                              Long64_t max_entries = -1,
                              const char* count_path = "")
{
    const std::vector<std::string> inputs = AnaSkim::split_inputs(input_path);
    if (inputs.empty()) AnaSkim::fail("no input file specified");

    TChain chain(tree_name);
    for (const std::string& input : inputs) {
        if (chain.Add(input.c_str()) == 0) {
            AnaSkim::fail(TString("cannot add input file/tree: ") + input.c_str());
        }
    }
    if (!chain.GetListOfBranches() || !chain.GetListOfBranches()->FindObject("pdg")) {
        AnaSkim::fail(TString("missing final-state pdg branch in ") + input_path);
    }

    gSystem->mkdir(gSystem->DirName(output_path), kTRUE);
    TFile output(output_path, "RECREATE");
    if (output.IsZombie()) AnaSkim::fail(TString("cannot create output file: ") + output_path);

    const char* selection = "Sum$(pdg==3122||pdg==3212||pdg==3322||pdg==3312||pdg==3334)>0";
    TTree* skim = max_entries >= 0
        ? chain.CopyTree(selection, "", max_entries, 0)
        : chain.CopyTree(selection);
    if (!skim) AnaSkim::fail(TString("CopyTree failed for ") + input_path);
    skim->Write("", TObject::kOverwrite);
    output.Write();

    if (count_path && TString(count_path) != "") {
        std::ofstream count(count_path);
        count << skim->GetEntries() << "\n";
    }

    std::cout << "skim_final_state_hyperon: " << skim->GetEntries()
              << " / " << chain.GetEntries()
              << " events -> " << output_path << std::endl;
}
