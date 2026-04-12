#include <iostream>

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
}

void skim_final_state_hyperon(const char* input_path,
                              const char* output_path,
                              const char* tree_name = "FlatTree_VARS")
{
    TFile input(input_path, "READ");
    if (input.IsZombie()) AnaSkim::fail(TString("cannot open input file: ") + input_path);

    TTree* tree = nullptr;
    input.GetObject(tree_name, tree);
    if (!tree) AnaSkim::fail(TString("missing tree ") + tree_name + " in " + input_path);
    if (!tree->GetBranch("pdg")) AnaSkim::fail(TString("missing final-state pdg branch in ") + input_path);

    gSystem->mkdir(gSystem->DirName(output_path), kTRUE);
    TFile output(output_path, "RECREATE");
    if (output.IsZombie()) AnaSkim::fail(TString("cannot create output file: ") + output_path);

    TTree* skim = tree->CopyTree("Sum$(pdg==3122||pdg==3212||pdg==3322||pdg==3312||pdg==3334)>0");
    if (!skim) AnaSkim::fail(TString("CopyTree failed for ") + input_path);
    skim->Write("", TObject::kOverwrite);
    output.Write();

    std::cout << "skim_final_state_hyperon: " << skim->GetEntries()
              << " / " << tree->GetEntries()
              << " events -> " << output_path << std::endl;
}
