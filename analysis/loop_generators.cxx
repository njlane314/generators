#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

#include <iostream>
#include <vector>

{
    std::vector<TString> sample_paths;
    std::vector<TString> sample_labels;

    gROOT->ProcessLine(".L analysis/analyser.cxx+");

    for (int i = 0; i < (int)sample_paths.size(); ++i) {
        if (gSystem->AccessPathName(sample_paths[i])) {
            std::cout << "Skipping missing sample: " << sample_paths[i] << std::endl;
            continue;
        }
        gROOT->ProcessLine("analyser(\"" + sample_paths[i] + "\", \"" + sample_labels[i] + "\").loop()");
    }

    gROOT->ProcessLine(".q");
}
