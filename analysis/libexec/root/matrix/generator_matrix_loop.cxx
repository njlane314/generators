#include <TROOT.h>

{
    gROOT->ProcessLine(".L analysis/libexec/root/matrix/generator_matrix.cxx+");
    gROOT->ProcessLine("run_generator_matrix(\"analysis/share/config/sample_matrix.tsv\", \"analysis/output/matrix\")");
    gROOT->ProcessLine(".q");
}
