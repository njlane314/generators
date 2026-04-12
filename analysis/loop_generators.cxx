#include <TROOT.h>

{
    gROOT->ProcessLine(".L analysis/run_generator_matrix.cxx+");
    gROOT->ProcessLine("run_generator_matrix(\"ana/config/generator_loop_plan.tsv\", \"analysis/output/matrix\")");
    gROOT->ProcessLine(".q");
}
