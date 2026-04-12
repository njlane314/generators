#include <TROOT.h>

void generator_overlay()
{
    gROOT->ProcessLine(".L analysis/plot_generator_matrix.cxx+");
    gROOT->ProcessLine("plot_generator_matrix()");
}
