#include <TROOT.h>

void plot_generator_matrix_overlay()
{
    gROOT->ProcessLine(".L analysis/libexec/root/matrix/plot_generator_matrix.cxx+");
    gROOT->ProcessLine("plot_generator_matrix()");
}
