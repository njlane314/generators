#!/usr/bin/env bash
set -euo pipefail

status_csv="${1:-analysis/output/matrix/generator_matrix_status.csv}"
input_dir="${2:-analysis/output/matrix}"
output_dir="${3:-analysis/output/matrix/plots}"

root -l -b -q "analysis/plot_generator_matrix.cxx+(\"${status_csv}\",\"${input_dir}\",\"${output_dir}\")"
