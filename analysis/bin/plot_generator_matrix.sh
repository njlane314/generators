#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${script_dir}/../lib/common.sh"
status_csv="${1:-analysis/output/matrix/generator_matrix_status.csv}"
input_dir="${2:-analysis/output/matrix}"
output_dir="${3:-analysis/output/matrix/plots}"

ana_cd_repo
root -l -b -q "${ANA_ROOT_DIR}/matrix/plot_generator_matrix.cxx+(\"${status_csv}\",\"${input_dir}\",\"${output_dir}\")"
