#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${script_dir}/../lib/common.sh"
if [[ $# -gt 7 ]]; then
  cat >&2 <<'EOF'
Usage: analysis/bin/plot_enu_q2_coverage.sh [INPUT_SOURCE] [OUTPUT_DIR] [LABELS_OR_GROUPINGS] [OUTPUT_STEM] [SELECTION] [PROTON_THRESHOLD] [PIMINUS_THRESHOLD]

This is the E_nu-vs-Q2 shortcut for analysis/bin/plot_kinematic_coverage.sh.
EOF
  exit 1
fi

input_source="${1:-analysis/output/matrix/generator_matrix_status.csv}"
output_dir="${2:-analysis/output/coverage/enu_q2}"
labels_or_groupings="${3:-variation beam generator}"
output_stem="${4:-enu_q2_coverage}"
selection="${5:-final_hyperon}"
proton_threshold="${6:-0.30}"
piminus_threshold="${7:-0.07}"

ana_cd_repo
root -l -b -q "${ANA_ROOT_DIR}/coverage/plot_enu_q2_coverage.cxx+(\"${input_source}\",\"${output_dir}\",\"${labels_or_groupings}\",\"${output_stem}\",\"${selection}\",${proton_threshold},${piminus_threshold})"
