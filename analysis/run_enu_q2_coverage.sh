#!/usr/bin/env bash
set -euo pipefail

if [[ $# -gt 10 ]]; then
  cat >&2 <<'EOF'
Usage: analysis/run_enu_q2_coverage.sh [INPUT_SOURCE] [OUTPUT_DIR] [LABELS_OR_GROUPINGS] [OUTPUT_STEM] [SELECTION] [FLUX_FILE] [PROPOSAL_EMIN] [PROPOSAL_EMAX] [PROTON_THRESHOLD] [PIMINUS_THRESHOLD]

This is the E_nu-vs-Q2 shortcut for analysis/run_kinematic_coverage.sh.
EOF
  exit 1
fi

input_source="${1:-analysis/output/matrix/generator_matrix_status.csv}"
output_dir="${2:-analysis/output/coverage/enu_q2}"
labels_or_groupings="${3:-variation beam generator}"
output_stem="${4:-enu_q2_coverage}"
selection="${5:-final_hyperon}"
flux_file="${6-analysis/flux/microboone_numi_flux_5mev.root}"
proposal_emin="${7:-0.0}"
proposal_emax="${8:-10.0}"
proton_threshold="${9:-0.30}"
piminus_threshold="${10:-0.07}"

root -l -b -q "analysis/plot_enu_q2_coverage.cxx+(\"${input_source}\",\"${output_dir}\",\"${labels_or_groupings}\",\"${output_stem}\",\"${selection}\",\"${flux_file}\",${proposal_emin},${proposal_emax},${proton_threshold},${piminus_threshold})"
