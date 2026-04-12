#!/usr/bin/env bash
set -euo pipefail

if [[ $# -gt 6 ]]; then
  cat >&2 <<'EOF'
Usage: analysis/run_beam_kinematic_envelope.sh [STATUS_CSV] [OUTPUT_DIR] [VARIABLES] [SELECTION] [PROTON_THRESHOLD] [PIMINUS_THRESHOLD]

STATUS_CSV defaults to analysis/output/matrix/generator_matrix_status.csv.
OUTPUT_DIR defaults to analysis/output/beam_kinematic_envelope.
VARIABLES defaults to "enu q2 w hyperon_p"; use "all" for the same default.
SELECTION defaults to detector_visible_lambda; alternatives include final_hyperon, detector_branching, and all.
EOF
  exit 1
fi

status_csv="${1:-analysis/output/matrix/generator_matrix_status.csv}"
output_dir="${2:-analysis/output/beam_kinematic_envelope}"
variables="${3:-enu q2 w hyperon_p}"
selection="${4:-detector_visible_lambda}"
proton_threshold="${5:-0.30}"
piminus_threshold="${6:-0.07}"

root -l -b -q "analysis/beam_kinematic_envelope.cxx+(\"${status_csv}\",\"${output_dir}\",\"${variables}\",\"${selection}\",${proton_threshold},${piminus_threshold})"
