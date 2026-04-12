#!/usr/bin/env bash
set -euo pipefail

if [[ $# -gt 7 ]]; then
  cat >&2 <<'EOF'
Usage: analysis/run_generator_kinematic_envelope.sh [STATUS_CSV] [OUTPUT_DIR] [VARIABLES] [SELECTION] [PROTON_THRESHOLD] [PIMINUS_THRESHOLD] [FIX_FSI]

STATUS_CSV defaults to analysis/output/matrix/generator_matrix_status.csv.
OUTPUT_DIR defaults to analysis/output/generator_kinematic_envelope.
VARIABLES defaults to "enu q2 w hyperon_p"; use "all" for the same default.
SELECTION defaults to detector_visible_lambda; alternatives include final_hyperon, detector_branching, and all.
FIX_FSI defaults to false, so FSI state is included in the generator/config envelope. Set true to compare generator/config variations at fixed FSI state.
EOF
  exit 1
fi

status_csv="${1:-analysis/output/matrix/generator_matrix_status.csv}"
output_dir="${2:-analysis/output/generator_kinematic_envelope}"
variables="${3:-enu q2 w hyperon_p}"
selection="${4:-detector_visible_lambda}"
proton_threshold="${5:-0.30}"
piminus_threshold="${6:-0.07}"
fix_fsi="${7:-false}"

case "${fix_fsi}" in
  1|true|TRUE|yes|YES|y|Y) include_fsi_state_in_group=true ;;
  *) include_fsi_state_in_group=false ;;
esac

root -l -b -q "analysis/generator_kinematic_envelope.cxx+(\"${status_csv}\",\"${output_dir}\",\"${variables}\",\"${selection}\",${proton_threshold},${piminus_threshold},${include_fsi_state_in_group})"
