#!/usr/bin/env bash
set -euo pipefail

if [[ $# -gt 11 ]]; then
  cat >&2 <<'EOF'
Usage: analysis/run_generator_kinematic_envelope.sh [STATUS_CSV] [OUTPUT_DIR] [VARIABLES] [SELECTION] [FLUX_FILE] [FLUX_FLOOR] [PROPOSAL_EMIN] [PROPOSAL_EMAX] [PROTON_THRESHOLD] [PIMINUS_THRESHOLD] [FIX_FSI]

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
flux_file="${5:-analysis/flux/microboone_numi_flux_5mev.root}"
flux_floor="${6:-0.0}"
proposal_emin="${7:-0.0}"
proposal_emax="${8:-10.0}"
proton_threshold="${9:-0.30}"
piminus_threshold="${10:-0.07}"
fix_fsi="${11:-false}"

case "${fix_fsi}" in
  1|true|TRUE|yes|YES|y|Y) include_fsi_state_in_group=true ;;
  *) include_fsi_state_in_group=false ;;
esac

root -l -b -q "analysis/generator_kinematic_envelope.cxx+(\"${status_csv}\",\"${output_dir}\",\"${variables}\",\"${selection}\",\"${flux_file}\",${flux_floor},${proposal_emin},${proposal_emax},${proton_threshold},${piminus_threshold},${include_fsi_state_in_group})"
