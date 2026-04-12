#!/usr/bin/env bash
set -euo pipefail

if [[ $# -gt 9 ]]; then
  cat >&2 <<'EOF'
Usage: analysis/run_expected_event_yields.sh [STATUS_CSV] [OUTPUT_DIR] [FLUX_FILE] [FLUX_FLOOR] [PROPOSAL_EMIN] [PROPOSAL_EMAX] [EXPOSURE_SCALE] [PROTON_THRESHOLD] [PIMINUS_THRESHOLD]

STATUS_CSV defaults to analysis/output/matrix/generator_matrix_status.csv.
OUTPUT_DIR defaults to analysis/output/expected_event_yields.
FLUX_FILE defaults to analysis/flux/microboone_numi_flux_5mev.root; pass "" to disable downstream flux reweighting for already flux-shaped samples.
PROPOSAL_EMIN/PROPOSAL_EMAX default to the 0-10 GeV analysis energy window.
EXPOSURE_SCALE defaults to 1.0; pass an exposure/POT/target scale for absolute event counts.
PROTON_THRESHOLD/PIMINUS_THRESHOLD default to 0.30/0.07 GeV for the Lambda decay daughters.
EOF
  exit 1
fi

status_csv="${1:-analysis/output/matrix/generator_matrix_status.csv}"
output_dir="${2:-analysis/output/expected_event_yields}"
flux_file="${3-analysis/flux/microboone_numi_flux_5mev.root}"
flux_floor="${4:-0.0}"
proposal_emin="${5:-0.0}"
proposal_emax="${6:-10.0}"
exposure_scale="${7:-1.0}"
proton_threshold="${8:-0.30}"
piminus_threshold="${9:-0.07}"

root -l -b -q "analysis/expected_event_yields.cxx+(\"${status_csv}\",\"${output_dir}\",\"${flux_file}\",${flux_floor},${proposal_emin},${proposal_emax},${exposure_scale},${proton_threshold},${piminus_threshold})"
