#!/usr/bin/env bash
set -euo pipefail

if [[ $# -gt 7 ]]; then
  cat >&2 <<'EOF'
Usage: analysis/run_expected_event_yields.sh [STATUS_CSV] [OUTPUT_DIR] [EXPOSURE_SCALE] [PROTON_THRESHOLD] [PIMINUS_THRESHOLD] [ENERGY_MIN] [ENERGY_MAX]

STATUS_CSV defaults to analysis/output/matrix/generator_matrix_status.csv.
OUTPUT_DIR defaults to analysis/output/expected_event_yields.
EXPOSURE_SCALE defaults to 1.0; pass an exposure/POT/target scale for absolute event counts.
PROTON_THRESHOLD/PIMINUS_THRESHOLD default to 0.30/0.07 GeV for the Lambda decay daughters.
ENERGY_MIN/ENERGY_MAX default to the 0-10 GeV analysis energy window.
EOF
  exit 1
fi

status_csv="${1:-analysis/output/matrix/generator_matrix_status.csv}"
output_dir="${2:-analysis/output/expected_event_yields}"
exposure_scale="${3:-1.0}"
proton_threshold="${4:-0.30}"
piminus_threshold="${5:-0.07}"
energy_min="${6:-0.0}"
energy_max="${7:-10.0}"

root -l -b -q "analysis/expected_event_yields.cxx+(\"${status_csv}\",\"${output_dir}\",${exposure_scale},${proton_threshold},${piminus_threshold},${energy_min},${energy_max})"
