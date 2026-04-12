#!/usr/bin/env bash
set -euo pipefail

if [[ $# -gt 11 ]]; then
  cat >&2 <<'EOF'
Usage: analysis/run_kinematic_coverage.sh [INPUT_SOURCE] [OUTPUT_DIR] [LABELS_OR_GROUPINGS] [OUTPUT_STEM] [PAIRS] [SELECTION] [FLUX_FILE] [PROPOSAL_EMIN] [PROPOSAL_EMAX] [PROTON_THRESHOLD] [PIMINUS_THRESHOLD]

INPUT_SOURCE defaults to analysis/output/matrix/generator_matrix_status.csv.
For a status CSV, LABELS_OR_GROUPINGS is a space/comma list from: variation beam generator all.
For explicit ROOT files, INPUT_SOURCE is a comma-separated file list and LABELS_OR_GROUPINGS is the matching label list.
SELECTION defaults to final_hyperon, using final-state PDG in {3122,3212,3322,3312,3334} with no momentum cut where FlatTree final-state PDGs are available. Use detector_branching or detector_visible_lambda to apply the detector-visibility envelope.
EOF
  exit 1
fi

input_source="${1:-analysis/output/matrix/generator_matrix_status.csv}"
output_dir="${2:-analysis/output/coverage}"
labels_or_groupings="${3:-variation beam generator}"
output_stem="${4:-coverage}"
pairs="${5:-enu_q2 enu_w q2_w enu_lambda_p w_lambda_p lambda_p_costheta}"
selection="${6:-final_hyperon}"
flux_file="${7:-analysis/flux/microboone_numi_flux_5mev.root}"
proposal_emin="${8:-0.0}"
proposal_emax="${9:-10.0}"
proton_threshold="${10:-0.30}"
piminus_threshold="${11:-0.07}"

root -l -b -q "analysis/plot_kinematic_coverage.cxx+(\"${input_source}\",\"${output_dir}\",\"${labels_or_groupings}\",\"${output_stem}\",\"${pairs}\",\"${selection}\",\"${flux_file}\",${proposal_emin},${proposal_emax},${proton_threshold},${piminus_threshold})"
