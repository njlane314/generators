#!/usr/bin/env bash
set -euo pipefail

if [[ $# -gt 6 ]]; then
  cat >&2 <<'EOF'
Usage: analysis/run_detector_visibility_coverage.sh [INPUT_SOURCE] [OUTPUT_DIR] [GROUPINGS] [OUTPUT_STEM] [PAIRS] [ASSUMPTIONS]

ASSUMPTIONS is a comma list of label:selection:p_proton:p_piminus rows.
Use "default" for truth, branching, loose, nominal, and tight.
EOF
  exit 1
fi

input="${1:-analysis/output/matrix/generator_matrix_status.csv}"
out="${2:-analysis/output/coverage/detector_visibility}"
groups="${3:-variation beam generator}"
stem="${4:-coverage}"
pairs="${5:-enu_q2 enu_w q2_w enu_lambda_p w_lambda_p lambda_p_costheta}"
assumptions="${6:-default}"

mkdir -p "${out}"
manifest="${out}/detector_visibility_coverage_assumptions.tsv"
printf 'label\tselection\tproton_threshold_gev\tpiminus_threshold_gev\toutput_dir\n' > "${manifest}"

safe_name() {
  printf '%s' "$1" | tr ' /:;|(),#{}^=+' '_' | sed 's/-/minus/g'
}

run_one() {
  local label="$1" selection="$2" proton="$3" piminus="$4"
  local safe_label
  safe_label="$(safe_name "${label}")"
  local dir="${out}/${safe_label}"
  printf '%s\t%s\t%s\t%s\t%s\n' "${label}" "${selection}" "${proton}" "${piminus}" "${dir}" >> "${manifest}"
  root -l -b -q "analysis/plot_kinematic_coverage.cxx+(\"${input}\",\"${dir}\",\"${groups}\",\"${stem}_${safe_label}\",\"${pairs}\",\"${selection}\",${proton},${piminus})"
}

case "${assumptions}" in
  ""|default|all)
    run_one truth final_hyperon 0.30 0.07
    run_one branching detector_branching 0.30 0.07
    run_one loose detector_visible_lambda 0.20 0.05
    run_one nominal detector_visible_lambda 0.30 0.07
    run_one tight detector_visible_lambda 0.45 0.12
    ;;
  *)
    IFS=',' read -r -a rows <<< "${assumptions}"
    for row in "${rows[@]}"; do
      IFS=':' read -r label selection proton piminus extra <<< "${row}"
      if [[ -n "${extra:-}" || -z "${label:-}" || -z "${selection:-}" || -z "${proton:-}" ]]; then
        printf 'bad assumption: %s\n' "${row}" >&2
        exit 2
      fi
      if [[ -z "${piminus:-}" ]]; then
        piminus="${proton}"
        proton="${selection}"
        selection="detector_visible_lambda"
      fi
      run_one "${label}" "${selection}" "${proton}" "${piminus}"
    done
    ;;
esac
