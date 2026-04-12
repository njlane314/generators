#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ "$#" -gt 1 ]; then
  cat >&2 <<'EOF'
Usage: analysis/run_pipeline.sh [CONFIG_ENV]

CONFIG_ENV is an optional shell environment file such as
analysis/config/milestone_10k.env. Values in the caller environment still win
over defaults set by the config file. With no config, this runs local 10,000
proposal-event generation for every active row in generator_loop_plan.tsv.
EOF
  exit 1
fi

config="${1:-}"
if [ -n "${config}" ]; then
  [ -r "${config}" ] || { printf 'ERROR: missing pipeline config: %s\n' "${config}" >&2; exit 1; }
  # shellcheck disable=SC1090
  . "${config}"
fi

events="${events:-10000}"
produce_samples="${produce_samples:-1}"
run_analysis="${run_analysis:-0}"
run_plots="${run_plots:-0}"
dry_run="${dry_run:-0}"
plan="${plan:-${script_dir}/config/generator_loop_plan.tsv}"
matrix_dir="${matrix_dir:-analysis/output/milestone_10k/matrix}"
yield_dir="${yield_dir:-analysis/output/milestone_10k/expected_event_yields}"
sensitivity_dir="${sensitivity_dir:-analysis/output/milestone_10k/phenomenology_sensitivity}"

printf 'Pipeline plan: %s\n' "${plan}"
printf '  proposal events per sample: %s\n' "${events}"
printf '  generation flux: NuMI inputs from the plan row\n'
printf '  run analysis: %s\n' "${run_analysis}"

if [ "${produce_samples}" = 1 ]; then
  dry_run="${dry_run}" events="${events}" "${script_dir}/run_sample_matrix.sh" "${plan}"
fi

if [ "${dry_run}" = 1 ] || [ "${run_analysis}" != 1 ]; then
  exit 0
fi

"${script_dir}/run_generator_matrix.sh" "${plan}" "${matrix_dir}" nominal
"${script_dir}/run_expected_event_yields.sh" "${matrix_dir}/generator_matrix_status.csv" "${yield_dir}"
"${script_dir}/run_phenomenology_sensitivity.sh" "${yield_dir}/expected_event_yields_by_variation.csv" "${sensitivity_dir}"

if [ "${run_plots}" = 1 ]; then
  "${script_dir}/run_beam_kinematic_envelope.sh" "${matrix_dir}/generator_matrix_status.csv" analysis/output/beam_kinematic_envelope "enu q2 w hyperon_p" final_hyperon
  "${script_dir}/run_generator_kinematic_envelope.sh" "${matrix_dir}/generator_matrix_status.csv" analysis/output/generator_kinematic_envelope "enu q2 w hyperon_p" final_hyperon
  "${script_dir}/run_kinematic_coverage.sh" "${matrix_dir}/generator_matrix_status.csv" analysis/output/coverage "variation beam generator" coverage "enu_q2 q2_w" final_hyperon
fi

printf '%s\n' \
  "${matrix_dir}/generator_matrix_status.csv" \
  "${yield_dir}/expected_event_yields_by_variation.csv" \
  "${sensitivity_dir}"
