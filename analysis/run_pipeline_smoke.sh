#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
events="${events:-100}"
produce_samples="${produce_samples:-1}"
run_plots="${run_plots:-0}"
dry_run="${dry_run:-0}"
plan="${plan:-${script_dir}/config/generator_smoke_plan.tsv}"
matrix_dir="${matrix_dir:-analysis/output/matrix}"
yield_dir="${yield_dir:-analysis/output/expected_event_yields}"
sensitivity_dir="${sensitivity_dir:-analysis/output/phenomenology_sensitivity}"
flux_file="${flux_file:-analysis/flux/microboone_numi_flux_5mev.root}"
generation_flux_mode="${generation_flux_mode:-numi}"
if [ -z "${analysis_flux_file+x}" ]; then
  case "${generation_flux_mode}" in
    numi) analysis_flux_file="" ;;
    *) analysis_flux_file="${flux_file}" ;;
  esac
fi

if [ "${produce_samples}" = 1 ]; then
  dry_run="${dry_run}" events="${events}" generation_flux_mode="${generation_flux_mode}" "${script_dir}/run_sample_matrix.sh" "${plan}"
fi

if [ "${dry_run}" = 1 ]; then
  exit 0
fi

"${script_dir}/run_generator_matrix.sh" "${plan}" "${matrix_dir}" nominal "${analysis_flux_file}"
"${script_dir}/run_expected_event_yields.sh" "${matrix_dir}/generator_matrix_status.csv" "${yield_dir}" "${analysis_flux_file}"
"${script_dir}/run_phenomenology_sensitivity.sh" "${yield_dir}/expected_event_yields_by_variation.csv" "${sensitivity_dir}"

if [ "${run_plots}" = 1 ]; then
  "${script_dir}/run_beam_kinematic_envelope.sh" "${matrix_dir}/generator_matrix_status.csv" analysis/output/beam_kinematic_envelope "enu q2 w hyperon_p" final_hyperon "${analysis_flux_file}"
  "${script_dir}/run_generator_kinematic_envelope.sh" "${matrix_dir}/generator_matrix_status.csv" analysis/output/generator_kinematic_envelope "enu q2 w hyperon_p" final_hyperon "${analysis_flux_file}"
  "${script_dir}/run_kinematic_coverage.sh" "${matrix_dir}/generator_matrix_status.csv" analysis/output/coverage "variation beam generator" coverage "enu_q2 q2_w" final_hyperon "${analysis_flux_file}"
fi

printf '%s\n' \
  "${matrix_dir}/generator_matrix_status.csv" \
  "${yield_dir}/expected_event_yields_by_variation.csv" \
  "${sensitivity_dir}"
