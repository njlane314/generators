#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${script_dir}/../lib/common.sh"

if [ "$#" -gt 1 ]; then
  cat >&2 <<'EOF'
Usage: analysis/bin/pipeline.sh [CONFIG_ENV]

CONFIG_ENV is an optional shell environment file such as
analysis/share/config/milestone_10k.env. Values in the caller environment still
win over defaults set by the config file. With no config, this runs local
10,000 proposal-event generation for every active local row in sample_matrix.tsv.
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
plan="${plan:-${ANA_SHARE_DIR}/config/sample_matrix.tsv}"
generator_filter="${generator_filter:-NuWro,GiBUU}"
beam_energy_min_gev="${beam_energy_min_gev:-${proposal_min_gev:-0}}"
beam_energy_max_gev="${beam_energy_max_gev:-${proposal_max_gev:-10}}"
analysis_energy_min_gev="${analysis_energy_min_gev:-${beam_energy_min_gev}}"
analysis_energy_max_gev="${analysis_energy_max_gev:-${beam_energy_max_gev}}"
matrix_dir="${matrix_dir:-analysis/output/milestone_10k/matrix}"
yield_dir="${yield_dir:-analysis/output/milestone_10k/expected_event_yields}"

printf 'Pipeline plan: %s\n' "${plan}"
printf '  proposal events per sample: %s\n' "${events}"
printf '  generation filter: %s\n' "${generator_filter:-all enabled rows}"
printf '  generation flux: NuMI inputs from the plan row\n'
printf '  beam energy envelope: %s-%s GeV\n' "${beam_energy_min_gev}" "${beam_energy_max_gev}"
printf '  analysis energy window: %s-%s GeV\n' "${analysis_energy_min_gev}" "${analysis_energy_max_gev}"
printf '  run analysis: %s\n' "${run_analysis}"

if [ "${produce_samples}" = 1 ]; then
  dry_run="${dry_run}" \
  events="${events}" \
  generator_filter="${generator_filter}" \
  beam_energy_min_gev="${beam_energy_min_gev}" \
  beam_energy_max_gev="${beam_energy_max_gev}" \
    "${script_dir}/generate_sample_matrix.sh" "${plan}"
fi

if [ "${dry_run}" = 1 ] || [ "${run_analysis}" != 1 ]; then
  exit 0
fi

"${script_dir}/analyse_generator_matrix.sh" "${plan}" "${matrix_dir}" nominal "${analysis_energy_min_gev}" "${analysis_energy_max_gev}"
"${script_dir}/compute_expected_yields.sh" "${matrix_dir}/generator_matrix_status.csv" "${yield_dir}" 1.0 0.30 0.07 "${analysis_energy_min_gev}" "${analysis_energy_max_gev}"

if [ "${run_plots}" = 1 ]; then
  "${script_dir}/plot_beam_kinematic_envelope.sh" "${matrix_dir}/generator_matrix_status.csv" analysis/output/beam_kinematic_envelope "enu q2 w hyperon_p" final_hyperon
  "${script_dir}/plot_generator_kinematic_envelope.sh" "${matrix_dir}/generator_matrix_status.csv" analysis/output/generator_kinematic_envelope "enu q2 w hyperon_p" final_hyperon
  "${script_dir}/plot_kinematic_coverage.sh" "${matrix_dir}/generator_matrix_status.csv" analysis/output/coverage "variation beam generator" coverage "enu_q2 q2_w" final_hyperon
fi

printf '%s\n' \
  "${matrix_dir}/generator_matrix_status.csv" \
  "${yield_dir}/expected_event_yields_by_variation.csv"
