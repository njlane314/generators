#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${script_dir}/../lib/common.sh"
repo_root="${ANA_REPO_ROOT}"
plan="${1:-${ANA_SHARE_DIR}/config/sample_matrix.tsv}"

jobs_per_sample="${jobs_per_sample:-10}"
events_per_job="${events_per_job:-200000}"
beam_energy_min_gev="${beam_energy_min_gev:-${proposal_min_gev:-0}}"
beam_energy_max_gev="${beam_energy_max_gev:-${proposal_max_gev:-10}}"
target="${target:-1000180400}"
grid_checkout="${grid_checkout:-R-3_06_02}"
grid_output_base="${grid_output_base:-/pnfs/uboone/persistent/users/$(whoami)/genie_grid}"
grid_resources_dir="${grid_resources_dir:-/pnfs/uboone/persistent/users/$(whoami)/grid}"
args_dir="${args_dir:-${repo_root}/analysis/output/grid/genie/args}"
manifest="${manifest:-${repo_root}/analysis/output/grid/genie/submissions.tsv}"
make_spline="${make_spline:-0}"
dry_run="${dry_run:-0}"
production_role="${production_role:-Analysis}"
submit_script="${submit_script:-${repo_root}/production/submit_genie.sh}"

case "${jobs_per_sample}" in ''|*[!0-9]*) ana_die "jobs_per_sample must be a positive integer" ;; esac
case "${events_per_job}" in ''|*[!0-9]*) ana_die "events_per_job must be a positive integer" ;; esac
[ "${jobs_per_sample}" -gt 0 ] || ana_die "jobs_per_sample must be > 0"
[ "${events_per_job}" -gt 0 ] || ana_die "events_per_job must be > 0"
[ -r "${plan}" ] || ana_die "missing plan: ${plan}"
[ -x "${submit_script}" ] || ana_die "missing GENIE submit script: ${submit_script}"

if [ "${dry_run}" != 1 ]; then
  ana_check_cmds jobsub_submit
  mkdir -p "${grid_resources_dir}"
  cp "${repo_root}/production/genie_grid.sh" "${grid_resources_dir}/genie_grid.sh"
fi
mkdir -p "${args_dir}" "$(dirname "${manifest}")"
printf 'sample\toutput_flat\tgrid_output_dir\tbeam_mode\tbeam_species\tinteraction\ttune\tspline\tflux_file\tflux_hist\targs_file\n' > "${manifest}"

split_csv() {
  printf '%s\n' "$1" | tr ',' '\n'
}

expand_template() {
  local text="$1"
  text="${text//\{version\}/${version}}"
  text="${text//\{knob\}/${knob}}"
  text="${text//\{beam_mode\}/${beam_mode}}"
  text="${text//\{generation_beam_mode\}/${generation_beam_mode}}"
  text="${text//\{beam_species\}/${species}}"
  text="${text//\{interaction\}/${interaction}}"
  text="${text//\{fsi_state\}/${fsi_state}}"
  printf '%s\n' "${text}"
}

safe_name() {
  printf '%s' "$1" | tr ' /:;|(),#{}^=+' '_' | sed 's/-/minus/g'
}

ensure_spline() {
  if [ -s "$1" ]; then
    return 0
  fi
  if [ "${dry_run}" = 1 ]; then
    printf '  dry run: spline is not present yet: %s\n' "$1"
    return 0
  fi
  case "${make_spline}" in
    1|true|TRUE|yes|YES|auto)
      ana_check_cmds gmkspl
      mkdir -p "$(dirname "$1")"
      gmkspl -p "$2" -t "${target}" -e "${beam_energy_max_gev}" \
        -o "$1" --tune "$3" --event-generator-list "$4"
      ;;
    *) ;;
  esac
  [ -s "$1" ] || ana_die "missing GENIE spline: $1"
}

submit_one() {
  local output_path sample_label safe_sample args_file output_dir probe version_tag spline flux_file flux_hist
  output_path="${repo_root}/$(expand_template "${input_template}")"
  sample_label="$(expand_template "${sample_template}")"
  safe_sample="$(safe_name "${sample_label}")"
  args_file="${args_dir}/${safe_sample}.args.txt"
  output_dir="${grid_output_base}/${safe_sample}"
  probe="$(ana_probe_pdg "${species}")"
  version_tag="$(ana_genie_version_tag "${version}")"
  spline="${ANA_SHARE_DIR}/splines/${probe}_${target}_${interaction}_${version_tag}_${knob}.xml"
  flux_file="$(ana_flux_root)"
  flux_hist="$(ana_flux_hist "${beam_mode}" "${species}")"

  ensure_spline "${spline}" "${probe}" "${knob}" "${interaction}"
  printf -- '-n %s -e %s,%s -p %s -t %s --tune %s --event-generator-list %s\n' \
    "${events_per_job}" "${beam_energy_min_gev}" "${beam_energy_max_gev}" \
    "${probe}" "${target}" "${knob}" "${interaction}" > "${args_file}"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${sample_label}" "${output_path}" "${output_dir}" "${beam_mode}" "${species}" \
    "${interaction}" "${knob}" "${spline}" "${flux_file}" "${flux_hist}" "${args_file}" >> "${manifest}"

  printf 'GENIE grid sample: %s\n' "${sample_label}"
  printf '  jobs/events: %s x %s\n' "${jobs_per_sample}" "${events_per_job}"
  printf '  beam energy envelope: %s-%s GeV\n' "${beam_energy_min_gev}" "${beam_energy_max_gev}"
  printf '  grid output: %s\n' "${output_dir}"
  printf '  local flat target: %s\n' "${output_path}"

  if [ "${dry_run}" = 1 ]; then
    printf '  dry run: %s %s %s %s %s %s %s %s %s\n' \
      "${submit_script}" "${jobs_per_sample}" "${safe_sample}" "${args_file}" \
      "${spline}" "${flux_file}" "${flux_hist}" "${output_dir}" "${grid_checkout}"
    return 0
  fi

  submit_args=()
  case "${production_role}" in
    Production|production|prod|PROD) submit_args+=(--prod) ;;
    Analysis|analysis) ;;
    *) ana_die "unsupported production_role: ${production_role}" ;;
  esac

  GRID_RESOURCES_DIR="${grid_resources_dir}" "${submit_script}" "${submit_args[@]}" \
    "${jobs_per_sample}" "${safe_sample}" "${args_file}" "${spline}" "${flux_file}" \
    "${flux_hist}" "${output_dir}" "${grid_checkout}"
}

count=0
while IFS=$'\t' read -r enabled generator versions knobs beam_modes beam_species interactions fsi_states input_template sample_template variation_template run_primary run_nuclear_exit notes; do
  case "${enabled}" in
    true|TRUE|1|yes|YES) ;;
    *) continue ;;
  esac
  [ "${generator}" = "GENIE" ] || continue

  for version in $(split_csv "${versions}"); do
    for knob in $(split_csv "${knobs}"); do
      for beam_mode in $(split_csv "${beam_modes}"); do
        generation_beam_mode="${beam_mode}"
        for species in $(split_csv "${beam_species}"); do
          for interaction in $(split_csv "${interactions}"); do
            for fsi_state in $(split_csv "${fsi_states}"); do
              count=$((count + 1))
              submit_one
            done
          done
        done
      done
    done
  done
done < <(tail -n +2 "${plan}")

printf 'GENIE grid samples considered: %d\n' "${count}"
printf 'submission manifest: %s\n' "${manifest}"
