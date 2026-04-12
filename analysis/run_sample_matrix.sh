#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
. "${script_dir}/sample_common.sh"
plan="${1:-${script_dir}/config/generator_loop_plan.tsv}"
events="${events:-100}"
dry_run="${dry_run:-0}"
skip_existing="${skip_existing:-0}"
generator_filter="${generator_filter:-}"
max_samples="${max_samples:-0}"
seed_base="${seed_base:-1000}"
gibuu_events_per_ensemble="${gibuu_events_per_ensemble:-5}"
proposal_events="${events}"

case "${proposal_events}" in
  ''|*[!0-9]*) printf 'ERROR: events must be a positive integer proposal-event count\n' >&2; exit 1 ;;
esac
[ "${proposal_events}" -gt 0 ] || { printf 'ERROR: events must be > 0\n' >&2; exit 1; }
case "${gibuu_events_per_ensemble}" in
  ''|*[!0-9]*) printf 'ERROR: gibuu_events_per_ensemble must be a positive integer\n' >&2; exit 1 ;;
esac
[ "${gibuu_events_per_ensemble}" -gt 0 ] || { printf 'ERROR: gibuu_events_per_ensemble must be > 0\n' >&2; exit 1; }

want_generator() {
  [ -z "${generator_filter}" ] && return 0
  case ",${generator_filter}," in
    *,"$1",*) return 0 ;;
    *) return 1 ;;
  esac
}

split_csv() {
  local value="$1"
  printf '%s\n' "${value}" | tr ',' '\n'
}

expand_template() {
  local text="$1"
  text="${text//\{version\}/${version}}"
  text="${text//\{knob\}/${knob}}"
  text="${text//\{beam_mode\}/${beam_mode}}"
  text="${text//\{beam_species\}/${species}}"
  text="${text//\{interaction\}/${interaction}}"
  text="${text//\{fsi_state\}/${fsi_state}}"
  printf '%s\n' "${text}"
}

run_generator_sample() {
  local n="$1"
  local sample_label="$2"
  local output_path="$3"
  local selected_count_path="$4"
  local outdir workdir

  outdir="$(dirname "${output_path}")"
  workdir="${repo_root}/analysis/output/work/${generator}/${sample_label}"

  case "${generator}" in
    GENIE)
      events="${proposal_events}" \
      version="${version}" \
      tune="${knob}" \
      beam_mode="${beam_mode}" \
      beam_species="${species}" \
      interaction="${interaction}" \
      run="$((seed_base + n))" \
      sample="${sample_label}" \
      outdir="${outdir}" \
      workdir="${workdir}" \
      skim_final_state=1 \
      skim_count_path="${selected_count_path}" \
      "${script_dir}/run_genie.sh"
      ;;
    NuWro)
      events="${proposal_events}" \
      seed="$((seed_base + n))" \
      version="${version}" \
      knob="${knob}" \
      beam_mode="${beam_mode}" \
      beam_species="${species}" \
      interaction="${interaction}" \
      fsi_state="${fsi_state}" \
      sample="${sample_label}" \
      outdir="${outdir}" \
      workdir="${workdir}" \
      skim_final_state=1 \
      skim_count_path="${selected_count_path}" \
      "${script_dir}/run_nuwro.sh"
      ;;
    GiBUU)
      events="${proposal_events}" \
      gibuu_events_per_ensemble="${gibuu_events_per_ensemble}" \
      num_runs_same_energy="${num_runs_same_energy:-1}" \
      version="${version}" \
      knob="${knob}" \
      beam_mode="${beam_mode}" \
      beam_species="${species}" \
      interaction="${interaction}" \
      fsi_state="${fsi_state}" \
      sample="${sample_label}" \
      outdir="${outdir}" \
      workdir="${workdir}" \
      skim_final_state=1 \
      skim_count_path="${selected_count_path}" \
      "${script_dir}/run_gibuu.sh"
      ;;
    *)
      printf 'ERROR: unsupported generator: %s\n' "${generator}" >&2
      exit 1
      ;;
  esac
}

run_sample() {
  local n="$1"
  local output_path sample_label selected_count_path selected note
  output_path="${repo_root}/$(expand_template "${input_template}")"
  sample_label="$(expand_template "${sample_template}")"
  selected_count_path="${repo_root}/analysis/output/work/${generator}/${sample_label}/${sample_label}.selected.count"

  if [ "${skip_existing}" = 1 ] && [ -s "${output_path}" ]; then
    printf 'skip existing [%04d] %s\n' "${n}" "${output_path}"
    return 0
  fi

  note=""
  if [ "${generator}" = "GiBUU" ]; then
    if [ -n "${num_ensembles:-}" ]; then
      note=" num_ensembles=${num_ensembles}"
    else
      note=" corrected_num_ensembles=$(( (proposal_events + gibuu_events_per_ensemble - 1) / gibuu_events_per_ensemble ))"
    fi
  fi
  printf 'sample [%04d] %s version=%s knob=%s beam=%s species=%s interaction=%s fsi=%s proposal_events=%s%s\n' \
    "${n}" "${generator}" "${version}" "${knob}" "${beam_mode}" "${species}" "${interaction}" "${fsi_state}" "${proposal_events}" "${note}"

  if [ "${dry_run}" = 1 ]; then
    return 0
  fi

  rm -rf "${repo_root}/analysis/output/work/${generator}/${sample_label}"
  run_generator_sample "${n}" "${sample_label}" "${output_path}" "${selected_count_path}"

  [ -s "${selected_count_path}" ] || { printf 'ERROR: missing selected-event count: %s\n' "${selected_count_path}" >&2; exit 1; }
  selected="$(cat "${selected_count_path}")"
  case "${selected}" in
    ''|*[!0-9]*) printf 'ERROR: invalid selected-event count in %s: %s\n' "${selected_count_path}" "${selected}" >&2; exit 1 ;;
  esac
  printf 'selected [%04d] %s: %s/%s proposal events -> %s\n' "${n}" "${sample_label}" "${selected}" "${proposal_events}" "${output_path}"
}

[ -r "${plan}" ] || { printf 'ERROR: missing plan: %s\n' "${plan}" >&2; exit 1; }

count=0
while IFS=$'\t' read -r enabled generator versions knobs beam_modes beam_species interactions fsi_states input_template sample_template variation_template run_primary run_nuclear_exit notes; do
  case "${enabled}" in
    true|TRUE|1|yes|YES) ;;
    *) continue ;;
  esac
  want_generator "${generator}" || continue

  for version in $(split_csv "${versions}"); do
    for knob in $(split_csv "${knobs}"); do
      for beam_mode in $(split_csv "${beam_modes}"); do
        for species in $(split_csv "${beam_species}"); do
          for interaction in $(split_csv "${interactions}"); do
            for fsi_state in $(split_csv "${fsi_states}"); do
              count=$((count + 1))
              if [ "${max_samples}" -gt 0 ] && [ "${count}" -gt "${max_samples}" ]; then
                printf 'stopped after max_samples=%s\n' "${max_samples}"
                exit 0
              fi
              run_sample "${count}"
            done
          done
        done
      done
    done
  done
done < <(tail -n +2 "${plan}")

printf 'matrix samples considered: %d\n' "${count}"
