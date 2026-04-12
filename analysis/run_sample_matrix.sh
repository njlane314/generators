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
batch_events_was_set=0
if [ -n "${batch_events+x}" ]; then
  batch_events_was_set=1
else
  batch_events="${events}"
fi
batch_events_min="${batch_events_min:-1000}"
max_batches="${max_batches:-200}"
gibuu_events_per_ensemble="${gibuu_events_per_ensemble:-5}"
target_events="${events}"

case "${target_events}" in
  ''|*[!0-9]*) printf 'ERROR: events must be a positive integer selected-event target\n' >&2; exit 1 ;;
esac
[ "${target_events}" -gt 0 ] || { printf 'ERROR: events must be > 0\n' >&2; exit 1; }
case "${batch_events_min}" in
  ''|*[!0-9]*) printf 'ERROR: batch_events_min must be a positive integer\n' >&2; exit 1 ;;
esac
[ "${batch_events_min}" -gt 0 ] || { printf 'ERROR: batch_events_min must be > 0\n' >&2; exit 1; }
case "${batch_events}" in
  ''|*[!0-9]*) printf 'ERROR: batch_events must be a positive integer\n' >&2; exit 1 ;;
esac
[ "${batch_events}" -gt 0 ] || { printf 'ERROR: batch_events must be > 0\n' >&2; exit 1; }
if [ "${batch_events_was_set}" = 0 ] && [ "${batch_events}" -lt "${batch_events_min}" ]; then
  batch_events="${batch_events_min}"
fi
case "${max_batches}" in
  ''|*[!0-9]*) printf 'ERROR: max_batches must be a positive integer\n' >&2; exit 1 ;;
esac
[ "${max_batches}" -gt 0 ] || { printf 'ERROR: max_batches must be > 0\n' >&2; exit 1; }
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

join_csv() {
  local IFS=,
  printf '%s' "$*"
}

run_generator_batch() {
  local n="$1"
  local batch="$2"
  local batch_sample="$3"
  local batch_outdir="$4"
  local batch_workdir="$5"
  local batch_count_file="$6"

  case "${generator}" in
    GENIE)
      events="${batch_events}" \
      version="${version}" \
      tune="${knob}" \
      beam_mode="${beam_mode}" \
      beam_species="${species}" \
      interaction="${interaction}" \
      run="$((seed_base + n * 1000 + batch))" \
      sample="${batch_sample}" \
      outdir="${batch_outdir}" \
      workdir="${batch_workdir}" \
      skim_final_state=1 \
      skim_count_path="${batch_count_file}" \
      "${script_dir}/run_genie.sh"
      ;;
    NuWro)
      events="${batch_events}" \
      seed="$((seed_base + n * 1000 + batch))" \
      version="${version}" \
      knob="${knob}" \
      beam_mode="${beam_mode}" \
      beam_species="${species}" \
      interaction="${interaction}" \
      fsi_state="${fsi_state}" \
      sample="${batch_sample}" \
      outdir="${batch_outdir}" \
      workdir="${batch_workdir}" \
      skim_final_state=1 \
      skim_count_path="${batch_count_file}" \
      "${script_dir}/run_nuwro.sh"
      ;;
    GiBUU)
      events="${batch_events}" \
      gibuu_events_per_ensemble="${gibuu_events_per_ensemble}" \
      num_runs_same_energy="${num_runs_same_energy:-1}" \
      version="${version}" \
      knob="${knob}" \
      beam_mode="${beam_mode}" \
      beam_species="${species}" \
      interaction="${interaction}" \
      fsi_state="${fsi_state}" \
      sample="${batch_sample}" \
      outdir="${batch_outdir}" \
      workdir="${batch_workdir}" \
      skim_final_state=1 \
      skim_count_path="${batch_count_file}" \
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
  local output_path sample_label reject_dir batch selected batch_sample batch_outdir batch_workdir batch_flat batch_count_file batch_selected count_file merged_inputs note
  local selected_files=()
  output_path="${repo_root}/$(expand_template "${input_template}")"
  sample_label="$(expand_template "${sample_template}")"

  if [ "${skip_existing}" = 1 ] && [ -s "${output_path}" ]; then
    printf 'skip existing [%04d] %s\n' "${n}" "${output_path}"
    return 0
  fi

  note=""
  if [ "${generator}" = "GiBUU" ]; then
    if [ -n "${num_ensembles:-}" ]; then
      note=" num_ensembles=${num_ensembles}"
    else
      note=" corrected_num_ensembles=$(( (batch_events + gibuu_events_per_ensemble - 1) / gibuu_events_per_ensemble ))"
    fi
  fi
  printf 'sample [%04d] %s version=%s knob=%s beam=%s species=%s interaction=%s fsi=%s target_selected=%s batch_events=%s max_batches=%s%s\n' \
    "${n}" "${generator}" "${version}" "${knob}" "${beam_mode}" "${species}" "${interaction}" "${fsi_state}" "${target_events}" "${batch_events}" "${max_batches}" "${note}"

  if [ "${dry_run}" = 1 ]; then
    return 0
  fi

  reject_dir="${repo_root}/analysis/output/work/rejection/${sample_label}"
  rm -rf "${reject_dir}"
  mkdir -p "${reject_dir}/selected" "$(dirname "${output_path}")"

  selected=0
  batch=0
  while [ "${selected}" -lt "${target_events}" ] && [ "${batch}" -lt "${max_batches}" ]; do
    batch=$((batch + 1))
    batch_sample="${sample_label}.batch_${batch}"
    batch_outdir="${reject_dir}/selected"
    batch_workdir="${reject_dir}/work_${batch}"
    batch_flat="${batch_outdir}/${batch_sample}.flat.root"
    batch_count_file="${batch_outdir}/${batch_sample}.count"

    run_generator_batch "${n}" "${batch}" "${batch_sample}" "${batch_outdir}" "${batch_workdir}" "${batch_count_file}"
    selected_files+=("${batch_flat}")

    [ -s "${batch_count_file}" ] || { printf 'ERROR: missing selected-event count: %s\n' "${batch_count_file}" >&2; exit 1; }
    batch_selected="$(cat "${batch_count_file}")"
    case "${batch_selected}" in
      ''|*[!0-9]*) printf 'ERROR: invalid selected-event count in %s: %s\n' "${batch_count_file}" "${batch_selected}" >&2; exit 1 ;;
    esac
    selected=$((selected + batch_selected))
    printf 'selected [%04d] %s: %s/%s after %s batch(es); last batch +%s\n' "${n}" "${sample_label}" "${selected}" "${target_events}" "${batch}" "${batch_selected}"
  done

  [ "${selected}" -ge "${target_events}" ] || {
    printf 'ERROR: only selected %s/%s events for %s after %s batches\n' "${selected}" "${target_events}" "${sample_label}" "${max_batches}" >&2
    exit 1
  }

  count_file="${reject_dir}/selected.count"
  merged_inputs="$(join_csv "${selected_files[@]}")"
  ana_skim_final_state_hyperon "${script_dir}/skim_final_state_hyperon.cxx" "${merged_inputs}" "${output_path}" "${target_events}" "${count_file}"
  selected="$(cat "${count_file}")"
  printf 'final selected [%04d] %s: %s/%s -> %s\n' "${n}" "${sample_label}" "${selected}" "${target_events}" "${output_path}"
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
