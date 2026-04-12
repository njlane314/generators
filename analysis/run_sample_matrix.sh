#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
plan="${1:-${script_dir}/config/generator_loop_plan.tsv}"
events="${events:-100}"
dry_run="${dry_run:-0}"
skip_existing="${skip_existing:-0}"
generator_filter="${generator_filter:-}"
max_samples="${max_samples:-0}"
seed_base="${seed_base:-1000}"

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

run_sample() {
  local n="$1"
  local output_path
  output_path="${repo_root}/$(expand_template "${input_template}")"

  if [ "${skip_existing}" = 1 ] && [ -s "${output_path}" ]; then
    printf 'skip existing [%04d] %s\n' "${n}" "${output_path}"
    return 0
  fi

  printf 'sample [%04d] %s version=%s knob=%s beam=%s species=%s interaction=%s fsi=%s\n' \
    "${n}" "${generator}" "${version}" "${knob}" "${beam_mode}" "${species}" "${interaction}" "${fsi_state}"

  if [ "${dry_run}" = 1 ]; then
    return 0
  fi

  case "${generator}" in
    GENIE)
      events="${events}" \
      version="${version}" \
      tune="${knob}" \
      beam_mode="${beam_mode}" \
      beam_species="${species}" \
      interaction="${interaction}" \
      run="$((seed_base + n))" \
      "${script_dir}/run_genie.sh"
      ;;
    NuWro)
      events="${events}" \
      seed="$((seed_base + n))" \
      version="${version}" \
      knob="${knob}" \
      beam_mode="${beam_mode}" \
      beam_species="${species}" \
      interaction="${interaction}" \
      fsi_state="${fsi_state}" \
      "${script_dir}/run_nuwro.sh"
      ;;
    GiBUU)
      num_ensembles="${num_ensembles:-${events}}" \
      num_runs_same_energy="${num_runs_same_energy:-1}" \
      version="${version}" \
      knob="${knob}" \
      beam_mode="${beam_mode}" \
      beam_species="${species}" \
      interaction="${interaction}" \
      fsi_state="${fsi_state}" \
      "${script_dir}/run_gibuu.sh"
      ;;
    *)
      printf 'ERROR: unsupported generator: %s\n' "${generator}" >&2
      exit 1
      ;;
  esac
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
