#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${script_dir}/../lib/common.sh"
repo_root="${ANA_REPO_ROOT}"

events="${events:-100}"
version="${version:-2025}"
knob="${knob:-DIS_only}"
beam_mode="${beam_mode:-FHC}"
beam_species="${beam_species:-numu}"
interaction="${interaction:-CC}"
fsi_state="${fsi_state:-fsi_on}"
num_runs_same_energy="${num_runs_same_energy:-1}"
gibuu_events_per_ensemble="${gibuu_events_per_ensemble:-5}"
case "${events}" in
  ''|*[!0-9]*) ana_die "events must be a positive integer" ;;
esac
case "${gibuu_events_per_ensemble}" in
  ''|*[!0-9]*) ana_die "gibuu_events_per_ensemble must be a positive integer" ;;
esac
[ "${events}" -gt 0 ] || ana_die "events must be > 0"
[ "${gibuu_events_per_ensemble}" -gt 0 ] || ana_die "gibuu_events_per_ensemble must be > 0"
if [ -z "${num_ensembles:-}" ]; then
  num_ensembles="$(( (events + gibuu_events_per_ensemble - 1) / gibuu_events_per_ensemble ))"
fi
case "${num_ensembles}" in
  ''|*[!0-9]*) ana_die "num_ensembles must be a positive integer" ;;
esac
[ "${num_ensembles}" -gt 0 ] || ana_die "num_ensembles must be > 0"
beam_energy_min_gev="${beam_energy_min_gev:-${proposal_min_gev:-0}}"
beam_energy_max_gev="${beam_energy_max_gev:-${proposal_max_gev:-10}}"
proposal_min_gev="${proposal_min_gev:-${beam_energy_min_gev}}"
proposal_max_gev="${proposal_max_gev:-${beam_energy_max_gev}}"
skim_final_state="${skim_final_state:-1}"

base_job="${base_job:-${ANA_SHARE_DIR}/cards/GiBUU2025_numu.job}"
if [ -z "${gibuu_bin:-}" ]; then
  for candidate in \
    "${repo_root}/GiBUU/release/testRun/GiBUU.x" \
    "${repo_root}/../GiBUU/release/testRun/GiBUU.x"
  do
    if [ -x "${candidate}" ]; then
      gibuu_bin="${candidate}"
      break
    fi
  done
  if [ -z "${gibuu_bin:-}" ]; then
    if command -v GiBUU.x >/dev/null 2>&1; then
      gibuu_bin="$(command -v GiBUU.x)"
    elif command -v gibuu >/dev/null 2>&1; then
      gibuu_bin="$(command -v gibuu)"
    else
      gibuu_bin="${repo_root}/GiBUU/release/testRun/GiBUU.x"
    fi
  fi
fi
if [ -z "${gibuu_input:-}" ]; then
  for candidate in \
    "${repo_root}/GiBUU/buuinput" \
    "${repo_root}/../GiBUU/buuinput"
  do
    if [ -d "${candidate}" ]; then
      gibuu_input="${candidate}"
      break
    fi
  done
  gibuu_input="${gibuu_input:-${repo_root}/GiBUU/buuinput}"
fi
generation_fluxfile="${generation_fluxfile:-$(ana_flux_dat "${beam_mode}" "${beam_species}")}"
preparation_fluxfile="${preparation_fluxfile:-$(ana_flux_root)}"
preparation_fluxhisto="${preparation_fluxhisto:-$(ana_flux_hist "${beam_mode}" "${beam_species}")}"
sample="${sample:-GiBUU${version}_NuMI_${beam_mode}_${beam_species}_${interaction}_all_strange_filter_${fsi_state}}"
outdir="${outdir:-${repo_root}/analysis/output/nuisance_flat/GiBUU}"
workdir="${workdir:-${repo_root}/analysis/output/work/GiBUU/${sample}}"
job_card="${job_card:-${workdir}/${sample}.job}"
gibuu_event_root="${gibuu_event_root:-EventOutput.Pert.00000001.root}"

ana_check_cmds PrepareGiBUU nuisflat awk
[ "${skim_final_state}" = 1 ] && ana_check_cmds root
ana_check_exe "${gibuu_bin}"
ana_check_files "${base_job}" "${preparation_fluxfile}"
[ -d "${gibuu_input}" ] || ana_die "missing GiBUU input directory: ${gibuu_input}"
mkdir -p "${outdir}" "${workdir}"

process_id="$(ana_gibuu_process_id "${interaction}" "${beam_species}")"
case "${fsi_state}" in
  fsi_on) num_time_steps="${num_time_steps:-150}" ;;
  fsi_off) num_time_steps="${num_time_steps:-0}" ;;
  *) ana_die "unsupported GiBUU fsi_state: ${fsi_state}" ;;
esac

case "${knob}" in
  DIS_only|dis_only)
    includeQE=F
    includeDELTA=F
    includeRES=F
    include1pi=F
    includeDIS=T
    include2p2hQE=F
    include2pi=F
    ;;
  all_strange|nominal|all)
    includeQE=T
    includeDELTA=T
    includeRES=T
    include1pi=T
    includeDIS=T
    include2p2hQE=T
    include2pi=T
    ;;
  *) ana_die "unsupported GiBUU knob: ${knob}" ;;
esac

proposal_fluxfile="${proposal_fluxfile:-${workdir}/numi_${beam_mode}_${beam_species}_${proposal_min_gev}_${proposal_max_gev}gev.dat}"
proposal_label="NuMI ${beam_mode} ${beam_species}: ${proposal_fluxfile}"
ana_check_files "${generation_fluxfile}"
ana_write_flux_dat_window "${proposal_fluxfile}" "${generation_fluxfile}" "${proposal_min_gev}" "${proposal_max_gev}"

awk \
  -v process_id="${process_id}" \
  -v flux="${proposal_fluxfile}" \
  -v input="${gibuu_input}" \
  -v version="${version}" \
  -v num_runs_same_energy="${num_runs_same_energy}" \
  -v num_ensembles="${num_ensembles}" \
  -v num_time_steps="${num_time_steps}" \
  -v includeQE="${includeQE}" \
  -v includeDELTA="${includeDELTA}" \
  -v includeRES="${includeRES}" \
  -v include1pi="${include1pi}" \
  -v includeDIS="${includeDIS}" \
  -v include2p2hQE="${include2p2hQE}" \
  -v include2pi="${include2pi}" '
    /^[[:space:]]*process_ID[[:space:]]*=/ { printf "      process_ID      =  %s\n", process_id; next }
    /^[[:space:]]*FileNameFlux[[:space:]]*=/ { printf "        FileNameFlux = \047%s\047\n", flux; next }
    /^[[:space:]]*includeQE[[:space:]]*=/ { printf "      includeQE       = %s\n", includeQE; next }
    /^[[:space:]]*includeDELTA[[:space:]]*=/ { printf "      includeDELTA    = %s\n", includeDELTA; next }
    /^[[:space:]]*includeRES[[:space:]]*=/ { printf "      includeRES      = %s\n", includeRES; next }
    /^[[:space:]]*include1pi[[:space:]]*=/ { printf "      include1pi      = %s\n", include1pi; next }
    /^[[:space:]]*includeDIS[[:space:]]*=/ { printf "      includeDIS      = %s\n", includeDIS; next }
    /^[[:space:]]*include2p2hQE[[:space:]]*=/ { printf "      include2p2hQE   = %s\n", include2p2hQE; next }
    /^[[:space:]]*include2pi[[:space:]]*=/ { printf "      include2pi      = %s\n", include2pi; next }
    /^[[:space:]]*num_runs_SameEnergy[[:space:]]*=/ { printf "      num_runs_SameEnergy=%s\n", num_runs_same_energy; next }
    /^[[:space:]]*numEnsembles[[:space:]]*=/ { printf "      numEnsembles=%s\n", num_ensembles; next }
    /^[[:space:]]*numTimeSteps[[:space:]]*=/ { printf "      numTimeSteps=%s\n", num_time_steps; next }
    /^[[:space:]]*path_to_input[[:space:]]*=/ { printf "      path_to_input=\047%s\047\n", input; next }
    /^[[:space:]]*version[[:space:]]*=/ { printf "      version = %s\n", version; next }
    { print }
  ' "${base_job}" > "${job_card}"

prep="${sample}.prep.root"
flat="${outdir}/${sample}.flat.root"
raw_flat="${flat}"
[ "${skim_final_state}" = 1 ] && raw_flat="${workdir}/${sample}.inclusive.flat.root"

printf 'GiBUU sample: %s\n' "${sample}"
printf '  requested proposal events: %s\n' "${events}"
printf '  GiBUU numEnsembles: %s (events/ensemble correction: %s)\n' "${num_ensembles}" "${gibuu_events_per_ensemble}"
printf '  generation flux: %s\n' "${proposal_label}"
printf '  preparation flux: %s,%s\n' "${preparation_fluxfile}" "${preparation_fluxhisto}"
printf '  beam energy envelope: %s-%s GeV\n' "${proposal_min_gev}" "${proposal_max_gev}"
printf '  output: %s\n' "${flat}"

(
  cd "${workdir}"
  "${gibuu_bin}" < "${job_card}"
  [ -s "${gibuu_event_root}" ] || ana_die "GiBUU did not write ${workdir}/${gibuu_event_root}"
  PrepareGiBUU -i "${gibuu_event_root}" -f "${preparation_fluxfile},${preparation_fluxhisto}" -o "${prep}"
  nuisflat -i "GiBUU:${prep}" -o "${raw_flat}"
)
[ "${skim_final_state}" = 1 ] && ana_build_analysis_ntuple "${ANA_ROOT_DIR}/ntuples/build_analysis_ntuple.cxx" "${raw_flat}" "${flat}" -1 "${skim_count_path:-}"
