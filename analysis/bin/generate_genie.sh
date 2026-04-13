#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${script_dir}/../lib/common.sh"
repo_root="${ANA_REPO_ROOT}"

events="${events:-100}"
version="${version:-3.6.2}"
version_label="$(ana_genie_version_label "${version}")"
genie_version_tag="$(ana_genie_version_tag "${version}")"
tune="${tune:-AR23_20i_00_000}"
beam_mode="${beam_mode:-FHC}"
beam_species="${beam_species:-numu}"
probe="${probe:-$(ana_probe_pdg "${beam_species}")}"
target="${target:-1000180400}"
interaction="${interaction:-CC}"
run="${run:-1}"
beam_energy_min_gev="${beam_energy_min_gev:-${proposal_min_gev:-0}}"
beam_energy_max_gev="${beam_energy_max_gev:-${proposal_max_gev:-10}}"
proposal_min_gev="${proposal_min_gev:-${beam_energy_min_gev}}"
proposal_max_gev="${proposal_max_gev:-${beam_energy_max_gev}}"
skim_final_state="${skim_final_state:-1}"
make_spline="${make_spline:-auto}"

generation_fluxfile="${generation_fluxfile:-$(ana_flux_root)}"
generation_fluxhisto="${generation_fluxhisto:-$(ana_flux_hist "${beam_mode}" "${beam_species}")}"
sample="${sample:-GENIE_${version_label}_NuMI_${beam_mode}_${beam_species}_${interaction}_all_strange_filter_${tune}}"
outdir="${outdir:-${repo_root}/analysis/output/nuisance_flat/GENIE}"
workdir="${workdir:-${repo_root}/analysis/output/work/GENIE/${sample}}"
spline="${spline:-${ANA_SHARE_DIR}/splines/${probe}_${target}_${interaction}_${genie_version_tag}_${tune}.xml}"

ana_check_cmds gevgen gntpc PrepareGENIE nuisflat
[ "${skim_final_state}" = 1 ] && ana_check_cmds root
mkdir -p "${outdir}" "${workdir}"

if [ ! -s "${spline}" ]; then
  case "${make_spline}" in
    1|true|TRUE|yes|YES|auto)
      ana_check_cmds gmkspl
      mkdir -p "$(dirname "${spline}")"
      gmkspl \
        -p "${probe}" \
        -t "${target}" \
        -e "${proposal_max_gev}" \
        -o "${spline}" \
        --tune "${tune}" \
        --event-generator-list "${interaction}"
      ;;
    *) ;;
  esac
fi

proposal_flux="${proposal_flux:-${generation_fluxfile},${generation_fluxhisto}}"
prepare_flux="${prepare_flux:-${generation_fluxfile},${generation_fluxhisto}}"
proposal_label="NuMI ${beam_mode} ${beam_species}: ${generation_fluxfile},${generation_fluxhisto}"
ana_check_files "${spline}" "${generation_fluxfile}"

ghep="${workdir}/${sample}.ghep.root"
gst="${workdir}/${sample}.gst.root"
prep="${workdir}/${sample}.gprep.root"
flat="${outdir}/${sample}.flat.root"
raw_flat="${flat}"
[ "${skim_final_state}" = 1 ] && raw_flat="${workdir}/${sample}.inclusive.flat.root"

printf 'GENIE sample: %s\n' "${sample}"
printf '  generation flux: %s\n' "${proposal_label}"
printf '  beam energy envelope: %s-%s GeV\n' "${proposal_min_gev}" "${proposal_max_gev}"
printf '  output: %s\n' "${flat}"

gevgen \
  -n "${events}" \
  -p "${probe}" \
  -t "${target}" \
  -e "${proposal_min_gev},${proposal_max_gev}" \
  -f "${proposal_flux}" \
  -r "${run}" \
  --event-generator-list "${interaction}" \
  --tune "${tune}" \
  --cross-sections "${spline}" \
  -o "${ghep}"

gntpc -f gst -i "${ghep}" -o "${gst}" --tune "${tune}"
PrepareGENIE -i "${ghep}" -f "${prepare_flux}" -t "${target}[1]" -o "${prep}"
nuisflat -i "GENIE:${prep}" -o "${raw_flat}"
[ "${skim_final_state}" = 1 ] && ana_build_analysis_ntuple "${ANA_ROOT_DIR}/ntuples/build_analysis_ntuple.cxx" "${raw_flat}" "${flat}" -1 "${skim_count_path:-}"
