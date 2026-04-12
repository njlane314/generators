#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
. "${script_dir}/sample_common.sh"

events="${events:-100}"
test_events="${test_events:-1000000}"
seed="${seed:-0}"
version="${version:-25.03.1}"
knob="${knob:-all_strange}"
beam_mode="${beam_mode:-FHC}"
beam_species="${beam_species:-numu}"
beam_particle="${beam_particle:-$(ana_probe_pdg "${beam_species}")}"
interaction="${interaction:-CC}"
fsi_state="${fsi_state:-fsi_on}"
proposal_min_gev="${proposal_min_gev:-0}"
proposal_max_gev="${proposal_max_gev:-10}"
proposal_bin_width_gev="${proposal_bin_width_gev:-0.005}"
skim_final_state="${skim_final_state:-1}"
sanitize_card="${sanitize_card:-1}"

base_card="${base_card:-${repo_root}/analysis/cards/NuWroCard_CC_Ar_numu.txt}"
generation_fluxfile="${generation_fluxfile:-$(ana_flux_dat "${repo_root}" "${beam_mode}" "${beam_species}")}"
sample="${sample:-NuWro_${version}_NuMI_${beam_mode}_${beam_species}_${interaction}_${knob}_filter_${fsi_state}}"
outdir="${outdir:-${repo_root}/analysis/output/proxy_flat/NuWro}"
workdir="${workdir:-${repo_root}/analysis/output/work/NuWro/${sample}}"
run_card="${run_card:-${workdir}/${sample}.card.txt}"

ana_check_cmds nuwro PrepareNuWroEvents nuisflat
[ "${skim_final_state}" = 1 ] && ana_check_cmds root
ana_check_files "${base_card}" "${generation_fluxfile}"
mkdir -p "${outdir}" "${workdir}"

beam_energy="$(ana_nuwro_beam_energy_from_flux "${generation_fluxfile}" "${proposal_min_gev}" "${proposal_max_gev}" "${proposal_bin_width_gev}")"
proposal_label="NuMI ${beam_mode} ${beam_species}: ${generation_fluxfile}"
generation_flux_integral="$(ana_flux_integral "${generation_fluxfile}" "${proposal_min_gev}" "${proposal_max_gev}")"
if [ "${sanitize_card}" = 1 ]; then
  awk '
    /^[[:space:]]*(sf_nuclearRecoil|sf_CoulombDistortion|sf_src|sf_pb)[[:space:]]*=/ {
      print "# disabled by analysis/run_nuwro.sh for installed NuWro compatibility: " $0
      next
    }
    { print }
  ' "${base_card}" > "${run_card}"
else
  cp "${base_card}" "${run_card}"
fi

{
  printf '\n################################################################################\n'
  printf '# ana Enu proposal and generator-variation override\n'
  printf '################################################################################\n\n'
  printf 'number_of_test_events = %s\n' "${test_events}"
  printf 'number_of_events = %s\n' "${events}"
  printf 'random_seed = %s\n' "${seed}"
  printf 'beam_type = 0\n'
  printf 'beam_particle = %s\n' "${beam_particle}"
  printf 'beam_energy = %s\n' "${beam_energy}"
  printf '# generation_flux_file = %s\n' "${generation_fluxfile}"
  printf '# generation_flux_integral_0_10gev = %s\n\n' "${generation_flux_integral}"

  case "${knob}" in
    all_strange|nominal)
      printf 'hyp_lambda = 1\nhyp_sigma_zero = 1\nhyp_sigma_minus = 1\n'
      ;;
    dis_only)
      printf 'dyn_qel_cc = 0\ndyn_res_cc = 0\ndyn_dis_cc = 1\ndyn_coh_cc = 0\ndyn_mec_cc = 0\ndyn_hyp_cc = 0\n'
      ;;
    hyp_all)
      printf 'dyn_qel_cc = 0\ndyn_res_cc = 0\ndyn_dis_cc = 0\ndyn_coh_cc = 0\ndyn_mec_cc = 0\ndyn_hyp_cc = 1\n'
      printf 'hyp_lambda = 1\nhyp_sigma_zero = 1\nhyp_sigma_minus = 1\n'
      ;;
    hyp_lambda_only)
      printf 'dyn_qel_cc = 0\ndyn_res_cc = 0\ndyn_dis_cc = 0\ndyn_coh_cc = 0\ndyn_mec_cc = 0\ndyn_hyp_cc = 1\n'
      printf 'hyp_lambda = 1\nhyp_sigma_zero = 0\nhyp_sigma_minus = 0\n'
      ;;
    hyp_sigma0_only)
      printf 'dyn_qel_cc = 0\ndyn_res_cc = 0\ndyn_dis_cc = 0\ndyn_coh_cc = 0\ndyn_mec_cc = 0\ndyn_hyp_cc = 1\n'
      printf 'hyp_lambda = 0\nhyp_sigma_zero = 1\nhyp_sigma_minus = 0\n'
      ;;
    hyp_sigmam_only)
      printf 'dyn_qel_cc = 0\ndyn_res_cc = 0\ndyn_dis_cc = 0\ndyn_coh_cc = 0\ndyn_mec_cc = 0\ndyn_hyp_cc = 1\n'
      printf 'hyp_lambda = 0\nhyp_sigma_zero = 0\nhyp_sigma_minus = 1\n'
      ;;
    *) ana_die "unsupported NuWro knob: ${knob}" ;;
  esac

  case "${fsi_state}" in
    fsi_on) printf 'FSI_on = 1\n' ;;
    fsi_off) printf 'FSI_on = 0\n' ;;
    *) ana_die "unsupported NuWro fsi_state: ${fsi_state}" ;;
  esac
} >> "${run_card}"

native="${sample}.native.root"
prep="${sample}.prep.root"
flat="${outdir}/${sample}.flat.root"
raw_flat="${flat}"
[ "${skim_final_state}" = 1 ] && raw_flat="${workdir}/${sample}.inclusive.flat.root"

printf 'NuWro sample: %s\n' "${sample}"
printf '  generation flux: %s\n' "${proposal_label}"
printf '  output: %s\n' "${flat}"

(
  cd "${workdir}"
  nuwro -i "${run_card}" -o "${native}"
  PrepareNuWroEvents -f "${native}" -o "${prep}"
  nuisflat -i "NuWro:${prep}" -o "${raw_flat}"
)
[ "${skim_final_state}" = 1 ] && ana_skim_final_state_hyperon "${script_dir}/skim_final_state_hyperon.cxx" "${raw_flat}" "${flat}" -1 "${skim_count_path:-}"
