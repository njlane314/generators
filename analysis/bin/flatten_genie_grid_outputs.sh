#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${script_dir}/../lib/common.sh"
repo_root="${ANA_REPO_ROOT}"

if [ "$#" -lt 2 ]; then
  cat >&2 <<'EOF'
Usage: analysis/bin/flatten_genie_grid_outputs.sh OUTPUT_FLAT INPUT_GHEP [INPUT_GHEP ...]

Convert GENIE GHEP files returned from the grid into one common NUISANCE
FlatTree analysis ntuple. Set beam_mode and beam_species for the sample being
flattened, for example:

  beam_mode=FHC beam_species=numu \
    analysis/bin/flatten_genie_grid_outputs.sh analysis/output/nuisance_flat/GENIE/sample.flat.root /pnfs/.../*.ghep.root
EOF
  exit 1
fi

output_flat="$1"
shift

beam_mode="${beam_mode:-FHC}"
beam_species="${beam_species:-numu}"
target="${target:-1000180400}"
prepare_flux="${prepare_flux:-$(ana_flux_root),$(ana_flux_hist "${beam_mode}" "${beam_species}")}"
prepare_flux_file="${prepare_flux%%,*}"
sample="${sample:-$(basename "${output_flat%.flat.root}")}"
workdir="${workdir:-${repo_root}/analysis/output/work/GENIE_GRID/${sample}}"
skim_final_state="${skim_final_state:-1}"
count_path="${count_path:-${output_flat%.root}.selected.count}"

inputs=("$@")
if [ "${#inputs[@]}" -eq 1 ] && [ -d "${inputs[0]}" ]; then
  inputs=("${inputs[0]}"/*.ghep.root)
fi
[ "${#inputs[@]}" -gt 0 ] || ana_die "no grid GHEP inputs supplied"
[ -e "${inputs[0]}" ] || ana_die "no grid GHEP inputs matched"

ana_check_cmds PrepareGENIE nuisflat
[ "${skim_final_state}" = 1 ] && ana_check_cmds root
ana_check_files "${prepare_flux_file}"
mkdir -p "${workdir}" "$(dirname "${output_flat}")"

flat_inputs=""
for ghep in "${inputs[@]}"; do
  [ -s "${ghep}" ] || ana_die "missing grid GHEP input: ${ghep}"
  name="$(basename "${ghep}")"
  name="${name%.ghep.root}"
  prep="${workdir}/${name}.gprep.root"
  raw_flat="${workdir}/${name}.flat.root"

  PrepareGENIE -i "${ghep}" -f "${prepare_flux}" -t "${target}[1]" -o "${prep}"
  nuisflat -i "GENIE:${prep}" -o "${raw_flat}"

  if [ -z "${flat_inputs}" ]; then
    flat_inputs="${raw_flat}"
  else
    flat_inputs="${flat_inputs},${raw_flat}"
  fi
done

if [ "${skim_final_state}" = 1 ]; then
  ana_build_analysis_ntuple "${ANA_ROOT_DIR}/ntuples/build_analysis_ntuple.cxx" "${flat_inputs}" "${output_flat}" -1 "${count_path}"
else
  ana_check_cmds hadd
  IFS=',' read -r -a flat_array <<< "${flat_inputs}"
  hadd -f "${output_flat}" "${flat_array[@]}"
fi
