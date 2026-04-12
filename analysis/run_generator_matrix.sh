#!/usr/bin/env bash
set -euo pipefail

plan="${1:-ana/config/generator_loop_plan.tsv}"
output_dir="${2:-analysis/output/matrix}"
working_point="${3:-nominal}"
flux_file="${4:-example/numi/flux/microboone_numi_flux_5mev.root}"
flux_floor="${5:-0.0}"
proposal_emin="${6:-0.0}"
proposal_emax="${7:-10.0}"

root -l -b -q "analysis/run_generator_matrix.cxx+(\"${plan}\",\"${output_dir}\",\"${working_point}\",\"${flux_file}\",${flux_floor},${proposal_emin},${proposal_emax})"
