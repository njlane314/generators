#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${script_dir}/../lib/common.sh"
plan="${1:-${ANA_SHARE_DIR}/config/sample_matrix.tsv}"
output_dir="${2:-analysis/output/matrix}"
working_point="${3:-nominal}"
energy_min="${4:-0.0}"
energy_max="${5:-10.0}"

ana_cd_repo
root -l -b -q "${ANA_ROOT_DIR}/matrix/generator_matrix.cxx+(\"${plan}\",\"${output_dir}\",\"${working_point}\",${energy_min},${energy_max})"
