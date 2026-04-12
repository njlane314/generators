#!/usr/bin/env bash
set -euo pipefail

plan="${1:-analysis/config/generator_loop_plan.tsv}"
output_dir="${2:-analysis/output/matrix}"
working_point="${3:-nominal}"
energy_min="${4:-0.0}"
energy_max="${5:-10.0}"

root -l -b -q "analysis/run_generator_matrix.cxx+(\"${plan}\",\"${output_dir}\",\"${working_point}\",${energy_min},${energy_max})"
