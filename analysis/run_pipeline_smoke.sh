#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
events="${events:-100}"
plan="${plan:-${script_dir}/config/generator_smoke_plan.tsv}"

export events plan
exec "${script_dir}/run_pipeline.sh"
