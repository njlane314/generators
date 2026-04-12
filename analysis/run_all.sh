#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
generators="${generators:-GENIE NuWro GiBUU}"

for generator in ${generators}; do
  case "${generator}" in
    GENIE|genie) "${script_dir}/run_genie.sh" ;;
    NuWro|nuwro) "${script_dir}/run_nuwro.sh" ;;
    GiBUU|gibuu) "${script_dir}/run_gibuu.sh" ;;
    *) printf 'ERROR: unknown generator: %s\n' "${generator}" >&2; exit 1 ;;
  esac
done
