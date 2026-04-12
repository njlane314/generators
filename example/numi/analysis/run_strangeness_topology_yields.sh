#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ $# -lt 1 || $# -gt 5 ]]; then
  echo "Usage: $0 INPUT.flat.root [OUTPUT.csv] [SAMPLE_LABEL] [GENERATOR] [KNOB]" >&2
  exit 1
fi

input="$1"
output="${2:-}"
label="${3:-}"
generator="${4:-}"
knob="${5:-}"

if [[ ! -f "${input}" ]]; then
  echo "ERROR: missing input flat tree: ${input}" >&2
  exit 1
fi

if [[ -z "${output}" && -z "${label}" && -z "${generator}" && -z "${knob}" ]]; then
  root -l -b -q "${script_dir}/strangeness_topology_yields.cxx+(\"${input}\")"
elif [[ -z "${label}" && -z "${generator}" && -z "${knob}" ]]; then
  root -l -b -q "${script_dir}/strangeness_topology_yields.cxx+(\"${input}\",\"${output}\")"
elif [[ -z "${generator}" && -z "${knob}" ]]; then
  root -l -b -q "${script_dir}/strangeness_topology_yields.cxx+(\"${input}\",\"${output}\",\"${label}\")"
elif [[ -z "${knob}" ]]; then
  root -l -b -q "${script_dir}/strangeness_topology_yields.cxx+(\"${input}\",\"${output}\",\"${label}\",\"${generator}\")"
else
  root -l -b -q "${script_dir}/strangeness_topology_yields.cxx+(\"${input}\",\"${output}\",\"${label}\",\"${generator}\",\"${knob}\")"
fi
