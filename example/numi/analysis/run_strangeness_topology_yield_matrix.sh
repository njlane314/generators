#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 MANIFEST.tsv OUTPUT.csv" >&2
  echo "Manifest columns: input_flat_root<TAB>generator<TAB>knob<TAB>sample_label" >&2
  echo "The generator, knob, and sample_label columns may be left blank." >&2
  exit 1
fi

manifest="$1"
output="$2"

if [[ ! -f "${manifest}" ]]; then
  echo "ERROR: missing manifest: ${manifest}" >&2
  exit 1
fi

mkdir -p "$(dirname "${output}")"
root -l -b -q "${script_dir}/yield_matrix.cxx+(\"${manifest}\",\"${output}\")"
