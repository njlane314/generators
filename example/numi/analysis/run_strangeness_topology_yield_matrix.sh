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

tmpdir="$(mktemp -d)"
trap 'rm -rf "${tmpdir}"' EXIT

mkdir -p "$(dirname "${output}")"
rm -f "${output}"

line_number=0
written_header=0

while IFS=$'\t' read -r input generator knob sample_label extra || [[ -n "${input:-}" ]]; do
  line_number=$((line_number + 1))

  [[ -z "${input// }" ]] && continue
  [[ "${input}" == \#* ]] && continue

  if [[ -n "${extra:-}" ]]; then
    echo "ERROR: too many columns in ${manifest}:${line_number}" >&2
    exit 1
  fi

  generator="${generator:-}"
  knob="${knob:-}"
  sample_label="${sample_label:-}"

  if [[ ! -f "${input}" ]]; then
    echo "ERROR: missing input flat tree in ${manifest}:${line_number}: ${input}" >&2
    exit 1
  fi

  tmp_csv="${tmpdir}/sample_${line_number}.csv"
  "${script_dir}/run_strangeness_topology_yields.sh" \
    "${input}" "${tmp_csv}" "${sample_label}" "${generator}" "${knob}"

  if [[ "${written_header}" == "0" ]]; then
    cat "${tmp_csv}" > "${output}"
    written_header=1
  else
    tail -n +2 "${tmp_csv}" >> "${output}"
  fi
done < "${manifest}"

if [[ "${written_header}" == "0" ]]; then
  echo "ERROR: manifest had no samples: ${manifest}" >&2
  exit 1
fi

echo "Wrote ${output}"
