#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${script_dir}/../lib/common.sh"

if [ "$#" -lt 2 ]; then
  cat >&2 <<'EOF'
Usage: analysis/bin/build_analysis_ntuple.sh OUTPUT_FLAT INPUT_FLAT [INPUT_FLAT ...]

Build the common analysis ntuple from one or more NUISANCE FlatTree files.
The current pass keeps events with a final-state Lambda, Sigma0, Xi0, Xi-, or
Omega-. Add custom per-analysis branches in
analysis/libexec/root/ntuples/build_analysis_ntuple.cxx when needed.
EOF
  exit 1
fi

output_flat="$1"
shift
count_path="${count_path:-${output_flat%.root}.selected.count}"
max_entries="${max_entries:--1}"

input_csv=""
for input in "$@"; do
  [ -s "${input}" ] || ana_die "missing input flat: ${input}"
  if [ -z "${input_csv}" ]; then
    input_csv="${input}"
  else
    input_csv="${input_csv},${input}"
  fi
done

ana_build_analysis_ntuple "${ANA_ROOT_DIR}/ntuples/build_analysis_ntuple.cxx" "${input_csv}" "${output_flat}" "${max_entries}" "${count_path}"
