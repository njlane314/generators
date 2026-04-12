#!/usr/bin/env bash
set -euo pipefail

if [[ $# -gt 4 ]]; then
  cat >&2 <<'EOF'
Usage: analysis/plot_phenomenology_sensitivity.sh [MEMBERS_TSV] [OUTPUT_SVG] [MAX_ROWS] [METRIC]

MEMBERS_TSV defaults to analysis/output/phenomenology_sensitivity/detector_visible_Lambda_to_p_piminus_beam_members.tsv.
OUTPUT_SVG defaults to the input path with _members.tsv replaced by _ranked_impacts.svg.
MAX_ROWS defaults to 24.
METRIC is one of frac, pull, or events; default frac.
EOF
  exit 1
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
input="${1:-analysis/output/phenomenology_sensitivity/detector_visible_Lambda_to_p_piminus_beam_members.tsv}"
output="${2:-}"
max_rows="${3:-24}"
metric="${4:-frac}"

if ! command -v python3 >/dev/null 2>&1; then
  printf 'missing command: python3\n' >&2
  exit 2
fi

args=(--input "${input}" --max-rows "${max_rows}" --metric "${metric}")
if [[ -n "${output}" ]]; then
  args+=(--output "${output}")
fi

python3 "${script_dir}/plot_phenomenology_sensitivity.py" "${args[@]}"
