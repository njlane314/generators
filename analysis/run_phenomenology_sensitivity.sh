#!/usr/bin/env bash
set -euo pipefail

if [[ $# -gt 6 ]]; then
  cat >&2 <<'EOF'
Usage: analysis/run_phenomenology_sensitivity.sh [YIELD_CSV] [OUTPUT_DIR] [CATEGORY] [GROUPING] [MIN_PULL] [MIN_FRAC]

YIELD_CSV defaults to analysis/output/expected_event_yields/expected_event_yields_by_variation.csv.
CATEGORY defaults to detector_visible_Lambda_to_p_piminus.
GROUPING is one of beam, beam_envelope, all_beam, generator; default beam.
MIN_PULL and MIN_FRAC define the "sensitive" flags; defaults are 2.0 and 0.20.
EOF
  exit 1
fi

input="${1:-analysis/output/expected_event_yields/expected_event_yields_by_variation.csv}"
out="${2:-analysis/output/phenomenology_sensitivity}"
category="${3:-detector_visible_Lambda_to_p_piminus}"
grouping="${4:-beam}"
min_pull="${5:-2.0}"
min_frac="${6:-0.20}"

if [[ ! -r "${input}" ]]; then
  printf 'missing input: %s\n' "${input}" >&2
  exit 2
fi

mkdir -p "${out}"
stem="$(printf '%s_%s' "${category}" "${grouping}" | tr ' /:;|(),#{}^=+' '_' | sed 's/-/minus/g')"

awk -f analysis/yield_sensitivity.awk \
  -v category="${category}" \
  -v grouping="${grouping}" \
  -v min_pull="${min_pull}" \
  -v min_frac="${min_frac}" \
  -v summary="${out}/${stem}_summary.tsv" \
  -v members="${out}/${stem}_members.tsv" \
  -v axes="${out}/${stem}_axes.tsv" \
  "${input}"

cat > "${out}/${stem}_README.txt" <<EOF
input	${input}
category	${category}
grouping	${grouping}
min_pull	${min_pull}
min_frac	${min_frac}

summary	${out}/${stem}_summary.tsv
members	${out}/${stem}_members.tsv
axes	${out}/${stem}_axes.tsv

members: favouring a row means shifting from the central row to that variation.
members: excluding a row is the leave-one-out envelope contraction.
axes: heuristic grouping of variations into primary_mechanism, transport,
generator_configuration, generator_version, or generator_model.
EOF

printf '%s\n' \
  "${out}/${stem}_summary.tsv" \
  "${out}/${stem}_members.tsv" \
  "${out}/${stem}_axes.tsv" \
  "${out}/${stem}_README.txt"
