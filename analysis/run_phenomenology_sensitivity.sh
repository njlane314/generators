#!/usr/bin/env bash
set -euo pipefail

if [[ $# -gt 6 ]]; then
  cat >&2 <<'EOF'
Usage: analysis/run_phenomenology_sensitivity.sh [YIELD_CSV] [OUTPUT_DIR] [CATEGORY] [GROUPING] [MIN_PULL] [MIN_FRAC]

YIELD_CSV defaults to analysis/output/expected_event_yields/expected_event_yields_by_variation.csv.
CATEGORY defaults to detector_visible_Lambda_to_p_piminus.
GROUPING is one of beam, beam_envelope, all_beam, generator; default beam.
MIN_PULL and MIN_FRAC define the "sensitive" flags; defaults are 2.0 and 0.20.
Set make_plot=0 to skip the ranked-impact SVG.
EOF
  exit 1
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
input="${1:-analysis/output/expected_event_yields/expected_event_yields_by_variation.csv}"
out="${2:-analysis/output/phenomenology_sensitivity}"
category="${3:-detector_visible_Lambda_to_p_piminus}"
grouping="${4:-beam}"
min_pull="${5:-2.0}"
min_frac="${6:-0.20}"
make_plot="${make_plot:-1}"

if [[ ! -r "${input}" ]]; then
  printf 'missing input: %s\n' "${input}" >&2
  exit 2
fi

mkdir -p "${out}"
stem="$(printf '%s_%s' "${category}" "${grouping}" | tr ' /:;|(),#{}^=+' '_' | sed 's/-/minus/g')"
summary_path="${out}/${stem}_summary.tsv"
members_path="${out}/${stem}_members.tsv"
axes_path="${out}/${stem}_axes.tsv"
readme_path="${out}/${stem}_README.txt"
plot_path="${out}/${stem}_ranked_impacts.svg"

awk -f "${script_dir}/yield_sensitivity.awk" \
  -v category="${category}" \
  -v grouping="${grouping}" \
  -v min_pull="${min_pull}" \
  -v min_frac="${min_frac}" \
  -v summary="${summary_path}" \
  -v members="${members_path}" \
  -v axes="${axes_path}" \
  "${input}"

if [[ "${make_plot}" = 1 ]]; then
  "${script_dir}/plot_phenomenology_sensitivity.sh" "${members_path}" "${plot_path}" >/dev/null
else
  plot_path=""
fi

cat > "${readme_path}" <<EOF
input	${input}
category	${category}
grouping	${grouping}
min_pull	${min_pull}
min_frac	${min_frac}

summary	${summary_path}
members	${members_path}
axes	${axes_path}
plot	${plot_path:-not_written}

members: favouring a row means shifting from the central row to that variation.
members: excluding a row is the leave-one-out envelope contraction.
axes: heuristic grouping of variations into primary_mechanism, transport,
generator_configuration, generator_version, or generator_model.
plot: ranked-impact SVG of the largest signed fractional yield shifts from the
central row, colored by model_axis and annotated with the statistical pull.
EOF

if [[ -n "${plot_path}" ]]; then
  printf '%s\n' "${summary_path}" "${members_path}" "${axes_path}" "${plot_path}" "${readme_path}"
else
  printf '%s\n' "${summary_path}" "${members_path}" "${axes_path}" "${readme_path}"
fi
