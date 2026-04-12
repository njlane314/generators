#!/usr/bin/env bash

ana_die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

ana_lower() {
  printf '%s' "$1" | tr '[:upper:]' '[:lower:]'
}

ana_check_cmds() {
  local missing
  missing=""
  while [ "$#" -gt 0 ]; do
    if ! command -v "$1" >/dev/null 2>&1; then
      missing="${missing} $1"
    fi
    shift
  done
  [ -z "${missing}" ] || ana_die "missing required commands:${missing}"
}

ana_check_exe() {
  case "$1" in
    */*) [ -x "$1" ] || ana_die "missing executable: $1" ;;
    *) command -v "$1" >/dev/null 2>&1 || ana_die "missing executable: $1" ;;
  esac
}

ana_check_files() {
  while [ "$#" -gt 0 ]; do
    [ -s "$1" ] || ana_die "missing file: $1"
    shift
  done
}

ana_probe_pdg() {
  case "$1" in
    numu) printf '%s\n' 14 ;;
    numubar) printf '%s\n' -14 ;;
    *) ana_die "unknown beam_species: $1" ;;
  esac
}

ana_gibuu_process_id() {
  local interaction species
  interaction="$(ana_lower "$1")"
  species="$2"
  case "${interaction}:${species}" in
    cc:numu) printf '%s\n' 2 ;;
    cc:numubar) printf '%s\n' -2 ;;
    nc:numu) printf '%s\n' 3 ;;
    nc:numubar) printf '%s\n' -3 ;;
    *) ana_die "unsupported GiBUU interaction/species: $1/$2" ;;
  esac
}

ana_genie_version_label() {
  case "$1" in
    v*) printf '%s' "$1" | sed 's/^v//; s/_/./g' ;;
    *) printf '%s' "$1" ;;
  esac
}

ana_genie_version_tag() {
  case "$1" in
    v*) printf '%s' "$1" ;;
    *) printf 'v%s' "$(printf '%s' "$1" | tr . _)" ;;
  esac
}

ana_flux_root() {
  printf '%s/example/numi/flux/microboone_numi_flux_5mev.root\n' "$1"
}

ana_flux_hist() {
  local mode species
  mode="$(ana_lower "$1")"
  species="$2"
  case "${mode}:${species}" in
    fhc:numu|fhc:numubar|rhc:numu|rhc:numubar)
      printf '%s/%s/Detsmear/%s_CV_AV_TPC_5MeV_bin\n' "${mode}" "${species}" "${species}"
      ;;
    *) ana_die "unsupported flux beam/species: $1/$2" ;;
  esac
}

ana_flux_dat() {
  local mode species
  mode="$(ana_lower "$2")"
  species="$3"
  case "${mode}:${species}" in
    fhc:numu|fhc:numubar|rhc:numu|rhc:numubar)
      printf '%s/example/numi/flux/gibuu_numi_%s_%s.dat\n' "$1" "${mode}" "${species}"
      ;;
    *) ana_die "unsupported flux beam/species: $2/$3" ;;
  esac
}

ana_nuwro_beam_energy() {
  awk -v min_gev="$1" -v max_gev="$2" -v bin_width_gev="$3" '
    BEGIN {
      emin = min_gev * 1000.0
      emax = max_gev * 1000.0
      step = bin_width_gev * 1000.0
      if (step <= 0 || emax <= emin) exit 1
      n = int(((emax - emin) / step) + 0.5)
      printf "%.12g %.12g", emin, emax
      for (i = 0; i < n; ++i) printf " 1"
      printf "\n"
    }'
}

ana_write_flat_flux_dat() {
  mkdir -p "$(dirname "$1")"
  awk -v min_gev="$2" -v max_gev="$3" -v bin_width_gev="$4" '
    BEGIN {
      print "# energy_GeV flat_weight"
      n = int(((max_gev - min_gev) / bin_width_gev) + 0.5)
      for (i = 0; i < n; ++i) {
        e = min_gev + (i + 0.5) * bin_width_gev
        printf "%.7f 1.000000000000e+00\n", e
      }
    }' > "$1"
}

ana_flux_integral() {
  awk -v min_gev="$2" -v max_gev="$3" '
    /^[[:space:]]*#/ || NF < 2 { next }
    $1 >= min_gev && $1 < max_gev { sum += $2 }
    END { printf "%.12g\n", sum }' "$1"
}

ana_skim_final_state_hyperon() {
  ana_check_cmds root
  root -l -b -q "$1(\"$2\",\"$3\")"
  [ -s "$3" ] || ana_die "skim did not write output: $3"
}
