#!/bin/bash
set -euo pipefail

#https://nuwro.github.io/user-guide/

export beam_particle="${beam_particle:-14}"
export proposal_min_gev="${proposal_min_gev:-0}"
export proposal_max_gev="${proposal_max_gev:-10}"
export proposal_bin_width_gev="${proposal_bin_width_gev:-0.005}"
export reweight_fluxfile="${reweight_fluxfile:-numi/flux/gibuu_numi_fhc_numu.dat}"
export outdir="${outdir:-samples}"
export run_card="${run_card:-NuWroCard_CC_Ar_numu.numi.tmp.txt}"

beam_energy="$(
  awk -v min_gev="${proposal_min_gev}" -v max_gev="${proposal_max_gev}" -v bin_width_gev="${proposal_bin_width_gev}" '
    BEGIN {
      emin = min_gev * 1000.0
      emax = max_gev * 1000.0
      step = bin_width_gev * 1000.0
      if (step <= 0 || emax <= emin) {
        exit 1
      }
      n = int(((emax - emin) / step) + 0.5)
      printf "%.12g %.12g", emin, emax
      for (i = 0; i < n; i++) {
        printf " 1"
      }
      printf "\n"
    }'
)"

if [[ ! -s "${reweight_fluxfile}" ]]; then
  echo "ERROR: missing NuMI reweight flux file: ${reweight_fluxfile}" >&2
  exit 1
fi

reweight_flux_integral="$(
  awk -v min_gev="${proposal_min_gev}" -v max_gev="${proposal_max_gev}" '
    /^[[:space:]]*#/ || NF < 2 { next }
    END {
      printf "%.12g\n", sum
    }
    $1 >= min_gev && $1 < max_gev { sum += $2 }
  ' "${reweight_fluxfile}"
)"

cp NuWroCard_CC_Ar_numu.txt "${run_card}"
cat >> "${run_card}" <<EOF

################################################################################
# Generated flat Enu proposal override
################################################################################

beam_type = 0
beam_particle = ${beam_particle}
beam_energy = ${beam_energy}

# NuMI flux reweight target, not the generation spectrum:
# reweight_flux_file = ${reweight_fluxfile}
# reweight_flux_integral_0_10gev = ${reweight_flux_integral}
EOF

mkdir -p "${outdir}"

# Generate NuWro events
nuwro -i "${run_card}" -o NuWroCard_CC_Ar_numu.prep.root

#convert to nuisance format
PrepareNuWroEvents -f NuWroCard_CC_Ar_numu.prep.root -o NuWro.prep.root

# Convert to Nuisance flat tree
nuisflat -i NuWro:NuWro.prep.root -o "${outdir}/NuWro.flat.root"

# Remove unnecessary files
rm -f NuWroCard_CC_Ar_numu.prep.root
rm -f NuWro.*
rm -f NuWroCard_CC_Ar_numu.prep.root.*
rm -f q*.txt
rm -f T.txt
rm -f random_seed
rm -f totals.txt
