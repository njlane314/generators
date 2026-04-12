#!/bin/bash
set -euo pipefail

export proposal_min_gev="${proposal_min_gev:-0}"
export proposal_max_gev="${proposal_max_gev:-10}"
export proposal_bin_width_gev="${proposal_bin_width_gev:-0.005}"
export proposal_fluxfile="${proposal_fluxfile:-./numi/flux/flat_0_10gev_5mev.dat}"
export reweight_fluxfile="./numi/flux/microboone_numi_flux_5mev.root"
export reweight_fluxhisto="fhc/numu/Detsmear/numu_CV_AV_TPC_5MeV_bin"

awk -v min_gev="${proposal_min_gev}" -v max_gev="${proposal_max_gev}" -v bin_width_gev="${proposal_bin_width_gev}" 'BEGIN {
  print "# energy_GeV flat_weight"
  n = int(((max_gev - min_gev) / bin_width_gev) + 0.5)
  for (i = 0; i < n; i++) {
    e = min_gev + (i + 0.5) * bin_width_gev
    printf "%.7f 1.000000000000e+00\n", e
  }
}' > "${proposal_fluxfile}"

# Generate GiBUU events
../GiBUU/release/testRun/./GiBUU.x < GiBUU2025_numu.job

# Convert to Nuisance format
#for i in {1..9}; do PrepareGiBUU -i EventOutput.Pert.0000000${i}.root -f MCC9_FluxHist_volTPCActive.root,hEnumu_cv -o GiBUU_${i}.prep.root; done
#for i in {10..99}; do PrepareGiBUU -i EventOutput.Pert.000000${i}.root -f MCC9_FluxHist_volTPCActive.root,hEnumu_cv -o GiBUU_${i}.prep.root; done
#for i in {100..300}; do PrepareGiBUU -i EventOutput.Pert.00000${i}.root -f MCC9_FluxHist_volTPCActive.root,hEnumu_cv -o GiBUU_${i}.prep.root; done

# Keep the NuMI target flux attached at NUISANCE preparation; the GiBUU
# event proposal itself is the generated flat file above.
PrepareGiBUU -i EventOutput.Pert.00000001.root -f "${reweight_fluxfile},${reweight_fluxhisto}" -o GiBUU.prep.root

# Convert to Nuisance flat tree format
#for i in {1..300}; do nuisflat -i GiBUU:GiBUU_${i}.prep.root -o samples/GiBUU_${i}.flat.root; done
nuisflat -i GiBUU:GiBUU.prep.root -o samples/GiBUU.flat.root

#cd samples
#hadd GiBUU2023.flat.root GiBUU*.flat.root
#mv GiBUU2023.flat.root /pnfs/uboone/persistent/users/apapadop/GiBUU_Samples/GiBUU2023/GiBUU2023_300runs.root
##rm *.root
#cd ..

# Remove unnecessary files
#rm *.prep.root
rm -f *.dat
rm -f GiBUU_database_decayChannels.txt
rm -f GiBUU_database.tex
rm -f main.run
rm -f PYR.RG
rm -f EventOutput.Pert.0000*.root
rm -f *.prep.root
