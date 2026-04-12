MicroBooNE NuMI Flux Files
==========================

These files were imported from `../old.generators/jobcards/numi/flux`.

Generate the common analysis samples directly from the NuMI flux by default.
Use `microboone_numi_flux_5mev.root` as the ROOT-format flux for GENIE
generation and GiBUU preparation.  The central-value histograms are:

  * `fhc/numu/Detsmear/numu_CV_AV_TPC_5MeV_bin`
  * `fhc/numubar/Detsmear/numubar_CV_AV_TPC_5MeV_bin`
  * `rhc/numu/Detsmear/numu_CV_AV_TPC_5MeV_bin`
  * `rhc/numubar/Detsmear/numubar_CV_AV_TPC_5MeV_bin`

Use the `gibuu_numi_*.dat` files when a text-format target NuMI flux is
needed for NuWro or GiBUU generation.  They are two-column tables:

  `energy_GeV flux_per_5MeV_bin_per_1e6POT`

The active generator-loop manifest is `analysis/config/generator_loop_plan.tsv`;
the flux path mapping used by the local sample runners is in
`analysis/sample_common.sh`.

The GiBUU and NuWro analysis runners now read these NuMI text fluxes directly.
