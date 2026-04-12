Strategy for flat analysis ntuples
==================================

Goal
----

Create one common flat ROOT ntuple shape for GENIE, NuWro, and GiBUU so
the `ana` fixed-envelope analysis can compare the same final-state
hyperon sample across generators.

The sample definition is:

  keep event if final-state PDG contains any of
  {3122, 3212, 3322, 3312, 3334}

corresponding to Lambda, Sigma0, Xi0, Xi-, or Omega-.  This selection has
no momentum, angle, visibility, or topology requirement.

Generate the common samples with a flat neutrino-energy proposal over
0-10 GeV.  Use the local NuMI flux bundle in `example/numi/flux` only as
the target reweighting flux.  The ROOT file
`microboone_numi_flux_5mev.root` is the target for ROOT-based reweighting
and analysis macros.  The four `gibuu_numi_*.dat` two-column files are the
text-format target fluxes.  The file and histogram mapping is recorded in
`config/numi_flux_manifest.tsv`; the flat proposal definition is recorded
in `config/flat_enu_proposal.tsv`.

Make two ntuple tiers
--------------------

Tier 1: NUISANCE proxy flats

  Purpose:
    immediate cross-generator validation using the tools already in this
    repository.

  Common tree:
    `FlatTree_VARS` from `nuisflat`.

  Required branches for the existing proxy macros:
    `Enu_true`, `nfsp`, `pdg`, `px`, `py`, `pz`, `E`,
    `nvertp`, `pdg_vert`, `Weight` or `weight`,
    `fScaleFactor` or `scale_factor`, and mode/lepton metadata where
    available.

  Limitation:
    this tier is not the final detector-envelope input.  It has final
    particles and some primary-vertex content, but it does not reliably
    preserve nuclear-exit species, Geant4 post-exit decay ancestry,
    Lambda decay coordinates, or detector-visible topology flags.

Tier 2: `ana` reduced flats

  Purpose:
    final fixed-envelope projection and conditional composition maps.

  Common tree:
    also use `FlatTree_VARS`, but add a documented branch contract rather
    than relying only on NUISANCE defaults.

  Required additions:
    Stage B nuclear-exit arrays such as `n_exitp`, `pdg_exit`,
    `px_exit`, `py_exit`, `pz_exit`, `E_exit`; Stage C visible-Lambda
    ancestry and decay coordinates; Stage D truth-visible topology flags;
    and analysis metadata for generator, tune, flux, target, seed, and
    original generator weight.

Generator production flow
-------------------------

GENIE:

  1. Generate GHEP with `gevgen`, using `-e 0,10` and no flux file in the
     generation command.  The MicroBooNE NuMI flux is applied later as a
     reweighting target.
  2. Keep the GHEP file as the authoritative native record.  It is the
     place to extract particle history, status codes, primary strange
     ancestry, and any pre-final transport labels before conversion loses
     information.
  3. Make the proxy flat with the existing path:

       `gntpc -f gst`
       `PrepareGENIE`
       `nuisflat -i GENIE:<gprep.root> -o <genie>.flat.root`

  4. Build the Tier 2 reduced flat from GHEP plus any detector/transport
     augmentation product.  Do not try to infer nuclear exit from the
     final-state `pdg` array.

NuWro:

  1. Generate native NuWro output from an argon card with the appended
     flat 0-10 GeV beam histogram.  The current
     `example/NuWroCard_CC_Ar_numu.txt` is a starting point, but its beam
     block must be overridden by the flat proposal.
  2. Keep the native NuWro output as the authoritative record for primary
     and FSI information.
  3. Make the proxy flat with the existing path:

       `nuwro -i <card> -o <nuwro>.prep.root`
       `PrepareNuWroEvents -f <nuwro>.prep.root -o <nuwro>.nuisance.root`
       `nuisflat -i NuWro:<nuwro>.nuisance.root -o <nuwro>.flat.root`

  4. Build the Tier 2 reduced flat from the native NuWro record plus any
     detector/transport augmentation product, preserving primary species,
     nuclear-exit species when available, and post-exit Lambda ancestry.

GiBUU:

  1. Generate GiBUU ROOT event output with `WritePerturbativeParticles=T`
     and `EventFormat=4`.
  2. Use a flat 0-10 GeV proposal consistently.  The current
     `example/GiBUU2025_numu.job` is configured as an example card, not a
     frozen production card; use `nuExp=99` with the generated
     `flat_0_10gev_5mev.dat` proposal file, then reweight to the NuMI
     target flux in analysis.
  3. Keep `EventOutput.Pert.*.root` as the authoritative native record.
  4. Make the proxy flat with the existing path:

       `PrepareGiBUU -i EventOutput.Pert.<run>.root -f <flux.root>,<hist> -o <gibuu>.prep.root`
       `nuisflat -i GiBUU:<gibuu>.prep.root -o <gibuu>.flat.root`

  5. Build the Tier 2 reduced flat from the GiBUU event output plus any
     detector/transport augmentation product.  GiBUU includes hadronic
     transport through the nucleus, but it is still not a substitute for
     detector-material Geant4 decay ancestry and topology flags.

Common reduced-ntuple contract
------------------------------

Keep the NUISANCE branch names where they already work, then add the
missing `ana` branches.  At minimum:

  Event identity:
    `run`, `subrun`, `event`, `generator`, `tune`, `sample_label`,
    `seed`, `beam_mode`, `target_pdg`, `target_a`, `target_z`.

  Weights and flux:
    `Weight`, `fScaleFactor`, `proposal_id`, `proposal_emin_gev`,
    `proposal_emax_gev`, `reweight_flux_file`, `reweight_flux_hist`, `Enu_true`,
    `pdg_nu`, `cc`, `Mode`, `PDGLep`, `ELep`.

  Final-state particles:
    `nfsp`, `pdg`, `px`, `py`, `pz`, `E`.

  Primary particles:
    `nvertp`, `pdg_vert`, and primary four-vectors if available.

  Final-state sample flags:
    `ana_hyperon_final_state_flag`,
    `ana_hyperon_final_state_pdg_mask`,
    `ana_hyperon_final_state_n`.

  Stage B nuclear exit:
    `n_exitp`, `pdg_exit`, `px_exit`, `py_exit`, `pz_exit`, `E_exit`,
    plus an `exit_source_code` or metadata string explaining how the exit
    table was produced.

  Stage C visible Lambda candidate:
    `post_exit_decay_source`, `exit_seed_pdg`, `lambda_decay_x`,
    `lambda_decay_y`, `lambda_decay_z`, `lambda_decay_length`,
    `daughter_proton_px`, `daughter_proton_py`, `daughter_proton_pz`,
    `daughter_proton_E`, `daughter_piminus_px`, `daughter_piminus_py`,
    `daughter_piminus_pz`, `daughter_piminus_E`,
    `daughter_opening_angle`, `extra_visible_energy`.

  Stage D topology flags:
    `fiducial_vertex_flag`, `lambda_decays_in_active_flag`,
    `proton_visible_flag`, `pion_visible_flag`, `detached_vertex_flag`,
    `visible_muon_flag`, `visible_kaon_flag`, `visible_gamma_flag`,
    `low_extra_activity_flag`, `topology_loose_flag`,
    `topology_nominal_flag`, `topology_tight_flag`.

Implementation order
--------------------

1. Freeze common run metadata:
   generator version, tune/model, target, beam polarity, flat proposal
   range, reweight flux file, reweight flux histogram, event count, random
   seed policy, and output naming.

2. Make small Tier 1 proxy flats for all three generators using the
   existing scripts after replacing their beam settings with the flat
   0-10 GeV proposal.

3. Validate the proxy flat structure with the existing topology-yield
   macro.  If a final-state hyperon flag has already been added, require
   it to agree with a direct scan of `pdg`; otherwise use the direct scan
   as the Tier 1 selection check.

4. Write one adapter per generator that reads the authoritative native
   record and writes the common Tier 2 `FlatTree_VARS` contract.

5. Add detector/transport augmentation before the adapter drops ancestry:
   either consume a MicroBooNE/LArSoft output product or a dedicated
   truth-augmentation product that records the Stage C/D Lambda fields.

6. Re-run the validation macros on the Tier 2 flats.  The
   `analysis/nuclear_exit.cxx` macro should find explicit `n_exitp` and
   `pdg_exit` branches.  If those branches are absent, the macro refuses
   to use final-state PDGs as an exit proxy, and the ntuple is not
   acceptable for the fixed-envelope analysis.

7. Only after the small Tier 2 files validate, scale the generator
   production and create final-state hyperon skims.

Suggested output layout
-----------------------

Use stable paths so manifests can be generated mechanically:

  `ana/output/proxy_flat/<generator>/<sample>.flat.root`
  `ana/output/reduced_flat/<generator>/<sample>.analysis.flat.root`
  `ana/output/manifests/proxy_flat_manifest.tsv`
  `ana/output/manifests/reduced_flat_manifest.tsv`

Manifest columns:

  `input_flat_root<TAB>generator<TAB>tune_or_knob<TAB>sample_label<TAB>beam_mode<TAB>tier`

Validation checks
-----------------

For every file:

  * `FlatTree_VARS` exists.
  * The event count, sum of weights, generator, tune, target, proposal
    range, target reweight flux file, and target reweight flux histogram
    are written to metadata.
  * The final-state hyperon selection is an OR over PDGs
    {3122, 3212, 3322, 3312, 3334}, with no momentum cut.
  * The `Enu_true` distribution is flat over 0-10 GeV before NuMI
    reweighting.
  * The final-state PDG spectrum is nonempty and includes the expected
    hyperon categories for the chosen sample.
  * Tier 2 only: `n_exitp/pdg_exit`, `post_exit_decay_source`, Lambda
    decay coordinates, daughter four-vectors, and topology flags exist.
  * Tier 2 only: `post_exit_decay_source` agrees with
    `config/post_exit_lambda_feeddown.tsv` for Lambda, Sigma0, Xi0, Xi-,
    and Omega- seeds.

Near-term practical path
------------------------

The fastest useful milestone is not a full detector-envelope production.
It is:

  1. make NuMI/argon Tier 1 proxy flats for GENIE, NuWro, and GiBUU;
  2. skim them by the final-state hyperon rule;
  3. run the existing topology-yield and primary-mechanism checks;
  4. inspect which native records preserve enough information for the
     Tier 2 adapter;
  5. implement the Tier 2 adapter first for the generator whose native
     record gives the cleanest primary-to-exit ancestry.

This gives immediate generator-comparison plots while making clear which
files are proxy inputs and which are valid fixed-envelope analysis inputs.
