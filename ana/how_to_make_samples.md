How to make samples for this analysis
=====================================

There are three different "sample" objects here.  Keep them separate.

1. Generator truth samples
--------------------------

These are ordinary neutrino-generator outputs.  In this repository, the
existing starting points are:

  * `example/run_genie.sh` for local GENIE production;
  * `example/run_nuwro.sh` for local NuWro production;
  * `production/submit_genie.sh` plus `production/genie_grid.sh` for grid
    GENIE production.

For `ana`, define the generator truth sample by final-state content:
keep an event if its generator final-state particle list contains at
least one of Lambda (3122), Sigma0 (3212), Xi0 (3322), Xi- (3312), or
Omega- (3334).  This selection has no momentum constraint: do not impose
minimum momentum, angular acceptance, visibility, or Lambda-topology cuts
while making the generator sample.

Those scripts can produce GHEP/GST/prepared/flat files, for example a
NUISANCE `FlatTree_VARS` file.  That is enough for a first proxy analysis,
but it is not enough for the fixed-envelope analysis unless the conversion
also preserves:

  * Stage A primary strange ancestry;
  * Stage B nuclear-exit species;
  * Stage C Geant4 post-exit decay ancestry and Lambda decay coordinates;
  * Stage D truth-visible topology flags.

So the practical chain is:

  generator event file
    -> detector/transport or truth-augmentation pass
    -> reduced Stage A/B/C/D ntuple
    -> fixed-envelope projection skim

Concrete repo entry points:

  local GENIE:

    cd example
    ./run_genie.sh

  local NuWro:

    cd example
    ./run_nuwro.sh

  grid GENIE:

    cd production
    ./submit_genie.sh NUM_JOBS JOB_NAME args.txt SPLINES_FILE FLUX_FILE FLUX_HIST_NAME OUTPUT_DIRECTORY GIT_CHECKOUT

The GENIE grid script currently calls `gevgen`, then `gntpc`.  The local
GENIE script also calls `PrepareGENIE` and `nuisflat`.  These are ordinary
generator-production steps.  The fixed-envelope analysis needs one more
stage after that: a reducer that records the detector-medium Lambda decay
candidate and its Stage A/B/C/D ancestry.

2. Topology-gun envelope samples
--------------------------------

These are not neutrino-generator samples.  They are detector-acceptance
throws used to measure A_det(y).

The minimal topology-gun sample throws final Lambda -> p pi- candidates
across the fixed y space:

  throw y = (z_decay, r_decay, p_p, p_pi, theta_p_pi,
             L_sep, E_extra_vis, mu/K/gamma flags)
    -> run detector geometry / truth-visible selection
    -> fill N_throw(y_bin) and N_pass(y_bin)
    -> write A_det(y_bin) and Omega_det

An augmented topology-gun sample can instead throw post-exit parents:

  throw post-exit parent species and kinematics
    -> Geant4 transport/decay
    -> find resulting Lambda -> p pi- candidate
    -> map final visible Lambda decay into the same y_bin
    -> fill N_throw and N_pass

Use the augmented parent-gun version when parent-decay companions, such as
photons or pions, affect the truth-visible topology or extra-visible-energy
axis.

3. Fixed-envelope generator skims
---------------------------------

These are derived analysis samples.  They are made only after A_det(y) and
Omega_det are frozen.

For each event in a reduced Stage A/B/C/D generator ntuple:

  * compute y_bin from the final visible Lambda candidate;
  * look up A_det(y_bin);
  * mark `in_omega_det`;
  * keep the original generator event weight;
  * copy ancestry labels:
      origin_class,
      mechanism_class,
      primary_species,
      exit_species,
      post_exit_decay_source,
      charm_ancestry_flag;
  * copy visible category flags:
      inclusive topology,
      CC muon,
      visible kaon,
      visible gamma,
      low extra activity,
      no visible muon.

This creates two useful tables:

  all_projectable_events
      every event with a valid y_bin and A_det lookup;

  accepted_envelope_events
      the subset with y_bin in Omega_det.

Use `accepted_envelope_events` for conditional maps such as
P(origin | y_bin, Omega_det).  Use `all_projectable_events` for acceptance
and lost-space plots.

What is possible in this repository now?
----------------------------------------

Possible now:

  * run or adapt the existing GENIE/NuWro scripts to create ordinary
    generator samples, then skim them with the `ana` final-state rule in
    `config/generator_final_state_selection.tsv`;
  * run the older proxy macros over NUISANCE FlatTree files;
  * run the manifest-driven generator matrix in
    `config/generator_loop_plan.tsv`, then make comparison plots with
    `../analysis/plot_generator_matrix.cxx`;
  * make the first expected-yield study with
    `../analysis/run_expected_event_yields.sh`, which reports stat
    uncertainties from `sqrt(sum weights^2)` and systematic uncertainties as a
    variation envelope over the generator matrix plus generator-envelope yield
    summaries, then adds fast detector-envelope rows corrected for
    Lambda -> p pi- branching and daughter momentum thresholds;
  * test phenomenological event-yield sensitivity with
    `../analysis/run_phenomenology_sensitivity.sh`, which reports yield shifts
    from favouring each generator variation and leave-one-out envelope
    contractions from excluding each variation;
  * make the FHC/RHC kinematic beam envelope with
    `../analysis/run_beam_kinematic_envelope.sh`, which compares the NuMI
    reweighted `enu`, `q2`, `w`, and leading-hyperon-momentum spectra for
    the same generator/variation/species/FSI group;
  * make the generator/configuration envelope with
    `../analysis/run_generator_kinematic_envelope.sh`, which varies generator,
    version, tune/configuration, knob, and optionally FSI state at fixed beam
    mode/species and interaction;
  * make old-jobcard-style kinematic coverage contours with
    `../analysis/plot_kinematic_coverage.cxx` or the E_nu-vs-Q2 shortcut
    `../analysis/plot_enu_q2_coverage.cxx`;
  * scan coverage through detector-visibility assumptions with
    `../analysis/run_detector_visibility_coverage.sh`, which loops the raw
    final-hyperon, branching-only, and loose/nominal/tight visible-daughter
    threshold assumptions for each requested beam/generator/variation envelope;
  * define the fixed y-space and envelope metadata.

Not yet implemented here:

  * a topology-gun producer;
  * a reduced Stage A/B/C/D ntuple maker that includes Geant4 post-exit
    decay ancestry;
  * the envelope builder that writes A_det(y);
  * the projection macro that creates fixed-envelope generator skims.

That means the next implementation target after proxy matrix plots should
be the reduced ntuple contract and projection code.  See
`flat_ntuple_strategy.md` for the staged GENIE, NuWro, and GiBUU
flat-ntuple production plan.

Geant4 status in this repository
--------------------------------

Geant4 is not currently being run by the sample scripts in this repository.

What exists:

  * `global_vars.sh` sets up a Geant4 product;
  * `build_genie.sh` configures GENIE with `--enable-geant4`.

What does not exist here yet:

  * a Geant4 application or LArSoft detector-simulation job;
  * a MicroBooNE geometry-driven topology-gun runner;
  * a post-exit transport/decay pass that writes Geant4 trajectory ancestry
    into the reduced ntuple.

The current production scripts run generator and conversion tools such as
`gevgen`, `gntpc`, `PrepareGENIE`, `nuwro`, `PrepareNuWroEvents`, and
`nuisflat`.  Those do not, by themselves, provide the Stage C Geant4
post-exit decay ancestry needed by the fixed-envelope analysis.

So the Geant4 step must be added as a new layer or supplied by an external
MicroBooNE/LArSoft production chain before the final fixed-envelope
projection is considered complete.

Recommended near-term production order
--------------------------------------

1. Make a small local generator sample with the existing scripts.
2. Inspect the resulting flat tree and confirm which ancestry fields are
   present.
3. Add or run a reducer that writes the Stage A/B/C/D columns listed
   below.  If Geant4 trajectory ancestry is not available in the flat tree,
   this reducer must operate before information is discarded.
4. Build a small topology-gun envelope over a coarse y grid.
5. Project the small reduced generator sample through the coarse envelope.
6. Only after the table shapes and labels validate, scale to grid
   production and targeted oversamples.

Minimal next implementation target
----------------------------------

Create a reduced ntuple with one row per visible Lambda candidate, with at
least these columns:

  run, subrun, event
  generator, sample_label, beam_weight
  y_bin, z_decay, r_decay, p_p, p_pi_minus, theta_p_pi, L_sep, E_extra_vis
  visible_muon_flag, visible_kaon_flag, visible_gamma_flag
  topology_loose_flag, topology_nominal_flag, topology_tight_flag
  origin_class, mechanism_class
  primary_species, exit_species, post_exit_decay_source
  charm_ancestry_flag

Then the envelope projection can be a deterministic table operation:

  reduced_ntuple + frozen_A_det(y_bin) -> accepted_envelope_events
