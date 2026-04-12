Fixed-envelope Lambda analysis
==============================

This directory is the restart point for the detector-facing Lambda study.
It is intentionally separate from `analysis/`, whose current macros are
useful references but are organized around mechanism summaries rather than
a frozen detector-visible coordinate system.

Core rule
---------

Define the detector-visible coordinate system first:

  y = (z_decay, r_decay, p_p, p_pi_minus, theta_p_pi,
       L_sep, E_extra_vis, visible_muon, visible_kaon, visible_gamma)

Then define a fixed detector envelope:

  A_det(y_bin) = N_pass_topology_gun(y_bin) / N_throw_topology_gun(y_bin)
  Omega_det   = accepted detector-visible bins in this fixed y space

Only after that should generator samples be projected through the same
envelope to measure conditional composition:

  f_origin(o | y_bin, Omega_det)
  f_mechanism(m | y_bin, Omega_det)
  M(primary_species -> exit_species | y_bin, Omega_det)
  G(exit_species -> post_exit_decay_source | y_bin, Omega_det)
  f_charm(y_bin | Omega_det)

The nuclear-exit state is not the last truth state.  Unstable particles
that leave the nucleus can still decay in detector material through
Geant4 and feed the visible Lambda topology.  Those post-exit decays must
be preserved as Stage C ancestry, then evaluated inside the same fixed
detector envelope.

Generator sample definition
---------------------------

The `ana` generator samples are final-state hyperon samples.  Select an
event if the generator final-state particle list contains at least one of:

  * Lambda, PDG 3122;
  * Sigma0, PDG 3212;
  * Xi0, PDG 3322;
  * Xi-, PDG 3312;
  * Omega-, PDG 3334.

This is an inclusive OR over final-state PDGs.  Do not apply a momentum
threshold, angular cut, detector visibility cut, or Lambda-topology
requirement when defining these generator samples.  The fixed detector
envelope is applied later as a projection step.

The machine-readable seed for this rule is
`config/generator_final_state_selection.tsv`.

What the detector envelope can create
-------------------------------------

The detector envelope can create detector-facing samples, not new physics.

It can be used to create:

  * topology-gun samples: flat or importance-sampled detached
    Lambda -> p pi- truth throws through detector geometry;
  * accepted generator skims: generator events tagged by y_bin and by
    whether they fall in Omega_det;
  * targeted oversample requests: a list of y bins or physics corners
    where the accepted generator statistics are weak;
  * plotting inputs: common-binning maps for acceptance, composition,
    origin migration, and companion-object fractions.

It should not be used to invent generator mechanisms or replace generator
weights.  A_det(y) is a geometry/topology acceptance object.  Direct Lambda,
Sigma0 feed-through, charged-Sigma exchange, associated K production, and
post-exit Geant4 decay feed-down remain conditional labels inside the
accepted envelope.

Recommended workflow
--------------------

1. Freeze `config/y_space.tsv`.
2. Produce a topology-gun ntuple with one row per thrown Lambda decay.
3. Fill `N_throw(y_bin)` and `N_pass(y_bin)` from that topology-gun ntuple.
4. Freeze an envelope artifact, for example
   `ana/output/envelope_microboone_nominal.root`.
5. Run each neutrino-generator reduced ntuple through the frozen envelope.
6. Write generator skims with event id, generator weight, y_bin, origin,
   primary species, exit species, post-exit decay source, visible category,
   and charm flag.
7. Make fixed-envelope plots by repainting the same accepted y bins with
   different conditional quantities.

Useful starting files
---------------------

  * `config/y_space.tsv`: canonical y variables, binning, and required
    reduced-ntuple sources.
  * `envelope_workflow.md`: implementation detail and formulas.
  * `sample_plan.md`: sample families and how to use the envelope without
    biasing generator interpretation.
  * `how_to_make_samples.md`: concrete sample-making path and what is
    already possible in this repository.
  * `flat_ntuple_strategy.md`: strategy for building common GENIE, NuWro,
    and GiBUU proxy and reduced flat ntuples.
  * `config/generator_loop_plan.tsv`: manifest-driven analysis loop over
    generator, version, knob/config, beam, interaction, and FSI variations.
  * `../analysis/expected_event_yields.cxx` and
    `../analysis/run_expected_event_yields.sh`: first-study macro for expected
    event yields with statistical uncertainties, variation-envelope
    systematic uncertainties, generator-envelope yield summaries, and a fast
    detector-appearance envelope using Lambda -> p pi- branching and daughter
    momentum thresholds.
  * `../analysis/yield_sensitivity.awk` and
    `../analysis/run_phenomenology_sensitivity.sh`: yield-sensitivity
    post-processing for favouring or excluding generator mechanism,
    configuration, and transport variations.
  * `../analysis/beam_kinematic_envelope.cxx` and
    `../analysis/run_beam_kinematic_envelope.sh`: FHC/RHC kinematic beam
    envelope plots and CSVs for the flat-sample-to-NuMI reweighting.
  * `../analysis/generator_kinematic_envelope.cxx` and
    `../analysis/run_generator_kinematic_envelope.sh`: generator, version,
    configuration, knob, and optional FSI kinematic envelope plots and CSVs.
  * `../analysis/plot_generator_matrix.cxx`: ROOT plot macro for the
    manifest-driven generator, variation, and beam comparison outputs.
  * `../analysis/coverage_plotter.cxx`: ROOT coverage contour class for
    old-jobcard-style E_nu/Q2/W/Lambda phase-space coverage comparisons,
    including detector-visibility weighting assumptions.
  * `../analysis/plot_kinematic_coverage.cxx` and
    `../analysis/plot_enu_q2_coverage.cxx`: compatibility entry points for
    that coverage plotter.
  * `../analysis/run_detector_visibility_coverage.sh`: shell driver for
    coverage scans across raw final-state, branching-only, and
    visible-daughter-threshold detector assumptions.
  * `config/flat_enu_proposal.tsv`: flat 0-10 GeV neutrino-energy
    proposal used before NuMI flux reweighting.
  * `config/numi_flux_manifest.tsv`: local MicroBooNE NuMI ROOT and
    two-column flux files used by the generator examples and analysis.
  * `geant4_decay_source_notes.md`: where the Geant4 decay sampler lives
    and why it is not a full detector-envelope replacement.
  * `feeddown_sampling_strategy.md`: how to use post-exit decay
    reachability to choose generator oversample regions.
  * `config/post_exit_lambda_feeddown.tsv`: compact post-exit
    species-to-Lambda reachability table.
