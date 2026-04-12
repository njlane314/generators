Sample plan for the fixed-envelope analysis
===========================================

Sample families
---------------

Layer 0: flux and weights

  Purpose:
    define FHC, RHC, and combined beam projections for the same truth
    samples.

  Output:
    flux weights or flux-envelope metadata, not detector acceptance.

Layer 1: topology-gun detector envelope

  Purpose:
    generate detached Lambda -> p pi- truth throws across y space, or
    post-exit parent throws that Geant4 decays into the same visible
    Lambda coordinates, and define A_det(y) with minimal
    neutrino-generator bias.

  Output:
    frozen detector-envelope artifacts, pass/fail maps, and Geant4
    decay-handling metadata.

Layer 2: final-state hyperon generator samples

  Purpose:
    provide realistic generator populations selected only by final-state
    content, then projected through the same Omega_det.  The sample rule is
    an inclusive OR over final-state Lambda, Sigma0, Xi0, Xi-, or Omega-
    with no momentum constraint.

  Output:
    conditional origin, mechanism, companion-object, and migration maps.

Layer 3: mechanism-control samples

  Purpose:
    improve interpretation of direct Lambda, Sigma0 feed-through,
    charged-Sigma exchange, associated K production, and charm feed-down.

  Output:
    ancestry-calibration and systematic comparison inputs, not the main
    detector-facing signal definition.

Layer 4: targeted oversamples

  Purpose:
    fill rare accepted corners after the first fixed-envelope projection
    exposes low statistics.

  Output:
    generator-specific oversamples with explicit weights and provenance.

Can the detector envelope be used to create samples?
---------------------------------------------------

Yes, with a precise meaning.

Good uses:

  * create the topology-gun sample used to measure A_det(y);
  * create an accepted generator skim by keeping events with y_bin in
    Omega_det and preserving their original generator weights;
  * create a prioritized request list for new generator oversamples in
    accepted bins where statistical uncertainty is too high;
  * choose importance-sampling regions for extra topology-gun throws.

Bad uses:

  * do not use A_det(y) to manufacture direct-Lambda or Sigma-specific
    event rates;
  * do not use the topology-gun distribution as a replacement for the
    generator distribution;
  * do not reweight origin fractions by hand unless the weights correspond
    to a documented generation or importance-sampling proposal;
  * do not define the signal as "events from a chosen mechanism that pass
    the detector envelope."

Operational recipe
------------------

1. Build the envelope from topology-gun throws.
2. Freeze the envelope file and metadata.
3. For each generator sample, first keep events with any final-state PDG
   in {3122, 3212, 3322, 3312, 3334}; do not impose a momentum cut.
4. Compute y_bin per kept event when a detector-visible Lambda candidate
   exists after the post-exit decay/transport step.
5. Add columns:

     envelope_version
     y_bin
     in_omega_det
     A_det_y_bin
     origin_class
     mechanism_class
     primary_species
     exit_species
     post_exit_decay_source
     visible_category
     charm_ancestry_flag

6. Make two generator-derived tables:

     all_projectable_events: every event with a valid y_bin;
     accepted_envelope_events: subset with in_omega_det == true.

7. Use accepted_envelope_events for conditional composition maps.
8. Use all_projectable_events and A_det_y_bin for acceptance and lost-space
   diagnostics.

Charm feed-down
---------------

Charm should not become a separate signal definition.  Treat it as an
ancestry flag inside the same accepted envelope:

  P(charm | y_bin, Omega_det)

Because charm is rare, plan a charm-enhanced or high-W oversample if the
first inclusive projection has weak accepted statistics.  Keep the charm
oversample as an additive weighted component with explicit provenance.

Post-exit Geant4 feed-down
--------------------------

The generator projection must not stop at nuclear exit.  An exit Xi,
Omega, Sigma0, or charm hadron can decay in the detector simulation and
feed the same visible Lambda -> p pi- topology.  Preserve that chain as:

  exit_species -> post_exit_decay_source -> visible Lambda candidate

Then make conditional maps such as:

  P(post_exit_decay_source | exit_species, y_bin, Omega_det)

This keeps the detector-visible topology fixed while still separating
direct exit Lambda survival from post-exit decay feed-down.
