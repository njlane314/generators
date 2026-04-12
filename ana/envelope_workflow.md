Detector envelope workflow
==========================

Definitions
-----------

Use one fixed detector-visible coordinate system:

  y = (z_decay_cm,
       r_decay_cm,
       p_proton_GeV,
       p_pi_minus_GeV,
       theta_p_pi_rad,
       lambda_sep_cm,
       extra_visible_energy_GeV,
       visible_muon_flag,
       visible_kaon_flag,
       visible_gamma_flag)

For each bin in this space, the topology-gun layer defines:

  A_det(y_bin) = N_pass(y_bin) / N_throw(y_bin)

The detector envelope is then:

  Omega_det = { y_bin with enough throws and an accepted topology definition }

The threshold for "enough throws" should be stored with the envelope
artifact.  Do not silently change this threshold between generator
comparisons.

Topology-gun layer
------------------

The topology-gun layer should define detector acceptance with minimal
neutrino-model assumptions.  The simplest version throws detached
Lambda -> p pi- truth decays through the detector geometry.  An augmented
version can throw post-exit parent particles, let Geant4 transport and
decay them, and then map any resulting Lambda -> p pi- candidate into the
same y space.

At minimum, each throw needs:

  * throw type, for example final_Lambda_gun or post_exit_parent_gun;
  * post-exit parent PDG if the throw starts before the final Lambda;
  * Lambda decay position;
  * proton and pi- four-vectors;
  * Lambda direction or pointing proxy;
  * visible companion flags or visible companion four-vectors if the
    topology selection depends on parent-decay products;
  * generated y_bin;
  * pass/fail flags for loose, nominal, and tight truth-visible topology;
  * the proposal density or sampling label if the throws are not flat.

If importance sampling is used, store the proposal weight separately from
the detector pass/fail flag.  The envelope numerator and denominator must
be formed with the same proposal correction.

If the first envelope is a final_Lambda_gun envelope, it is a core
Lambda-topology acceptance map.  It can still be used to project generator
feed-down events because the y_bin is defined from the final visible
Lambda decay.  If parent-decay companions such as photons or pions affect
the pass/fail topology rule, build an augmented post_exit_parent_gun
envelope for those axes rather than silently folding the effect into a
mechanism label.

Fast first-study envelope
-------------------------

The first expected-yield macro implements a narrower truth-level
detector-appearance envelope:

  final-state hyperon -> Lambda reachability
                      -> Lambda -> p pi- branching ratio
                      -> proton and pi- momentum thresholds

This is useful for the first yield table because it corrects the
generator final-state hyperon yield toward the visible Lambda channel
without changing the generator sample definition.  It is not the full
fixed detector envelope.  It has no geometry, fiducial volume,
containment, detached-vertex, reconstruction, or y-bin acceptance map.
Replace or multiply it by the frozen A_det(y_bin) object once the
topology-gun envelope and reduced Stage C/D ntuples exist.

The coverage scan in `../analysis/run_detector_visibility_coverage.sh`
uses the same fast first-study model to show how phase-space coverage
changes across detector assumptions.  The default scan compares:

  * raw generator final-state hyperon coverage;
  * Lambda reachability times B(Lambda -> p pi-) only;
  * loose, nominal, and tight visible daughter momentum thresholds.

These coverage plots are diagnostic plots for the current generator
matrix.  They are not a replacement for the frozen A_det(y_bin) envelope,
but they make the detector-visible assumption dependence explicit before
the topology-gun and reduced Stage C/D products exist.

The phenomenological yield-sensitivity pass is similarly a diagnostic:
given the current detector-visible and beam-reweighted yield table, it
asks whether favouring or excluding a generator variation moves the total
event yield by more than the configured statistical or fractional scale.
It should be interpreted as a generator-model sensitivity screen, not as a
claim about reconstruction efficiency.

Generator projection layer
--------------------------

The input generator samples are defined by final-state content before the
detector-envelope projection.  Keep an event if the generator final-state
particle list contains at least one of Lambda (3122), Sigma0 (3212), Xi0
(3322), Xi- (3312), or Omega- (3334).  This selection has no momentum
threshold or angular requirement.

After the envelope is frozen, each generator event should be mapped into
the same y_bin if it contains a post-transport Lambda decay candidate.
For accepted bins, fill conditional quantities with generator weights:

  f_o(y_bin) =
      sum_events w * I(origin == o) * I(y in y_bin) * I(y_bin in Omega_det)
      --------------------------------------------------------------------
      sum_events w * I(y in y_bin) * I(y_bin in Omega_det)

  f_m(y_bin) =
      sum_events w * I(mechanism == m) * I(y in y_bin) * I(y_bin in Omega_det)
      -----------------------------------------------------------------------
      sum_events w * I(y in y_bin) * I(y_bin in Omega_det)

  M_a_to_b(y_bin) =
      sum_events w * I(primary == a) * I(exit == b) * I(y in y_bin) * I(y_bin in Omega_det)
      ------------------------------------------------------------------------------------
      sum_events w * I(primary == a) * I(y in y_bin) * I(y_bin in Omega_det)

The post-exit Geant4 decay step gets its own migration term.  Do not fold
it into Stage B:

  G_b_to_c(y_bin) =
      sum_events w * I(exit == b) * I(post_exit_decay_source == c) * I(y in y_bin) * I(y_bin in Omega_det)
      --------------------------------------------------------------------------------------------------
      sum_events w * I(exit == b) * I(y in y_bin) * I(y_bin in Omega_det)

Charm feed-down is handled as another ancestry flag:

  f_charm(y_bin) =
      sum_events w * I(charm_ancestry) * I(y in y_bin) * I(y_bin in Omega_det)
      -----------------------------------------------------------------------
      sum_events w * I(y in y_bin) * I(y_bin in Omega_det)

Required reduced-ntuple content
-------------------------------

The fixed-envelope analysis needs more than bare final-state PDGs.
Do not treat final-state PDG as a substitute for Stage B nuclear exit.

Stage A, primary:

  * primary strange baryon species before prompt decay;
  * primary kaon content;
  * origin class O0..O7;
  * interaction mode, E_nu, W, Q2 when available.

Stage B, nuclear exit:

  * strange baryon species at nuclear exit;
  * kaon species at nuclear exit;
  * Lambda four-vector at exit or before decay.

Stage C, post-decay detector medium:

  * Geant4 particle trajectory or ancestry id for the visible Lambda
    candidate;
  * exit particle that seeded the visible Lambda candidate;
  * post-exit decay source class, for example exit_Lambda_survival,
    Sigma0_to_Lambda_gamma, Xi_to_Lambda_pi, Omega_to_Lambda_X,
    charm_to_Lambda_X, other_post_exit_feeddown;
  * decay process name or decay-chain code when available;
  * Lambda decay vertex x, y, z;
  * Lambda decay length or separation from primary vertex;
  * daughter proton four-vector;
  * daughter pi- four-vector;
  * daughter opening angle;
  * extra visible energy estimate.

Stage D, truth-visible tags:

  * fiducial_vertex_flag;
  * lambda_decays_in_active_flag;
  * proton_visible_flag;
  * pion_visible_flag;
  * detached_vertex_flag;
  * visible_muon_flag;
  * visible_kaon_flag;
  * visible_gamma_flag;
  * low_extra_activity_flag;
  * topology_loose_flag;
  * topology_nominal_flag;
  * topology_tight_flag.

Post-exit decay handling
------------------------

The visible Lambda topology can be seeded by a particle that is present at
nuclear exit but is not itself the Lambda that decays to p pi-.  Geant4 may
propagate and decay that particle in the detector medium before the final
Lambda decay.  The reduced ntuple therefore needs a bridge from Stage B to
Stage C:

  primary species -> nuclear-exit species -> post-exit decay source
                  -> visible Lambda candidate -> p pi- topology

For example, keep these as different cases:

  * primary Lambda -> exit Lambda -> Lambda -> p pi-;
  * primary Sigma0 -> exit Sigma0 or prompt Lambda gamma -> Lambda -> p pi-;
  * primary Xi -> exit Xi -> Lambda pi -> Lambda -> p pi-;
  * charm hadron -> Lambda X -> Lambda -> p pi-.

This is still a fixed-envelope analysis.  The y_bin is defined from the
detector-visible Lambda decay and companions.  The post-exit decay source
is a conditional label painted onto that same bin, not a new signal
definition.

Envelope artifact metadata
--------------------------

Every frozen envelope artifact should record:

  * y-space config path and checksum;
  * MicroBooNE geometry tag;
  * Geant4 physics-list and decay-handling tag;
  * truth-visible working point;
  * topology-gun generation range;
  * number of throws per populated bin;
  * minimum-bin-statistics rule;
  * any bin merging or smoothing rule;
  * software commit or local git status summary.
