Feed-down sampling strategy
===========================

Goal
----

Use decay reachability to decide which generator phase-space regions need
extra statistics for the fixed-envelope Lambda topology analysis.

This is narrower than a full detector simulation.  The first question is:

  Which particles present at nuclear exit can later feed a visible
  Lambda -> p pi- topology?

Once that is known, the second question is:

  Which generator regions produce those nuclear-exit parent particles?

Three maps
----------

Use three separate maps:

  R_A_to_B:
    generator transport map from primary species to nuclear-exit species.

  D_B_to_C:
    decay reachability map from nuclear-exit species to a final Lambda
    candidate after post-exit decays.

  A_det:
    detector-visible acceptance map in fixed y bins.

The oversampling target is not just "make more Lambda".  It is:

  generator final-state samples that contain at least one Lambda, Sigma0,
  Xi0, Xi-, or Omega- with no momentum constraint, especially regions that
  populate B species with nonzero D_B_to_C and weak accepted-envelope
  statistics.

Post-exit feed-down classes
---------------------------

Start with these classes:

  exit_Lambda_survival
    Nuclear-exit Lambda is the visible Lambda candidate.

  exit_Sigma0_to_Lambda_gamma
    Nuclear-exit Sigma0 feeds Lambda gamma.  This is a direct post-exit
    feed-down class and the gamma can also affect visible companion axes.

  exit_Xi_to_Lambda_pi
    Nuclear-exit Xi- or Xi0 feeds Lambda pi.  This can populate the same
    p pi- Lambda topology but with an extra pion companion.

  exit_Omega_to_Lambda_X
    Nuclear-exit Omega can feed Lambda directly or through Xi chains.

  charm_to_Lambda_X
    Charmed baryon or charm event ancestry feeds Lambda plus hadrons.
    This should be treated as an ancestry flag and usually needs a
    charm-enhanced or high-W oversample.

  other_hyperon_resonance_to_Lambda_X
    Higher hyperon resonances, if present in the generator/Geant4 decay
    table, that decay to Lambda plus mesons.

  no_true_Lambda_decay_feeddown
    Stable or weakly decaying particles that do not produce a true Lambda
    by decay.  Keep these out of the true-Lambda signal; later they may be
    studied as fake-topology backgrounds.

Important separation
--------------------

Charged Sigma production is still important, but usually as a Stage A to
Stage B transport/exchange question:

  primary charged Sigma -> nuclear-exit Lambda

Do not automatically treat a nuclear-exit charged Sigma as a post-exit
Lambda feed-down source.  The ordinary weak decays of charged Sigmas do
not make Lambda -> p pi-.

How to use Geant4 source here
-----------------------------

Use Geant4 decay tables to build D_B_to_C:

  parent PDG -> allowed decay channels -> daughter PDGs -> Lambda reachability

The relevant classes are documented in `geant4_decay_source_notes.md`.
For each candidate nuclear-exit parent:

  1. Select decay channels from the Geant4/PDG table.
  2. Walk daughter chains until a Lambda is found or the chain terminates.
  3. Record direct Lambda feed-down, indirect Lambda feed-down, and visible
     companion particles.
  4. Store a class label such as `exit_Xi_to_Lambda_pi`.

This is enough to decide which parent species and companion axes matter.
It is not enough to finalize A_det(y) without geometry validation.

Generator oversampling rule
---------------------------

For each generator sample, measure or estimate:

  N_gen(x_gen, primary_species)
  R_A_to_B(primary_species -> exit_species, x_gen)
  D_B_to_C(exit_species -> post_exit_decay_source)
  A_det(y_bin)

Then prioritize new generator production where:

  * exit_species has nonzero Lambda reachability;
  * accepted y bins have large statistical uncertainty;
  * the component is physically important or could dominate a small
    accepted corner;
  * inclusive production under-samples the region.

Concrete regions likely needing attention
-----------------------------------------

Direct Lambda:
  Sample ordinary Lambda-bearing events, including near-threshold and
  forward Lambda regions.

Sigma0 feed-through:
  Sample primary or exit Sigma0 and retain the gamma companion.  Do not
  use a no-shower topology definition as the headline signal.
  If the reducer records Sigma0 as already prompt-decayed before the
  nuclear-exit table, keep the Stage A Sigma0 origin label and the
  Lambda-gamma companion record rather than requiring an exit Sigma0 row.

Charged-Sigma exchange:
  Sample primary charged Sigma production and inspect whether transport
  produces nuclear-exit Lambda.  This is mostly R_A_to_B, not D_B_to_C.

Xi/Omega feed-down:
  Sample higher-strangeness and higher-multiplicity regions, especially
  where exit Xi or Omega appears.  Keep extra pion/kaon companions.

Associated K/Y:
  Sample events with primary kaon plus hyperon content and retain visible
  kaon companions.

Charm feed-down:
  Use charm-enhanced or high-W DIS samples.  Keep charm ancestry separate
  from the signal definition and combine with documented weights.

Output table
------------

The first useful output is a compact table:

  exit_species
  post_exit_decay_source
  lambda_reachable
  direct_or_indirect
  companion_species
  recommended_generator_region
  oversample_priority

This table becomes the bridge from Geant4 decay-source knowledge to the
generator production plan.

The current seed version of this table is:

  config/post_exit_lambda_feeddown.tsv
