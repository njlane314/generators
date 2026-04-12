Geant4 decay-source notes
=========================

Question
--------

Could we find the Geant4 source that samples particle decays and use it
directly?

Short answer
------------

Yes for decay-channel and decay-kinematics understanding.  No as a full
replacement for the detector-envelope layer.

The relevant Geant4 source is concentrated in a few classes:

  * `G4Decay`
      process wrapper that acts on a `G4Track`, selects a decay channel,
      calls the channel's `DecayIt`, boosts products to the lab frame, and
      creates secondaries.

  * `G4DecayTable`
      chooses the decay channel using branching ratios and parent-mass
      threshold checks.

  * `G4PhaseSpaceDecayChannel`
      samples one-body, two-body, three-body, and many-body phase-space
      decay kinematics.

For the local setup request `geant4 v4_11_2_p02`, the matching upstream
tag is approximately Geant4 `v11.2.2`.  Source pointers:

  * https://github.com/Geant4/geant4/blob/v11.2.2/source/processes/decay/src/G4Decay.cc
  * https://github.com/Geant4/geant4/blob/v11.2.2/source/particles/management/src/G4DecayTable.cc
  * https://github.com/Geant4/geant4/blob/v11.2.2/source/particles/management/src/G4PhaseSpaceDecayChannel.cc

What this gives us
------------------

This is enough to reproduce or audit simple decay sampling such as:

  * choosing a branching channel;
  * sampling daughter directions and phase space in the parent rest frame;
  * boosting daughter four-vectors to the lab frame.

For example, a two-body decay such as Lambda -> p pi- is simple: fixed
daughter momentum in the parent rest frame, isotropic direction, then a
boost by the parent four-vector.

What this does not give us
--------------------------

The detector envelope needs more than just a decay sampler.  It needs:

  * decay position from lifetime and transported parent trajectory;
  * geometry boundary checks;
  * material and active-volume definitions;
  * particle tracking until decay or escape;
  * secondary ancestry records;
  * visible companion activity from parent-decay products;
  * consistency with the MicroBooNE/LArSoft physics list and geometry.

So a standalone sampler based on these source classes is useful for a
fast topology-gun prototype, but it should be labeled as a toy or
truth-level envelope until it is validated against a real Geant4/LArSoft
transport pass.

Recommended use here
--------------------

Use the Geant4 decay source in two ways:

  * as a reference implementation for a fast final-Lambda topology-gun
    prototype;
  * as a validation target for Lambda -> p pi- and simple parent feed-down
    kinematics.

Do not use it alone as the final A_det(y) unless the analysis explicitly
chooses a geometry-only or truth-only envelope.  For a MicroBooNE detector
envelope, use the actual Geant4/LArSoft chain or validate the standalone
sampler against it.

