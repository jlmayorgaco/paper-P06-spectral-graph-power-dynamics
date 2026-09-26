# Model reconciliation

## Frozen canonical model

The frozen custom model is a SynchronousMachine with first-order AVR and a
two-state IEEEST PSS, zero damping, constant mechanical power, and no
governor. Its source, network, and configuration snapshots are in
`raw/same_model/canonical_source/` and are hash-listed in
`raw/same_model/canonical_source_manifest.json`.

## Python/Julia evidence

The GFL11 current-source device and transfer harness were implemented in both
languages. The canonical sweep has 64 cases and the transfer audit has 602
frequency points. The maximum state relative error is approximately
5.14e-16; the network residual is approximately 2.09e-13 in the canonical
Python run and 4.22e-9 in the Julia harness, below the declared 1e-6 residual
limit.

This establishes device/port parity and a documented network harness. It does
not establish exact full IEEE-39 parity for the frozen custom synchronous
machine model.

## Alternative models

The PowerDynamics census uses official machines/AVR/governors. The
SimpleGFLDC census uses an official inverter model. Both are preserved as
alternative-model evidence and explicitly cannot be substituted for the
frozen custom model.

## Gate outcome

The exact Julia implementation required for device derivatives, initialization,
Jacobian, frequency response, equilibrium reconciliation, transverse
projection, and mode MAC was not completed. The same-model gate is therefore
`STOPPED_BY_GATE`; no downstream mechanism or robustness claim is promoted.
