# True same-model Julia/Python reconciliation

status: STOPPED_BY_GATE
result_label: STOPPED_BY_GATE
evidence_class: STOPPED_BY_GATE

The canonical frozen model is the custom SynchronousMachine with first-order AVR and two-state IEEEST PSS, D=0, constant mechanical power, and no governor. Its source and network snapshots are frozen under `raw/true_same_model/canonical_source/`.

The fresh P2 implementation uses the official PowerDynamics IEEE-39 machines/AVR/governors around the GFL replacement, and P5 uses official SimpleGFLDC. Those are valid alternative-model negative results but are not the same-model cross-code gate. The exact Julia port required for device derivatives/current/initialization/Jacobian/frequency response, full equilibrium reconciliation, A-perp, and mode MAC was not completed; those checks remain NOT_RUN. No retuning or model substitution is used to upgrade this gate.

Manifest: `raw/true_same_model/canonical_source_manifest.json`.
