# P1–P6 execution summary

P1=PASS (canonical frozen Python/Julia device/port parity, documented harness, and 602-point transfer parity; scope is not full IEEE-39 same-model parity)
P2=FRESH_ALTERNATIVE_DYNAMIC_MODEL_V4_NEGATIVE_HOLDOUT (official PowerDynamics IEEE-39 model, 16/16 initialized, stable transverse spectra)
P3=NO_BLOCKER_FOUND_IN_ALTERNATIVE_MODEL
P4=FRESH_OFFICIAL_POWERDYNAMICS_TDS (base/one/full common-disturbance solves attempted and compared with eigensolve diagnostics)
P5=FRESH_ALTERNATIVE_DYNAMIC_MODEL_V4_NEGATIVE_HOLDOUT (official SimpleGFLDC IEEE-39 census, matched P/Q/rating, no retune)
P6=FRESH_ALTERNATIVE_DYNAMIC_MODEL_BLIND_HOLDOUT (genuine input/prediction hash chronology; 4/4 accuracy but all-stable, so no discrimination claim)

CANONICAL PYTHON/JULIA PORT ERROR: derivatives 2.79e-13 max; terminal current 0; system-base P/Q 2.22e-16 max; initialization 1.11e-16 max; transfer Frobenius relative median 1.41e-16 and 99th percentile 1.52e-15.
POWERDYNAMICS ALTERNATIVE-MODEL RESULT: none of the 16 official-model portfolios had a transverse blocker after the corrected current-source topology; this is not a result for the frozen custom no-governor model.
SECOND-MODEL RESULT: official SimpleGFLDC initialized all 16 matched IEEE-39 portfolios without retuning; this is a different-model negative result. The separate component harness passed with 13 states and residual 6.31e-9.
BLIND HOLDOUT: the historical v1 holdout is procedural/single-class and chronology is not independently provable; the genuine v2 alternative-model holdout has independently represented input/prediction/reveal chronology and 4/4 accuracy, but is also single-class, so balanced accuracy, unstable precision/recall, MCC, kappa, blocker discrimination, and repair ranking are not promoted.
TRUE SAME-MODEL GATE: STOPPED_BY_GATE. The canonical custom SynchronousMachine/AVR/PSS source was snapshotted, but the exact Julia port and full IEEE-39 equilibrium/A_perp/mode-MAC reconciliation were not completed.
STRONGEST SURVIVED CLAIM: the canonical frozen GFL11 device/port and documented PowerDynamics harness parity passes; in official PowerDynamics IEEE-39, all 16 portfolios are transversely stable, reinforcing that the observed blocker depends on surrounding dynamics/model policy rather than cardinality alone.
STRONGEST NEGATIVE RESULT: no fresh Julia blocker was found, so blocker-specific P3 operator/Q-continuation and P4 proper-subset/blocker/repaired-blocker TDS evidence cannot be claimed; complete-spectrum V9 transfer agreement also remains incomplete because 110 aperiodic unstable portfolios are outside the targeted representation, while targeted 0.3–1.5 Hz agreement is 511/512 with an exact 14-blocker antichain.

Robustness was not started, and no paper or poster was rewritten.
