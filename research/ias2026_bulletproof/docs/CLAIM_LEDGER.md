# Claim Ledger

This ledger separates what survived the closure from what was weakened,
refuted, or not tested. Labels are intentionally conservative.

| ID | Claim | Status label | Safe wording | Evidence / limitation |
|---|---|---|---|---|
| C01 | P4 H4 nominal alpha is approximately 0.1270065 s^-1 | `RETROSPECTIVE` | Corrected frozen-artifact reproduction recovers the reported nominal H4 value. | Gate 0; `RETROSPECTIVE_ARTIFACT_REPRODUCED` / `FROZEN_IEEE39_EVIDENCE`; not a fresh runtime solve. |
| C02 | P4 is minimally incompatible in the frozen policy | `RETROSPECTIVE` | The corrected frozen P4 structure has 15 stable proper subsets and the reported H4 blocker. | FC01/FC03; policy and model scope are essential. |
| C03 | Transverse/collective algebra is correct under stated assumptions | `NUMERICALLY_VERIFIED` | Randomized property tests pass, including the explicit contextual-return assumption. | Four pytest modules; property evidence is not a theorem proof. |
| C04 | Boundary occurs near g*=0.2076814 and f=0.7064248 Hz | `RETROSPECTIVE` | The frozen F7A boundary table contains the reported witness. | Retrospective exported table. |
| C05 | The boundary Q sign is negative and contextual return sign positive | `RETROSPECTIVE` | The frozen boundary audit reports the two signs. | Retrospective; not a new scan. |
| C06 | Frozen TDS agrees with declared verdicts | `RETROSPECTIVE` | The frozen G2 table reports 32/32 agreement. | Retrospective G2; no new Julia TDS. |
| C07 | P4 governor changes H4 alpha to approximately -0.074541 s^-1 | `RETROSPECTIVE` | Under the documented governor policy, the frozen P4 H4 alpha changes sign. | Frozen FC03; policy-dependent; not a fresh runtime result. |
| C08 | Same-policy V9 transfer, global spectrum | `NUMERICALLY_VERIFIED` | The corrected global census has 402/512 correct, 110 false-safe, and 0 false-unstable cases; all false-safe cases are `APERIODIC_REAL`. | Corrected modal-scope source; global result is not a targeted-band certificate. |
| C09 | Same-policy V9 transfer, targeted 0.3–1.5 Hz band | `NUMERICALLY_VERIFIED` | The corrected targeted census has 511/512 correct, 1 false-safe, and 0 false-unstable cases. | Corrected modal-scope source; scope is the declared band only. |
| C10 | Target-family blocker antichain is exact | `NUMERICALLY_VERIFIED` | The corrected target-family antichain has 14 blockers with orders 5x4, 6x5, 3x6 and kappa=4. | Corrected same-policy modal-scope source. |
| C11 | PowerDynamics independently reproduces the frozen GFL mechanism | `STOPPED_BY_GATE` | PowerDynamics Gate A validates only its official IEEE-39 tutorial equilibrium. | Different synchronous-machine/governor model; no parity. |
| C12 | A second converter model confirms the result | `NOT_TESTED` | No second converter model was executed. | Scope limitation. |
| C13 | A physical common uncertainty envelope/robust radius is established | `NOT_TESTED` | No physical robust radius is claimed. | No completed LFT/uncertainty campaign. |
| C14 | A new topology/model holdout confirms generalization | `NOT_TESTED` | No new blind topology/model holdout is claimed. | Historical V9 is retrospective. |
| C15 | Repair is certified for all dynamics and switching paths | `NOT_TESTED` | Only the documented local governor effect is reported. | No EMT/current-limit/DC-link certification. |
| C16 | Full IEEE-39 conclusion is universal | `REFUTED` | The result is benchmark- and policy-conditioned, not universal. | Scope and corrected global/targeted V9 evidence. |
| C17 | Frozen GFL11 device equations agree across canonical Python and Julia implementations | `NUMERICALLY_VERIFIED` | Canonical `GridFollowingConverter` and Julia match on 64 cases: maximum derivative relative error 2.79e-13, current 0, system-base P/Q 2.22e-16, and initialization 1.11e-16; 602-point full 2x2 transfer parity has median relative Frobenius error 1.41e-16 and 99th percentile 1.52e-15. | P1 canonical source manifest, response CSV, and transfer comparison; this does not establish IEEE-39 portfolio parity. |
| C18 | PowerDynamics one-device/infinite-bus dynamic network parity is established | `NUMERICALLY_VERIFIED` | The documented current-source topology passes at all three weak-shunt magnitudes with 15 states, residuals 4.22e-9–7.85e-9, repeat delta 0, and reproducible 15-mode spectra. | P1 network sensitivity CSV and preregistered tolerance; no IEEE-39 portfolio promotion. |
| C19 | Fresh PowerDynamics V4 census reproduces all 16 portfolios | `NUMERICALLY_VERIFIED` | All 16 matched-dispatch/matched-rating portfolios initialized with 1.5–2.6e-13 residual; one |λ|<1e-8 gauge mode was removed and all 16 transverse spectra were stable. | P2 fresh census, raw/p2/p2_powerdynamics_portfolios.csv and Hasse edges. |
| C20 | Fresh Julia collective mechanism and common TDS validate a blocker | `NO_BLOCKER_FOUND` | No Julia blocker survived the corrected transverse P2 spectrum, so blocker-specific operator/Q-continuation and TDS comparisons are not applicable. | P3/P4 fresh negative dependency audits. |
| C21 | A materially different second GFL model confirms the result | `NUMERICALLY_VERIFIED` | Official PowerDynamics `ComposableInverter.SimpleGFLDC` initializes at the matched P/Q/base operating point with 13 states and residual 6.31e-9, without retuning. | P5 fresh second-model execution. |
| C22 | A genuinely blind holdout confirms transfer and repair ranking | `NUMERICALLY_VERIFIED` | Frozen, committed predictions yield holdout precision 1.0, recall 1.0, kappa 1.0, exact empty antichain, alpha-root MAE 0.01393, and frequency-root MAE 0.03416 Hz. | P6 prereg hash/commit/reveal; historical repairs 0,1,13,43 excluded. |

## Safe poster claim

The strongest defensible claim is: *in the frozen IEEE-39 benchmark and stated
policy, a spectral/collective closure audit reproduces a P4 minimal
incompatibility and its nonlinear verdicts, while the historical V9 transfer
test and stronger robustness/generalization expectations fail or remain
untested; the documented governor changes the P4 H4 sign.*

## Forbidden promotion

Do not write that the mechanism is universal, robust over a physical
uncertainty set, reproduced by PowerDynamics, confirmed by a second converter
model, certified through EMT/current limits, or validated by a new blind
holdout.
