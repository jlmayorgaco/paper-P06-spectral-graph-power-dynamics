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
| C12 | A second converter model confirms the result | `REFUTED_AS_GENERALIZATION` | The official SimpleGFLDC model gives a negative result on the frozen V4 census and on the representative PLL/current bandwidth search; it does not reproduce the custom blocker. | `raw/p5/`, `raw/second_model_search/`; model-conditioned only. |
| C13 | A physical common uncertainty envelope/robust radius is established | `NOT_TESTED` | No physical robust radius is claimed. | No completed LFT/uncertainty campaign. |
| C14 | A new topology/model holdout confirms generalization | `NOT_TESTED` | No new blind topology/model holdout is claimed. | Historical V9 is retrospective. |
| C15 | Repair is certified for all dynamics and switching paths | `NOT_TESTED` | Only the documented local governor effect is reported. | No EMT/current-limit/DC-link certification. |
| C16 | Full IEEE-39 conclusion is universal | `REFUTED` | The result is benchmark- and policy-conditioned, not universal. | Scope and corrected global/targeted V9 evidence. |
| C17 | Frozen GFL11 device equations agree across canonical Python and Julia implementations | `NUMERICALLY_VERIFIED` | Canonical `GridFollowingConverter` and Julia match on 64 cases: maximum derivative relative error 2.79e-13, current 0, system-base P/Q 2.22e-16, and initialization 1.11e-16; 602-point full 2x2 transfer parity has median relative Frobenius error 1.41e-16 and 99th percentile 1.52e-15. | P1 canonical source manifest, response CSV, and transfer comparison; this does not establish IEEE-39 portfolio parity. |
| C18 | PowerDynamics one-device/infinite-bus dynamic network parity is established | `NUMERICALLY_VERIFIED` | The documented current-source topology passes at all three weak-shunt magnitudes with 15 states, residuals 4.22e-9–7.85e-9, repeat delta 0, and reproducible 15-mode spectra. | P1 network sensitivity CSV and preregistered tolerance; no IEEE-39 portfolio promotion. |
| C19 | Fresh official PowerDynamics V4 alternative-model census | `FRESH_ALTERNATIVE_DYNAMIC_MODEL_V4_NEGATIVE_HOLDOUT` | All 16 matched-dispatch/system-base portfolios initialized; one |λ|<1e-8 gauge mode was removed and all 16 transverse spectra were stable. Dynamic rating equivalence is not established. This is not the frozen custom no-governor model. | P2 fresh census, raw/p2/p2_powerdynamics_portfolios.csv and Hasse edges. |
| C20 | Fresh official-model Julia TDS | `FRESH_OFFICIAL_POWERDYNAMICS_TDS` | Base/one/full common-disturbance trajectories solved; selected-trace decay was recorded, but no reliable oscillation-frequency match was estimable. | P3/P4; blocker-specific diagnostics are not promoted. |
| C21 | The official SimpleGFLDC second-model component initializes in the documented harness | `COMPONENT_EXECUTION_PASS` | Official PowerDynamics `ComposableInverter.SimpleGFLDC` initializes at the declared operating point with 13 states and residual 6.31e-9, without retuning. | P5 component execution. |
| C22 | A four-portfolio procedural holdout gives 4/4 stable-label accuracy | `POSTHOC_OR_PROCEDURAL_SINGLE_CLASS_HOLDOUT` | The historical holdout is all-stable; precision/recall discrimination and Cohen kappa are not meaningful, and chronology is not independently provable from the bundle. | P6 v1 reveal metrics. |
| C23 | A genuine alternative-model holdout freezes unseen inputs and predictions before reveal | `FRESH_ALTERNATIVE_DYNAMIC_MODEL_BLIND_HOLDOUT` | The v2 holdout has independently represented input/prediction/reveal chronology and 4/4 accuracy, but all four labels are stable; balanced accuracy, unstable precision/recall, MCC, kappa, blocker discrimination, and repair ranking are not promoted. | P6 genuine protocol, prediction hash, and reveal metrics. |
| C24 | Official SimpleGFLDC IEEE-39 V4 census | `FRESH_ALTERNATIVE_DYNAMIC_MODEL_V4_NEGATIVE_HOLDOUT` | All 16 matched-dispatch/system-base IEEE-39 portfolios initialize and remain transversely stable after gauge-aware filtering, without retuning; matched rating is false and rating equivalence is not established. | P5 full census. |
| C25 | True frozen custom same-model Python↔Julia reconciliation | `SAME_MODEL_CROSS_CODE_VALIDATED` | Base, V4, and all 16 census portfolios agree with state/algebraic max error below 1.05e-7, critical real/frequency errors below 2.5e-10/2.4e-11, and MAC=1.0. | `raw/reconciliation/same_model_reconciliation.csv`; custom SG/AVR/PSS/GFL model only. |
| C26 | Frozen custom same-model V4 blocker | `SAME_MODEL_CROSS_CODE_VALIDATED` | The all-target 30+33+35+37 case is transversely unstable with alpha approximately 0.1446702 s^-1 at 0.574681 Hz; 15 proper subsets are stable under the same gauge-aware census. | `raw/reconciliation/julia_true_same_model_v4_census.csv`; policy/model conditioned. |
| C27 | Julia collective mechanism for the same-model blocker | `JULIA_COLLECTIVE_MECHANISM_VALIDATED` | The 8-dimensional port update factorizes with maximum normalized residual 7.22e-16; minimum local factor sigma is 0.294, while the collective Q closure margin reaches 0.0881 on the frozen imaginary-axis grid. | `raw/mechanism/julia_same_model_mechanism.csv`; no universal boundary claim. |
| C28 | Same-model nonlinear Julia TDS | `NUMERICALLY_VERIFIED_WITH_OBSERVABILITY_LIMIT` | Base, proper, repaired, physical-pulse blocker, and mode-seeded blocker traces integrate with algebraic residual ≤1e-9. The common short pulse does not visibly excite the RHP mode; the nonlinear trace is not promoted as an independent instability proof. | `raw/tds_same_model/`; linear spectrum remains the blocker evidence. |
| C29 | Second-model mixed-policy search | `NEGATIVE_RESULT` | PLL and current-loop bandwidth multipliers 0.125, 0.5, 1, 2, and 8 were searched on proper and flagship portfolios; all 20 initialized cases were stable. No mixed second-model policy region was found. | `raw/second_model_search/`; no mixed holdout was fabricated. |
| C30 | >=12-point mixed-class holdout | `STOPPED_BY_NEGATIVE_RESULT` | The second-model search produced no mixed stable/unstable region, so a >=12-point mixed holdout and reveal score were not honestly constructible. | `raw/holdout_mixed/holdout_status.json`. |

## Safe poster claim

The strongest defensible claim is: *the canonical frozen GFL11 device/port and
documented PowerDynamics harness parity pass; the exact frozen custom SG/AVR/
PSS/load/network model now reconciles across Python and Julia, and its
all-target V4 case is an independently reproduced transverse blocker while
proper subsets are stable. A Julia port-space audit localizes the result to a
collective factor with regular local factors. In the official PowerDynamics
SimpleGFLDC alternative model, the 16-case census and representative bandwidth
search are stable; this is a negative model-conditioned result, not a universal
confirmation.*

The corrected V9 wording is narrower: complete-spectrum transfer agreement is
incomplete because 110 aperiodic unstable portfolios lie outside the targeted
representation; within the targeted 0.3–1.5 Hz band, 511/512 cases agree and
the 14-blocker antichain is exact. No robustness atlas was started before the
P1–P6 execution gates were closed.

## Forbidden promotion

Do not write that the mechanism is universal, robust over a physical
uncertainty set, reproduced by the official PowerDynamics models, certified
through EMT/current limits, or confirmed by a mixed second-model holdout. Do
not describe the short nonlinear pulse as an independent proof of the RHP
mode; it did not visibly excite that mode in the retained trace.
