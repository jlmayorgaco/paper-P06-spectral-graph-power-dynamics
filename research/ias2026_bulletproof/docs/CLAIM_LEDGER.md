# Claim Ledger

This ledger separates what survived the closure from what was weakened,
refuted, or not tested. Labels are intentionally conservative.

| ID | Claim | Status label | Safe wording | Evidence / limitation |
|---|---|---|---|---|
| C01 | P4 H4 nominal alpha is approximately 0.1270065 s^-1 | `IEEE39_VALIDATED` | Frozen P4 reproduces the reported nominal H4 value. | Gate 0, frozen FC03; artifact-level reproduction. |
| C02 | P4 is minimally incompatible in the frozen policy | `IEEE39_VALIDATED` | The frozen P4 structure has 15 stable proper subsets and the reported H4 blocker. | FC01/FC03; policy and model scope are essential. |
| C03 | Transverse/collective algebra is correct under stated assumptions | `NUMERICALLY_VERIFIED` | Randomized property tests pass, including the explicit contextual-return assumption. | Four pytest modules; property evidence is not a theorem proof. |
| C04 | Boundary occurs near g*=0.2076814 and f=0.7064248 Hz | `RETROSPECTIVE` | The frozen F7A boundary table contains the reported witness. | Retrospective exported table. |
| C05 | The boundary Q sign is negative and contextual return sign positive | `RETROSPECTIVE` | The frozen boundary audit reports the two signs. | Retrospective; not a new scan. |
| C06 | Frozen TDS agrees with declared verdicts | `NONLINEAR_TDS_VALIDATED` | The frozen G2 table reports 32/32 agreement. | Retrospective G2; no new Julia TDS. |
| C07 | P4 governor changes H4 alpha to approximately -0.074541 s^-1 | `IEEE39_VALIDATED` | Under the documented governor policy, the P4 H4 alpha changes sign. | Frozen FC03; policy-dependent. |
| C08 | Same-policy V9 transfer, global spectrum | `NUMERICALLY_VERIFIED` | The corrected global census has 402/512 correct, 110 false-safe, and 0 false-unstable cases; all false-safe cases are `APERIODIC_REAL`. | Corrected modal-scope source; global result is not a targeted-band certificate. |
| C09 | Same-policy V9 transfer, targeted 0.3–1.5 Hz band | `NUMERICALLY_VERIFIED` | The corrected targeted census has 511/512 correct, 1 false-safe, and 0 false-unstable cases. | Corrected modal-scope source; scope is the declared band only. |
| C10 | Target-family blocker antichain is exact | `NUMERICALLY_VERIFIED` | The corrected target-family antichain has 14 blockers with orders 5x4, 6x5, 3x6 and kappa=4. | Corrected same-policy modal-scope source. |
| C11 | PowerDynamics independently reproduces the frozen GFL mechanism | `STOPPED_BY_GATE` | PowerDynamics Gate A validates only its official IEEE-39 tutorial equilibrium. | Different synchronous-machine/governor model; no parity. |
| C12 | A second converter model confirms the result | `NOT_TESTED` | No second converter model was executed. | Scope limitation. |
| C13 | A physical common uncertainty envelope/robust radius is established | `NOT_TESTED` | No physical robust radius is claimed. | No completed LFT/uncertainty campaign. |
| C14 | A new topology/model holdout confirms generalization | `NOT_TESTED` | No new blind topology/model holdout is claimed. | Historical V9 is retrospective. |
| C15 | Repair is certified for all dynamics and switching paths | `NOT_TESTED` | Only the documented local governor effect is reported. | No EMT/current-limit/DC-link certification. |
| C16 | Full IEEE-39 conclusion is universal | `REFUTED` | The result is benchmark- and policy-conditioned, not universal. | Scope and corrected global/targeted V9 evidence. |

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
