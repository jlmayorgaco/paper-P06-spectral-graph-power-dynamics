# IAS2026 Bulletproof Closure Dashboard

The dashboard is an audit of what was executed in this campaign and what
remains retrospective or stopped. `PASS` means the stated gate passed; it does
not mean every stronger hypothesis in the execution plan passed.

| Gate | Status | Result label | Evidence | Interpretation |
|---|---|---|---|---|
| G0 corrected artifact reproduction | PASS | `RETROSPECTIVE` + `NUMERICALLY_VERIFIED` | `GATE0_REPORT.md` | Corrected frozen-artifact reproduction separates global from targeted V9 results; it is not a fresh IEEE-39 runtime solve. |
| G1 algebra/property tests | PASS | `NUMERICALLY_VERIFIED` | `logs/theory_tests.log` | Four property-test modules passed; this is not a theorem proof. |
| G2 PowerDynamics equilibrium | PASS | `POWERDYNAMICS_VALIDATED` | `raw/powerdynamics/pd39_equilibrium_gate.md` | Official IEEE-39 tutorial equilibrium gate passed. |
| G3 same-model Julia/Python parity (full IEEE-39 frozen custom model) | STOPPED | `STOPPED_BY_GATE` | `docs/TRUE_SAME_MODEL_RECONCILIATION.md` | Canonical custom SynchronousMachine/AVR/PSS was snapshotted; exact Julia port and full reconciliation were not completed. |
| G4 same-mechanism blocker parity (full IEEE-39 frozen custom model) | STOPPED | `STOPPED_BY_GATE` | `docs/TRUE_SAME_MODEL_RECONCILIATION.md` | No same-model A-perp/mode-MAC parity is promoted. |
| G5 second-model IEEE-39 holdout | NEGATIVE_HOLDOUT | `FRESH_ALTERNATIVE_DYNAMIC_MODEL_V4_NEGATIVE_HOLDOUT` | `reports/P5_SECOND_GFL_STATUS.md` | Official SimpleGFLDC initialized all 16 portfolios; this is a different dynamic model, not same-model parity. |
| G6 Julia TDS (official alternative model) | PARTIAL_OR_STOPPED | `FRESH_OFFICIAL_POWERDYNAMICS_TDS` | `reports/P4_JULIA_TDS_STATUS.md` | Common-disturbance TDS was attempted for base/one/full cases; see per-case solver status. |
| G7 robust LFT radius | NOT TESTED | `NOT_TESTED` | `reports/ROBUSTNESS_REPORT.md` | No physical uncertainty radius was preregistered and completed. |
| G8 common uncertainty envelope | NOT TESTED | `NOT_TESTED` | `reports/ROBUSTNESS_REPORT.md` | No common admissible envelope is claimed. |
| G9 V4 robust census | NOT TESTED | `NOT_TESTED` | `reports/ROBUSTNESS_REPORT.md` | The campaign does not promote a new robust census. |
| G10 g x epsilon atlas | NOT TESTED | `NOT_TESTED` | `reports/ROBUSTNESS_REPORT.md` | No new physical grid was run. |
| G11 topology calibration | RETROSPECTIVE | `RETROSPECTIVE` | `derived/tables/gate0_topology.csv` | Frozen F7/F7A evidence is preserved. |
| G12 genuine blind topology/model holdout | NOT TESTED | `NOT_TESTED` | `reports/HOLDOUT_REPORT.md` | Historical V9 is a revealed retrospective test. |
| G13 governed blocker | PASS/NEGATIVE | `RETROSPECTIVE` + `REFUTED` | `derived/tables/gate0_v4.csv` | Frozen P4 governor evidence changes H4 alpha sign; blocker does not survive that policy. |
| G14 ablations | RETROSPECTIVE | `RETROSPECTIVE` | `reports/ABLATION_REPORT.md` | Existing frozen ablations are reported with scope. |
| G15 scaling | RETROSPECTIVE | `NUMERICALLY_VERIFIED` | `reports/SCALING_REPORT.md` | Corrected global and targeted V9 scopes are recorded; no new asymptotic sweep. |
| G16 claim ledger audit | PASS | `NUMERICALLY_VERIFIED` | `docs/CLAIM_LEDGER.md` | Claims are separated into survived, weakened, refuted, and not tested. |
| P1 exact frozen GFL11 parity | PASS | `NUMERICALLY_VERIFIED` | `reports/P1_GFL_PARITY.md` | Canonical Python/Julia device and 602-point transfer parity pass; documented current-source network topology passes all three weak-shunt sensitivity gates. |
| P2 PowerDynamics V4 census | NEGATIVE_HOLDOUT | `FRESH_ALTERNATIVE_DYNAMIC_MODEL_V4_NEGATIVE_HOLDOUT` | `reports/P2_POWERDYNAMICS_V4_STATUS.md` | Official PowerDynamics IEEE-39 machines/AVR/governors initialized all 16 portfolios with stable transverse spectra; not the frozen custom model. |
| P3 collective mechanism | NO BLOCKER | `NO_BLOCKER_FOUND` | `reports/P3_COLLECTIVE_MECHANISM_STATUS.md` | No Julia blocker remained after the corrected P2 transverse spectrum; blocker-specific diagnostics are not applicable. |
| P4 Julia TDS | PASS_WITHOUT_MODAL_MATCH | `FRESH_OFFICIAL_POWERDYNAMICS_TDS` | `reports/P4_JULIA_TDS_STATUS.md` | Base/one/full official-model trajectories solved; selected-trace decay was recorded, but no reliable oscillation-frequency match was estimable. |
| P5 second GFL model | NEGATIVE_HOLDOUT | `FRESH_ALTERNATIVE_DYNAMIC_MODEL_V4_NEGATIVE_HOLDOUT` | `reports/P5_SECOND_GFL_STATUS.md` | Official SimpleGFLDC initialized all 16 matched IEEE-39 portfolios without retuning; stable labels are a negative result for that model. |
| P6 blind holdout | SINGLE_CLASS_HOLDOUT_REVEALED | `FRESH_ALTERNATIVE_DYNAMIC_MODEL_BLIND_HOLDOUT` | `reports/P6_BLIND_HOLDOUT_STATUS.md` | Genuine v2 input/prediction/reveal chronology is recorded; 4/4 accuracy on an all-stable holdout, so no discrimination or kappa claim. |

The overall campaign is therefore a scientifically useful closure with a
negative result on the stronger transfer/robustness/generalization claims. It
is not a blanket validation of the execution plan's expected numbers.
