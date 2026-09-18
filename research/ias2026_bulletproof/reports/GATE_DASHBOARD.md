# IAS2026 Bulletproof Closure Dashboard

The dashboard is an audit of what was executed in this campaign and what
remains retrospective or stopped. `PASS` means the stated gate passed; it does
not mean every stronger hypothesis in the execution plan passed.

| Gate | Status | Result label | Evidence | Interpretation |
|---|---|---|---|---|
| G0 frozen artifact reproduction | FAIL | `IEEE39_VALIDATED` + `REFUTED` | `GATE0_REPORT.md` | Nominal P4 reproduces; stronger V9 expectation is refuted. |
| G1 theorem/unit tests | PASS | `THEOREM_VALIDATED` | `logs/theory_tests.log` | Four property-test modules passed. |
| G2 PowerDynamics equilibrium | PASS | `POWERDYNAMICS_VALIDATED` | `raw/powerdynamics/pd39_equilibrium_gate.md` | Official IEEE-39 tutorial equilibrium gate passed. |
| G3 same-model Julia/Python parity | STOPPED | `STOPPED_BY_GATE` | `docs/POWERDYNAMICS_RECONCILIATION.md` | Available tutorial model is not the frozen L0 GFL model. |
| G4 same-mechanism blocker parity | STOPPED | `STOPPED_BY_GATE` | `docs/POWERDYNAMICS_RECONCILIATION.md` | No converter-port/mode parity established. |
| G5 second-model holdout | NOT TESTED | `NOT_TESTED` | `docs/SECOND_MODEL_SPEC.md` | No independent converter model was available in scope. |
| G6 Julia TDS | NOT TESTED | `NOT_TESTED` | `reports/TDS_REPORT.md` | Frozen G2 TDS is retrospective, not a new Julia TDS run. |
| G7 robust LFT radius | NOT TESTED | `NOT_TESTED` | `reports/ROBUSTNESS_REPORT.md` | No physical uncertainty radius was preregistered and completed. |
| G8 common uncertainty envelope | NOT TESTED | `NOT_TESTED` | `reports/ROBUSTNESS_REPORT.md` | No common admissible envelope is claimed. |
| G9 V4 robust census | NOT TESTED | `NOT_TESTED` | `reports/ROBUSTNESS_REPORT.md` | The campaign does not promote a new robust census. |
| G10 g x epsilon atlas | NOT TESTED | `NOT_TESTED` | `reports/ROBUSTNESS_REPORT.md` | No new physical grid was run. |
| G11 topology calibration | RETROSPECTIVE | `RETROSPECTIVE` | `derived/tables/gate0_topology.csv` | Frozen F7/F7A evidence is preserved. |
| G12 genuine blind topology/model holdout | NOT TESTED | `NOT_TESTED` | `reports/HOLDOUT_REPORT.md` | Historical V9 is a revealed retrospective test. |
| G13 governed blocker | PASS/NEGATIVE | `IEEE39_VALIDATED` + `REFUTED` | `derived/tables/gate0_v4.csv` | Documented P4 governor changes H4 alpha sign; blocker does not survive that policy. |
| G14 ablations | RETROSPECTIVE | `RETROSPECTIVE` | `reports/ABLATION_REPORT.md` | Existing frozen ablations are reported with scope. |
| G15 scaling | RETROSPECTIVE | `REFUTED` | `reports/SCALING_REPORT.md` | V9 transfer is weaker than the stronger plan expectation. |
| G16 final claim audit | PASS | `CLAIM_AUDITED` | `docs/CLAIM_LEDGER.md` | Claims are separated into survived, weakened, refuted, and not tested. |

The overall campaign is therefore a scientifically useful closure with a
negative result on the stronger transfer/robustness/generalization claims. It
is not a blanket validation of the execution plan's expected numbers.
