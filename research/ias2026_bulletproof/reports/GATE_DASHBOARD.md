# IAS2026 Bulletproof Closure Dashboard

The dashboard is an audit of what was executed in this campaign and what
remains retrospective or stopped. `PASS` means the stated gate passed; it does
not mean every stronger hypothesis in the execution plan passed.

| Gate | Status | Result label | Evidence | Interpretation |
|---|---|---|---|---|
| G0 corrected artifact reproduction | PASS | `RETROSPECTIVE` + `NUMERICALLY_VERIFIED` | `GATE0_REPORT.md` | Corrected frozen-artifact reproduction separates global from targeted V9 results; it is not a fresh IEEE-39 runtime solve. |
| G1 algebra/property tests | PASS | `NUMERICALLY_VERIFIED` | `logs/theory_tests.log` | Four property-test modules passed; this is not a theorem proof. |
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
| G13 governed blocker | PASS/NEGATIVE | `RETROSPECTIVE` + `REFUTED` | `derived/tables/gate0_v4.csv` | Frozen P4 governor evidence changes H4 alpha sign; blocker does not survive that policy. |
| G14 ablations | RETROSPECTIVE | `RETROSPECTIVE` | `reports/ABLATION_REPORT.md` | Existing frozen ablations are reported with scope. |
| G15 scaling | RETROSPECTIVE | `NUMERICALLY_VERIFIED` | `reports/SCALING_REPORT.md` | Corrected global and targeted V9 scopes are recorded; no new asymptotic sweep. |
| G16 claim ledger audit | PASS | `NUMERICALLY_VERIFIED` | `docs/CLAIM_LEDGER.md` | Claims are separated into survived, weakened, refuted, and not tested. |
| P1 exact frozen GFL11 parity | STOPPED | `STOPPED_BY_GATE` | `reports/P1_GFL_PARITY.md` | Fresh Julia/Python device parity and PowerDynamics compilation pass; the minimal current-source infinite-bus network retained zero states. |
| P2 PowerDynamics V4 census | STOPPED | `STOPPED_BY_GATE` | `reports/P2_POWERDYNAMICS_V4_STATUS.md` | Blocked by the complete P1 gate; no fresh 16-portfolio census is claimed. |
| P3 collective mechanism | STOPPED | `STOPPED_BY_GATE` | `reports/P3_COLLECTIVE_MECHANISM_STATUS.md` | No fresh Julia blocker was available. |
| P4 Julia TDS | STOPPED | `STOPPED_BY_GATE` | `reports/P4_JULIA_TDS_STATUS.md` | No fresh proper-subset/blocker/repaired-blocker TDS was run. |
| P5 second GFL model | STOPPED | `STOPPED_BY_GATE` | `reports/P5_SECOND_GFL_STATUS.md` | No materially different model was frozen and executed. |
| P6 blind holdout | STOPPED | `STOPPED_BY_GATE` | `reports/P6_BLIND_HOLDOUT_STATUS.md` | No fresh prediction hash/reveal exists. |

The overall campaign is therefore a scientifically useful closure with a
negative result on the stronger transfer/robustness/generalization claims. It
is not a blanket validation of the execution plan's expected numbers.
