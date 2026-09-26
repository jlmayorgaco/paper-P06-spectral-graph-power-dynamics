# Deviations from the IAS2026 master execution prompt

| Item | Required | Observed | Adjudication |
|---|---|---|---|
| Branch and root | New final branch and `research/ias2026_final_closure/` | Created | Satisfied |
| Frozen source | Same frozen Python model and exact Julia port | Python source and manifest frozen; exact Julia device/network port incomplete | `STOPPED_BY_GATE` |
| P1 | Canonical Python/Julia GFL11 parity | 64 cases and 602 transfer points verified; network harness residual within declared limit | `NUMERICALLY_VERIFIED` |
| P2 | Same-model SG/AVR/PSS parity | Not established | `STOPPED_BY_GATE` |
| P3 | Collective mechanism | Dependent gate stopped | `NOT_TESTED` |
| P4 | Julia TDS | Official alternative model base/one/full solves completed; no reliable modal crossing | `NUMERICALLY_VERIFIED` with scope limit |
| P5 | Second dynamic model | PowerDynamics SimpleGFLDC 16/16 census completed | `SECOND_MODEL_VALIDATED` / negative result |
| P6 | Genuine mixed holdout | Executed, but all four labels stable | `BLIND_HOLDOUT` with single-class limitation |
| Robustness | Physical common uncertainty envelope and robust radius | Not executed | `NOT_TESTED` |
| Scaling | New model-size/asymptotic sweep | Not executed | `NOT_TESTED` |
| Submission artifacts | Rewrite paper/poster after gates | Prohibited until gates pass | Not performed |

The deviation is deliberate: dependent claims were not promoted after the
same-model gate stopped. No retuning or model substitution was used to close
that gate.
