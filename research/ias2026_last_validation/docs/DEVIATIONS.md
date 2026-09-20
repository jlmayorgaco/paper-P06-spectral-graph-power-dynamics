# Deviations from the IAS2026 master execution prompt

| Item | Required | Observed | Adjudication |
|---|---|---|---|
| Branch and root | New final branch and `research/ias2026_final_closure/` | Created | Satisfied |
| Frozen source | Same frozen Python model and exact Julia port | Python source and frozen custom Julia device/network implementation retained | `SAME_MODEL_CROSS_CODE_VALIDATED` |
| P1 | Canonical Python/Julia GFL11 parity | 64 cases and 602 transfer points verified; network harness residual within declared limit | `NUMERICALLY_VERIFIED` |
| P2 | Same-model SG/AVR/PSS parity | Base, V4, and all 16 frozen custom portfolios reconciled across Python/Julia | `SAME_MODEL_CROSS_CODE_VALIDATED` |
| P3 | Collective mechanism | Julia terminal-operator factorization and local/collective diagnostics completed for V4 | `JULIA_COLLECTIVE_MECHANISM_VALIDATED` |
| P4 | Julia TDS | Same-model base/proper/repaired/blocker traces solved; common pulse did not excite the RHP mode | `NUMERICALLY_VERIFIED` with observability limit |
| P5 | Second dynamic model | PowerDynamics SimpleGFLDC 16/16 census plus 20-case representative bandwidth search completed | `SECOND_MODEL_VALIDATED` / negative result |
| P6 | Genuine mixed holdout | Second-model search found no mixed region; >=12-point mixed holdout not fabricated | `STOPPED_BY_NEGATIVE_RESULT` |
| Robustness | Physical common uncertainty envelope and robust radius | Not executed | `NOT_TESTED` |
| Scaling | New model-size/asymptotic sweep | Not executed | `NOT_TESTED` |
| Submission artifacts | Rewrite paper/poster after gates | Prohibited until gates pass | Not performed |

The remaining limitation is deliberate: no mixed holdout was fabricated after
the second-model policy search returned a negative result. No retuning or
model substitution was used to close that gate.
