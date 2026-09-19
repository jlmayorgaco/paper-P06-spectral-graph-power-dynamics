# IAS2026 final closure dashboard

The machine-readable adjudication is in `GATE_DASHBOARD.csv`. The campaign
stops dependent gates when the true same-model cross-code requirement is not
met. See `docs/NEGATIVE_RESULTS.md` and `docs/DEVIATIONS.md` for the full
negative-result register.

| Gate family | Status | Meaning |
|---|---|---|
| Canonical GFL11 parity | `NUMERICALLY_VERIFIED` | Device and documented harness evidence only |
| Frozen SG/AVR/PSS and IEEE-39 parity | `STOPPED_BY_GATE` | Exact Julia port not completed |
| Official TDS | `NUMERICALLY_VERIFIED` | Alternative-model decay diagnostics; no modal match |
| SimpleGFLDC second model | `SECOND_MODEL_VALIDATED` | 16/16 stable; negative and model-conditioned |
| Genuine holdout | `BLIND_HOLDOUT` | 4/4 stable single-class; not discriminative |
| Robustness and scaling | `NOT_TESTED` / `RETROSPECTIVE` | No new promoted claim |
