# TX3 Final Manuscript and Review Package

This package closes the TX3 materiality search at E05D and contains the final
manuscript, its reproducible source, the frozen claim envelope, the descriptive
critical-mode baseline audit, and two independent adversarial reviewer audits.

## Terminal interpretation

- E05D C3c remains **UNRESOLVED** because all eight empty-coalition baselines
  violate the immutable 5% damping eligibility threshold at kappa = 1.
- The terminal E05D status is **BASELINE_INELIGIBLE**. No E05D coalition outcome
  was inspected after that gate failed.
- The TX3 materiality search is closed. This package does not introduce E06,
  change any threshold, stress branch, tracking rule, or direction, and does not
  add a ParaEMT campaign.
- The kappa = 1 critical-mode work is a descriptive post-mortem, not a new gate
  or a reopening of materiality testing.

## Recommended reading order

1. `manuscript/TX3_FINITE_CONNECTED_EXTERNALITIES_FINAL_2026-08-28.pdf`
2. `claim_freeze/TX3_FINAL_CLAIM_FREEZE.md`
3. `critical_mode_audit/report/TX3_CRITICAL_MODE_BASELINE_AUDIT.md`
4. `reviews/TX3_FINAL_REVIEWER1_AUDIT.md`
5. `reviews/TX3_FINAL_REVIEWER2_AUDIT.md`
6. `VALIDATION_SUMMARY.md`

The complete numerical audit, source code, unit test, and terminal decision
records are included for independent inspection. `MANIFEST_SHA256.txt` records
the SHA-256 digest of every packaged file other than the manifest itself.

## Manuscript scope

The paper establishes a connected finite-intervention calculus for AC
re-equilibrated converter-rich power-system models. Its central boundary is
that the existence and exact reconstructability of a connected externality do
not imply operational materiality. The finite-pole mechanism is supported; the
immutable damping-margin materiality claim is not established because E05D is
baseline-ineligible.

