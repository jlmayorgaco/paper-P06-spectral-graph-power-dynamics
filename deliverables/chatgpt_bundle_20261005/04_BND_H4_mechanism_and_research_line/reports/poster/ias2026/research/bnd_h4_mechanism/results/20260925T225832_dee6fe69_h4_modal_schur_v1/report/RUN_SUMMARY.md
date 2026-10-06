# H4 modal-Schur run

Status: **BLOCKED_M1_STRICT**

Baseline: `C:\Users\walla\Documents\Github\paper-P06-spectral-graph-power-dynamics\reports\poster\ias2026\research\bnd_h4_mechanism\results\20260925T204017_549c3c07_h4_crossmode_v3`

The run stopped before M2/M3 because the retained port operator did not meet the preregistered strict scalar-Schur residual. No figures were generated.

- full DAE pole: `[0.12700646782830968, 3.9098984763819242]`
- `abs(h_k(lambda_c))`: `9.453013584752564e-08`
- strict threshold: `1e-08`
- independent Schur root: `[0.12700684036014512, 3.9098984759249666]`
- root distance to full pole: `3.7253211569341967e-07`
- complement condition: `4.551910971502867`

Interpretation: this is a retained-operator numerical gate failure, not evidence against the collective mechanism. The next safe action is to reconcile the port linearization with the full DAE before running eta continuation.
