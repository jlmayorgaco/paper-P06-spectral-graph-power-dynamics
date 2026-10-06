# Claim ledger — Experiment A

| ID | Claim | Mathematical basis | Evidence | Scope and limitation | Status |
|---|---|---|---|---|---|
| A-01 | The bus-33 PowerDynamics DAE equilibrium and state map are reproducible. | PowerDynamics PF + initialize_from_pf!; E/A state mapping | TABLE_A01–A03 | Installed PowerDynamics 5.0.0 IEEE-39 data and this operating point only. | SUPPORTED |
| A-02 | Algebraic elimination reproduces finite descriptor poles. | Schur elimination of G_y and generalized QZ | descriptor matrices; RESULTS JSON | Numerical conditioning and DAE model only. | SUPPORTED |
| A-03 | The declared retained synchronization operator satisfies the exact Schur identity. | S_r=T_rr-T_rc T_cc^{-1}T_cr | TABLE_A06; FIG_A03 | Conditioned sample points and declared state partition. | SUPPORTED |
| A-04 | Retained-visible full-system poles satisfy the Schur zero residual. | σ_min(S_r(λ)) certificate | TABLE_A05; FIG_A02 | Certificate at full poles; no independent rational root solver. | SUPPORTED |
| A-05 | The retained angle dynamics admit the exact generalized second-order bridge; global Σ(s) remains undefined. | q̇=Cv substitution and exact elimination | TABLE_A06; FIG_A04; RESULTS JSON | Bus 33 and 0.1–5 Hz bridge validation; Π_q(0) is singular, so Σ(s) and graph-modal analysis are not claimed. | SUPPORTED |
| A-06 | A +0.1% local pulse agrees between nonlinear and linearized dynamics within tolerance. | Small-signal linearization comparison | TABLE_A08; FIG_A10 | One bus/load and one local pulse; not transient-stability or EMT validation. | SUPPORTED |
| A-07 | The extraction pipeline reproduces for buses 30, 33, 35 and 37. | Same frozen GFL template and Schur/pole protocol | TABLE_A09; FIG_A11 | Single replacements only, nominal controller settings. | SUPPORTED |
