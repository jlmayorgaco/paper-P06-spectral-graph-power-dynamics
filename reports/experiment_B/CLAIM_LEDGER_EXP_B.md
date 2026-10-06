# CLAIM_LEDGER_EXP_B

| ID | Statement | Type | Derivation / evidence | Assumptions and limitations | Status |
|---|---|---|---|---|---|
| B-C01 | Graph congruence preserves the exact nonlinear eigenvalue problem zeros. | EXACT | `T̂=ΦᴴTΦ`; Φ is nonsingular and `det(T̂)=det(Φᴴ)det(T)det(Φ)`. | M positive definite and Φ square/nonsingular. | SUPPORTED |
| B-C02 | Γₖ is exactly the Schur self-energy of complementary graph modes. | EXACT | Block Gaussian elimination and the Schur residual/determinant diagnostics. | T_rr nonsingular at the evaluation point. | SUPPORTED |
| B-C03 | Off-diagonal coupling O(ε) gives Γₖ=O(ε²). | ASYMPTOTIC | Two outer cross-block factors are O(ε); synthetic log-log fits are in TABLE_B06. | Complement remains nonsingular; local expansion only. | SUPPORTED |
| B-C04 | The pole predictor is accurate in a weak-coupling interval. | EMPIRICAL | Tracked augmented-system poles versus `s₀−Γ/D` in TABLE_B05. | Simple root, nonzero derivative, no branch collision; case-specific. | SUPPORTED |
| B-C05 | Near-resonance amplifies collective self-energy. | EMPIRICAL | Detuning and complementary-stiffness sweep in TABLE_B08. | Synthetic pair, not a universal threshold. | SUPPORTED |
| B-C06 | Positive diagonal graph dissipative entries do not ensure a positive semidefinite total operator. | EXACT COUNTEREXAMPLE | Optional S5 has positive diagonal and a negative minimum eigenvalue. | This does not by itself imply instability. | SUPPORTED |
| B-C07 | χ_comm measures graph/self-energy noncommutativity, not stability margin. | EXACT DEFINITION | Commutator identity and S0/S1-S4 diagnostics. | Interpretation depends on graph eigenvalue degeneracy. | SUPPORTED |
| B-C08 | χ_G is a sufficient stability certificate only under full analytic/small-gain hypotheses. | CONJECTURE / PARTIAL | No proof of all right-half-plane and infinity assumptions is claimed here. | Finite sampled jω values do not establish the theorem. | PARTIAL (numerical gate: INCONCLUSIVE) |
