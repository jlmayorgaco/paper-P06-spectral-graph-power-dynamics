# TX3 theory freeze for C1/C2

**Freeze ID:** `TX3-TF-C1-C2-1.0`  
**Pre-pilot freeze date:** 2026-08-27  
**E02C aggregate model hash:** `0f40a8b33e7b93e16a7d8a9f5d1b7ed9f065a64f52c074591890fae9cea6e072`

For normalized physical actions `alpha_i`, each operator evaluation applies the physical parameters, solves the AC power flow, initializes the complete ANDES DAE, and rebuilds `T(s,alpha)=sE(alpha)-J(alpha)`. Frozen-linearization results are diagnostic only.

The primary functional is

`Phi(alpha)=integral w(omega)[log|det T(s,alpha)|-log|det T(s,0)|]domega`,

with the primary/secondary bands and action-independent contour chosen only by the deterministic conditioning rule in `FREQUENCY_BAND_FREEZE.json`. The weight is `1/[omega log(omega_max/omega_min)]`.

Direct pair/triple externalities are the finite Möbius differences over independently re-equilibrated vertices. Connected derivatives use the validated set-partition formula from TX3-TF-0.2 and E04B. Pairs are split into `Re tr(GT_ij)` curvature and `-Re tr(GT_iGT_j)` feedback. Triples are split into intrinsic `T_ijk`, mixed pair-vertex/singleton, and pure singleton-feedback classes. Sparse factorized solves are mandatory; `T^-1` is never formed.

Mixed operator derivatives are evaluated by tensor polynomial differentiation on interior Gauss-Legendre nodes. Three central-difference radii (0.04, 0.02, 0.01) and action-grid orders 2/3 are compared in the pilot. The first production rule satisfying the predeclared convergence target (connected q2/q3 discrepancy no more than 20% of the pilot direct significance floor and derivative refinement materially decreasing) is frozen in `NUMERICAL_METHOD_FREEZE.json`. No setting is selected from effect desirability.

Uncertainty is the conservative sum of frequency, derivative/action-grid, hypercube, sparse-logdet, and PF/Jacobian components. C1 uses `|Delta Phi| >= 5 epsilon_direct`; smaller results are UNRESOLVED. C2 uses `|Delta_connected-Delta_direct| <= 5 epsilon_combined`, with the population rules frozen in `C1_C2_GATE_FREEZE.json`. Correlation is descriptive only.

The finite-cube bridge is an exact identity. E04 therefore verifies numerical closure, while physical evidence additionally requires finite AC re-equilibration, externalities above uncertainty, mechanistic contribution classes, and holdout replication. C3–C6, cycle localization, surgery, causality, and EMT transfer are not assessed.
