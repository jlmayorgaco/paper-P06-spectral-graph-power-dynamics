# Claim ledger for a possible poster

## C1 — Reproduced baseline

- **CLAIM:** The frozen 87.5% GFL seed reproduces the equilibrium, 203-pole physical spectrum, and all five zero-delay events.
- **STATUS:** NUMERICALLY_VALIDATED.
- **MODEL SCOPE:** Repository IEEE-39 fixed-interior PowerDynamics model, declared frozen limits.
- **ASSUMPTIONS:** Original generator dispatch defines the MW denominator.
- **EVIDENCE FILE:** `M0_BASELINE_REPRODUCTION.csv`, `M0_EVENT_REPRODUCTION.csv`.
- **NUMERICAL VALUE:** 4727.415954 MW GFL; 675.345136 MW retained SG; critical real part −0.0808091 s⁻¹.
- **WHAT IT DOES NOT PROVE:** Maximum replacement or delayed security.

## C2 — Best point on frozen zero-delay search

- **CLAIM:** A preregistered interpolation search found a co-designed point with 88.4551408% GFL passing all five frozen zero-delay events.
- **STATUS:** NUMERICALLY_VALIDATED, BEST_FOUND_ONLY.
- **MODEL SCOPE:** Same nonlinear ODE/DAE event model and limits as C1, with zero PLL delay.
- **ASSUMPTIONS:** Frozen candidate path, bounds and event metrics.
- **EVIDENCE FILE:** `M1_CANDIDATE_AUDIT.csv`, `M1_ZERO_DELAY_EVENTS.csv`, `M1_ZERO_DELAY_DESIGN.toml`.
- **NUMERICAL VALUE:** 4779.019928 MW GFL; 623.741162 MW retained SG; +51.603975 MW over seed; worst frequency deviation 0.485513 Hz.
- **WHAT IT DOES NOT PROVE:** Local or global optimality, positive-delay feasibility, hard current-limiter safety.

## C3 — Exact action-space upper rank

- **CLAIM:** At a fixed equilibrium and characteristic frequency, variations in `rho,Kp,Ki` change the 282-dimensional descriptor through an algebraic factor of rank at most 30.
- **STATUS:** EXACT_IDENTITY for the declared descriptor factorization; NUMERICALLY_VALIDATED on 21 model-accepted designs times three heterogeneous delay fields.
- **MODEL SCOPE:** Fixed model architecture and identical network topology; equilibrium is recomputed with design changes.
- **ASSUMPTIONS:** Model descriptor and regular reference pencil at tested `s`.
- **EVIDENCE FILE:** `TABLE_Q06_MASTER_ACTION_SPACE.csv`, `M2_ACTION_SPACE_RANK.csv`, `M2_RECONSTRUCTION_VALIDATION.csv`, `M2_DETERMINANT_IDENTITY.csv`.
- **NUMERICAL VALUE:** Max reconstruction relative error `3.23e-16`; max determinant log-magnitude error `1.20e-12`; numerical difference rank 30.
- **WHAT IT DOES NOT PROVE:** 30-dimensional globally solvable optimization, performance improvement, or exact minimal rank at every parameter/frequency point. The 20 out-of-box draws were rejected before descriptor assembly.

## C4 — Uniform-delay PLL instability cluster

- **CLAIM:** For the two fixed zero-delay designs, numerical contour counting of the exact delayed characteristic gives zero roots to the right of `−0.05 s⁻¹` at 20 and 30 ms; at 40 ms it gives 6 roots for the seed and 12 for the higher-replacement design.
- **STATUS:** NUMERICALLY_VALIDATED for the stated points; not interval-certified.
- **MODEL SCOPE:** Linear DDE about each solved equilibrium with uniform exogenous PLL measurement delay.
- **ASSUMPTIONS:** Simple bounded contour, Float64 analytic characteristic, no root on contour.
- **EVIDENCE FILE:** `M3_DDE_ORACLE.csv`, `M3_ALL_REFINED_ROOTS.csv`, `M3_ORACLE_METHOD.md`.
- **NUMERICAL VALUE:** Nine distinct positive-imaginary roots refined at 40 ms across both designs, with conjugates matching counts; frequencies approximately 5.18–5.40 Hz.
- **WHAT IT DOES NOT PROVE:** Delayed nonlinear event feasibility, a delay-specific maximum replacement frontier, heterogeneous-delay effect, or a certified global root count.

## C5 — Conditional PI target-root law

- **CLAIM:** A fixed-bus prescribed characteristic root admits a 2 by 2 real conditional PI gain equation derived from the exact rank-one determinant lemma.
- **STATUS:** EXACT_IDENTITY under regularity assumptions; NUMERICALLY_VALIDATED target residuals.
- **MODEL SCOPE:** Single chosen PLL gain pair varied about a frozen design/target root.
- **ASSUMPTIONS:** Reference pencil regular and 2 by 2 gain matrix nonsingular.
- **EVIDENCE FILE:** `M5_CLOSED_FORM_PI_VALIDATION.csv`, `THEORY_ANALYTICAL_CODESIGN.md` §5.
- **NUMERICAL VALUE:** 45/45 target-root residuals ≤`1.02e-16`; only 5/45 gain pairs satisfy the frozen bounds.
- **WHAT IT DOES NOT PROVE:** A target is the rightmost pole, globally stabilizing gains, or a usable full co-design law.

## C6 — Singular zero-frequency Schur route

- **CLAIM:** The proposed zero-frequency hidden-state Schur damping coefficient is not defined for the tested gauge-deflated IEEE-39 model.
- **STATUS:** NEGATIVE_RESULT.
- **MODEL SCOPE:** Tested zero-frequency retained/hidden partition of the current model.
- **ASSUMPTIONS:** Current partition and gauge deflation.
- **EVIDENCE FILE:** `M0_ZERO_FREQUENCY_SCHUR_REPRODUCTION.csv`.
- **NUMERICAL VALUE:** Hidden block condition `3.095842e18`.
- **WHAT IT DOES NOT PROVE:** Impossibility of a different finite-frequency or port-preserving graph reduction.

## C7 — Reduced LP partial-retention theorem

- **CLAIM:** An extreme optimum of the reduced one-mode/one-RoCoF retention LP has at most two partially retained SGs.
- **STATUS:** PROVED_REDUCED_MODEL.
- **MODEL SCOPE:** The LP in `THEORY_ANALYTICAL_CODESIGN.md` §9–10 only.
- **ASSUMPTIONS:** Two independent active aggregate inequalities plus box constraints.
- **EVIDENCE FILE:** `THEORY_ANALYTICAL_CODESIGN.md` §10.
- **NUMERICAL VALUE:** At most two fractional components in an extreme LP solution.
- **WHAT IT DOES NOT PROVE:** The structure of a full nonlinear delayed multi-event optimum; no physical LP coefficients were validated here.

**Poster decision:** The chain from a delay-aware analytical law to a fully validated replacement frontier is incomplete. C3 and C4 are real research leads, but the proposed award poster headline is not yet supported. Do not present M1 as `P_GFL^max(tau)` or the M3 fixed-design counts as a frontier.
