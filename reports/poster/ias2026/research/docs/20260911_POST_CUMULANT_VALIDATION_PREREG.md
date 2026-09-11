# Preregistration: final post-cumulant validation campaign (2026-09-11)

This document is committed **alone**, before any new numerical result is
computed or inspected. The scientific state is fixed at `cdfcefcf`, with the
dossier at `a0aa831b`.

- Frozen results are only read, never rewritten. New outputs go to new files:
  `results/PCV/…`, plus the named deliverables.
- There is no search for a stronger example and no retuning.
- Thresholds are not changed after results are seen. Any deviation is recorded
  in a "Deviations" section of the report.
- Connected cumulants are a secondary diagnostic only.
- Environment rule: `OPENBLAS_NUM_THREADS = OMP_NUM_THREADS = MKL_NUM_THREADS = 1`
  is set in the launching shell before any multiprocessing run.
- Unrelated processes, such as the pmu_hybrid job, are not touched.

## 0. Thesis under test

Local device changes lead to exact network-mediated collective closure. That
closure produces policy-dependent minimal failing coalitions, which make
replacement portfolios non-composable. The claimed payoff is exact diagnosis,
planning and targeted design.

The campaign tests whether this survives:

- fair baselines;
- a clean one-coordinate policy counterfactual;
- a lower-order hierarchy;
- targeted remediation;
- targeted robustness;
- an independent reproduction of the custom GFL.

## 1. Hypotheses

| id | hypothesis | decided by |
|---|---|---|
| H1 | At P4 the complete portfolio H4 = {30,33,35,37} is unstable while every proper subset is stable (κ = 4, 𝓗 = {H4}). | B0 recomputation (§4) |
| H2 | Lower-order portfolio screening (B3 modal sensitivity, B4 additive, B5 pairwise) approves H4 at P4. | Task B/C at P4 |
| H3 | Aggregate MW/MVA and static grid strength (SCR, gSCR, MIIF, electrical distance) cannot distinguish two policy points with identical network, locations, ratings and scheduled P/Q whose dynamic verdicts differ. | Task D, §6 |
| H4 | Along the clean converter-only line (k = 1.425, t = 1.5, h = 1), only g changes, and the verdict of H4 changes at g* ≈ 0.20768. | §6 |
| H5 | The exact port representation predicts the local direction and the rank of how physical/control coordinates move the actual boundary. | Task E, §8 |
| H6 | Non-composability persists under the preregistered envelopes (R1); the exact witness H4 persists (R2); κ = 4 persists (R3); the g* boundary direction persists (R4). | §9 |
| H7 | An equation-equivalent second-simulator implementation of the custom GFL reproduces the 32 frozen verdicts, 𝓗 and κ. | §10 |

## 2. Frozen policy points and portfolio sets

| name | θ = (g, k, t, h) | role |
|---|---|---|
| P4 | (0.03625, 1.425, 1.5, 1) | canonical unstable point |
| G_S (primary clean stable point) | (0.25, 1.425, 1.5, 1) | frozen FC03 path point on the stable side: 𝓗 = ∅, frozen α(H4) = −0.017385 |
| G_S2 (secondary clean stable point) | (1.0, 1.425, 1.5, 1) | frozen FC03 path point: 𝓗 = ∅, frozen α(H4) = −0.101629 |
| P_inf (secondary, multi-coordinate) | (1, 0.5, 1.5, 1) | changes g **and** k; never used as a one-variable causal comparison |
| E12 census policy | solve_case defaults (fixed Q, k = t = 1) | census of 9 candidates, 512 portfolios |
| F10 points | the 52 frozen points of `results/F10/F10_baselines.csv` (rows not in the 36-point line) | baseline recomputation |
| F10 line | the 36 frozen F7B points t = 0.8515625, g = 0 … 0.35 | secondary (Task B/D) |

**Portfolio sets.**

- The 16 subsets of H4 at P4, G_S, G_S2 and P_inf.
- The 16 subsets of H4 at each of the 52 + 36 F10 points.
- Census portfolios of size 4 (126, primary for Task A), and sizes 5 (126) and
  6 (84) as secondary.

## 3. Truth labels (B0 oracle)

- **Core lattices (P4, G_S, G_S2, P_inf, the F10 points and line).** The FC01
  part-3 direct path, recomputed:
  - `_f7_common.solve_subset(members, Theta(g,k,t,h))`;
  - `certification.physical.physical_matrices` (Jacobians at steps h and 2h; the
    error matrix is their difference);
  - `rotation_generator`, `frequency_partner`, `transverse_operator`;
  - `classify_spectrum(A_perp, Zᵀ dA Z, SAFETY)`.
  - Labels are STABLE, UNSTABLE or BOUNDARY_OR_UNRESOLVED. 𝓗 and κ come from
    FC01 `h_of`.
  - α_perp is the largest real part of the transverse spectrum. The critical
    frequency is the imaginary part (Hz) of the rightmost 0.3–1.5 Hz transverse
    eigenvalue.
- **Census (E12 policy).** The frozen transverse labels and α_perp of
  `FC10_census_transverse.csv` (post re-audit). They are not recomputed.
- Unresolved points are reported and excluded from scoring; their count is
  reported.

## 4. Baselines (fixed definitions; no tuning)

| id | method | inputs | needs critical frequency? | needs the full portfolio evaluated? | output |
|---|---|---|---|---|---|
| B0 | full transverse eigenanalysis of every portfolio (oracle, not a competitor) | full DAE per portfolio | no | yes (all 2^m) | α_perp, verdict, 𝓗 |
| B1 | displaced active dispatch Pg [MW]; replaced rating Sn [MVA]; penetration = Pg / total system Pg | dispatch, ratings | no | no | scalar score |
| B2 | min nodal SCR and gSCR (repository `diagnosis.baselines`: short-circuit network with machines behind x″d, loads omitted, replaced machines excluded; gSCR = min eigenvalue of the rating-normalized Kron-reduced susceptance, labelled an approximation of the gSCR literature); max off-diagonal MIIF = abs(Z_ij)/abs(Z_jj) | network, ratings, PF voltages | no | no | scalar score |
| B3 | first-order modal sensitivity (F10 B5): dλ/dρ_a of the base least-damped 0.3–1.5 Hz mode by a 1 % partial replacement at each candidate; λ̂(S) = λ0 + Σ_{a∈S} dλ_a; verdict unstable iff Re λ̂ > 0 and 0.3 ≤ Im λ̂/2π ≤ 1.5 Hz. Correction relative to F10: the full θ including h is used (F10 ignored h). | 1 + m full solves | no | no | α̂, verdict |
| B4 | additive single-device extrapolation: order-1 Möbius truncation of the exact transverse α_perp | exact α of ∅ and singles | no | no | α̂, verdict |
| B5 | pairwise truncation: order-2 Möbius truncation of α_perp; exact for \|S\| ≤ 2 | exact α of subsets of size ≤ 2 | no | no | α̂, verdict |
| B5b | third-order truncation: order-3 Möbius truncation; exact for \|S\| ≤ 3 | exact α of subsets of size ≤ 3 (15 of 16 core members) | no | no | α̂, verdict |
| B6 | participation factors of the base least-damped band mode, summed over the (δ, ω) states of each candidate SG | 1 base eigen-analysis | no | no | ranking (localization only, no portfolio classifier) |
| B7 | electrical compactness: mean over i < j in S of d_ij = abs(Z_ii + Z_jj − Z_ij − Z_ji), with Z the inverse of the B2 short-circuit Ybus excluding S; score = −mean d | network | no | no | scalar score |
| B8a | closure distance, **non-oracle**: min over 121 uniform frequencies in 0.3–1.5 Hz of min_i abs(λ_i(Q_SS(jω)) + 1) | port operator of S (C2 realization) | no | no (no eigen-analysis of S) | scalar |
| B8b | closure distance, **oracle-assisted**: the same at jω*, with ω* the imaginary part of S's band-critical eigenvalue | port operator plus S's eigenvalue | **yes** | yes | scalar |
| B9a / B9b | abs(chi_S) at the B8a minimizing frequency (non-oracle) / at jω* (oracle-assisted); secondary diagnostic | port operator | no / yes | no / yes | scalar |
| FW | our framework: exact closure on every principal sub-loop (GN) plus the hypergraph | port operator plus all principal minors | no | no eigen-analysis per portfolio; one base and one all-replaced operator per point (F11) | exact labels (proved equal to B0; F11 0 count errors); not re-run as a competitor |

**Orientations (higher = riskier).**

- Scores used as they are: Pg, Sn, penetration, B3/B4/B5/B5b α̂, max MIIF, B7
  score, abs(chi).
- Scores negated: gSCR, min SCR, closure distance.

## 5. Tasks and metrics

- **Task A: risk ranking.**
  - Metrics: ROC-AUC (unstable = truth UNSTABLE), average precision, and
    Spearman with α_perp, plus permutation p (10 000 permutations, seed
    20260916).
  - Primary: census size 4 (n = 126).
  - Secondary: census sizes 5 and 6; F10 52 points pooled over the 15 non-empty
    subsets.
- **Task B: minimum blocking-set recovery.** Applies to B3, B4, B5 and B5b (B0 is
  the reference).
  - Metrics:
    - exact 𝓗;
    - exact κ;
    - cardinality error abs(κ̂ − κ), with a finite-versus-∞ mismatch counted
      separately;
    - false blocking units (units in the union of predicted hyperedges but in
      no true one);
    - missed blocking units.
  - Datasets: P4, G_S, G_S2, P_inf, the F10 52 points, the F10 line, and the
    census.
  - B1, B2, B7, B8 and B9 are N/A (no stability threshold defined). B6 is N/A as
    a classifier.
- **Task C: stability/margin prediction** (B3, B4, B5, B5b).
  - Metrics: sign accuracy, MAE and RMSE of α̂ − α_perp over the non-empty
    subsets.
  - Datasets: P4, G_S, P_inf lattices; F10 points; census.
- **Task D: policy discrimination.**
  - (i) Clean pair P4 versus G_S for H4: does the method's verdict or score
    change, and is the change correct?
  - (ii) F10 52 points: within each subset S, the AUC of the score over policy
    points against the truth label, averaged over the subsets that have both
    classes. Static methods are constant, so their within-portfolio AUC is 0.5
    by construction; this is reported, not scored as failure.
- **Task E: intervention ranking at P4** (complete portfolio).
  - **Lines.** The 12 frozen holdout lines (seed 20260911). The actual effect is
    exact Δα_perp(H4) at γ = 1.5, the primary finite intervention; γ = 1.25 and
    2.0 are also reported. The C2 realization is re-solved with the scaled
    network.
    - Predictors:
      - (a) port derivative ds*/dγ at the P4 critical eigenvalue × 0.5;
      - (b) conventional full-DAE eigenvalue sensitivity (central difference at
        γ = 1 ± 0.002) × 0.5;
      - (c) static: ΔgSCR(H4) at γ = 1.5;
      - (d) static: 1/x of the line;
      - (e) static: x-weighted betweenness (frozen F2 file).
    - Metrics: Spearman and Kendall between the predicted stabilization and the
      actual stabilization −Δα; top-1 agreement.
  - **Policy coordinates g, k, t, h.** Finite step 0.10 × SCALE_a (SCALE = 1,
    1.8, 2.5, 2), both directions, exact Δα. The actual stabilizing effect is
    the more negative Δα. The prediction is |port Re ds*/da| × step. Spearman and
    Kendall over the 4 coordinates are reported and flagged as n = 4.

## 6. Clean policy counterfactual

- Compare P4 against G_S (primary) and G_S2 (secondary). Only g changes.
- Quantities that must be identical: network Ybus (hash), candidate locations,
  SG/GFL ratings, scheduled P and Q (PF solution, hash), k, t, h, and every B1,
  B2 and B7 value.
- Quantities compared: α_perp(H4), critical frequency, 𝓗, κ, RHP count, and B3,
  B4, B5, B8 and B9.
- The direct bisection of g* on the line (H4 only; bracket [0.2025, 0.2256],
  40 steps on the B0 label) must reproduce the frozen 0.20768 within 1e-4
  (stopping rule S2).
- P4 versus P_inf is reported as a **multi-coordinate** comparison only.

## 7. Lower-order hierarchy (P4 and G_S; complete portfolio)

**A. α space.**

- Levels: B3 first-order modal sensitivity, B4 additive, B5 pairwise, B5b
  third-order, exact.
- Reported per level: predicted α̂(H4), sign, residual to the exact value, and
  whether H4 is flagged.

**B. Characteristic-function space.** F(s) = Π_i det(I + M_ii(s)) · T_r(s), with
T_r the truncation of det(I + Q_SS) to Boolean orders ≤ r:

- r = 1 is local factors only;
- r = 2, 3;
- r = 4 is exact.

For each r, the predicted critical eigenvalue is the zero of F_r nearest the
exact critical eigenvalue:

- Newton from the exact λ, complex-step slope h = 1e-5, radius 1.0 s⁻¹;
- "no zero" if the search leaves the radius.

Reported: sign and Re of that zero; the residual abs(T_r(λ_exact)); whether
instability is predicted.

The report must state that κ = 4 is a minimum destabilizing cardinality, not an
irreducible four-device connected interaction. CC03 found the flagship composite
(ν = 0.069; deletion displacement 0.024 s⁻¹).

## 8. Remediation (P4 complete portfolio; documented families only)

| intervention | actual finite result | prediction | boundary |
|---|---|---|---|
| Q/V gain g | frozen continuation g* = 0.20768 (FC18); new α at G_S | port derivative at P4 (frozen FC18: one step 0.0518, Newton 2e-10) | direct bisection (§6) |
| AVR gain scale k (at g = 0.03625, t = 1.5) | frozen F7A bracket k ∈ [1.303125, 1.30625] | port derivative dRe s/dk at P4 | port-guided Newton on k (≤ 8 steps, tolerance abs(Re s) ≤ 1e-7) against a direct bisection in [1.30, 1.31] (40 steps) |
| damped condenser | frozen threshold 2.48 % (106.0 MVA); frozen TDS 2 % unstable / 3 % stable | N/A (not in the port family) | frozen |
| documented governors | frozen FC03: α −0.0745 | N/A (model change) | — |
| lines (12 frozen) | exact α at γ = 1.25, 1.5, 2.0 | port derivative and full-DAE sensitivity | whether any single line up to 2× restores stability |

No common economic cost is defined, and abs(chi) suppression is not an
objective.

## 9. Targeted robustness (deterministic stress-test envelopes)

The fractions below are **coverage fractions over declared envelopes, not
probabilities**.

| envelope | source | parameters (multiplicative factors) | sampling |
|---|---|---|---|
| EM-f (fleet-wide machines) | the documented E37 project envelope | inertia M [0.80, 1.20]; x′d and x′q (common) [0.90, 1.10]; K_A [0.85, 1.15]; T_E [0.85, 1.15]; PSS gain [0.80, 1.20]; PSS time constants T4, T5, T6 (common) [0.80, 1.20] | one factor per group per draw, applied to all SGs (E37 style) |
| EM-u (per-unit machines) | same bounds | same groups | an independent factor per SG per group |
| EC (converters) | **stress test; no project source** | PLL (kp, ki common) [0.80, 1.20]; outer loops (kp_p, ki_p, kp_q, ki_q common) [0.80, 1.20]; current loop (kp_i, ki_i common) [0.80, 1.20]; τ_p [0.80, 1.20]; x_f [0.90, 1.10]; Q/V gains excluded (they are the policy coordinate) | an independent factor per converter per group |
| EMC | EM-u together with EC | — | combined |

**Draws and evaluation.**

- Latin hypercube, N = 100 per envelope, seed 20260917 (+ envelope index).
- One-at-a-time extremes: every group at its lower and upper bound, others
  nominal: 6 × 2 machine and 5 × 2 converter cases.
- Policy fixed at P4. Network, loads, dispatch, D = 0 and "no governor" are fixed.
  Governor dependence is already frozen in FC03 and not resampled.
- Loads and the network are excluded: they change the operating point and the
  matched-dispatch semantics. This is stated as a limitation.
- Per draw:
  - all 16 subsets at P4 through the direct path;
  - FC01 logic: sign count if min abs(Re) ≥ 2e-3, otherwise the classifier with
    the h/2h error matrix;
  - 𝓗, κ, α(H4) and its critical frequency;
  - g* by bisection on g ∈ [0.03625, 1.0] at k = 1.425 (25 steps, H4 only), only
    if H4 is unstable at P4 and stable at g = 1; otherwise "no crossing".

**Metrics.**

- R1: fraction of draws with non-composability, i.e. a stable base, all singles
  stable, and some unstable portfolio (κ ∈ {2, 3, 4}).
- R2: fraction with 𝓗 = {H4}.
- R3: the κ distribution (1, 2, 3, 4, ∞, base-unstable).
- R4: g* median and IQR; "direction robust" if α(H4, g = 1) < α(H4, P4) in
  ≥ 90 % of draws.

**Decision per envelope.** PASS if the fraction is ≥ 0.80, PARTIAL for 0.50–0.80,
FAIL below 0.50. Applied separately to R1 and R2.

## 10. Independent custom-GFL reproduction (Phase 8)

- **Environment.** A separate ANDES working copy (from vendor tag v2.0.0,
  `eda5163c`) and a separate venv `.venv/xtool-andes-gfl`, with an isolated
  generated-code cache. `.venv/tx3-andes` and `vendor/andes` are not modified.
- **Specification first.** `docs/20260911_GFL_REPRODUCTION_SPEC.md` maps every
  equation, base, sign convention and initialization. It is committed before any
  comparison. If exact equivalence of any block cannot be achieved, Phase 8
  stops and the reason is documented. No library GFL is substituted.
- **Gate 0, device level.** For 20 random states and terminal voltages per
  device, max abs(f_ANDES − f_internal) / max abs(f) ≤ 1e-9 for the GFL (and the
  PSS if custom). If this fails, STOP.
- **Validation set.** The 16 subsets of H4 at P4 and at G_S (32 cases). P_inf is
  optional and secondary.
- **ANDES transverse spectrum.** Remove the two eigenvalues of smallest modulus.
  They must have modulus < 1e-3, and no other eigenvalue may have modulus below
  1e-2 (checked, reported).
- **Primary tolerances.**
  - equilibrium residual ≤ 1e-6;
  - abs(Δα_perp) ≤ 1e-3 s⁻¹;
  - abs(Δ critical frequency) ≤ 1e-3 Hz;
  - RHP counts equal in 32/32;
  - verdicts equal in 32/32;
  - 𝓗 and κ equal at P4 and at G_S.
- **Secondary tolerance.** Band nearest-mode ≤ 1e-4.
- **Stopping.** Any verdict or 𝓗 mismatch, or abs(Δα) > 1e-2, stops the stronger
  claims. The first mismatching equation or convention is then identified.
  Nothing is retuned.

## 11. Stopping rules

- **S1.** If the recomputed B0 α_perp of any P4, P_inf or G_S subset differs from
  the frozen FC03 value by more than 1e-6, STOP and diagnose before any baseline
  scoring.
- **S2.** If the recomputed g* differs from 0.20768 by more than 1e-4, STOP the
  Phase 4 and Phase 6 claims.
- **S3.** More than 10 % unresolved F10 points is reported; it does not stop the
  campaign.
- **S4.** A Phase 8 Gate 0 failure stops Phase 8.
- **S5.** No threshold, envelope or dataset changes after results. Deviations are
  disclosed.

## 12. Scientific gates (decided at the end)

- **GATE 1: the P4 case survives a fair comparison.** PASS iff:
  - H1 is confirmed by B0;
  - B3, B4 and B5 approve H4 at P4;
  - B1, B2 and B7 are unchanged across the clean counterfactual while the
    verdict changes.
- **GATE 2: additional information beyond enumerating B0.** PASS iff at least
  three are demonstrated quantitatively:
  - (a) exact reduced representation (identity ≤ 1e-10 at the tested points);
  - (b) analytic boundary sensitivity correct in sign for every intervention
    coordinate, with Newton reaching the direct boundary within 1e-6;
  - (c) Spearman(port-predicted, actual) ≥ 0.9 for the lines;
  - (d) policy mapping and any-order planning derived from 𝓗 (proved and
    frozen);
  - (e) structured attribution (principal-minor anatomy, frozen FC18).

  Never claimed: higher stability accuracy than eigen-analysis, or polynomial
  complexity.
- **GATES 3 and 4.** R1 and R2 per the §9 rule.
- **GATE 5.** Phase 8 primary tolerances met. Otherwise FAIL, or NOT RUN with the
  reason.
- **GATE 6.** The strongest defensible practical claim (a narrative decision).

## 13. Deliverables and commits

**Deliverables.**

- `docs/20260911_PORTFOLIO_DECISION_BENCHMARK.md`
- `results/20260911_BASELINE_COMPARISON.csv`
- `results/20260911_POLICY_COUNTERFACTUAL.csv`
- `results/20260911_ROBUSTNESS.csv`
- `docs/20260911_GFL_REPRODUCTION_SPEC.md`, `docs/20260911_GFL_REPRODUCTION.md`
  and `results/20260911_GFL_REPRODUCTION.csv`
- `docs/20260911_FINAL_VALIDATION_LEDGER.md` and
  `results/20260911_FINAL_VALIDATION_MATRIX.csv`
- the author decision memo
- `figures/20260911_portfolio_decision_case.*` (4 panels, no cumulant panel)

**Commits** are local, in this order:

1. preregistration
2. baseline recomputation
3. counterfactual and hierarchy
4. remediation
5. robustness
6. GFL specification
7. GFL reproduction
8. ledger, figure and memo

There is no push. Deterministic stages are rerun and byte identity is checked.
