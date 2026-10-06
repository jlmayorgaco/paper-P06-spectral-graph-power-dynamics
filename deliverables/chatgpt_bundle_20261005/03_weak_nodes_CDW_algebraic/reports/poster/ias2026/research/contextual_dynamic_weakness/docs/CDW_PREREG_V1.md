# CDW preregistration V1 — frozen before any CDW numerical result

Date: 2026-09-12. Branch `research/contextual-dynamic-weakness`, based on tag
`TX4_FINAL_MANUSCRIPT_FREEZE` (69f200df). This file is committed **before** any
CDW scientific computation.

**Only pre-commit computation.** A timing-only pilot ran three `solve_case`
calls and printed wall times, about 0.35 s per portfolio evaluation. It printed
no α or verdict.

**Change rule.** After commit, no hypothesis, threshold, seed, test set, baseline
or definition may change. Any deviation is recorded in
`docs/CDW_PREREG_V1_DEVIATIONS.md` with a timestamp and a reason, and is never
applied silently.

Definitions are in `docs/CDW_THEORY_V1.md`; the model is in
`theory/CDW_GRAPH_DAE_V1.md`.

## 0. Frozen inputs

| input | value |
|---|---|
| model | TX4 frozen IEEE-39 (`configs/ias2026/ieee39_network.json`, source sha256 `9c2048dc…`); TX4 `solve_case`, device models, classifier and transverse operator, unchanged |
| converter | TX4 GFL with voltage_control=True, leak 0.05 rad/s, policy gain `g` |
| census candidates | V9 = {30,31,32,33,34,35,36,37,38}: all 512 portfolios |
| core | V4 = {30,33,35,37}: 16 portfolios |
| holdout policies | `results/prereg_inputs/holdout_policies.json`, sha256 `4ac8c95d…023b349`, 24 points H01–H24 (LHS, seed 20260920) |
| CDW envelope draws | `results/prereg_inputs/cdw_envelope_draws.json`, sha256 `864a5969…0a681abf`: EM-f, EM-u, EC, EMC with N = 40 each, seed 20260921. TX4 group bounds; converter groups at all nine candidates |
| TX4 envelope draws (discovery) | regenerated with the frozen PCV05 generator (seed 20260917, N = 100 per envelope). They must reproduce the frozen `results/PCV/PCV05/PCV05_draws.csv` factors exactly, or they are not used. Converters outside V4 are nominal |
| numerics | α_⊥ from the TX4 transverse operator; status from the TX4 classifier; BLAS pinned to one thread |

### 0.1 Discovery policies D01–D15, as (g, k, t, h)

| id | point | role |
|---|---|---|
| D01 | (0.03625, 1.425, 1.5, 1) | P4 |
| D02 | (1.0, 0.5, 1.5, 1) | P_inf (multi-coordinate reference) |
| D03 | (0.25, 1.425, 1.5, 1) | G_S (clean g-only stable point) |
| D04 | (0.0, 1.425, 1.5, 1) | g-only line |
| D05 | (0.1, 1.425, 1.5, 1) | g-only line |
| D06 | (0.5, 1.425, 1.5, 1) | g-only line |
| D07 | (1.0, 1.425, 1.5, 1) | g-only line (G_S2) |
| D08 | (0.1, 0.75, 1.5, 1) | F7A |
| D09 | (0.5, 1.75, 1.5, 1) | F7A |
| D10 | (0.9, 2.2, 1.5, 1) | F7A |
| D11 | (0.2, 1.0, 0.8515625, 1) | F7B, frozen path t |
| D12 | (0.6, 1.0, 0.8515625, 1) | F7B, frozen path t |
| D13 | (0.3, 1.0, 2.5, 1) | F7B |
| D14 | (0.2, 1.0, 1.0, 0.5) | F7C |
| D15 | (0.7, 1.0, 1.0, 1.5) | F7C |

**Holdout policies H01–H24** (g, k, t, h), frozen in the JSON:
H01 (0.18945, 1.4747, 0.6878, 1.1064); H02 (0.02346, 0.5976, 1.3129, 0.2088);
H03 (0.01530, 1.9434, 0.5290, 0.9935); H04 (0.23126, 1.3305, 2.9979, 1.3093);
H05 (0.00096, 0.7405, 2.4917, 1.5652); H06 (0.80414, 1.7659, 1.2132, 0.4136);
H07 (0.46847, 0.6713, 0.9122, 0.7570); H08 (0.34716, 0.5183, 1.7816, 0.0308);
H09 (0.93736, 0.9551, 2.4472, 1.7281); H10 (0.65284, 1.2584, 1.6821, 1.8266);
H11 (0.76549, 1.4806, 2.8248, 1.4581); H12 (0.03915, 1.0797, 2.0554, 1.6502);
H13 (0.00433, 0.8740, 1.8800, 1.2194); H14 (0.28964, 0.9011, 1.4124, 1.0652);
H15 (0.07751, 1.6134, 1.1162, 1.9559); H16 (0.11904, 2.2984, 0.7348, 0.6345);
H17 (0.89553, 1.8045, 0.9595, 0.2815); H18 (0.51711, 1.6789, 2.5846, 0.6913);
H19 (0.16423, 2.1013, 2.7419, 0.1185); H20 (0.31479, 2.1759, 1.5332, 0.5580);
H21 (0.62141, 1.1387, 2.3380, 1.9161); H22 (0.39196, 2.0478, 1.6299, 1.3697);
H23 (0.10462, 1.2460, 2.2078, 0.8480); H24 (0.05137, 1.8697, 2.0668, 0.4597).

**Split rule.**
- Discovery points may be used to explore, and to select instances such as
  corridor candidates or the baseline to carry forward.
- Every GOLD gate is evaluated on the **holdout** policies and the **CDW-seed**
  envelope draws only.
- A policy with an unstable base (α(∅) ≥ 0) is kept in the data and excluded
  from the gate denominators. The count is reported.

## 1. Thresholds (with justification; not copied from TX4)

| symbol | value | justification |
|---|---|---|
| τ_mat | 0.01 s⁻¹ | A 0.01 s⁻¹ change of the real part is about 0.2 % damping ratio at 0.7 Hz, and about 8 % of the P4 margin of 0.127. It is well above the numerical error of α (≤ 1e-7). It is the smallest change a planner would act on |
| τ_res | 1e-3 s⁻¹ | Resolution only. It is ten times the classifier's own ambiguity band, 1e-4. It is used for sensitivity analyses and violation counts |
| MAC_min, Δf_max | 0.80, 0.15 Hz | Mode matching on bus-voltage shapes. MAC 0.8 is the conventional "same mode" threshold. 0.15 Hz is less than half the spacing (≥ 0.26 Hz) between the inter-area families seen in the frozen modal library |
| stable context | α(S) < 0 and status STABLE | The planner's actual decision situation |
| coverage bands | PASS ≥ 0.75; PARTIAL 0.50–0.75; FAIL < 0.50 | "Typical" means present in at least three quarters of the envelope; "majority" means at least one half. These are coverage fractions, **not probabilities** |
| agreement ceiling | 0.90 | A node-only ranking that orders ≥ 90 % of resolved pairs correctly would be operationally adequate. Insufficiency must be shown below that ceiling |
| decision-failure rate | ≥ 10 % of stable contexts | One in ten next-replacement decisions materially wrong (regret ≥ τ_mat) is a planning-relevant failure rate |
| ranking gain | Δρ ≥ 0.20 | The difference between a moderate and a strong rank correlation. Smaller gains are within the variation across conditions |
| strong-rank bar | ρ ≥ 0.70 | Conventional "strong" monotone association |

## 2. Hypotheses, experiments, metrics and gates

### E1 — contextual sign-reversal census (H1, H1b)

- **Inputs.** Policies D01–D15 and H01–H24. At each policy, all 512 V9
  portfolios are solved (SPR, the TX4 `solve_case`).
- **Recorded per portfolio.**
  - status, α, the rightmost eigenvalue and its frequency;
  - the bus-voltage mode shapes of the transverse modes in 0.1–2.0 Hz;
  - infeasibility, if any.
- **Marginals.** `Δ_i α(S)` for all 2304 pairs (S, i ∉ S). Tracked marginals
  follow C2 (MAC_min, Δf_max).
- **Per intervention and policy.**
  - STAB, DESTAB and NEUTRAL fractions;
  - min, median and max;
  - variance, sign entropy (three classes, log₂);
  - the number of contexts in each class;
  - the global reversal indicator, the stable-context reversal indicator and
    the mode-tracked reversal indicator;
  - the strongest reversal pair, i.e. the pair maximizing
    `min(|Δ(S1)|, |Δ(S2)|)`.
- **Node-only ranking test (H1b).** For each base-stable policy:
  - Resolved triples are (S stable, i, j ∉ S) with `|Δ_i − Δ_j| ≥ τ_mat`.
  - The optimal fixed ranking maximizes pairwise agreement over the policy's
    **own** triples, found exactly by dynamic programming over subsets of V9.
    This is the most favourable case for node-only ranking. Its agreement is
    `p*`.
  - The best next unit under that ranking is compared with the oracle
    `argmin Δ_i α(S)` among the available units. The regret is the difference
    in α.
- **Gates** (holdout, base-stable policies):
  - **A1 (recurrence).** The fraction of policies with ≥ 1 stable-context global
    reversal. PASS if ≥ 0.75.
  - **A2 (mode-tracked).** The fraction of policies with ≥ 1 mode-tracked
    reversal. PASS if ≥ 0.50. Reported; it does not block GOLD-A but it
    qualifies it.
  - **A3 (insufficiency).** PASS if median `p*` ≤ 0.90 **and** the median
    fraction of stable contexts (with ≥ 2 available units) that have regret
    ≥ τ_mat is ≥ 0.10.
- **Negative stopping rule.** If A1 < 0.25, H1 is REFUTED in the model class.
  No weakness index is constructed. E2, E10 and E12 still run, descriptively.

### E1b — submodularity (HS)

- **Data.** From the E1 data (no new solves), for each policy: all one-step
  second differences
  `d_ij(S) = α(S∪i∪j) − α(S∪i) − α(S∪j) + α(S)`.
- **Violations.**
  - Submodularity fails where `d_ij > τ_res`.
  - Supermodularity fails where `d_ij < −τ_res`.
- **Reported:**
  - counts and magnitudes;
  - cardinality `|S|`;
  - the modal family (tracked or switch);
  - explicit minimal counterexamples (the largest violation per cardinality).
- **Tracked-mode margin.** The same test is repeated with
  `α̃(S) = Re λ` of the mode tracked from the base critical mode along chains.
  This is secondary and does not change the primary claim.
- **Claim form.** "α_⊥ is or is not submodular or supermodular on this class."
  Never "greedy fails".

### E2 — robust contextuality (H1c)

- **Data.**
  - Policies P4 (D01) and G_S (D03); the core lattice (16 portfolios) for every
    draw:
    - TX4 discovery draws: 4 × 100;
    - CDW holdout draws: 4 × 40.
  - Additionally, all 512 V9 portfolios for the first 10 CDW draws per
    envelope, at P4.
- **Per draw:**
  - H, κ, the witness change relative to nominal;
  - the persistence of each nominal Δ_i sign class;
  - reversal presence;
  - Kendall τ, Spearman ρ and top-3 overlap of the Δ-ranking per context
    against nominal.
- **Gate A4** (CDW holdout draws, core lattice, P4).
  - For each envelope, the coverage of "≥ 1 stable-context global reversal" is
    measured.
  - PASS if ≥ 0.50 in ≥ 3 of the 4 envelopes.
  - Also reported: reversal coverage **conditioned on the witness having
    changed**. This is the robust-phenomenon vs robust-identity test.
- **GOLD-A** = A1 ∧ A3 ∧ A4.

### E3 — total re-equilibrated node sensitivity (H2, nodes)

- **Conditions.**
  - Discovery: D01, D03, D08, D11 and D14.
  - Holdout: H01–H12 (base-stable ones), plus the first 5 CDW draws per
    envelope at P4 (20 conditions, H4 target only).
  - Targets: T_a = H4 = {30,33,35,37} and T_b = V9.
- **Parameters.** Ranges are for normalization only.

  | code | parameter | range |
  |---|---|---|
  | NG_i | g_i for i ∈ T, absolute | [0, 1] |
  | NP_i | PLL scale c_i (kp_pll and ki_pll × c), i ∈ T | [0.8, 1.2] |
  | NK_i | K_A scale of each surviving SG, including 39 | [0.85, 1.15] |
  | NV_i | V_set at all 10 generator buses | ±0.02 p.u. |
  | NL_j | load scale (P and Q) at each of the 19 load buses | ±10 % |

- **Semantics.** SPR (primary) and RP (secondary). Under RP, NV is the device
  voltage reference.
- **Quantities computed:**
  - the frozen partial;
  - the total via IFT (`dw*/da = −R_w^{-1} R_a` on the semantic's residual,
    then `D_w A [dw*]` by a directional central difference of the Jacobians);
  - `dα/da = Re(v^H dA u / v^H u)` for the rightmost eigenvalue. A near-tie,
    with the second-rightmost within 1e-3, is flagged and excluded from sign
    metrics.
- **Validation.** Finite re-solves under the same semantics:
  - small central step: relative 1e-3 (NP, NK, NL), absolute 1e-3 (NG) or
    1e-4 p.u. (NV);
  - large step: ±10 % (NP, NK, NL), ±0.1 clipped to [0, 1] (NG), ±0.01 p.u.
    (NV).
- **Implementation-validity gate IV.** The total derivative must match the
  small-step finite difference within 5 % relative, or within 1e-4 absolute in
  per-range units, for ≥ 95 % of the (parameter, condition) pairs.
  - If IV fails, every derivative-based E3–E9 claim is FAIL and only finite
    differences are reported.
- **Structural checks.** Remarks R1 and R2 are checked numerically: frozen =
  total (NG, NP, NK under SPR) within 1e-6 relative, and frozen NV = 0.
- **H2-node (frozen vs total material).** Among the parameters with
  `|total| × range ≥ τ_mat`, material in a condition if either:
  - the sign disagreement rate between frozen and total is ≥ 10 %; or
  - the Kendall τ between frozen and total is ≤ 0.6.

  H2-node is TRUE if material in ≥ 25 % of holdout conditions. It is reported
  separately for SPR and RP.
- **Prediction.**
  - Spearman and Kendall between the derivative-predicted large-step Δα and
    the finite large-step Δα;
  - sign accuracy on parameters with `|finite| ≥ τ_mat`.

### E4 — total re-equilibrated link sensitivity (H2, links; GOLD-B)

- **Conditions and targets.** As in E3.
- **Parameters.** `γ_e` for all 46 branches, with the whole two-port scaled and
  the tap unchanged.
- **Steps.**
  - small: ±1e-3 relative;
  - large, and the truth for GOLD-B: finite ×1.5 under SPR (Δα);
  - also ×2 (SPR), and ×1.5 under RP (discovery and holdout policies only).
- **Predictors.** Oriented so that a higher score means a larger stabilizing
  effect of strengthening, i.e. correlation with −Δα.
  - **Static:**
    - S1: |P_e|;
    - S2: |S_e|;
    - S3: |z_e|;
    - S4: endpoint electrical distance
      `|Z_ff + Z_tt − 2 Z_ft|`, where Z is the inverse of the Ybus with SG
      subtransient reactances;
    - S5: effective resistance of the edge on L_B;
    - S6: Fiedler edge score `(u₂_f − u₂_t)²` on L_B;
    - S7: weighted edge betweenness (weight |x_e|);
    - S8: the maximum endpoint `dV/dQ` from the PF Jacobian;
    - S9: the gSCR change of T under ×1.5.
  - **Dynamic:**
    - Dconv: the conventional eigenvalue sensitivity of the **base** (empty
      portfolio) rightmost mode;
    - Dfrozen: the frozen partial for T;
    - Dtotal: the SPR total for T;
    - DtotalRP: the RP total for T;
    - Dport: the TX4 port derivative, a cross-check at D01/H4 only.
- **Metrics per condition:**
  - Spearman and Kendall against the finite ×1.5 (and ×2);
  - top-5 precision;
  - sign accuracy on branches with `|finite| ≥ τ_mat`.
- **Selection.** The static baseline carried forward is the one with the best
  median Spearman on **discovery**. The per-condition oracle-best static is also
  reported, and is favourable to static.
- **GOLD-B (links).** On holdout conditions with target H4, PASS if all hold:
  - median (ρ_Dtotal − ρ_static,selected) ≥ 0.20;
  - median ρ_Dtotal ≥ 0.70;
  - median top-5 precision of Dtotal ≥ 0.6.

  V9 is reported alongside.
- **H2-link (frozen vs total).** The same rule as H2-node, over the branches.

### E5 — static vs dynamic weakness table (descriptive)

This uses the E1, E2, E3 and E4 data; there are no new solves.

- **Per candidate unit:**
  - static: SCR, Thevenin |Z|, dV/dQ, Fiedler entry;
  - dynamic: NG and NP totals;
  - contextual: sign fractions, entropy, reversal frequency, Shapley value;
  - robustness: E2 rank persistence.
- **Per branch:** the static scores against the Dfrozen, Dtotal and finite
  effects.
- **Descriptive categories** (rules fixed now).
  - Static-weak is the top 3 by lowest SCR (units) or by S4 (branches).
  - Dynamic-critical is the top 3 by `|W̄_total|`.
  - The categories are:
    - STATIC + DYNAMIC;
    - STATIC ONLY;
    - DYNAMIC ONLY;
    - CONTEXT-DEPENDENT: reversal in ≥ 50 % of holdout policies;
    - ROBUSTLY NEUTRAL: NEUTRAL in ≥ 80 % of contexts across holdout policies.
- **Divergence.** The Kendall τ between the static and dynamic rankings.
- No combined score is built.

### E6 — dynamic weak corridors (H3)

- **Corridor constructions.** Fixed now and computed from L_B only, before any
  dynamics:
  - K2, K3 and K4: spectral clustering with k-means on the first k non-trivial
    eigenvectors of the normalized L_B (seed 20260922, 50 initializations).
    Every inter-community cutset is a corridor.
  - TXcore: branches {34, 40, 42, 44} (the core step-up transformers).
  - TXother: the other 8 transformers.
  - TXall: all 12 transformers.
- **Coordinate.** ρ_C is a common scale on all branches of C.
  - The derivative is `Σ_{e ∈ C} d/dγ_e`, exact by the chain rule.
  - The finite truth is ×1.5 on C under SPR.
- **Non-additivity index.**
  `NA = |Δ_C α − Σ_e Δ_e α| / max(|Σ_e Δ_e α|, τ_res)`.
- **Robust weak corridor.** A corridor in the top 3 by finite stabilizing effect
  in ≥ 75 % of holdout conditions (H4 target) that improves α by ≥ τ_mat in
  ≥ 75 % of them.
  - H3 is TRUE if at least one exists.
  - Also reported: the corridor-ranking Spearman of Dtotal on holdout (≥ 0.70 is
    "predictive").

### E7 — topology reconfiguration (H4-topology)

- **Operational actions:** single-branch outages that keep the network
  connected and the PF converged, with all bus voltages in [0.90, 1.10] p.u.
  - The [0.95, 1.05] compliance and the generator Q against its limits are
    reported but not enforced.
  - No thermal ratings exist in the data, so thermal limits are **not
    enforced**. Flows are reported.
- **Planning actions:** doubling each branch (×2, a parallel identical circuit),
  labelled planning reinforcement. There is no cost model.
- **Scheduled P and Q** are fixed (SPR).
- **Evaluations:**
  - the core lattice for every admissible action at D01, D03, D11 and D14 and
    at H01–H06;
  - all 512 V9 portfolios at D01 for every admissible action.
- **Reported:** H, κ, α(H4), the V9 hypergraph at P4, and flow and voltage
  changes.
- **Questions:**
  - (i) Does any single action make H4 stable at P4? (topology alone removes
    an incompatibility)
  - (ii) Does any action create a new hyperedge of size ≤ 3 at a policy where
    none existed?
  - (iii) Static prediction: the Spearman between Δα(H4) and each of ΔFiedler
    (λ₂ of L_B), Δ Kirchhoff index, ΔgSCR(H4) and Δ min SCR over the core.
    "Predicts" means |ρ| ≥ 0.6 in ≥ 75 % of the evaluated policies.

### E8 — local stability-equivalent exchange rate

- **Point.** The TX4 H4 boundary on the P4 line, θ* = (0.20768, 1.425, 1.5, 1),
  where α(H4) ≈ 0. Discovery P4 is also reported.
- **Controls.** g_i (i ∈ H4) and K_A scales of the six surviving SGs. The links
  are the top 5 by `|dα/dγ_e|` total at θ*; the selection rule is fixed here.
- **Rate.** `r = −(dα/dθ_i)/(dα/dγ_e)` (SPR totals) for the 50 pairs.
- **Validation.**
  - Paired finite intervention: `θ_i += δ` and `γ_e += r δ`, with
    δ_g = 0.02 (0.1 as the large step) and δ_K = 0.02 relative.
  - Compensation ratio `|Δα_pair| / |Δα_control only|`.
  - PASS if ≤ 0.2 in ≥ 80 % of the pairs at the small δ.

### E9 — plan-level multi-constraint design (H5; GOLD-D)

- **Targets**, run in order; a target is used only if unsafe (Φ_T ≥ 0) at
  nominal:
  - T1 = (H4, D01);
  - T2 = (V9, D01);
  - T3 = (H4, D11);
  - T4 = (V9, D03).
- **Control families:**
  - (a) control only: g_i ∈ [0, 1] for i ∈ T, and K_A scale ∈ [0.7, 1.5] per
    surviving SG;
  - (b) topology only: γ_e ∈ [1, 2], reinforcement only;
  - (c) joint.
- **Design problem.**
  - Weights `W = diag(1/range²)`; design margin ε = 0.02 s⁻¹.
  - The active set is the subsets with `α ≥ max(Φ_T − 0.05, −0.05)`.
- **Single-boundary tuning:** a QP with only the constraint on T.
- **Plan-level tuning:** a QP with every active constraint.
- **Iteration.** Sequential QP (up to 6 iterations), with totals (SPR), then an
  exact re-evaluation of every subset of T after every step. The QP solver is
  cvxpy with CLARABEL.
- **Infeasible QP.** Solve the Gordan LP `min ‖Σ λ_j g_j‖²` over the simplex.
  A value ≤ 1e-8 relative is a local conflict witness.
- **GOLD-D PASS** if, for at least one target and family, both hold:
  - single-boundary tuning ends with α(T) < 0 but Φ_T ≥ 0;
  - plan-level tuning ends with Φ_T < 0 (all subsets resolved stable).

  If single-boundary tuning already reaches Φ_T < 0, the result is recorded as
  "single-boundary suffices" (negative for H5).

### E10 — Shapley and context diagnostics (descriptive)

- **Data.** From E1.
- **Exact Shapley value** of each unit in the game `v(S) = α(S) − α(∅)` per
  policy; the marginal variance; the sign entropy.
- **Flag "hidden reversal".** `|φ_i| ≤ τ_mat` together with
  `std(Δ_i) ≥ 3 τ_mat` and a reversal present.
- No causal claim is made.

### E11 — spectral graph baselines (H6, static part)

- **Bases:**
  - (A) L_B (exact tap decomposition, C8);
  - (B) the Kron reduction of B_bus onto the 10 generator buses (Laplacian part
    of the reduced matrix);
  - (C) the dynamic-port SVD of K_o(jω_c) for H4 at D01.
- **Reported:** eigenvalues, eigenvectors, the Fiedler vector, effective
  resistances, and the graph-Fourier content of the critical modes.
- **Predictive tests** (holdout):
  - Node scores (Fiedler entry, effective resistance to bus 39) against:
    - the E1 DESTAB fraction;
    - the Shapley value;
    - the NG totals.
  - Edge scores (S5, S6) against the E4 Dtotal and the finite effect.
  - "Predicts" means median |ρ| ≥ 0.6.

### E12 — controller-weighted modal mixing (H6)

- **Computed** at D01–D15 and H01–H24, for the 16 core portfolios:
  - `T̂(jω_c)` with U from L_B;
  - `μ_mix = ‖R_mix‖_F / ‖T̂‖_F`;
  - the network-only mixing `μ_mix(Y)` as a reference.
- **Q1.** The relative range (max − min)/median of μ_mix(H4) over policies
  (with the network fixed). "Materially policy-dependent" if ≥ 0.20.
- **Q2.** The Spearman across base-stable holdout policies between the reversal
  count (E1) and μ_mix(H4). "Correlates" if |ρ| ≥ 0.5 and the permutation p
  (10⁴ permutations, seed 20260924) is ≤ 0.01.
- **Q3.** The ten frozen FC18 events on the F7B line (E1–E9 and F).
  - The critical mode of the witness is computed at g* ± 1 % of g*.
  - The support is the smallest set of L_B modes carrying ≥ 80 % of the
    graph-Fourier energy of the angle mode shape.
  - "Reproducible modal-support transition" if the support changes at ≥ 7 of
    the 9 witness-change events **and** at ≤ 20 % of the 9 control midpoints.
- **Q4.** The Jaccard overlap between the top-3 dynamically critical graph modes
  (by energy) and the 3 lowest nonzero L_B modes.

### E13 — low-dimensional and reduced models (H6 spectral closure; H7; GOLD-C)

- **Families:**
  - GM(r): Galerkin projection of the bus-voltage coordinates onto the r
    leading L_B modes, with r ∈ {2, 4, 8, 12, 16, 24, 39}. Only the network is
    reduced; the device states are exact.
  - POD(r): the same, with a basis from the SVD of the discovery critical-mode
    voltage shapes (D01–D15, core 16), r ∈ {2, 4, 8, 12, 16}.
  - TS(a/b/c): quasi-steady elimination of the fast GFL states:
    - (a) the current loop (i_d, i_q, x_id, x_iq);
    - (b) (a) plus the power filters;
    - (c) (b) plus the PLL.
- **Holdout:** the core 16 at H01–H24, plus 64 V9 portfolios (seed 20260923) at
  H01–H06.
- **Metrics:**
  - verdict accuracy, exact match of H, match of κ, α error;
  - the Kendall τ of the core Δ_i α(S) against the full model;
  - the state ratio;
  - the eigen-stage speedup and the end-to-end speedup (timed).
- **GOLD-C PASS** if a family has all of:
  - state ratio ≤ 0.70;
  - end-to-end speedup ≥ 3×;
  - verdict accuracy ≥ 0.99 on non-abstained cases;
  - abstention ≤ 20 %;
  - zero false certifications;
  - H exact on ≥ 90 % of the holdout policies.

### E14 — decision-preserving certificate (H7)

- **Families:** TS (reduced D̃, same K) and GM/POD (reduced K̃).
- **Contour:** the rectangle `Γ = {Re s ∈ [0.002, 20], |Im s| ≤ 2π·5}`,
  sampled with 400 points per side.
- **Abstain** if either:
  - a pole of D, D̃, K or K̃ lies inside Γ (checked from the device-alone and
    base spectra);
  - the full model has an eigenvalue with Re > 20 or |Im| > 2π·5 in the right
    half plane (scope check).
- **Certificate.** `max_Γ ‖R̃^{-1}(R − R̃)‖₂ < 0.9` gives SAMPLED-CERTIFIED
  (labelled non-rigorous). Otherwise the answer is ABSTAIN.
- **Metrics:** coverage, false certifications (must be 0), and speedup.
- **Stated limitation.** The certificate uses the full error term, so it cannot
  deliver online speedup by itself.

### E16 — modal energy vs stability-limiting mode (H8; secondary)

- **Scope:** stable core portfolios at D01–D15 and H01–H24.
- **Disturbances:**
  - a 0.2-s load pulse of 0.1 p.u. at bus 16;
  - the same at bus 20.
- **Channels:** the SG speeds and all 39 bus voltage magnitudes.
- **Modes:** the transverse modes in 0.1–2.0 Hz. Modal residues are
  `r_ck = (C u_k)(v_k^H x0)`.
- **Energy** over [2, 30] s is computed in closed form.
- **Measure.** How often `argmax energy ≠ argmax Re λ`.
- **"Systematic"** if they differ in ≥ 25 % of stable cases in ≥ 50 % of
  holdout policies.
- This is a phasor-domain study, **not** an EMT result.

### E17 — Africano / PV material gate

- **Search:** the repository and local project files, for the Africano source
  (thesis or paper), the feeder data, the weak-node definition (PVR or other),
  the hosting-capacity method and the reported results.
- **PASS** only if all of these are present. Otherwise:
  - write `docs/CDW_AFRICANO_MISSING_INPUTS.md`;
  - E18–E21 (Phases 18–21, GOLD-E and GOLD-F) are BLOCKED;
  - no generic feeder is substituted.

### E22 — nonlinear recovery (optional)

- Runs only if E1–E16 are complete and the E9 outcome makes it relevant.
- Uses the TX4 FC05 threshold machinery for H4 before and after the E9 design.
- Otherwise the status is DEFERRED.

### E23 — cross-model holdout (H9)

- **Alternative model.** The TX4 `StaticInjection` (constant-power converter
  limit, `device="static"`). It is clearly versioned, not retuned, and has no
  `g`.
- **Repeated:**
  - E1 at D01, D03, D11, D14 and H01–H06, using the (k, t, h) coordinates;
  - E4 Dtotal and finite ×1.5 for H4 at the same points;
  - the E6 corridors.
- **"Transfers"** if all hold:
  - reversal presence in ≥ 50 % of these policies;
  - median link Kendall(GFL, static) ≥ 0.5;
  - corridor top-3 overlap ≥ 2/3.

  Otherwise the result is "model-specific".

## 3. Stopping, blocking and negative-result policy

- **Independence.** A FAIL blocks only its logical dependents; independent
  experiments continue.
  - E1b and E10 depend on E1.
  - E5 depends on E1–E4.
  - E8 and E9 depend on the IV gate for their derivative parts; their finite
    evaluations are always reported.
- **No weakness index** is built, whatever the outcome.
- **GOLD labels** are applied only by the rules above. Results are ranked:
  - STRONG POSITIVE: the gate passes with a margin of ≥ 20 % of the threshold;
  - POSITIVE;
  - NEGATIVE BUT INFORMATIVE;
  - INCONCLUSIVE: numerical or unresolved cases > 20 %;
  - BLOCKED.
- **Negative results** are reported with the same prominence as positive ones.
- **Determinism.** The E1 census at D01, D11 and H01, E4 at D01/H4 and E9 T1 are
  rerun. The summary outputs must be byte-identical, excluding wall-time
  fields.
