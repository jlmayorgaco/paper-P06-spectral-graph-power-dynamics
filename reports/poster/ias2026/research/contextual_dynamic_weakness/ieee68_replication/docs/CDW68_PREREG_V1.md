# CDW68 preregistration V1: confirmatory replication on the IEEE 68-bus system

Branch `research/cdw-ieee68-replication`, from `e03a5312`. Committed **before any portfolio containing a converter is evaluated**.

## Status at freeze

- The model audit (`docs/CDW68_MODEL_AUDIT.md`) evaluated the base case only: S = ∅, excitation scale k = 1, all three variants.
- No converter portfolio, branch action, draw or TDS run has been computed.
- The frozen inputs are listed in `inputs/cdw68_inputs_manifest.json`, with sha256 for each file.

## 0. Question

Does the contextual-weakness phenomenon, the portfolio-conditioned dynamic reinforcement ranking, or both, replicate on a different network with documented primary-frequency dynamics? And what transfers across two frozen converter models?

**Correction of the campaign brief.** On IEEE-39, the new-policy median Spearman of the total re-equilibrated sensitivity was 0.996. The value 0.818 in the brief is the post-hoc fixed branch list learned from other clusters, and 0.712 is the frozen derivative.

## 1. Frozen model (R1)

- **Network and machines:** `configs/ieee68/ieee68_network.json` (Singh & Pal v3.3; Gate 3) and `src/ibr_cycles/models/ieee68_devices.py`, used read-only.
- **Primary-frequency data:** `inputs/pst_primary_frequency_v1.json`.
- **Variants:**
  - **REAL** (confirmatory): PST governors and PST rotor damping on every surviving machine.
  - **NOGOV** (ablation B): governors removed.
  - **SP33** (ablation C): the published benchmark, D = 0 and no governors.
- **Rules:**
  - Every headline result uses REAL.
  - No parameter of the benchmark, the governors, the damping or either converter may be changed after this freeze.
  - D = 0 is never imposed on REAL.

## 2. Converter models (R2)

- **Model A:** the TX4 GFL. Rating |S_gen|/0.8. Policy coordinates: the Q/V gain g and the excitation scale k.
- **Model B:** the WECC library chain TX3-GFL-0.1, frozen, constant-Q. Policy coordinate: k only, taken from the Model A policies.
- **Rules:**
  - Neither converter may be tuned, and neither may be changed to resemble the other.
  - No grid-forming model is run (none is validated).

## 3. Transverse classification (unchanged TX4 classifier)

- **Structural centre:**
  - span{R_x} for REAL and NOGOV (a single structural zero, verified in R1);
  - span{R_x, w} for SP33 (w derived as in Gate 3).
- **Verdict:**
  - STABLE if every transverse eigenvalue has Re ≤ −2e-3;
  - UNSTABLE if any has Re ≥ 2e-3;
  - otherwise the frozen four-state classifier with the two-scale Jacobian error (`classify_spectrum`, SAFETY).
- Dead states (identically zero rows) are deleted first.
- α⊥ is the largest transverse real part. No magnitude cut-off shortcut is used.

## 4. Candidate set (R4)

**V68 = {G3, G4, G6, G9, G11, G12}.** These are the six physical plants (G1–G12) with the largest documented scheduled active power:

| unit | P (pu) |
|---|---|
| G12 | 13.50 |
| G11 | 10.00 |
| G9 | 8.00 |
| G6 | 7.00 |
| G3 | 6.50 |
| G4 | 6.32 |

- **Exclusions:** the area equivalents G13–G16, including the slack G16.
- **Why this rule:**
  - It extends the repository's Gate 3 candidate rule (G9, G6, G3, G4) to both physical areas: NETS gets 4 units, NYPS gets 2.
  - No tie occurs.
  - It uses no dynamic quantity.
- **Alternative not used:** PST machine ratings would give {G2, G6, G9, G10, G11, G12}. They were not used because the model's own source (Singh & Pal) documents only a per-unit base, and the repository's rule uses dispatch.
- **Size:** 2^6 = 64 portfolios.
- **Replacement:** full (rho = 1), matched dispatch, one fixed operating point for every portfolio.
- The candidate set may not be replaced, whatever the outcome.

## 5. Model B composition and qualification (R2, R10)

**Composition.** Everything is built in ANDES 2.0.0 (`.venv/xtool-andes-gfl`):
- The network is rebuilt from the same JSON: Bus, Line with taps, PV, Slack, and PQ converted to constant impedance.
- The Singh & Pal machine, DC4B, ST1A, manual excitation, speed PSS, PST governor and rotor damping are transcribed as a custom ANDES model **SG68**, the way SG2AX was for IEEE-39.
- Converted buses get the library chain exactly as in H17 (REGCP1 Sn = |S_gen|/0.8 × 100 MVA, forced-PQ equilibrium).
- The spectrum uses the H17 descriptor QZ solve and H17 refinements 1, 1b and 2.

**Qualification (before any B portfolio with a converter).**

| check | requirement |
|---|---|
| Q1 | The ANDES all-SG REAL base reproduces the internal α⊥ within 1e-5 s⁻¹ at all 16 k values, with the same status. |
| Q2 | Initialization residual ≤ 1e-6 at every B case. |
| Q3 | Descriptor QZ agrees with the reduced-state spectrum within 1e-8 relative at one case. |
| Q4 | No limiter active outside the H17 exemption list. |

If Q1 fails, Model B is **BLOCKED** on IEEE-68. It is reported as such and never adjusted.

## 6. Policies (R5)

**Model A: 16 policies P68_01–P68_16.**
- Maximin Latin hypercube over the validated Gate 3 window: u ∈ [0, 1] with g = u², k ∈ [0.5, 2.0].
- Best of 2000 designs, seed 20260913, minimum pairwise distance 0.1748 in the unit square.
- Values are in `inputs/cdw68_design_v1.json`.
- Split: P68_01–04 are **discovery** (selection of the static comparator only); P68_05–16 are **holdout**.

**Model B: 16 conditions B68_01–16**, carrying the k of the corresponding A policy.

- Policies whose all-SG base is unstable are retained. They are excluded only from tests that require a stable base.
- No policy search is allowed.

## 7. Exhaustive census (R6)

- **Scope:** every variant, model, policy and S ⊆ V68 in scope.
  - REAL, Model A: 16 × 64 portfolios.
  - Model B: 16 × 64.
  - Ablations: 6 × 64 for each of NOGOV and SP33.
  - Draws: see R15.
- **Recorded:**
  - status and α⊥;
  - critical eigenvalue and frequency;
  - RHP count; gap to the next distinct eigenvalue;
  - bus-voltage mode shapes (136 real components) of the band modes (0.1–2 Hz, Re ≥ −1) and of the critical mode;
  - the rightmost EM-band mode;
  - rotor-angle/speed participation of the critical and EM-top modes;
  - feasibility; equilibrium residual; cond(g_z).
- **Per condition:** H_RHP (minimal unstable portfolios) and κ = the smallest hyperedge size.
- **Output:** `results/CDW68_PORTFOLIOS.parquet`.

## 8. Primary replication: nested same-mode contextual reversal (R7)

- **Marginal:** Δ_i α(S) = α(S ∪ {i}) − α(S), with τ = 0.01 s⁻¹ (unchanged). The sign class is −1 if Δ ≤ −τ, +1 if Δ ≥ +τ, otherwise 0.
- **Levels** (identical to the hardening campaign):
  - **A:** S stable; marginal resolved (S ∪ {i} not unresolved).
  - **B:** the critical mode of S has a match in S ∪ {i} (MAC ≥ 0.8 on the bus-voltage shapes among the recorded modes, |Δf| ≤ 0.15 Hz), and that match is the critical mode of S ∪ {i}.
  - **C:** B, with both critical frequencies in 0.1–2.0 Hz.
  - **D:** for a nested pair S1 ⊊ S2 (i ∉ S2, both stable), both marginals are level C and the critical modes of S1 and S2 match each other (MAC ≥ 0.8, |Δf| ≤ 0.15 Hz).
- **Level-D nested reversal:** the two marginals are in opposite non-zero sign classes. Fast-mode jumps, mode switches and unresolved transitions never count.

**PRIMARY GATE (Model A, REAL).**
- At least 50 % of the eligible (all-SG base STABLE) policies contain at least one level-D nested reversal, AND the median material reversal magnitude, min(|Δ(S1)|, |Δ(S2)|) over level-D pairs, is ≥ 0.01 s⁻¹.
- The second clause holds by construction whenever a level-D pair exists; it is recorded for completeness.

Model B uses the same gate and is reported independently. A and B are never pooled or averaged.

**PREREGISTERED SECONDARY: EM-tracked reversal (level T).**
1. For a stable context S, take the rightmost EM-band mode (0.1–2 Hz, excluding structural zeros).
2. Match it into S ∪ {i} among EM-band modes (MAC ≥ 0.8, |Δf| ≤ 0.15 Hz). The match must be the rightmost EM-band mode of S ∪ {i}.
3. The marginal is Δ^T_i(S) = Re(matched) − Re(tracked in S).
4. A level-T nested reversal additionally needs the tracked modes of S1 and S2 to match (MAC ≥ 0.8, |Δf| ≤ 0.15 Hz), and τ is unchanged.

Level T does not require the EM mode to be the global rightmost. It is reported next to the primary and never substitutes for it in the gate. Rationale: in REAL the base rightmost mode is a slow real mode (R1); level T asks the electromechanical question directly.

Reported descriptively:
- the direction split (s→d, d→s, both);
- magnitude quartiles;
- coverage against τ ∈ {0.01, 0.0125, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05};
- gap-robust coverage (pairs whose rightmost gap ≥ τ at all four portfolios).

## 9. Discrete curvature (R8)

- For every nested reversal (levels A–D and T), verify the chain identity Δ_i f(S_m) = Δ_i f(S_0) + Σ_r d_{i j_r}(S_{r−1}) on the canonical chain (added units in increasing bus order).
- **Acceptance:** residual ≤ 1e-12.
- **Reported:** direction, cumulative sum, first context where the partial sum crosses zero, and positive and negative contributions.
- The identity is classical and is not claimed as new.

## 10. Fixed node ranking (R9)

- **Construction:** per REAL policy, the exact oracle fixed ordering of the six units that maximizes pairwise agreement (all 720 orders enumerated).
- **Measures:**
  - p*;
  - material next-action regret (the chosen unit's Δα exceeds the best available by ≥ τ) over stable contexts with at least two resolved candidates;
  - the regret-optimal order (enumerated);
  - the stability-screened order (first ranked unit whose S ∪ {i} is stable, on contexts with a stable option).
- **Strata:**
  - **FULL** (all resolved decisions);
  - **EM** (EM→EM transitions);
  - **LEVEL-D-eligible** (only candidates whose marginal is level C);
  - **TRACKED** (the level-T marginal as objective).
- **Transfer:** an order fitted at policy a is applied at policy b (base-stable policies).
- **Insufficiency bar** (unchanged): p* ≤ 0.90 and regret ≥ 0.10.
- **Headline rule:** a global ranking failure may be headlined only if the LEVEL-D-eligible or TRACKED stratum also passes the bar. A stratum with fewer than 20 decisions per policy in the median is "not evaluable".

## 11. Branch reinforcement (R10–R12)

- **Target:** the portfolio V68 (all six candidates converted) at each condition.
- **Actions:** every in-service branch (83). Whole two-port admittance scaled by γ ∈ {1.10, 1.25, 1.50}.
- **Truth:** finite re-equilibrated effect −[α⊥(γ) − α⊥(1)] under the scheduled semantics SPR-68:
  - every unit keeps its scheduled P and V;
  - the slack G16 keeps its voltage and angle and absorbs losses;
  - every load keeps its scheduled P and Q at the new equilibrium (the benchmark's load conversion);
  - converter references are re-solved.
- **Also recorded:** power-flow convergence, bus-voltage magnitudes, max |ΔV|, and buses outside [0.8, 1.2] pu.

**Model A predictors:**
- **D_tot:** the implicit-function total derivative of the critical eigenvalue under SPR-68, dλ/dγ = v^H(∂_γA + ∂_wA[ẇ])u / v^H u, ẇ = −R_w⁻¹ R_γ.
- **D_fro:** the frozen term alone.
- **D_conv:** the base-case (S = ∅) total derivative.

**Model B predictor:** no analytical engine exists for the ANDES library model. D_tot^B is the central re-equilibrated difference at h = 1e-3, labelled "numerical total derivative".

**Static baselines:**
- |P_e|, |S_e|, |z_e|;
- electrical distance (endpoint distance on the sub-transient bus impedance matrix);
- effective resistance and Fiedler edge score (coupling Laplacian);
- weighted edge betweenness;
- endpoint dV/dQ;
- ΔgSCR of the target portfolio at γ = 1.5.

**Selection of the static comparator.** Per model, the static index with the highest median Spearman on the four discovery conditions. It is frozen for the holdout.

**R11 derivative validation (before any ranking claim).**
- Scope: 20 seeded branches (`r11_branches`) × the 4 discovery policies, for Model A.
- Test: D_tot against the central re-equilibrated difference at h = 1e-4.
- **Acceptance:** sign agreement ≥ 95 %, median relative error ≤ 5 %, and Spearman ≥ 0.95.
- If this fails: stop the Model A ranking claim, diagnose, do not retune.
- Model B has no analytical derivative, so its check is numerical consistency only: h = 1e-3 against h = 1e-4, with the same three acceptance numbers. This is a documented limitation of the B arm, not a validation of an engine.

**R12 gate, per model independently, on the 12 holdout conditions at γ = 1.5.**
- Median Spearman(D_tot, finite) ≥ 0.70;
- AND median paired advantage over the discovery-selected static comparator ≥ 0.20;
- AND median top-5 precision ≥ 0.60.

Also reported: Kendall, sign accuracy, NDCG@5, all 16 conditions, and degradation across γ = 1.10/1.25/1.50. Identical line rankings under A and B are **not** required.

## 12. Governor and damping ablation (R13)

- **Scope:** Model A, policies P68_01–06, all 64 portfolios, for NOGOV and SP33. REAL comes from the main census.
- **Reported per variant:**
  - level-D coverage (and level T);
  - magnitude quartiles;
  - the base rightmost mode type (EM or slow real);
  - H and κ.
- The headline uses REAL only.
- The question it answers: is reversal created by the no-governor/no-damping simplification, or does it survive documented primary-frequency dynamics?

## 13. Nonlinear phasor TDS holdout (R14)

- **When:** after the linear results are frozen.
- **Cases per model**, selected by rule before any TDS:
  - the 3 strongest level-D reversals (largest min |Δ|, distinct units first; if fewer than 3 exist, use level T instead and label it);
  - the 2 nearest-threshold reversals (smallest min |Δ| ≥ τ);
  - 2 negative controls: nested pairs meeting every level-D condition, but with the same non-zero sign in both contexts, largest min |Δ|.
- **Each case** simulates the four portfolios S1, S1 ∪ {i}, S2, S2 ∪ {i}.
- **Common disturbance:** a 50 MVAr (0.5 pu) shunt reactor at bus 3, on at t = 1 s and off at t = 11 s, run to t = 31 s. This is the documented benchmark case (b).
- **Integrator:** the validated `ibr_cycles.nonlinear.tds.simulate` (BDF, network solved at each step) with generic guards (bus voltage 0.8–1.2 pu, speed and PLL frequency −3/+1.8 Hz).
- **Decay estimate:** the tracked mode's decay rate σ̂ comes from the modal coordinate q(t) = v^H (x(t) − x*), v the left eigenvector of the tracked mode, by least-squares slope of log|q(t)| over 12–30 s.
- **Corroboration of a reversal pair:** sign(σ̂(S1 ∪ {i}) − σ̂(S1)) and sign(σ̂(S2 ∪ {i}) − σ̂(S2)) are opposite, and equal to the linear signs. If TDS contradicts the linear sign in two or more of the three strongest cases, the reversal claim is weakened.
- Model B runs TDS only if its ANDES TDS reproduces the linear decay of the base case within 10 %; otherwise B TDS is "not run".
- No EMT claim.

## 14. Uncertainty holdout (R15)

- **Envelopes** (deterministic stress envelopes, not probability distributions): E05 = ±5 % and E10 = ±10 %.
  - Fleet multipliers on the surviving machines: H; x'd and x'q (together); regulator gain; PSS gain.
  - Model A converters: PLL gains (together); current-loop gains (together); outer-loop gains (together).
- **Draws:** Latin hypercube with fixed seeds. Model A: 20 per envelope (40). Model B: 10 per envelope (20), machine multipliers only.
- **Reference policy:** the base-stable policy nearest the window centroid (u = 0.5, k = 1.25) in normalized coordinates; ties go to the lower index.
- **Per draw:**
  - all 64 portfolios;
  - level-D and level-T reversal presence;
  - the strongest witness (i, S1, S2);
  - branch-ranking Spearman on the first 5 draws of each A envelope.
- **Phenomenon robustness:** the fraction of draws with at least one level-D (or T) reversal.
- **Witness identity:** the fraction of draws whose strongest witness equals the nominal one.

## 15. Negative-result policy (R16)

- If the REAL Model A primary gate fails:
  - no other network is searched, no policy is tuned, no candidate is replaced;
  - the conclusion is "IEEE-39 contextual sign reversal is benchmark/model dependent under the present evidence".
- A pass on A with a fail on B means model dependence.
- A ranking failure is never rescued with a line subset.

## 16. Statistics (R17)

See `docs/CDW68_STATISTICAL_PLAN.md`.
- Coverages are deterministic fractions of designed conditions, never probabilities.
- Paired bootstrap over conditions is primary; paired sign-flip and Wilcoxon are secondary.
- Holm family F68 = {T1 (p* < 0.90 sign test, REAL A), T2A (ranking advantage, A holdout), T2B (ranking advantage, B holdout)}.

## 17. Paper decision (R19), exhaustive

**Definitions.**
- **REV68A** = the primary gate on REAL Model A.
- **REV68B** is one of:
  - PASS, if the same gate passes on B;
  - PARTIAL, if level-D coverage on B is > 0 but below 50 %, or if level-T coverage on B is ≥ 50 % while level D is below 50 %;
  - FAIL, otherwise, or if B is blocked.
- **RANK68A / RANK68B** = the R12 gate on each model.

**Decision rule.** Exactly one case applies.

| case | condition | consequence |
|---|---|---|
| A | REV68A = PASS, REV68B ∈ {PASS, PARTIAL}, RANK68A = PASS and RANK68B = PASS | CDW is a serious TPWRS candidate; revise into a cross-benchmark manuscript |
| B | RANK68A = PASS and RANK68B = PASS, and not case A | pivot the paper to portfolio-conditioned dynamic reinforcement ranking; contextual reversal is reported as benchmark-specific evidence |
| C | every other outcome | do not send CDW as a TPWRS claims paper; keep it as a methodological case study or preprint |

- The final report states each component (REV68A, REV68B, RANK68A, RANK68B).
- If Model B is blocked, RANK68B = FAIL, so neither case A nor case B applies.
- No other interpretation may be added after results.

## 18. Paper update (R23)

The CDW manuscript is not edited before the R19 decision. Afterwards it is edited only as the decided case prescribes.

## 19. Stop rules

- A failed preregistered gate is reported and never retuned.
- A failed R11 on Model A stops the A ranking claim.
- A failed Q1 blocks Model B.
- Execution fixes that change no definition are logged in `docs/CDW68_DEVIATIONS.md`, with system-clock times and commit hashes.
