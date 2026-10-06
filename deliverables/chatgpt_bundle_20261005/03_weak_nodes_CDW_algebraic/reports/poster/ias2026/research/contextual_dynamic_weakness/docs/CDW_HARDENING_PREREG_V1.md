# CDW hardening preregistration V1 — frozen before any hardening numerics

Date: 2026-09-12. Branch `research/contextual-dynamic-weakness-hardening`, created
from `e73dd355f20020f508238e7cd680e92bf665e272` (the completed CDW campaign).
TX4 (tag `TX4_FINAL_MANUSCRIPT_FREEZE`, commit `69f200df`) is read only. No push.

This file, `CDW_HARDENING_EXPERIMENT_MATRIX.md`, `CDW_HARDENING_CLAIM_MATRIX_V1.csv`
and `CDW_HARDENING_STATISTICAL_PLAN.md` are committed together, **before any new
numerical experiment**. The only computation before this commit is
`experiments/cdw_hardening/H00_prereg_inputs.py`, which generates design points,
envelope draws and null groups from seeds. It evaluates no model.

**What was known when this was written.** Every result of the first CDW campaign
(`docs/CDW_FULL_CAMPAIGN_REPORT.pdf`) was known. The hardening tests were designed
with that knowledge, so they are confirmatory only on data that did not exist then:
- the 24 new policies HARDENING_H01–H24;
- the 40 fresh envelope draws;
- the equal-budget corridor and null evaluations;
- the fixed-frequency mixing;
- the electromechanical (EM) topology tracking;
- the alternative converter model.

Re-analyses of old data (nested reversals on the old census, curvature witnesses,
H_EM on the old E7 census) are labelled **re-analysis** and never count as
independent confirmation.

**Change rule.** After commit, no hypothesis, threshold, seed, test set, baseline,
definition, figure rule or stopping rule may change. A deviation is recorded in
`docs/CDW_HARDENING_DEVIATIONS.md`, with a timestamp and reason, before the affected
result is read. Pure execution fixes (a crash, a path) are recorded there too,
together with a statement that no definition changed.

Definitions from `docs/CDW_THEORY_V1.md` and `docs/CDW_PREREG_V1.md` hold unless
redefined here:
- α = α_⊥ is the transverse spectral abscissa;
- status comes from the TX4 classifier;
- τ_mat = 0.01 s⁻¹ and τ_res = 1e-3 s⁻¹;
- MAC_min = 0.80 and Δf_max = 0.15 Hz;
- the EM band is 0.1–2.0 Hz;
- V9 = {30,…,38} and V4 = H4 = {30,33,35,37}.

## 0. Frozen inputs

| input | value |
|---|---|
| model | the CDW model library `experiments/cdw/_cdw.py` and `_sens.py` (TX4 `solve_case`, devices, classifier, transverse operator; custom GFL, voltage_control, leak 0.05), unchanged |
| new policy holdout | `results/hardening/prereg_inputs/hardening_policies.json`, sha256 `53640af6…4dc7614`: 24 points HARDENING_H01–H24. Deterministic maximin LHS: seed 20260929, 2000 candidate designs, best min-distance 0.3233 at candidate 1137. Same domain as the old holdout: u∈[0,1] with g = u², k∈[0.5,2.3], t∈[0.5,3.0], h∈[0,2]. No point was chosen, replaced or screened for stability |
| fresh uncertainty holdout | `hardening_envelope_draws.json`, sha256 `3d033647…c304a64`: 10 draws per envelope EM-f, EM-u, EC, EMC. The frozen CDW00 generator and bounds, new seed 20260930 (+ envelope index), all evaluated at P4 = D01 |
| corridor null groups | `corridor_null_groups.json`, sha256 `fd53204e…9a560d00`, seed 20260931. Family A (arbitrary): 500 per size for k ∈ {2,3,4,8,12}; all 46 branches for k = 1. Family B (connected edge subgraphs): the exact complete population for k ≤ 4 (46, 82, 168, 358 groups, all smaller than 500); 500 sampled by random growth for k ∈ {8,12} |
| frozen corridors | `results/CDW_E6_corridor_definitions.json` (old prereg), unchanged: K2_01 {0,4,18}; K3_01 = K4_03 {24,25}; K3_02 = K4_02 {4,14}; K3_12 = K4_23 {18}; K4_13 {21}; TXcore {34,40,42,44}; TXother {35,36,37,38,39,41,43,45}; TXall (12) |
| old data re-used | the E1 census raw (`raw/E01`, 39 policies × 512 portfolios with mode shapes), E4/E6/E7 outputs, E2 draws |
| discovery-selected static baseline | S1_absP, selected on discovery in the old E4 and never reselected |
| ω_ref | 2π × 0.6222796695779355 rad/s: the critical frequency of H4 at P4 = D01 (old `CDW_E12_mixing.csv`, the TX4 flagship mode) |
| θ_ref | P4 = D01 = (0.03625, 1.425, 1.5, 1) |
| seeds | bootstrap 20260932; permutations 20260933 |

The policy IDs are HARDENING_H01 … HARDENING_H24; the coordinates are in the JSON.
"Old holdout" means H01–H24 of `docs/CDW_PREREG_V1.md`. "Discovery" means D01–D15.

**Base-unstable rule.** A policy whose empty portfolio has α(∅) ≥ 0 is labelled
BASE_UNSTABLE.
- It is kept in every dataset and counted.
- It is excluded only from tests whose definition requires a stable base: the
  reversal gates, fixed ranking, and the H13 correlations.
- It is never replaced.

## 1. Theory (H2) — deliverable, no numerics

`theory/CDW_CONTEXTUAL_CURVATURE_THEOREM.md` will contain:
- the chain (telescoping) identity for Δ_i f(S_m) in terms of the one-step second
  differences d_{i,j_r}(S_{r−1});
- Propositions A (submodular ⇒ no − → + move along a growing chain),
  B (supermodular ⇒ no + → − move) and C (both directions observed ⇒ neither);
- a quantitative witness bound: some curvature term on the chain has the required
  sign and magnitude ≥ |Δ_i f(S_m) − Δ_i f(S_0)|/m;
- a meet lemma: every arbitrary-context sign reversal implies a nested sign change
  through S₁∩S₂.

Classical ingredients are credited:
- Topkis (increasing differences);
- Lovász (1983);
- Fujishige (2005);
- Bach (2013).

The algebraic identity is not claimed as new.

## 2. Hypotheses, measurements and frozen gates

### H1 — the new policy holdout census
- All 512 V9 portfolios at every HARDENING policy, with the E1 recording and mode
  shapes (the E1 `run_task` code path).
- Also the core lattice (16 portfolios) at P4 for each of the 40 fresh draws (L7).

### H3 — nested reversal (primary GOLD-A hardening)

A marginal is (S, i), i ∉ S, with Δ_i α(S) = α(S∪i) − α(S). Its sign class uses
τ_mat. Marginal eligibility has four levels:

| level | condition on (S, i) |
|---|---|
| A global | S STABLE |
| B mode-tracked | A, and the critical mode of S has a MAC ≥ 0.8, abs(Δf) ≤ 0.15 Hz match in S∪i that is the critical mode of S∪i (the old A2 "same_mode") |
| C EM same-mode | B, and the critical frequencies of S and S∪i both lie in 0.1–2.0 Hz |
| D same-family (secondary) | C at both ends, and crit(S₁), crit(S₂) matched to each other (MAC ≥ 0.8, abs(Δf) ≤ 0.15 Hz) |

- **Nested reversal at level X.** For one intervention i there are contexts
  S₁ ⊊ S₂ with i ∉ S₂. Both (S₁, i) and (S₂, i) are level-X eligible, and their
  sign classes are opposite (−1 and +1).
- **Direction.**
  - "stabilizing → destabilizing" (class −1 at S₁, +1 at S₂) violates
    submodularity.
  - "destabilizing → stabilizing" violates supermodularity.
- **Curvature witness.** For each nested reversal, take the canonical chain from
  S₁ to S₂ (added elements in increasing bus order) and compute every
  d_{i,j_r}(S_{r−1}) from census α.
  - The largest term with the required sign is recorded; by the identity it
    exists.
  - Also recorded: whether a witness exists whose four portfolios are all STABLE
    with EM-band critical modes ("EM-clean witness").
- **Computed on:** discovery, old holdout (re-analysis of `raw/E01`), and new
  holdout.
- **Reported:**
  - fraction of eligible policies;
  - counts;
  - |S₁|, |S₂| and |S₂∖S₁|;
  - reversal magnitude min(|Δ(S₁)|, |Δ(S₂)|);
  - witness magnitude;
  - direction counts.
- **PRIMARY GATE H3.** Among base-stable NEW holdout policies, the fraction with at
  least one level-C nested reversal must be ≥ 0.75.
  - If it is lower, the old GOLD-A remains as reported, but "nested electromechanical
    same-mode reversal" is NOT_SUPPORTED.
- Level D is reported with the same 0.75 bar, as a secondary, non-blocking
  qualifier.

### H4 — fixed-ranking insufficiency and transfer

- **Per policy.** The exact DP optimal fixed ranking over the policy's own stable
  contexts, as in old E1. Report:
  - p*;
  - material regret rate: the fraction of stable contexts with ≥ 2 available units
    whose ranked choice is ≥ τ_mat worse than the context-aware argmin;
  - median and max regret.
- **Stratification** (H26):
  - FULL: global Δα.
  - EM: contexts with an EM-band crit(S), and units with an EM-band crit(S∪i).
  - same-mode: level-C marginals.
- **Gate H4-new** (the old A3 rule, new holdout, FULL): median p* ≤ 0.90 **and**
  median material regret rate ≥ 0.10.
- **Transfer matrix.** Base-stable policies of all three sets. The ranking fitted at
  θ_a is evaluated at θ_b for pairwise accuracy and material regret rate.
  - **Distance.** Euclidean in (u, (k−0.5)/1.8, (t−0.5)/2.5, h/2).
  - **Reported:**
    - diagonal (= p*) vs off-diagonal medians;
    - Spearman(distance, degradation), where degradation = p*(b) − acc(a→b),
      with a Mantel permutation p (descriptive).
- **Layer L6 positive** if the median over test policies of the off-diagonal
  material regret rate is ≥ 0.10 **and** exceeds the median diagonal regret rate.

### H5 — GOLD-A statistics

This follows `CDW_HARDENING_STATISTICAL_PLAN.md`. The effect sizes are:
- p*;
- regret rate;
- nested (C) coverage;
- level-C reversal magnitude.

Each is reported for new, old and pooled policies, with median, IQR and a cluster
bootstrap over policies, and the old-vs-new shift with its bootstrap interval.

### H6 — GOLD-B new link holdout

- **Conditions.**
  - All 24 HARDENING policies, targets H4 and V9.
  - The 40 fresh draws at P4, target H4.
  - A condition is **eligible** if the target solves and the derivative engine
    returns (R0 ≤ 1e-8). Every eligible condition is included. Near-ties
    (gap2 < 1e-3) are flagged, and a sensitivity analysis excludes them.
- **Truth.** The finite re-equilibrated Δα of the target under γ_e ∈ {1.10, 1.25,
  1.50} for all 46 branches: the SPR re-solve with the TX4 `solve_case`. Recorded
  for each step:
  - convergence;
  - min and max bus voltage;
  - admissibility, with every bus in [0.90, 1.10] p.u.
- **Primary truth:** γ = 1.50, all converged steps.
- **Predictors.**
  - Static: S1–S9 as frozen in old E4, plus NDCG-compatible orientation.
  - Dynamic:
    - Dconv (base portfolio);
    - Dfrozen;
    - Dtotal (SPR);
    - Dport (H4 targets at the 24 policies).
- **Metrics per condition:**
  - Spearman;
  - Kendall;
  - pairwise concordance over pairs with abs(truth difference) ≥ τ_res;
  - top-3;
  - top-5;
  - NDCG@5, with relevance = truth − min(truth) over the 46 branches;
  - sign accuracy on abs(truth) ≥ τ_mat.

  Orientation: a higher score means a more stabilizing strengthening (−Δα).
- **Stratification** (H26):
  - FULL: truth = α.
  - EM: conditions whose target critical frequency is in 0.1–2.0 Hz.
  - same-mode: truth = the change of the tracked eigenvalue (nearest to λ₀ after
    the step).

### H7 — GOLD-B paired test

- **Primary set.** H4 targets over 24 HARDENING policies plus the 40 fresh draws.
- **Paired difference** per condition: Δρ = ρ_Dtotal − ρ_S1_absP, on the same
  branches and outcomes.
- **Clusters** are the 24 policies and the 4 envelopes.
- **PRIMARY GATE H7.** All three must hold:
  - the lower 95 % cluster-bootstrap bound of the **median** Δρ is > 0.20;
  - median ρ_Dtotal ≥ 0.80;
  - median top-5 precision of Dtotal ≥ 0.70.

  If any fails, GOLD-B is WEAKENED.
- **Secondary inferential summary.** A one-sided cluster sign-flip permutation test
  of the median of (Δρ − 0.20) > 0 (10⁴ flips). It is in the Holm family F_conf.
- **Also reported:**
  - the median Δρ against the per-condition oracle-best static (favourable to
    static);
  - the V9 targets;
  - γ = 1.10 and 1.25.

### H8 — why total beats frozen

- **Decomposition.** Per coordinate and condition: total = frozen + indirect.
  Reported:
  - ratio abs(indirect)/abs(frozen);
  - sign cancellation: sign(frozen) ≠ sign(indirect);
  - frozen-vs-total sign flip;
  - rank displacement abs(rank_total − rank_frozen);
  - Kendall(frozen, total).
- **Coordinates.**
  - Branches: from H6.
  - Node coordinates at the 24 HARDENING policies (H4, SPR): g_i, PLL_i (i ∈ H4),
    K_A of the surviving SGs, V_set at the 10 generator buses, load scale at the
    19 load buses. These are the old E3 specs.
- **Structural expectation** (old R1). Under SPR, frozen = total for controller
  coordinates (g, PLL, K_A) within 1e-6 relative; network and operating-point
  coordinates can be dominated by the indirect term.
- The report selects examples by a fixed rule:
  - "agree": the branch with minimal rank displacement among the 5 largest abs(total);
  - "differ": maximal rank displacement;
  - "flip": the largest abs(total) among sign flips.

### H9 — equal-budget corridors

- **Frozen corridors** are listed in §0.
- **Budgets:**
  - primary L1: Σ_{e∈C} abs(Δγ_e) = 0.50, uniform, Δγ_e = 0.50/abs(C);
  - secondary L1 = 1.00;
  - L2 sensitivity: ‖Δγ_C‖₂ = 0.50, uniform, Δγ_e = 0.50/√abs(C).
- **Truth.** Finite re-equilibrated Δα(H4), SPR.
- **Conditions:**
  - discovery (D01, D03, D08, D11, D14);
  - old holdout (H01–H12 plus the old 20 draws, as in old E6);
  - new holdout (the 24 HARDENING policies plus the 40 fresh draws).
- The old unnormalized ×1.5 corridor ranking is **not** final evidence.

### H10 — size-matched corridor null

- **Coverage.** Every corridor size k ∈ {1,2,3,4,8,12}, both families, under the
  primary budget, at every old- and new-holdout condition (96 conditions). The
  k = 1 null equals the single-branch ×1.5 results: H6 γ = 1.5 for new conditions,
  an explicit re-solve for old ones.
- **Percentile** of corridor C: the fraction of null groups whose stabilizing
  effect −Δα is strictly below C's, plus half the ties.
- **Primary null** is family A. Family B is the matched null for connected
  corridors ({24,25} and the size-1 corridors) and is secondary for all others.
- **PRIMARY GATE H10** ("robust weak corridor"). A corridor qualifies if:
  - its family-A percentile is ≥ 0.95 in ≥ 75 % of the NEW holdout conditions
    (24 policies + 40 draws);
  - the old-holdout percentile coverage is reported as replication.

  Tested for TXother, TXall and every spectral cutset, with TXcore reported.
  Otherwise the phrase "robust weak corridor" is prohibited.

### H11 — corridor cross-validation

Reported per set (discovery, old holdout, new holdout policies, fresh draws):
- the stabilizing fraction (Δα ≤ −τ_mat);
- the top-3 fraction among the 11 frozen corridor names;
- the median equal-budget effect;
- the rank under the fresh draws;
- the size-null percentile.

This table replaces the unnormalized E6 headline.

### H12/H13 — modal mixing without the frequency confound

- **Definition.** μ(θ, ω) = ‖T̂ − blockdiag(T̂)‖_F / ‖T̂‖_F, with
  T̂ = (U⊗I)ᵀ T(jω) (U⊗I), T(s) = g_z + g_x (sI − f_x)⁻¹ f_z at the H4 equilibrium
  of policy θ, and U from L_B. This is the old E12 construction.
- **Frequencies.** ω_c(θ) is H4's critical frequency, with 0.7 Hz if it is below
  1e-3 Hz (old rule).
- **Policies:** all 63 (15 discovery, 24 old, 24 new). Computed:
  - μ(θ, ω_ref);
  - μ(θ, 2π·{0.60, 0.70, 0.80} Hz);
  - μ(θ, ω_c(θ));
  - μ(θ_ref, ω_c(θ)).
- **Two-factor decomposition.** With μ00 = μ(θ_ref, ω_ref), μ10 = μ(θ, ω_ref),
  μ01 = μ(θ_ref, ω_c(θ)) and μ11 = μ(θ, ω_c(θ)):
  - control C = ½[(μ10−μ00) + (μ11−μ01)];
  - frequency F = ½[(μ01−μ00) + (μ11−μ10)];
  - C + F = μ11 − μ00 exactly.

  This is a decomposition only, with no causal reading.
- **Correlation targets** (base-stable policies):
  - n_rev: the number of interventions with a stable-context global reversal;
  - n_nested: the number of interventions with a level-A nested reversal;
  - regret: the material regret rate.
- **PRIMARY confirmatory test H13.** Spearman(μ(θ, ω_ref), n_rev) over base-stable
  NEW holdout policies. It "correlates" if abs(ρ) ≥ 0.5 **and** the Holm-adjusted
  two-sided permutation p ≤ 0.01. The permutation p is in F_conf.
- **Secondary** (Holm within family F_mix, labelled secondary):
  - the 3 fixed frequencies × 3 targets;
  - μ(θ, ω_ref) against n_nested and regret;
  - the control contribution C against the 3 targets;
  - the same tests on the old holdout (18 policies, a re-test of old Q2 at fixed
    frequency).
- If the primary test fails while the old μ(θ, ω_c(θ)) correlation replicates, the
  old E12 mechanistic reading is WEAKENED: the correlation is then attributed to
  frequency movement or unresolved.

### H14 — topology, EM-only

- **Actions.** The old admissible set (35 outages, 46 doublings), frozen in
  `results/CDW_E7_admissibility.csv`.
- **Policies:** the old 10 (D01, D03, D11, D14, H01–H06) and HARDENING_H01–H06 (by
  ID).
- **Evaluation.** Core lattice (16) with mode shapes, nominal topology and every
  action.
- **Tracking.**
  - For each (policy, action, S), m₀ is the rightmost transverse mode of S at
    nominal topology with frequency in 0.1–2.0 Hz (the EM-band rightmost, not
    necessarily the global one).
  - m₁ is its best MAC match among the post-action EM-band modes, with MAC ≥ 0.8
    and abs(Δf) ≤ 0.15 Hz.
  - Δα_EM = Re m₁ − Re m₀. Δα_full = α_post − α_pre.
- **EM class:**
  - UNRESOLVED: a solve fails, the status is BOUNDARY_OR_UNRESOLVED on either side,
    or there is no EM-band m₀;
  - MODE-SWITCH: no match;
  - EM-STABILIZING: Δα_EM ≤ −τ_mat;
  - EM-DESTABILIZING: Δα_EM ≥ +τ_mat;
  - EM-NEUTRAL: otherwise.
- **Fast flag.** FAST-MODE DOMINATED if the rightmost transverse mode before or
  after the action lies outside 0.1–2.0 Hz **and** abs(Δα_full − Δα_EM) ≥ τ_mat.
- **Combined label** (the five prompt categories): UNRESOLVED > MODE-SWITCH >
  FAST-MODE DOMINATED > EM-STABILIZING / EM-DESTABILIZING (EM-NEUTRAL separate).
- Fast-flagged events never count as EM evidence.

### H15 — H_EM

- **EM-failing.** A portfolio is EM-failing iff its status is UNSTABLE (resolved)
  and its rightmost transverse mode is in 0.1–2.0 Hz.
- H_EM is the set of minimal EM-failing portfolios, and κ_EM is the minimal size.
  It is only defined on a STABLE base. Canonical H and κ are unchanged.
- **Computed** for the H14 core lattices, and for the old E7 P4 V9 census
  (re-analysis, using the stored rightmost frequency).
- **EM removal.** An action makes H4 at P4 STABLE, starting from an EM-failing H4.
- **EM creation.** An action creates a hyperedge of H_EM of size ≤ 3 where the
  nominal κ_EM > 3 (or H_EM was empty). Fast-mode events are reported separately.
- **Question 13 answered YES** only if both EM removal and EM creation occur on the
  core lattices (any policy).

### H16 — static topology scores against EM truth

- **Scores:** ΔFiedler(L_B), ΔKirchhoff index, ΔgSCR(H4), Δmin SCR(core), Δ(max
  branch abs(S) flow), Δ(total active losses), Δ(mean effective resistance among
  {30,33,35,37,39}).
- **Target:** Δα_EM of H4 (matched, not fast-flagged), per policy.
- **Rule** (old E7). "Predicts" if abs(Spearman) ≥ 0.6 in ≥ 75 % of the evaluated
  policies (policies with ≥ 10 matched actions).
- Also against Δα_full (old truth) on the same policies.

### H17 — genuine dynamic alternative converter (decided before running)

- **Inventory result** (repository search, 2026-09-12). One existing, validated
  alternative dynamic converter exists: the **TX3 WECC library GFL chain**
  PLL2 + REGCP1 + REECB1 + REPCA1 (+ BusFreq).
  - ANDES 2.0.0 library equations, vendor tag v2.0.0.
  - Frozen parameter set `configs/tx3_gfl_common.yaml`, freeze_id TX3-GFL-0.1,
    sha256 `f07a6a40…54623d5`.
  - Validated on IEEE-39 against ParaEMT (TX3 E02C gates G1–G7 PASS; QZ spectra
    agree with ANDES EIG to 3e-13).
- **No validated GFM exists in the repository.** The REGF1 cases in `temp/` and in
  `artifacts/tx3/E02_andes_mix60.json` are uncalibrated or failed. GFL11 in
  `.venv/xtool-andes-gfl` is a transcription of the custom GFL, not an alternative.
- **Frozen composition** (model ALT-WECC):
  - **Network and machines.** The TX4 network, and the SG2AX transcription of the TX4
    machine, AVR and PSS in ANDES (`.venv/xtool-andes-gfl`). This was validated in
    PCV06 against the internal DAE: abs(Δα) ≤ 1.3e-6.
  - **Converted buses.** The replaced SG is removed and the TX3 chain is attached with
    **every TX3-GFL-0.1 value unchanged**. REGCP1 Sn is the StaticGen Sn of the
    replaced unit, following the TX3 rule that the rating comes from the static
    generator record.
  - **Converted-bus equilibrium.** The TX3 forced-PQ formulation
    (`dicgrid.adapters.andes.enforce_gfl_pq_equilibrium`): the unit's base-PF P and Q
    are held.
  - **Loads.** Constant power.
  - **Policies.** The alternative model has no g. Policies enter only through the
    machine coordinates (k, t, h): KA and TE of each SG2AX, exported from the
    internal model at that policy (exact per-bus values). This is the old E23
    precedent.
- **Nothing is retuned.** If any limit is active at the initial point, the case is
  flagged LIMIT_ACTIVE and excluded from linear verdicts (never adjusted).
- **Implementation-validity gates**, run first and in this order:
  - **Q1.** At each of the 8 H18 policies, the ANDES all-SG base reproduces the
    internal α_⊥ within 1e-4 s⁻¹ and the same status.
  - **Q2.** The ALT-WECC core-lattice cases initialize with a residual ≤ 1e-6 in
    ≥ 80 % of cases.
  - **Q3.** The descriptor spectrum (ANDES Jacobians, QZ) agrees with ANDES EIG.mu
    to ≤ 1e-6 relative on the P4 H4 case.
  - **Q4.** The structural treatment is fixed:
    - remove eigenvalues of identically zero state rows (dead states);
    - then remove the structural pair (two smallest abs(λ) < 1e-3), with all others
      ≥ 1e-2;
    - otherwise the case is BOUNDARY_OR_UNRESOLVED.
- **Gate consequences.**
  - Q1 fails ⇒ H17/H18 BLOCKED (the machine side is not reproduced).
  - Q2 or Q3 fails ⇒ H18 INCONCLUSIVE.

  Nothing is re-tuned to pass a gate.
- **Status verdicts** use the TX4 classifier bands on the transverse spectrum: min
  abs(Re) ≥ 2e-3 is resolved (PCV06 rule).

### H18 — small cross-model matrix (ALT-WECC vs custom GFL)

- **Policies.**
  - Old holdout: H02, H04, H05, H07 (the first four base-stable old-holdout IDs).
  - New holdout: HARDENING_H01–H04 (by ID).
- **Evaluated on both models:**
  - the core lattice (16);
  - the 12 frozen lines [3, 9, 15, 17, 22, 26, 31, 32, 35, 39, 42, 44] (cross-tool
    seed 20260911) at γ = 1.5 and γ = 1 ± 1e-3 (finite total derivative);
  - corridors TXother, TXall, K2_01 (the old top-3) at primary budget L1 = 0.5;
  - 6 topology actions.

  All on H4.
- **Topology actions**, fixed from old custom-GFL **discovery** data only (D01, D03,
  D11, D14). The median Δα(H4) is taken among admissible actions whose post-action
  α(H4) < 1 s⁻¹ at all four (an EM-scale filter).
  - Stabilizing: dbl45, dbl41, out2.
  - Destabilizing: out33, out26, out20.

**Tests.** "Tested policies" means those whose ALT base is STABLE.

| id | test | TRANSFERS | PARTIAL | MODEL-SPECIFIC |
|---|---|---|---|---|
| A | stable-context reversal on the core lattice | ≥ 0.50 of tested policies | > 0 but < 0.50 | 0 |
| B | same-mode reversal (only if ALT has an EM-band critical mode in ≥ 50 % of stable core contexts; else NOT_MEANINGFUL) | ≥ 0.50 | > 0 but < 0.50 | 0 |
| C | 12-line finite ranking, median over policies of Kendall(custom, ALT) Δα(H4) | ≥ 0.50 | 0.20–0.50 | < 0.20 |
| C2 | within-ALT GOLD-B analogue, median ρ(ALT FD-total, ALT finite) − ρ(S1_absP, ALT finite) (descriptive, 12 lines) | ≥ 0.20 | 0–0.20 | < 0 |
| D | corridors: top-1 of the 3 identical in the fraction of policies | ≥ 0.75 | ≥ 0.50 | otherwise |
| E | topology: sign agreement of Δα(H4) over pairs with abs(custom) ≥ τ_mat | ≥ 0.75 | 0.50–0.75 | < 0.50 |

Exact numerical equality is never required.

### H19 — master evidence table: "does a context-independent node ranking suffice?"

Eight layers, each positive or not by a frozen rule on the new data:

| layer | question | positive if (new holdout unless stated) |
|---|---|---|
| L1 | arbitrary-context reversal | ≥ 0.75 of base-stable policies have a stable-context global reversal (old A1 rule) |
| L2 | nested-context reversal | ≥ 0.75 have a level-A nested reversal |
| L3 | same-mode reversal | ≥ 0.50 have a level-B nested reversal (old A2 bar) |
| L4 | EM-only reversal | gate H3 (level C ≥ 0.75) |
| L5 | optimal fixed ranking failure | gate H4-new |
| L6 | ranking transfer failure | the H4 transfer rule |
| L7 | envelope robustness | fresh draws, core lattice at P4: stable-context reversal coverage ≥ 0.50 in ≥ 3 of 4 envelopes (old A4 rule) |
| L8 | dynamic-converter holdout | H18-A = TRANSFERS (BLOCKED or INCONCLUSIVE counts as not positive) |

**STRONG claim "weakness is contextual"** requires ≥ 5 of 8 positive, **including L2
and L5**. Otherwise the claim is stated at the strength the positive layers allow.

### H20–H22 — novelty boundary, literature, reviewer attack

- Documents, no gates.
- **Literature** (H21). Only sources whose DOI or publisher page was retrieved
  during the review are listed. arXiv only when no peer-reviewed version exists.
- "First" appears only if the gap matrix supports it, and then only with the scope
  of the search stated.

### H23–H30 — paper

- **Headline contributions** are limited to results that meet their gate.
- **Figure 8** (topology) goes in the main text iff question 13 is YES (H15 EM removal
  and EM creation both occur). Otherwise it goes to the supplement.
- **Corridors** become a contribution only if gate H10 passes.
- **"Spectral closure"** is never a contribution.

## 3. Stopping, blocking and reporting

- **Execution failure.** A phase whose tasks fail beyond 2 % is FAILED and blocks
  only its dependents. A scientific gate failure blocks nothing.
- **No rerun with changed definitions.** Execution reruns are deterministic and are
  logged.
- **Statuses** for the final claim matrix: PROVED (theory only), STRONGLY_SUPPORTED
  (gate passed with ≥ 20 % margin, and it holds in the EM stratum), SUPPORTED,
  NEGATIVE, INCONCLUSIVE and BLOCKED.
- **Wording.**
  - A deterministic policy or envelope fraction is "coverage across the tested …";
    never a probability.
  - Static metrics are never called universally useless.
  - Total sensitivity and topology switching are never called novel in isolation.
  - No EMT claims.
  - No Africano/PV result.
  - No simple-cycle or connected-cumulant causal story.
- **Negative results** are reported with the same prominence as positive ones.
- **Determinism.** The H1 census at HARDENING_H01 and the H6 link condition at
  HARDENING_H01/H4 are rerun into separate stores and must be identical on every
  summary field except wall-clock.
