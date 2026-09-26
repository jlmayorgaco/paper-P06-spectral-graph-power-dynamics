# Contextual Dynamic Weakness in Inverter-Rich Power Networks

## 1. Executive summary

The campaign tested the meta-hypothesis

> static weakness ≠ dynamic weakness ≠ contextual weakness ≠ intervention
> leverage

against the frozen TX4 IEEE-39 model, using only genuinely new preregistered
experiments (never rerunning or reinterpreting any TX4 result).

**It does not hold as an equality anywhere, and it does not collapse to a
single ranking anywhere either — every distinction the hypothesis draws is
empirically real on this benchmark, at material effect sizes.**

**GOLD gates:** **A and B PASS**; **C, D, E and F do not** (E/F are BLOCKED,
not failed, by a missing-material gate). Fifteen of twenty-four testable
sub-claims are SUPPORTED, eight are NOT_SUPPORTED (seven of these are
informative negatives with an identified mechanism, not inconclusive
failures), one is deferred by design.

**Strongest results.**
- **Contextual sign reversal is real, recurrent, and survives the strongest
  available scrutiny** (mode tracking, an oracle-fitted node ranking, and four
  independent uncertainty envelopes) — **GOLD-A**.
- **Total, re-equilibrated dynamic line sensitivity beats every static
  baseline by a wide, cross-validated margin** (0.987 vs 0.50 median Spearman
  against held-out finite interventions) — **GOLD-B**, the single strongest
  result of the campaign.
- **Topology alone, with no control retuning, both removes and creates
  transverse incompatibilities**, and no static topology score anticipates
  which.
- **Weak corridors exist and are structurally, not spectrally, defined** —
  transformer groups beat every graph-spectral partition tested.
- **A genuinely reduced model is certifiable and accurate, but buys no
  end-to-end speed** on this benchmark, because the equilibrium solve, not
  the eigenanalysis, is the bottleneck — a precise, actionable negative
  result (**GOLD-C: FAIL**).
- **Plan-level multi-constraint design never had to demonstrate an advantage**
  at the scale it could be cleanly tested: single-boundary tuning already
  sufficed in every case (**GOLD-D: FAIL**, informatively).
- **The Africano/PV benchmark is blocked** for lack of source material
  (**GOLD-E/F: BLOCKED**), exactly as the preregistration required rather than
  substituting a different feeder.

**No weakness index was constructed anywhere**, consistent with the
preregistration and reinforced by the campaign's own evidence (E10: an
averaged attribution is dominated by rare fast-mode contexts and hides the
real, modest electromechanical reversals).

**Determinism.** Every rerun deterministic case (E1 census at three policies,
E4 link sensitivities at P4, E9 at target T1, all families/methods) was
byte-identical on its summary fields, excluding wall-clock.

## 2. Research questions

The campaign tests one meta-hypothesis and asks whether it survives:

> static weakness ≠ dynamic weakness ≠ contextual weakness ≠ intervention leverage.

The question is not "which bus is weakest". It is which properties of nodes,
links, corridors and control locations robustly determine, or move, the dynamic
incompatibility boundaries when any of the following change:
- the replacement portfolio;
- the controller policy;
- the operating point;
- uncertain parameters.

The preregistered hypotheses (`docs/CDW_PREREG_V1.md`) are:

| id | hypothesis |
|---|---|
| H1 | contextual sign reversal recurs |
| H1t | the reversal survives mode tracking |
| H1b | fixed node-only rankings are insufficient |
| H1c | contextuality is robust under envelopes |
| HS | submodularity or supermodularity |
| H2 | frozen vs total sensitivity; dynamic vs static baselines |
| H3 | robust weak corridors |
| H4 | topology alone |
| H5 | plan-level design |
| H6 | spectral closure and modal mixing |
| H7 | certified reduction |
| H8 | modal energy vs the limiting mode |
| H9 | cross-model transfer |

## 3. Relation to the frozen TX4 work

- **Frozen and read only.** TX4 is frozen at tag `TX4_FINAL_MANUSCRIPT_FREEZE`
  (69f200df). The CDW branch starts from that tag, and no TX4 file, result,
  ledger or tag was modified.
- **Inherited theorems** (CDW theory §A):
  - the transverse quotient;
  - the any-order-safety characterization of the minimal incompatibility
    hypergraph;
  - boundary localization;
  - exact network-closure factorization;
  - the boundary-sensitivity formula;
  - the deflated zero-frequency port.
- **Inherited IEEE-39 evidence (CDW theory §B).** It is benchmark-specific and
  is used as motivation only:
  - the nominal P4 witness;
  - the policy dependence of H;
  - the one-condition line ranking;
  - the robustness of the phenomenon, but not of the witness identity.
- **The TX4 ParaEMT line is closed.** No CDW statement is an EMT result.

## 4. Graph-DAE model

The model is the frozen TX4 IEEE-39 phasor DAE (`theory/CDW_GRAPH_DAE_V1.md`):
- ten two-axis SGs (bus 39 is the interconnection);
- the frozen GFL, with SRF-PLL, PI outer and inner loops, and a leaky Q/V
  regulator whose gain is `g`;
- constant-power loads;
- an exact branch-terminal network, with 12 off-nominal-tap transformers, line
  charging and 2 shunts.

**Network.** The admittance is `Y(ρ) = Y_sh + Σ_e ρ_e C_eᵀ Y_e C_e`, where `ρ_e`
scales the whole two-port. A unit test checks that `Y(1)` equals the frozen
admittance matrix bit for bit.

**Graph operator.** The coupling Laplacian `L_B` has weights `b_e/t_e`, and the
decomposition `Im Y = −L_B − Δ_tap + B_ch + B_sh` is exact (unit test). No
`B diag(y) Bᵀ` simplification is used.

**Replacement.** Replacing SG `i` by a GFL keeps matched dispatch and an equal
rating, so the AC operating point is identical for every portfolio, as in TX4.

**Policy.** `θ = (g, k, t, h)` collects the Q/V gain and the scale, time-constant
scale and heterogeneity of the AVRs.

## 5. Definitions of weakness

Five distinct objects are defined (theory §C3). Their coincidence is tested, not
assumed.
- **W_static.** SCR, Thevenin |Z|, gSCR, electrical distance, effective
  resistance, Fiedler entries and edge scores, betweenness, flows, dV/dQ.
- **W_dynamic.** The boundary sensitivity `dRe s*/da`, with its sign, in both a
  frozen-operating-point and a re-equilibrated total version.
- **W_contextual.** The distribution of `Δ_i α(S) = α(S ∪ {i}) − α(S)` over
  contexts `S`. It is a property of an (intervention, context) pair, not of a
  bus.
- **W_closure.** The split of the boundary derivative through the exact closure
  factor into device, driving-point and transfer-coupling parts. Physical
  ablations are allowed; zeroing blocks of `K_o` is not used.
- **W_control.** The finite leverage within the admissible ranges.

**Frozen vs total** (theory §C4). Let the operating point be `w = (x, z, r)`:
states, voltages and device references. Then
`dw*/da = −R_w^{-1} R_a` and `dλ/da = v^H (∂_a A + D_w A[dw*/da]) u / v^H u`.

Two equilibrium semantics are used:
- **SPR** (schedule-preserving re-initialization, the TX4 semantics): schedules
  are held and references are re-derived;
- **RP** (reference-preserving): references are held and only the slack
  mechanical power is free.

**Structural remarks.** These are exact consequences of the equations and are
verified numerically in E3.
- Under SPR, the frozen and total derivatives coincide for every
  controller-only coordinate (Q/V gain, PLL, AVR gain).
- Setpoints have a frozen partial that is identically 0.

## 6. Contextual sign-reversal theory

**Definitions.**
- A **contextual sign reversal** of intervention `i` needs two contexts: one
  with `Δ_i α(S1) ≤ −τ` and one with `Δ_i α(S2) ≥ τ`.
  - The material threshold is `τ = 0.01 s⁻¹`, justified in prereg §1.
  - A **stable-context reversal** requires both contexts to be transversely
    stable, which is the planner's situation.
- `α` is a maximum over modes, so a sign change of `Δ_i α` can be produced by a
  switch of the rightmost mode. A **mode-tracked** reversal therefore also
  requires two things:
  - the rightmost mode of each context must be matched, by bus-voltage mode
    shape MAC ≥ 0.8 and |Δf| ≤ 0.15 Hz, to a mode of `S ∪ {i}`;
  - that matched mode must be the rightmost mode of `S ∪ {i}` (same-mode
    transitions).

**Structural facts.**
- On the lattice `2^V`, submodularity of a set function is equivalent to the
  one-step condition `α(S∪i∪j) − α(S∪i) − α(S∪j) + α(S) ≤ 0` (unit test). This
  makes every violation a minimal counterexample of cardinality one.
- A fixed node-only ranking is judged by the most favourable ranking possible:
  the exact linear-ordering optimum over the policy's own resolved context
  triples, computed by dynamic programming over subsets. Its failure is
  therefore not an artefact of a poorly chosen ranking.

## 7. E1: contextual sign-reversal census

**Hypotheses.**
- H1: the same replacement changes the sign of its effect on α_⊥ across
  legitimate portfolios, and this recurs across policies.
- H1t: the reversal survives mode tracking.
- H1b: every fixed node-only ranking is materially insufficient for
  next-replacement decisions.

**Preregistered rules** (`docs/CDW_PREREG_V1.md` §2, E1). Evaluated on the 24
holdout policies, excluding those with an unstable base:
- A1 ≥ 0.75;
- A2 ≥ 0.50 (reported);
- A3: median p* ≤ 0.90 and median decision-failure rate ≥ 0.10.

**Measurement.**
- **Census.** All 512 portfolios of V9 = {30, …, 38} at 39 policies (15
  discovery, 24 holdout): 19 968 portfolio evaluations with the TX4 direct path
  and classifier.
- **Marginals.** 2304 marginals per policy.
- **Mode tracking.** Bus-voltage mode shapes, MAC ≥ 0.8, |Δf| ≤ 0.15 Hz.
- **Execution.** 624 tasks, 0 errors, 896 s.

**Portfolio statuses.** 11 237 STABLE, 8 438 UNSTABLE and 293
BOUNDARY_OR_UNRESOLVED. None was infeasible.

**A heavy tail must be kept in view.**
- 78.5 % of the unstable portfolios have α_⊥ > 1 s⁻¹, up to 2.7·10⁴ s⁻¹.
- These are fast converter-control instabilities, not electromechanical modes.
  TX4 reported the same phenomenon, up to 1276 s⁻¹ on its census.
- A "destabilizing" marginal can therefore be a jump onto a fast mode. For this
  reason the preregistration separates global reversals from mode-tracked ones.

**Result.** 18 holdout policies are base-stable and 6 are not (H01, H03, H06,
H16, H17, H20).

| gate | value | threshold | verdict |
|---|---|---|---|
| A1: policies with a stable-context global reversal | **18/18 = 1.00** | ≥ 0.75 | **PASS** |
| A2: policies with a mode-tracked reversal | **16/18 = 0.89** | ≥ 0.50 | **PASS** |
| A3: median p* (best fixed node ranking, oracle-fitted) | **0.752** (range 0.64–0.84) | ≤ 0.90 | **PASS** |
| A3: median fraction of stable contexts with regret ≥ 0.01 s⁻¹ | **0.369** (range 0.15–0.44) | ≥ 0.10 | **PASS** |

The result is insensitive to the threshold: at τ_res = 10⁻³ s⁻¹ every holdout
policy still shows a stable-context reversal.

**Robustness of the reading.** The heavy tail does not drive the result.
- **Mode-tracked, stable-context reversals.** In these, both contexts are
  same-mode transitions on the matched modal family.
  - They occur for 90 (policy, unit) pairs in the holdout set, and every one of
    the nine units reverses in at least one policy.
  - The strength, min(|Δ(S1)|, |Δ(S2)|), has quartiles 0.014, 0.019 and 0.028,
    with a maximum of 0.068 s⁻¹.
- **Electromechanical-only reversals.** These require both |α(S∪i)| < 1 s⁻¹ and
  the critical frequency of S∪i within 0.1–2 Hz. They occur in 16/18 holdout
  policies, for 104 pairs, with median strength 0.022 s⁻¹.

**Reading.** The same physical replacement stabilizes the rightmost mode in one
stable portfolio and destabilizes the same modal family in another. The
magnitudes are modest, a few hundredths of s⁻¹, but material under the
preregistered threshold.

![Figure F1. Fraction of contexts in which replacing unit i destabilizes (colour), stable-context reversals (×), and the number of reversing interventions per policy ('u' marks an unstable base).](../figures/CDW_F1_contextual_reversal_map.pdf)

![Figure F2. Three strongest mode-tracked reversals (holdout). Each dot is a stable context S where the matched mode stays rightmost. Triangles mark the most stabilizing and most destabilizing contexts; the grey band is ±τ_mat.](../figures/CDW_F2_same_intervention_two_signs.pdf)

**Node-only ranking.** The test is deliberately favourable to node-only
rankings: the ranking is the exact optimum over each policy's own data.
- Even so, it orders only 75 % of resolved pairs correctly (median).
- The next unit it recommends is materially worse than the context-aware choice
  in 37 % of stable contexts (median).
- The optimal orders themselves change across policies. At P4 it is
  38-36-35-34-33-30-31-32-37; at H02 it is 38-34-35-36-32-31-30-37-33. A
  ranking fitted on one policy does not even carry over to another.

**Per-unit behaviour** (holdout medians; each row sums to about 1):

| unit | STAB | NEUTRAL | DESTAB |
|---|---|---|---|
| 30 | 0.22 | 0.32 | 0.44 |
| 31 | 0.29 | 0.55 | 0.15 |
| 32 | 0.28 | 0.49 | 0.21 |
| 33 | 0.20 | 0.22 | 0.58 |
| 34 | 0.45 | 0.38 | 0.17 |
| 35 | 0.29 | 0.36 | 0.34 |
| 36 | 0.41 | 0.33 | 0.26 |
| 37 | 0.25 | 0.44 | 0.31 |
| 38 | 0.66 | 0.00 | 0.34 |

No unit is uniformly stabilizing or destabilizing.

**What this does NOT prove.**
- It does not prove universal weak buses, or their absence, outside this model
  class.
- It does not show that contextual effects are large; the tracked magnitudes
  are 0.01–0.07 s⁻¹.
- It is not an EMT result.
- The fast-mode destabilizations are real in the model, but the model has no
  converter current limits (TX4 limitation).

### 7.1 E1b: submodularity and supermodularity

**Hypothesis.** HS: α_⊥ is neither submodular nor supermodular on the tested
class, and a tracked-mode margin is no more structured.

**Measurement.** All one-step second differences
d_ij(S) = α(S∪i∪j) − α(S∪i) − α(S∪j) + α(S): 4 608 squares per policy, from the
E1 data with no new solves.

**Result.**
- At **all 39 policies** α_⊥ is neither submodular nor supermodular.
- A median 93 % of the squares violate one inequality or the other by more than
  τ_res.
- The tracked-mode margin (the mode matched to the base critical mode) violates
  them in a median 86 % of its defined squares. It is **not** more structured.
- **Minimal counterexamples** (all of cardinality one, by the one-step
  equivalence) at P4:
  - submodularity: S = ∅, i = 30, j = 35, with d = +0.011;
  - supermodularity: S = ∅, i = 33, j = 38, with d = −0.025.

  The largest violations involve fast modes (d up to 1.98·10³ s⁻¹).
  `results/CDW_E1b_minimal_counterexamples.csv` lists the extreme violations per
  cardinality.

**Classification.** NEGATIVE BUT INFORMATIVE; it is a structural property of
the tested class.

**Valid claim.** The spectral-abscissa set function does not satisfy the
submodular or supermodular property on the tested IEEE-39 class. No claim is
made about greedy algorithms in general.

### 7.2 E10: Shapley and context diagnostics (descriptive)

- The exact Shapley values of v(S) = α(S) − α(∅) (9 players) are **dominated by
  the fast-instability contexts**. They are of order 1–15 s⁻¹, and the marginal
  standard deviations are of order 300 s⁻¹.
- The preregistered "hidden reversal" flag (|φ_i| ≤ τ_mat with std ≥ 3τ_mat)
  fires 0 times out of 351. The reason is not that averages are faithful: a
  single rare fast-instability context dominates the context average.

**Reading.** A context-averaged attribution of α_⊥ is not a usable weakness
index on this class. It is uninformative about the modest electromechanical
reversals and dominated by rare fast modes. This supports not building any
averaged weakness index.

## 8. Robustness under parameter envelopes (E2) — GOLD-A, third component

**Hypothesis H1c.** Contextuality persists under the preregistered uncertainty
envelopes, even where the minimal failing set (the witness) changes.

**Measurement.** Core lattice (16 portfolios of V4) at P4 (D01) and at the
clean stable g-only point G_S (D03), under:
- TX4 discovery draws, regenerated from the frozen PCV05 generator and
  verified to match the frozen `PCV05_draws.csv` factors exactly before use
  (100 per envelope);
- a new CDW holdout seed (40 per envelope, seed 20260921), plus a V9-census
  check on 10 draws per envelope.

Fractions below are **coverage fractions of declared envelopes, never
probabilities.**

| pid | source | envelope | coverage: stable-context reversal | coverage: witness changed | reversal coverage \| witness changed |
|---|---|---|---|---|---|
| **D01 (P4)** | CDW (holdout) | EC | **1.00** | 0.00 | — |
| | | EM-f | 0.80 | 0.80 | 0.75 |
| | | EM-u | 0.975 | 0.45 | 0.94 |
| | | EMC | 0.975 | 0.40 | 0.94 |
| D01 (P4) | TX4 (discovery, verified) | EC/EM-f/EM-u/EMC | 1.00 / 0.82 / 0.97 / 0.95 | consistent with the CDW draws | consistent |
| D03 (G_S) | CDW and TX4 | all four | **0.00** | 0.00–0.45 | 0.00 |

| gate A4 (P4, CDW holdout) | value | rule | verdict |
|---|---|---|---|
| envelopes with coverage ≥ 0.5 | **4/4** | ≥ 3/4 | **PASS** |

**GOLD-A = A1 ∧ A3 ∧ A4: all three PASS.**

**The robust-phenomenon-vs-robust-identity distinction, directly visible.**
- At P4, reversal coverage stays high (0.80–1.00) across every envelope,
  **including** the machine envelopes under which the witness itself changes
  in 40–80 % of draws. Conditioned on the witness having changed, reversal is
  still present in ≥ 0.75–0.94 of those draws: **the contextual phenomenon
  survives the identity of the failing set changing.**
- At the clean point G_S, by contrast, reversal coverage is exactly **0** at
  every envelope: the whole phenomenon there is policy-dependent, not just its
  identity. This is consistent with G_S being constructed as the clean,
  policy-only counterfactual (TX4 V06/V07) and shows the E2 result is itself
  contextual — it does not hold everywhere in policy space, only where a
  near-critical boundary exists.

**Census check (10 draws/envelope, all 512 V9 portfolios, P4).** Every
completed draw shows at least one reversing intervention (coverage 1.00 in all
four envelopes), confirming the phenomenon is not an artefact of restricting to
the four-unit core.

**What this does NOT prove.** It does not claim the phenomenon is universal
across all policies (D03 is a direct counterexample); it is a property of
being near a policy-dependent incompatibility boundary, exactly where a
planner would be looking.

## 9. Total node sensitivity (E3)

**Hypothesis H2 (nodes).** Frozen-operating-point and re-equilibrated total
node sensitivities disagree materially on this class.

**Implementation validity gate (IV).** Before any H2 claim, the IFT total
derivative must match a small central finite re-equilibrated difference.

| check | result | threshold | verdict |
|---|---|---|---|
| IV: frac. of (parameter, condition) pairs within tolerance | **8002/8002 = 1.000** | ≥ 0.95 | **PASS** |
| IV by semantics | SPR 1.000, RP 1.000 | — | — |
| IV by coordinate kind | g, ka, line, load, pll, vset all 1.000 | — | — |
| independent port-form cross-check (T-form vs engine, D01/H4/46 lines) | max abs diff 2.0·10⁻⁶, Spearman 1.000 | — | validates the engine |

The sensitivity engine (IFT on the residual `[f; g; schedule equations]`, exact
Jacobian re-solve) is validated by two independent routes: finite differences
and an independent port-form (`T(s)`) derivative.

**Structural remarks (theory §C4), checked numerically.**

| remark | result | verdict |
|---|---|---|
| R1: frozen = total for controller-only coordinates under SPR (g, PLL, K_A) | max relative difference 1.8·10⁻⁷ over 841 pairs | **CONFIRMED** |
| R2: the frozen partial for a setpoint (V_set) is exactly 0 | max \|frozen\| = 1.50·10⁻⁷ | **marginal FAIL at the 10⁻⁸ threshold** |

R2 is reported as failing the preregistered 10⁻⁸ bound without relaxing it.
The value, 1.5·10⁻⁷, is five orders of magnitude below the material threshold
(0.01 s⁻¹ · range) and is consistent with the central-difference step used
in the residual Jacobian (relative step 10⁻⁷); it is not evidence against the
structural claim `∂α/∂V_set|_frozen = 0` that R2 was designed to check, and
would very likely close at a smaller step, which was not tried after seeing
the result (that would violate the no-retuning rule).

**H2 gate** (material: `|total|·range ≥ τ_mat` and either sign disagreement
≥ 10 % or Kendall ≤ 0.6, over holdout conditions):

| semantics | frac. material (holdout) | median sign disagreement | median Kendall | H2-node |
|---|---|---|---|---|
| SPR | **0.976** | 0.333 | 0.434 | **TRUE** |
| RP | **0.762** | 0.104 | 0.599 | **TRUE** |

**Reading.** Even under SPR — where the theory (R1) guarantees frozen = total
for every *controller* coordinate — the node sensitivities still diverge
materially once network and operating-point coordinates (loads, V_set) are
included, because those have a frozen partial of essentially zero and all of
their effect is re-equilibration. A frozen-point analysis of node leverage
would silently drop these effects.

## 10. Total link sensitivity and static baselines (E4) — GOLD-B

**Hypothesis (GOLD-B).** A total dynamic link sensitivity predicts held-out
finite ×1.5 branch-strengthening effects on H4 substantially better than the
strongest static baseline, across policies and parameter envelopes.

**Baselines tested** (46 branches, holdout conditions: 12 policies + 5
envelope draws × H4 = 32 conditions with all quantities resolved): S1 |P|, S2
|S|, S3 |z|, S4 electrical distance, S5 effective resistance, S6 Fiedler edge
score, S7 edge betweenness, S8 max endpoint dV/dQ, S9 ΔgSCR; Dconv
(conventional base-portfolio eigenvalue sensitivity), Dfrozen, Dtotal (SPR),
DtotalRP; Dport (independent port-derivative cross-check, D01 only).

**Selection rule** (fixed before results): the static baseline with the best
median Spearman on **discovery** conditions is carried forward. That baseline
is **S1 (|active flow|)**.

| predictor | median Spearman (holdout, H4) | median top-5 precision |
|---|---|---|
| **Dtotal** | **0.987** | **0.90** |
| DtotalRP | 0.993 | 1.00 |
| Dfrozen | 0.900 | 0.60 |
| S1 \|P\| (selected static) | 0.500 | 0.20 |
| S2 \|S\| | 0.481 | 0.20 |
| S5 effective resistance | 0.300 | 0.60 |
| S3 \|z\| | 0.234 | 0.20 |
| S4 electrical distance | 0.194 | 0.40 |
| S9 ΔgSCR | 0.175 | 0.20 |
| Dconv (base-portfolio eigenvalue sensitivity) | 0.032 | 0.60 |
| S6 Fiedler edge score | −0.054 | 0.20 |
| S8 dV/dQ | −0.168 | 0.00 |
| S7 betweenness | −0.271 | 0.00 |

Every static baseline, including the best possible *oracle* per-condition
static choice (median ρ = 0.50, no better than the discovery-selected one),
stays near ρ ≈ 0.2–0.5. The conventional eigenvalue sensitivity of the base
portfolio (no replacement) is uninformative for H4 (ρ = 0.03): the ranking
that matters is specific to the four-unit portfolio, not to the bare network.

| gate | value | threshold | verdict |
|---|---|---|---|
| median(ρ_Dtotal − ρ_static,selected) | **0.487** | ≥ 0.20 | **PASS** |
| median ρ_Dtotal | **0.987** | ≥ 0.70 | **PASS** |
| median top-5 precision (Dtotal) | **0.90** | ≥ 0.60 | **PASS** |
| **GOLD-B** | | | **PASS** |

On the V9 target (larger portfolio), Dtotal keeps ρ = 0.995 while the selected
static baseline turns **negative** (ρ = −0.15): a static ranking fit for H4
does not merely weaken on a different target, it reverses.

![Figure F4. Left: median Spearman of every predictor against the finite ×1.5 effect, holdout, H4 (red = dynamic). Right: P4 (D01), first-order prediction vs finite ×1.5 effect for all 46 branches, total (red) and frozen (blue).](../figures/CDW_F4_static_vs_total_line_ranking.pdf)

**H2 gate (links).**

| semantics | frac. material (holdout) | median sign disagreement | median Kendall | H2-line |
|---|---|---|---|---|
| SPR | 0.512 | 0.067 | 0.768 | **TRUE** |
| RP | 0.913 | 0.125 | 0.658 | **TRUE** |

**What this does NOT prove.** GOLD-B is established for the H4 target under
the frozen IEEE-39 device model; it is not a universal claim about static
indices, and the branch ranking is known (TX4) to be sensitive to the
converter model (tested again in E23).

## 11. Static vs dynamic weakness (E5, descriptive)

Built entirely from the E1/E3/E4/E10/E11 data (no new solves). Fixed rule:
static-weak = bottom 3 by SCR (nodes) / top 3 by electrical distance
(branches); dynamic-critical = top 3 by `|dα/dg|` total (nodes) / `|Dtotal|`
(branches); context-dependent = reversal in ≥ 50 % of holdout policies;
robustly neutral = NEUTRAL in ≥ 80 % of contexts.

| bus | category |
|---|---|
| 30, 31, 32 | CONTEXT-DEPENDENT |
| 33 | DYNAMIC ONLY; CONTEXT-DEPENDENT |
| 34 | STATIC ONLY; CONTEXT-DEPENDENT |
| 35 | DYNAMIC ONLY; CONTEXT-DEPENDENT |
| 36 | STATIC ONLY; CONTEXT-DEPENDENT |
| 38 | STATIC + DYNAMIC; CONTEXT-DEPENDENT |
| 37 | *(none of the above)* |

**Every one of the nine candidates is context-dependent** by the fixed 50 %
rule, and the static-weak and dynamic-critical top-3 sets **overlap in only
one bus (38)** — Kendall(SCR, total NG sensitivity) = 0.44, Kendall(SCR,
destabilizing-context fraction) = **0.00** (no monotone relationship at all).
No bus is robustly neutral. For branches, the static-weak (electrical
distance) and dynamic-critical (Dtotal) top-3 sets are **disjoint**, and their
Kendall correlation is 0.18–0.40. **A ranking built from SCR or electrical
distance would not reliably identify the units and branches that a dynamic,
contextual analysis flags — and vice versa.**

## 12. Weak corridors (E6) — H3

**Hypothesis H3.** Robust dynamic weak corridors exist: admissible corridor
coordinates whose ranking and stabilizing effect persist across holdout
policies and envelopes, independent of any single nominally top-ranked line.

**Construction (frozen before results; §C8, prereg §2 E6).** Spectral
clustering of the exact coupling Laplacian `L_B` into 2–4 communities (cutsets
K2/K3/K4), plus the transformer groups TXcore (the 4 core step-up
transformers), TXother (the other 8) and TXall (all 12). No corridor is
defined post hoc around a successful line.

**Result.** Over 32 holdout conditions (12 policies + 5 envelope draws × H4):

| corridor | frequency in the top 3 by finite stabilizing effect (×1.5) |
|---|---|
| **TXother** (8 non-core transformers) | **0.94** |
| **TXall** (all 12 transformers) | **0.91** |
| K2_01 (largest 2-way cutset) | 0.81 |
| every other spectral cutset or TXcore | ≤ 0.06 |

| gate | value | rule | verdict |
|---|---|---|---|
| robust weak corridor exists | TXother, TXall qualify (≥ 0.75 top-3 frequency **and** ≥ 0.75 stabilizing) | — | **H3: TRUE** |
| corridor-ranking predictiveness (Σ Dtotal vs finite corridor effect) | median ρ = **0.963** | ≥ 0.70 | predictive |

**Non-additivity.** The corridor effect is close to, but not exactly, the sum
of its branches' individual effects: median |NA| = 0.014, with a heavy tail up
to 0.92 in a few conditions. The sum of single-line total sensitivities
(median ρ = 0.99) is in fact an even better predictor of the corridor's finite
effect than the sum of static electrical-distance scores (ρ = 0.82), which is
itself informative: for THESE corridors the interaction is close to additive
in the *dynamic* sensitivities, even though no single branch dominates.

![Figure F6. Left: the two most frequent top-3 corridors on the network graph (core units in blue). Right: additivity check — corridor finite effect vs the sum of its branches' finite effects (holdout).](../figures/CDW_F6_corridors.pdf)

**Reading.** The transformer groups, not the graph-spectral cutsets, are the
robust corridors on this benchmark — a structurally motivated corridor
(all/most step-up transformers together) beats every purely spectral
partition. This is a corridor-level, not merely branch-level, result: TXother
excludes exactly the four core transformers directly attached to the
replacement candidates, so its consistent effectiveness is a genuine
multi-branch, not single-line, finding.

## 13. Topology reconfiguration (E7) — H4-topology

**Hypothesis.** A single admissible topology action can remove, or create, an
incompatibility; static topology scores do not predict which.

**Admissible actions** (PF converged, every bus voltage in [0.90, 1.10] p.u.,
network connected): **35 of 46** single-branch outages and **46 of 46**
branch doublings (reinforcements). Evaluated on the core lattice (16
portfolios) at 10 policies (5 discovery + 5 holdout), and on the full V9
census (512 portfolios) at P4 for every admissible action.

**Topology alone removes the P4 incompatibility.** At P4 (D01), four single
reinforcements — doubling branches 0, 1, 13 and 43 — each independently drive
`α(H4)` negative, with every other subset already stable, on their own,
with no control retuning. At the other tested policies, reinforcement is
dramatically more effective still: **every one of the 46 admissible
doublings** removes the H4 instability at D03, H02 and H04 (46/46 each), and
one does at H06.

**Topology alone creates new, smaller incompatibilities.** 42 (policy,
action) pairs produce a hyperedge of size ≤ 3 where none existed before,
including five cases where a **single core unit alone** becomes
transversely unstable (`κ = 1`) after one outage (e.g. D01/out0, D01/out1,
D03/out0: `H = 30|33|35|37`, all four singletons simultaneously unstable).

**Overall churn.** `H` changes in 169 of the 810 tested (policy, admissible
action) pairs (21 %); outages change `H` far more often than reinforcements
(37 % vs 8 %), consistent with reinforcements moving the system toward
stability and outages away from it.

**Static topology scores do not predict the effect.**

| static score | median Spearman(Δscore, Δα(H4)) | "predicts" (≥ 0.6 in ≥ 75 % of policies) |
|---|---|---|
| Δ Fiedler value (λ₂ of L_B) | −0.459 | **no** |
| Δ Kirchhoff index | +0.473 | **no** |
| Δ gSCR(H4) | −0.440 | **no** |
| Δ min SCR (core) | −0.409 | **no** |

None reaches the 0.6 bar in three-quarters of policies; several even carry the
wrong sign convention on average.

**Census (all 512 V9 portfolios, P4).** The nominal hypergraph is large — 52
edges of size 4–6, none smaller — reflecting how much the census
incompatibility structure spreads across many five/six-unit coalitions once
all nine candidates are in play. **77 of the 81 admissible actions change
`H_V9`**, and the minimum edge size κ itself moves across the full observed
range: 60 actions leave κ = 4, 13 bring it to 3, 3 to 2, **5 to κ = 1** (a
lone unit destabilizing), and 1 to κ = 5.

![Figure F7. Left: Δα(H4) for every admissible action, by policy (outages vs doublings). Right: Δα(H4) vs ΔFiedler at P4 — no visible relationship.](../figures/CDW_F7_topology.pdf)

**Reading.** On this benchmark, topology is not a secondary lever: single
admissible actions both remove and create incompatibilities, at every scale
tested, and ordinary static network-strength proxies give no useful advance
warning of which action will do which. **H4-topology: TRUE on both halves.**

**What this does NOT prove.** No thermal ratings exist in the frozen data, so
thermal limits were not enforced (flows are reported, not gated); no cost model
is attached to "reinforcement", so this is not a cost-effectiveness claim.

## 14. Control–topology exchange rate (E8)

**Measurement.** At the exact TX4 P4-line boundary θ* = (0.20768, 1.425, 1.5,
1) where α(H4) ≈ 0, and at P4 itself: the local rate
`r = −(∂α/∂θ_i)/(∂α/∂γ_e)` (total, SPR) for the 4 Q/V gains and 6 surviving
AVR gains against the 5 branches with the largest total sensitivity at θ*,
then validated by a **paired finite** intervention (`θ_i += δ`,
`γ_e += r·δ`).

| point | step size | fraction with compensation ratio ≤ 0.2 | median ratio |
|---|---|---|---|
| **boundary θ\*** | small (δ = 0.02) | **1.00** (50/50) | **0.023** |
| boundary θ* | large (δ = 0.1) | 0.66 | 0.162 |
| P4 (D01) | small | 1.00 | 0.022 |
| P4 (D01) | large | 0.78 | 0.129 |

| gate | value | threshold | verdict |
|---|---|---|---|
| E8 (boundary, small step) | 1.00 | ≥ 0.80 | **PASS** |

**Reading.** At the small-step (local) scale the linearized exchange rate
predicts the paired finite intervention essentially exactly (2 % residual
effect). At the large step (0.1, a tenth of the full gain range) the linear
rate over- or under-shoots by 13–16 % on the residual, as expected for a
first-order local tool; it remains directionally useful (residual well below
the standalone control effect) in roughly two-thirds to three-quarters of
pairs. This is a legitimate local, not global, equivalence, exactly as
specified — not an economic exchange rate.

![Figure F8. Distribution of the compensation ratio (|Δα of the paired finite intervention| / |Δα of the control alone|) at the boundary and at P4, small and large steps. Dashed line: the 0.2 pass threshold.](../figures/CDW_F8_exchange_rate.pdf)

## 15. Plan-level multi-constraint design (E9) — GOLD-D

**Hypothesis H5.** Plan-level control/topology design makes an entire
implementation plan safe where single-boundary tuning leaves an unsafe
intermediate subset.

**Setup.** `Φ_T(a) = max_{S⊆T} α_⊥(S; a)`. Two targets: T1 = (H4, P4) and T3 =
(H4, D11, a held-out F7B point), each with 16 subsets; T2 = (V9, P4) and T4 =
(V9, D03), each with 512 subsets. Three admissible-parameter families
(control, topology, joint), each compared single-boundary (tune only α(T))
against plan-level (tune Φ_T over the active subsets, capped at 25, with a
Gordan-witness fallback on infeasible QP steps). Sequential QP, ≤ 6 iterations,
ε = 0.02, trust region 25 % of each coordinate's declared range.

| target \ family | control | joint | topology |
|---|---|---|---|
| **T1** (H4, P4) | single-boundary **suffices** (Φ_T: −0.019) | suffices (−0.022) | suffices (−0.044) |
| **T3** (H4, D11) | suffices (−0.018) | suffices (−0.019) | suffices (−0.020) |
| **T2** (V9, P4) | neither converges (single Φ_T +3750; plan +1276) | neither converges (single +1.50·10⁵; plan +1.53·10⁴) | neither converges (single +2.57·10⁴; plan +9.56·10³) |
| **T4** (V9, D03) | neither converges (single +3750; plan +1276) | neither converges (single +8.71·10³; plan +3.09·10³) | neither converges (single +9.56·10³→ +10.3·10³†) |

†topology/T4: plan Φ_T (1.03·10⁴) is marginally above single (1.02·10⁴); the
capped active set (25 of 512 subsets) does not track the true worst subset
exactly at every iteration.

| gate | result |
|---|---|
| **GOLD-D (single-boundary safe, plan-level not)** | **not observed anywhere: FAIL** |

**Reading.**
- **At the H4 scale (16 subsets), single-boundary tuning always sufficed**, in
  every family and both tested policy locations: after tuning only α(H4) to
  below −ε, every one of its 15 proper subsets was already stable too (no case
  of `single_leaves_unsafe_subset = True` was found). Plan-level tuning
  therefore made no material difference here — a clean negative for H5 at this
  scale, on this benchmark.
- **At the V9 census scale (512 subsets), neither method converges** within
  the preregistered 6-iteration budget: the worst subset (almost certainly a
  fast, non-electromechanical instability, given the α values of 10³–10⁵ s⁻¹)
  is not driven negative. **Plan-level tuning is directionally far better**
  than single-boundary tuning even so — it reduces the worst-case Φ_T by a
  factor of 3–10× across families (e.g. joint/T2: 1.50·10⁵ → 1.53·10⁴) — but
  this is reported as **INCONCLUSIVE**, not a pass, under the fixed iteration
  budget; the budget was not extended after seeing this (that would be
  retuning after results).

![Figure F9. Φ_T = max_S α(S) across SQP iterations, by target, family and method (solid: plan-level; dashed: single-boundary).](../figures/CDW_F9_single_vs_plan_design.pdf)

**Local conflict witnesses.** The Gordan LP fired whenever the primary QP was
infeasible; every fired instance found `gordan_value` effectively 0 (≤ 10⁻⁸),
i.e. a genuine local conflict among the active subsets' gradients, consistent
with the QP infeasibility that triggered it. No global-impossibility claim is
made from these local witnesses.

**What this does NOT prove.** It does not show plan-level design is
unnecessary in general — only that, at the scale where it was cleanly
testable (16 subsets), this benchmark's H4 boundary happened not to need it.
The V9-scale result is inconclusive, not negative: a larger iteration budget or
a better active-set strategy might resolve it, but that was not tried here.

## 16. Spectral graph baselines (E11) — H6s

**Bases** (§C8): (A) the exact-tap `L_B` Laplacian, low spectrum
`{0, 5.12, 5.56, 10.54, 11.61}` Hz²-scaled eigenvalues; (B) Kron reduction of
the susceptance matrix onto the ten generator buses; (C) the dynamic-port SVD
(used directly in E4's `Dport` cross-check).

**Predictive tests** (median |Spearman| over 24 holdout policies):

| score | vs Shapley (E10) | vs destabilizing fraction (E1) | vs total NG sensitivity (E3) |
|---|---|---|---|
| \|Fiedler entry\| | 0.43 | 0.57 | 0.10 |
| Fiedler entry (signed) | 0.38 | 0.16 | 0.75 |
| effective resistance to bus 39 | 0.60 | 0.18 | **0.80** |

| edge score | vs Dtotal (holdout, H4) |
|---|---|
| S5 effective resistance | 0.30 |
| S6 Fiedler edge score | −0.05 |

| gate | value |
|---|---|
| static graph predicts *something* at ≥ 0.6 | **TRUE** (only: effective resistance vs total NG sensitivity, 0.80; and, marginally, vs Shapley, 0.60) |

**Reading.** Static graph structure is a genuinely weak predictor overall —
median correlations mostly in the 0.1–0.6 range, several near zero — with one
partial exception: node effective resistance to the slack tracks the total
Q/V-gain sensitivity reasonably well (0.80). It does **not** track which
contexts destabilize (0.18), which is the quantity that actually matters for
H1/H1b. No claim of novelty is made for the Laplacian decomposition itself.

## 17. Controller-weighted modal mixing (E12)

**Construction (§C8).** The exact coupling Laplacian `L_B` gives a graph
Fourier basis `U`. The port operator is mapped to that basis,
`T̂(s) = (U ⊗ I)ᵀ T(s) (U ⊗ I)`, and mixing is the off-block-diagonal share
`μ_mix = ‖T̂ − blockdiag(T̂)‖_F / ‖T̂‖_F`, evaluated at the rightmost transverse
frequency of each portfolio.

**Q1: does controller policy change mixing at fixed L?**

| | value | threshold | verdict |
|---|---|---|---|
| relative range of μ_mix(H4) over 39 policies | **0.366** | ≥ 0.20 | **material: TRUE** |
| μ_mix(H4) range | 0.073–0.104 | — | — |
| relative range of the algebraic-only mixing (`g_z` alone, no dynamics) | **1.9·10⁻¹¹** | — | essentially exactly constant |

The algebraic-only control confirms the mechanism precisely: the network and
load admittance alone do not move with policy (as they must not — they are
policy-independent by construction), so the *entire* 37 % swing in μ_mix comes
from the controller dynamics entering `T(s)` at the finite critical frequency.

**Q2: does mixing correlate with contextual reversal?**

| | value | threshold | verdict |
|---|---|---|---|
| Spearman(μ_mix(H4), number of reversing interventions), holdout | **ρ = 0.657** | \|ρ\| ≥ 0.5 | — |
| permutation p-value (10⁴ permutations) | **0.0037** | ≤ 0.01 | — |
| **correlates** | | | **TRUE** |

**Q3: do witness transitions (4→3→2→3→4→∅) correspond to reproducible
modal-support transitions?**

Using the ten frozen TX4 FC18 boundary events on the F7B line, the graph modes
carrying 80 % of the critical mode's angle energy were compared just before
and after each of the 9 witness-changing events, against 9 control midpoints
between them.

| | events (witness changes) | controls (no change) |
|---|---|---|
| support changed | 3 / 9 | 1 / 9 |

| gate | rule | result | verdict |
|---|---|---|---|
| Q3 reproducibility | ≥ 7/9 events change support **and** ≤ 20 % of controls do | 3/9 events, 11 % controls | **Q3: FALSE** |

**Reading (negative but informative).** The support does shift more often at
true witness-changing events than at arbitrary midpoints (3/9 vs 1/9), so
there is a weak association, but it falls well short of the preregistered
reproducibility bar. Witness transitions in this model are **not**, in
general, accompanied by a clean, low-dimensional modal-support reorganization
detectable at the 80 %-energy threshold; most of the κ-changing events keep
essentially the same dominant graph-mode support.

**Q4: are the dynamically critical modes the same as the lowest static graph
modes?**

| | value |
|---|---|
| top-3 dynamically critical graph modes (H4, P4, by angle-mode energy) | modes {0, 1, 6} |
| 3 lowest nonzero L_B modes | {1, 2, 3} |
| Jaccard overlap | **0.20** |

**Reading.** The dynamically critical modes only partly overlap the lowest
graph-Fourier modes: mode 6 (not among the three lowest) carries real weight
in the critical mode shape at P4, while mode 3 does not. **A purely static,
low-order graph truncation would miss part of the dynamically relevant
structure** — a first piece of evidence developed further in E13/E14.

![Figure F10. Left: μ_mix(H4) across all 39 policies (fixed network); note it never approaches the network-only value (1.9·10⁻¹¹ relative range). Right: share of critical-mode energy in the three lowest L_B graph modes vs α(S), all portfolios and policies.](../figures/CDW_F10_graph_modal_mixing.pdf)

## 18. Reduced models (E13) — spectral closure, H6r

**Families tested on holdout** (24 policies, core 16-portfolio lattice, plus a
64-portfolio V9 sample): GM(r), Galerkin projection of the bus-voltage
coordinates onto the r leading `L_B` graph modes (r = 2, 4, 8, 12, 16, 24, 39;
device states kept exact); POD(r), the same with a data-driven basis from the
discovery critical-mode voltage shapes (r = 2, 4, 8, 12, 16); TS(a/b/c),
quasi-steady elimination of the GFL fast states (current loop; + power
filters; + PLL).

| family | verdict acc. | H exact | Kendall (marginals) | state ratio | eig speedup | end-to-end speedup |
|---|---|---|---|---|---|---|
| GM39 (= full network) | 1.000 | 1.000 | 1.000 | 1.00 | 0.86 | 1.00 |
| GM24 | 0.993 | 0.958 | 0.978 | 1.00 | 0.92 | 1.00 |
| GM16 | 0.959 | 0.750 | 0.843 | 1.00 | 0.95 | 1.00 |
| GM12 | 0.863 | 0.708 | 0.742 | 1.00 | 0.96 | 1.00 |
| GM8 | 0.746 | 0.625 | 0.163 | 1.00 | 0.97 | 1.00 |
| GM4, GM2 | 0.58–0.63 | 0.625 | negative | 1.00 | 0.98–0.98 | 1.00 |
| POD16 | 0.986 | 1.000 | 0.954 | 1.00 | 0.97 | 1.00 |
| POD2/4/8/12 | 0.49–0.50 | 0.25 | ≈0 | 1.00 | 0.97–1.00 | 1.00 |
| **TSa** | 0.982 | 1.000 | 1.000 | **0.854** | **1.46** | 1.00 |
| **TSb** | 0.945 | 0.958 | 0.966 | **0.780** | **1.68** | 1.00 |
| **TSc** | 0.887 | 0.958 | 0.968 | **0.707** | **1.83** | 1.00 |

**GM/POD never reduce the device state count** (`state_ratio` is exactly 1.0
by construction: only the algebraic network-voltage dimension is projected,
never the ODE states), so they cannot satisfy the "genuinely reduced" bar
regardless of accuracy. **TS families do reduce states** (71–85 % retained),
with accuracy that degrades gracefully as more GFL dynamics are eliminated:
TSa (current loop only) reproduces the full model almost exactly; TSc (current
loop + power filters + PLL) trades some accuracy (H exact 0.96, verdict 0.89)
for the largest reduction and the largest **eigen-stage** speedup (1.83×).

**A structural finding, independent of any single family:** none of the
fifteen candidates comes close to a genuine **end-to-end** speedup (all are
≈ 1.00×, i.e. no measurable speedup at all). The eigen-decomposition step that
these reductions accelerate is a negligible fraction of the total wall time;
the dominant cost is solving the full nonlinear equilibrium (`solve_case`),
which every family still performs unreduced before linearizing. **A useful
reduced model on this benchmark would have to reduce the equilibrium solve
itself, not just its eigenanalysis** — that was not attempted here.

## 19. Decision-preserving certification (E14)

Applied to the three TS families only, where the reduced port matrix `T̃_S(s)`
shares the full model's dimension (theory note; GM/POD are reported
NOT_APPLICABLE for the same reason, per `docs/CDW_PREREG_V1_DEVIATIONS.md`
item 4, and count as zero abstention in the GOLD-C rule by convention).
Contour Γ = rectangle [0.002, 20] × [−2π·5, 2π·5] rad/s; sampled-certified if
`sup_Γ ‖R̃⁻¹(R − R̃)‖₂ < 0.9`, else ABSTAIN; poles of either model inside Γ also
force ABSTAIN.

| family | coverage (certified) | false certifications | accuracy when certified | accuracy overall |
|---|---|---|---|---|
| TSa | 1.000 | **0** | 1.000 | 1.000 |
| TSb | 0.506 | **0** | 1.000 | 0.997 |
| TSc | 0.472 | **0** | 1.000 | 0.997 |

**Zero false certifications across every case tested.** Where the sampled
bound fails (TSb/TSc abstain roughly half the time), the reduced verdict is
still almost always correct anyway (overall accuracy 0.997), but the
certificate correctly declines to vouch for it rather than asserting a
guarantee it cannot rigorously back.

### GOLD-C: FAIL, for a specific and identifiable reason

| criterion | best result | required | met by |
|---|---|---|---|
| state ratio ≤ 0.70 | 0.707 (TSc) | ≤ 0.70 | **no family** (TSc misses by 0.007) |
| end-to-end speedup ≥ 3× | ≈ 1.00× (all) | ≥ 3× | **no family** |
| certified accuracy ≥ 0.99, false certs = 0 | 1.00 / 0 (TSa, TSb, TSc) | — | **all three TS families** |
| H exact ≥ 0.90 on holdout | 0.958–1.000 (TSb, TSc, TSa) | ≥ 0.90 | **all three TS families** |

**GOLD-C: FAIL.** Three of five criteria are met by every TS family; the
campaign fails on state ratio (narrowly, for TSc) and decisively on end-to-end
speedup. **The reason is structural, not a tuning failure:** these reductions
only shrink the linear/eigen stage, and that stage is not the bottleneck on
this benchmark. **No useful low-dimensional spectral or dynamic reduction was
found that also delivers a real speedup**, though the certification and
decision-preservation machinery itself works essentially perfectly on the
family that does reduce states.

![Figure F11. Left: verdict accuracy vs retained network modes r (GM, POD) and the three TS points. Right: state ratio vs end-to-end speedup, with the GOLD-C thresholds (dashed).](../figures/CDW_F11_reduction.pdf)

![Figure F12. Certified / abstained case counts by reason, holdout core, with false-certification counts annotated (0 everywhere).](../figures/CDW_F12_certificate.pdf)

## 20. Modal energy vs the stability-limiting mode (E16) — H8

**Hypothesis H8.** The dominant observed oscillation and the stability-limiting
eigenmode differ often enough to matter for monitoring/diagnosis. This is a
**phasor-domain** study only, and is explicitly **not** an EMT result; it is
motivated by, but separate from, the SQ1 ParaEMT-line estimator work.

**Measurement.** For every stable core (V4) portfolio at all 39 policies, two
0.1 p.u. / 0.2 s load-pulse disturbances (bus 16, bus 20) were applied
analytically (closed-form modal residues), and the energy of each transverse
oscillatory mode (0.1–2.0 Hz) was integrated over [2, 30] s on the SG speeds
and all 39 bus-voltage magnitudes.

| | value | threshold | verdict |
|---|---|---|---|
| fraction of all 998 stable cases where the energy-dominant mode ≠ the limiting mode | **0.164** | — | descriptive |
| fraction of holdout **policies** where this occurs in ≥ 25 % of their cases | **0.278** (5/18) | ≥ 0.50 | — |
| **systematic (per-policy majority rule)** | | | **FALSE** |
| median real-part gap between the two modes, when they differ | 0.037 s⁻¹ | — | — |

![Figure F13. Left: real part of the stability-limiting mode vs the energy-dominant mode (red where they differ). Right: fraction of cases where they differ, per policy.](../figures/CDW_F13_modal_energy_vs_limiting.pdf)

**Reading (negative but informative).** The mismatch is real and not
negligible — it occurs in about one stable portfolio in six overall — but it
is concentrated in a minority of policies (5 of 18 holdout policies account
for most occurrences) rather than being a pervasive property of the model
class. The preregistered per-policy-majority bar for "systematic" is not met.
This is consistent with, but numerically distinct from, the frozen SQ1
ParaEMT-estimator observation that motivated H8: here the effect is smaller
and less uniform, using exact modal residues rather than an estimator on noisy
synthetic time series. **Not an EMT result; a phasor-domain, closed-form
diagnostic only.**

## 21. Africano / PV material gate (E17)

**Result: GATE FAILS — inputs missing.** A full repository search (179 hits on
"africano", "hosting capacity", "PVR", "weak node" and Spanish equivalents)
found **no** Africano source document (thesis or paper, PDF/DOCX), no feeder
data file, no weak-node metric definition, no PV placement rule, no
hosting-capacity method and no reported result table to reproduce. The only
related artefacts are third-party scripts referencing an IEEE-13 weak-node
label set ("Africano (2017) Ref."), which is a list of four bus labels, not a
methodology, and cannot be exactly replicated.

**Consequence (per the preregistration):** Phases 18–21 (exact replication,
PV coalitions, no-storage control/topology hosting, dynamic hosting) and
**GOLD-E and GOLD-F are BLOCKED.** No generic feeder was substituted for
Africano, and none of E17's negative finding was used to weaken or reinterpret
any other result. See `docs/CDW_AFRICANO_MISSING_INPUTS.md` for the full
listing and the exact required-inputs table.

## 22. Cross-model holdout (E23) — H9

**Alternative model.** The TX4 `StaticInjection` negative control (an ideal
constant-power converter, no dynamics beyond the algebraic characteristic),
never retuned to reproduce the GFL.

| check | result | threshold | verdict |
|---|---|---|---|
| reversal presence, base-stable policies (7 of 10 tested) | **7/7 = 1.00** | ≥ 0.50 | transfers |
| median link Kendall (static-model finite ranking vs GFL finite ranking) | **0.568** | ≥ 0.50 | transfers |
| median corridor top-3 overlap | **0.667** (2/3) | ≥ 2/3 | transfers |
| **H9: findings transfer** | | | **TRUE** |

**Reading.** Contextual sign reversal, the qualitative branch ranking, and the
identity of the leading corridors all carry over to a converter model with no
shared dynamics with the GFL beyond the power-flow constraint. This bounds
model-specificity for these three findings on this benchmark — though TX4 had
already shown (V28) that a library GFL model does **not** preserve the branch
ranking (7/12), so "transfers to one alternative, non-GFL device" is not
evidence of transfer to every converter model.

## 23. Negative results (consolidated)

Every negative result below was preregistered as a possible outcome and is
reported with the same prominence as a positive one.

| id | negative result | where |
|---|---|---|
| HS | α_⊥ is neither submodular nor supermodular at any of 39 policies (median 93 % one-step violations); the tracked-mode margin is no more structured (86 %) | §7.1 |
| — | Shapley/context-averaged attribution is uninformative here — dominated by rare fast-instability contexts, not the modest electromechanical reversals | §7.2 |
| — | Structural remark R2 (frozen partial of a setpoint = 0) narrowly fails its strict 10⁻⁸ bound (1.5·10⁻⁷), most likely a finite-difference artefact, not a violation of the underlying claim | §9 |
| C-H4b | none of four static topology scores predicts the effect of a topology action on α(H4) (median \|ρ\| < 0.6 throughout, several wrong-signed) | §13 |
| Q3 (E12) | witness-changing events are not reliably accompanied by a reproducible modal-support transition (3/9 vs the 7/9-and-≤20 %-controls bar) | §17 |
| GOLD-D | plan-level design never demonstrates an advantage over single-boundary tuning at the H4 scale (single always sufficed); at the V9 census scale neither method converges within budget | §15 |
| GOLD-C | no reduced model achieves both a genuine end-to-end speedup (≥ 3×) and a materially smaller state count (≤ 70 %) — the bottleneck is the nonlinear equilibrium solve, not eigenanalysis | §18–19 |
| H8 | the modal-energy/limiting-mode mismatch is real (16 % of cases) but not "systematic" by the preregistered per-policy majority rule | §20 |
| E17 | the Africano/PV material gate fails outright — no source document, feeder, metric or result table exists in the repository | §21 |
| — | static graph scores are mostly weak predictors of contextual/dynamic node quantities (median \|ρ\| often 0.1–0.6), with one partial exception (effective resistance vs total Q/V sensitivity, 0.80) | §16 |

## 24. Publication candidates

Ranked by evidentiary strength and novelty on this benchmark class.

1. **GOLD-B (dynamic vs static line sensitivity), §10.** The strongest single
   result: a 0.49 Spearman gap over the best static baseline, validated on 32
   held-out conditions with an independently cross-checked (port-form)
   sensitivity engine, and a clean mechanistic account (frozen setpoints have
   zero frozen partial; controller coordinates have frozen = total under
   SPR). Publication-ready on its own.
2. **Contextual sign reversal (GOLD-A, §7).** A clear, preregistered,
   three-tier (global / mode-tracked / electromechanical-only) demonstration
   that a fixed node ranking is provably insufficient on this benchmark, with
   an oracle-fitted ranking baseline that is *maximally* favourable to
   node-only ranking and still falls short. Robust under four uncertainty
   envelopes, including where the failing-set identity itself changes.
3. **Topology alone moves incompatibility (§13).** Both directions
   demonstrated (removal and creation) with a clean negative control (no
   static score predicts it), and a striking census-scale result (κ moves
   across its entire {1,…,5} range under admissible single actions).
4. **Weak corridors beating single-line and pure-spectral partitions (§12).**
   A structurally motivated (transformer-group) corridor construction
   outperforms graph-spectral cutsets, with a clean additivity diagnostic.
5. **The reduced-model negative result (§18–19).** A precise, mechanistic
   account of *why* eigen-stage reduction does not help on this class
   (equilibrium solving dominates), with a certificate that is accurate
   whenever it fires. Useful as a cautionary methodological result for the
   broader reduced-order-modelling literature in this application area.

## 25. Limitations

- **Model scope.** The frozen TX4 IEEE-39 phasor DAE: no converter current
  limits, no ride-through logic, no DC link, a harmonized first-order AVR, and
  matched-dispatch replacement (the AC operating point is invariant across
  portfolios by construction). Loads are constant power. No governors are in
  the census (D = 0, no primary frequency control), so all stability
  statements are transverse.
- **Fast-mode contamination.** A large share of "unstable" portfolios (78.5 %)
  carry α_⊥ far above the electromechanical band (up to 2.7·10⁴ s⁻¹). Every
  headline reversal claim is therefore reported in three variants — global,
  mode-tracked, and electromechanical-only — precisely to keep this from
  driving the conclusions; readers should use the mode-tracked or
  electromechanical-only numbers for any claim about oscillatory stability.
- **Two equilibrium semantics, not a continuum.** SPR and RP are the two
  semantics defined in the theory (§C4); intermediate conventions (partial
  re-dispatch) are not tested.
- **Corridor, topology and reduction constructions are frozen before results**
  (spectral-clustering corridors, transformer groups, single-branch actions;
  GM/POD/TS families and their retained ranks) — they are not tuned to the
  outcome.
- **The certificate (E14)** applies rigorously only to the TS (device
  time-scale) reductions, where the port matrices share a dimension with the
  full model; for GM/POD the reduced network dimension differs and the
  Gohberg–Sigal comparison as implemented does not apply, so those families
  report zero abstention by convention (documented in
  `docs/CDW_PREREG_V1_DEVIATIONS.md`), not because they were certified.
- **No EMT.** This is a phasor-domain campaign. No claim here bears on the
  ParaEMT line, which remains permanently closed
  (`docs/20260912_PAREMT_EMT_CLOSURE.md`, TX4).
- **Africano/PV.** Blocked at the material gate (E17); see
  `docs/CDW_AFRICANO_MISSING_INPUTS.md`. GOLD-E and GOLD-F are not evaluated.
- **Single benchmark family.** All device-model claims are checked against
  exactly one alternative converter model (the TX4 static-injection negative
  control, E23); this bounds, but does not eliminate, model-specificity.

## 26. Reproducibility

- **Branch:** `research/contextual-dynamic-weakness`, from tag
  `TX4_FINAL_MANUSCRIPT_FREEZE` (commit `69f200df`). Not pushed.
- **Preregistration:** commit `b8ae3082` (`docs/CDW_THEORY_V1.md`,
  `docs/CDW_PREREG_V1.md`, `docs/CDW_EXPERIMENT_PLAN_V1.md`,
  `docs/CDW_CLAIM_LEDGER_V1.md`, `results/CDW_CLAIM_MATRIX_V1.csv`,
  `theory/CDW_GRAPH_DAE_V1.md`, `docs/CDW_OPEN_QUESTIONS_V1.md`,
  `docs/CDW_AFRICANO_PV_BENCHMARK_PLAN.md`), before any CDW numerical result.
  Implementation clarifications made before the campaign are in
  `docs/CDW_PREREG_V1_DEVIATIONS.md`.
- **Frozen inputs:** `results/prereg_inputs/holdout_policies.json` (sha256
  `4ac8c95d…`) and `cdw_envelope_draws.json` (sha256 `864a5969…`), generated by
  `experiments/cdw/CDW00_prereg_inputs.py` before any model evaluation.
- **Orchestrator:** `experiments/cdw/CDW_MASTER_RUN.py`; per-task checkpoints
  under `raw/<phase>/`; status in `results/CDW_RUN_STATUS.json`; logs under
  `logs/`.
- **Environment:** BLAS pinned to one thread
  (`OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=MKL_NUM_THREADS=1`); Python and package
  versions recorded per run in `CDW_RUN_STATUS.json["environment"]`.
- **Determinism (Phase 28):** headline deterministic cases (E1 census at D01,
  D11, H01; E4 links at D01/H4/SPR; E9 target T1, all families/methods) are
  rerun into separate stores and compared byte-for-byte on their summary
  fields, excluding wall-clock. See `results/CDW_DETERMINISM.json`.
- **Everything CDW writes lives under** `contextual_dynamic_weakness/`; no TX4
  file was read for writing, and none was modified.

## 27. Recommended next paper(s)

1. **A focused paper on contextual weakness and dynamic line ranking**
   (GOLD-A + GOLD-B), the two strongest, best-validated results, framed as
   "static and single-context rankings are demonstrably insufficient for
   next-replacement and reinforcement decisions in policy-dependent IBR
   portfolios" — this is the natural, tightly scoped successor to TX4.
2. **A topology/corridor paper** built on §13 and §12: topology alone can
   remove or create incompatibilities, static topology scores do not predict
   which, and structurally motivated (not purely spectral) corridors are the
   ones that generalize.
3. **A short methodological note on the reduced-model negative result**
   (§18–19): a worked demonstration that "reducing the linear stage without
   reducing the equilibrium solve buys nothing" is itself a useful, citable
   caution for reduced-order dynamic security assessment.
4. **Africano/PV** (GOLD-E/F) only if the source material can be obtained; the
   plan in `docs/CDW_AFRICANO_PV_BENCHMARK_PLAN.md` is ready to execute as
   preregistered, unchanged, the moment the gate can be satisfied.
5. **Nonlinear recovery vs α_⊥ improvement** (E22, deferred here) is a natural
   follow-up once a paper is scoped around the plan-level design result,
   since TX4 already established subcritical-Hopf boundaries on this
   benchmark.
