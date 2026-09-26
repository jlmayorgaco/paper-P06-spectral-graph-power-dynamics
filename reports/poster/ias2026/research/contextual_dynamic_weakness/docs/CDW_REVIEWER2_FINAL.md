# Reviewer 2 report: IEEE Transactions on Power Systems

**Manuscript:** "Portfolio-Conditioned Dynamic Weakness and Reinforcement Ranking in Inverter-Rich Power Networks"
**Source reviewed:** `reports/papers/cdw_contextual_dynamic_weakness/main.tex` (compiled `main.pdf`, 10 pp.), `cdw_numbers.tex`, `tab_*.tex`
**Supporting material read:** `CDW_HARDENING_PREREG_V1.md`, `CDW_HARDENING_STATISTICAL_PLAN.md`, `CDW_HARDENING_EXPERIMENT_MATRIX.md`, `CDW_HARDENING_DEVIATIONS.md`; `results/hardening/*` (H03, H04, H05, H06, H07, H13, H14/H15/H16, H18, H19, `CDWH_RUN_STATUS.json`, `CDWH_DETERMINISM.json`, `alt/*`); code in `experiments/cdw_hardening/*.py` and `experiments/cdw/{_cdw,_sens,E01_census,E34_sens}.py`; raw census `raw/H_H01/*.json`; git history of branch `research/contextual-dynamic-weakness-hardening` (commits `05b507e3`, `4de412d2`, `eb57d907`, and the uncommitted working tree).
**Date:** 2026-09-12
**Recommendation:** **Major revision.** C1 (nested same-mode reversal) and C3 (portfolio-conditioned branch ranking) can be supported with the changes below. As worded, C2 (fixed-ranking insufficiency "carried by fast converter-control modes"), the cross-model verdicts, and the flagship illustration are not supported by the result files. Several of the requested checks need only existing data.

---

## 1. Summary of the submission

The paper evaluates all 512 SG-to-GFL replacement portfolios of nine IEEE 39-bus generators at 63 controller policies. Of these, 24 policies are a new maximin-LHS holdout. The quantity studied is the marginal change of the transverse spectral abscissa, Δ_i α_⊥(S). There are four claims:
- **C1.** Nested sign reversals of the same replacement on a tracked EM mode occur in 17/19 base-stable new-holdout policies. A chain identity certifies discrete curvature.
- **C2.** An oracle fixed unit ranking has pairwise agreement p* = 0.75 and material regret 0.38. The failure is attributed to "fast converter-control modes": in the EM stratum, p* = 0.97 and regret = 0.03.
- **C3.** The re-equilibrated (IFT) eigenvalue sensitivity ranks finite branch reinforcements with median Spearman 0.985 (8 eligible conditions) or 0.987 (64 conditions), against about 0.5 for static indices.
- **C4.** Single branch actions both create and remove EM incompatibilities, and static topology scores do not predict them.

A cross-model holdout with a frozen WECC library GFL chain (ALT-WECC) is reported as reproducing "the ranking method but not the reversals or the ranking itself".

The work is careful in several respects:
- a preregistration with frozen seeds;
- a deviation log;
- deterministic reruns;
- numbers generated as LaTeX macros from result files;
- honest reporting of the failed mixing hypothesis and of the EM-stratum weakening of C2.

Most numbers in the text match the files (Section 6). The problems are in:
- what the dominant "fast" events physically are, and whether they are within model validity;
- how the fixed-ranking and static baselines are constructed;
- one mislabelled flagship example;
- the cross-model metric;
- the integrity of the audit trail.

---

## 2. Scores (1-10)

| criterion | score | justification |
|---|---|---|
| Novelty | 5 | Some ingredients are classical: the chain identity, feasible eigenvalue sensitivity, non-submodularity of spectral objectives (Olshevsky 2018), and Braess-type edge effects. The authors credit them. What is new is the systematic census with nested and mode-graded reversals, and the evaluation of the IFT sensitivity as a ranking against finite re-solves. This is a useful empirical contribution, not a conceptual one. C4 is largely expected physics: all 45 EM creations are outages and all EM removals are doublings (H15). |
| Correctness | 5 | The algebra is correct and most macros match the files. However: (i) the Fig. 1 / Sec. VI-A flagship pair is **not** a same-mode reversal (0.744 Hz vs 1.310 Hz; level D false). (ii) The abstract's direction-specific sentence holds in 9/19 policies, not 17/19. (iii) The "fast converter-control modes" are zero-frequency real eigenvalues with Re up to 2.7×10⁴ s⁻¹; their mechanism is not identified. (iv) The ALT "every replacement leaves α_⊥ unchanged or increases it" is an artifact of a pinned pole. (v) Test D uses a base-unstable policy that the prereg excludes. |
| Model adequacy | 3 | The model has D = 0, no governors, constant-power loads, no converter current limits, delays or DC link, a phasor network and matched dispatch. The FULL-stratum C2 result is driven by aperiodic eigenvalues of +10 to +26 563 s⁻¹, far outside phasor validity, whose magnitude co-varies with cond(g_z). The EM reversals (median 0.013 s⁻¹) are smaller than the damping a retired machine with realistic D would remove. There is only one benchmark and one custom GFL. |
| Statistics | 5 | Good framing: designed points, "coverage" not probability, cluster bootstrap, Holm family, effect-size gates. But: (i) T2 (GOLD-B sign-flip) fails under Holm (p = 0.061) and is not reported. (ii) The literal GOLD-B primary set is 8 conditions in 6 clusters, 6 of them at the discovery flagship P4. (iii) The H3 point estimate passes, but its own bootstrap interval (0.74-1.00) crosses the 0.75 gate. (iv) 3263 heavily overlapping pairs are presented as counts, and they include base-unstable policies. (v) "Rankings do not transfer" is inferred from a regret difference of 0.026. |
| Evidence | 5 | The volume is large (32 256 portfolio solves, 4048 IFT/FD checks). But the decisive comparisons are against baselines that cannot, by construction, respond to the conditions. The cross-model test uses a metric that hides the EM mode. The corridor phase (H10) is still RUNNING and is omitted. |
| Reproducibility | 6 | Deterministic, macro-generated numbers and determinism reruns (16+4 tasks identical). Against this: the repository branch is local (no URL or DOI); Table I is incomplete; deviation-log timestamps contradict git and file times; a post-hoc analysis (regret anatomy) is uncommitted and unlogged. |
| Industrial relevance | 3 | A planner choosing the next unit already evaluates each candidate in context: nine eigenanalyses per step, which is the "context-aware" choice itself. The cost argument for the branch sensitivity is not demonstrated (no timings) and misdescribes the implementation. There is no large system and no multi-step plan. |
| Writing | 7 | Concise and precise, with a good crediting paragraph and clear limitations. Weaknesses: the abstract is overloaded with numbers from two different condition sets; some jargon is undefined ("the P4 incompatibility", "harmonized benchmark"); some claims are stronger than the data. |
| Figures | 6 | Figs. 2, 3, 6 and 8 are informative. Fig. 1 carries an incorrect "(same tracked EM mode)" label and its fonts are illegible at column width. Fig. 4 (right) is a tautology. The Fig. 5 heatmap has no policy labels. |

---

## 3. Answers to the five mandatory questions

### Q1. "Isn't it caused by fast unstable modes?"

**For C2: yes, overwhelmingly, and the modes are not what the paper says they are.**

I reproduced the anatomy from `H03_marginals.parquet` (new holdout, level A):
- **All 14 130 out-of-band critical modes are real.** Every stable-context marginal whose post-replacement critical mode lies outside 0.1-2 Hz has |f| ≤ 0.0011 Hz, i.e. a real eigenvalue. None lies above 2 Hz.
- **Two distinct groups.** 4167 are slow stable real modes, with α between −0.24 and −0.05 s⁻¹. These are consistent with the PSS washout (1/T₆ = 0.238) and the Q/V leak (w = 0.05). The other 9886 are real eigenvalues with Re between +10 and +26 563 s⁻¹. The interval (0, 10] s⁻¹ is empty.
- **A penetration cliff.** In the new census the share of portfolios whose rightmost eigenvalue is such a large real one depends almost only on portfolio size: 0 % for |S| ≤ 3, 1 % at |S| = 4, 36 % at |S| = 5, and 93-95 % at |S| = 6-7.
- **Link to the algebraic Jacobian.** Among these cases, log α and log cond(g_z) correlate at 0.81. The median cond(g_z) is 2.1×10³ when the critical mode is oscillatory. It rises to 5.8×10⁴ when α > 10⁴ s⁻¹.
- **Link to the regret anatomy.** In `H04_gate.json`, 2032 of the 2089 "out-of-band" material regrets have Δα ≥ 1 s⁻¹.

This pattern is consistent with a singularity-induced-bifurcation-type phenomenon of the DAE, or with an ill-posed fast loop of the idealized GFL. Either one must be characterized before being called a "fast converter-control mode". It is not established by any participation analysis in the code: a repository search finds no participation-factor computation.

**For C1: largely no, with qualifications.**
- Level-C endpoints are EM-band. The EM-clean witness share is 94 % overall, but 88 % for the s→d direction that violates submodularity (97.5 % for d→s).
- Converter-parameter draws alone (EC) remove every EM same-mode reversal on the core lattice at P4 (0/10).
- In 36 % of level-C pairs, the real-part gap between the rightmost and the next eigenvalue at one of the four endpoint portfolios is smaller than the reversal magnitude itself. In 29 % it is smaller than τ_mat.

So the EM reversals are real but sit at a knife-edge of rightmost-mode identity.

### Q2. "Is the best fixed ranking baseline fair?"

It is favorable for p*, which is what the DP maximizes. It is **not** the most favorable baseline for the headline regret metric, and the regret metric itself is inflated.
- **A regret-optimal ranking does better.** I enumerated all 9! orders per new-holdout policy. The regret-optimal fixed ranking has median regret 0.345, against 0.376 for the DP ranking.
- **Many decisions have no feasible option.** In a median 21.8 % of the "stable contexts" counted, no candidate replacement leaves the portfolio STABLE. Regret there compares two infeasible portfolios, typically real eigenvalues of 10²-10⁴ s⁻¹.
- **A stability screen removes most of the failure.** The screened list takes the first unit whose S∪{i} is STABLE. On contexts that have at least one stable option, it has median material regret **0.065** (range 0.005-0.139). Without the screen the figure is 0.267.

The fixed-ranking failure as defined is therefore mostly "the list does not know which replacement crosses the penetration cliff". A one-eigenanalysis-per-candidate feasibility check fixes most of it, and that check is exactly the in-context evaluation the paper recommends. The comparison should be reported against these baselines (Major M2).

### Q3. "Are policies/envelope draws independent samples?"

No, and the authors mostly say so.
- **Policies.** They are designed LHS points in an author-chosen box: g = u², k ∈ [0.5, 2.3], t ∈ [0.5, 3]. This box contains AVR tunings that destabilize the all-SG system (5/24 new, 6/24 old). Coverage numbers are therefore relative to this box and this parameterization.
- **Envelope draws.** All 40 are at a single policy, P4, which is the discovery flagship. The 4 envelopes are the clusters.
- **Portfolios.** Within a policy, the 512 portfolios are a complete census, not a sample.
- **Pairs.** The 3263 level-C pairs share contexts heavily: only 604 distinct (policy, unit, S₁) marginals among the base-stable pairs. The pair count also includes 196 pairs from three base-unstable policies.

The percentile cluster bootstrap is a reasonable design-stability diagnostic. The p-values (T1-T3) are conditional-exchangeability summaries and should be labelled descriptive. Nothing in the paper supports inference to "inverter-rich networks" in general.

### Q4. "Doesn't the alternative converter invalidate the line ranking?"

Yes for the specific ranking: the median Kendall between models on 12 lines is 0.015. The paper concedes this. The claim that "the method transfers" rests on test C2:
- It is a descriptive test per the prereg.
- It uses 12 lines.
- The median advantage is 0.206 against a 0.20 threshold, and 3 of the 8 policies are below 0.20.
- In ALT the derivative is a central finite difference of the same α(γ) function that defines the truth.

So C2 shows mainly that α_⊥(γ) is smooth in ALT too. The line ranking is model-specific. Because the frozen chain has no voltage control, the cross-model test is not like-for-like for the mechanism the paper attributes reversals to (Major M6).

### Q5. "Why does a planner need contextual ranking instead of exhaustive eigenanalysis?"

The manuscript does not answer this. It should.
- **Unit replacement.** The "context-aware choice" is exhaustive in-context eigenanalysis: at most 9 portfolios per step, trivially cheap on this system. The paper offers no method that beats it. Its finding is diagnostic: do not use a fixed list.
- **Branch reinforcement.** The claimed saving is "one factorization per portfolio and one back-substitution per candidate branch" instead of a re-solve and eigenanalysis per candidate. But `_sens.py::Engine.derivatives` evaluates, per candidate:
  - two residuals;
  - one LU back-solve;
  - **four** full finite-difference A-matrix reconstructions (`A_of` calls `central_difference_jacobians`).
- **No timings.** The finite truth uses the same Newton solver plus one dense eigen-decomposition. No wall-clock comparison is reported. On a 39-bus system with 46 branches, exhaustive finite evaluation costs seconds.

The value of first-order screening appears only at scale:
- thousands of branches;
- multi-step plans;
- contingency sets;
- sparse eigen-triplets with analytic Hessian-vector products.

None of this is demonstrated (Major M8).

---

## 4. MAJOR comments

### M1. The "fast converter-control mode" mechanism behind C2 is mischaracterized and lies outside the model's validity

- **Quoted text.**
  - Abstract: "a failure carried mostly by transitions into fast converter-control modes".
  - Sec. VI-C: "In 2089 of the 2265 materially wrong decisions ... the unit chosen by the ranking moves the rightmost eigenvalue out of the EM band, almost always onto a fast converter-control mode."
  - Limitations: "the fast control-mode instabilities that drive the fixed-ranking failure are those of an idealized converter."
- **What is wrong.**
  - Every out-of-band critical mode in the new census is a real eigenvalue (|f| ≤ 1.1×10⁻³ Hz). Of the 14 130 such stable-context marginals, 9886 have Re(λ) ∈ [+10, +26 563] s⁻¹ and none lies in (0, 10] s⁻¹.
  - These events are a deterministic function of penetration: 0 % of portfolios at |S| ≤ 3, 36 % at |S| = 5, 93-95 % at |S| = 6-7.
  - Their magnitude co-varies with cond(g_z): log-log correlation 0.81, and the median cond(g_z) rises from 2.1×10³ to 5.8×10⁴.
  - No participation factors are computed anywhere in `experiments/`. The "converter-control" attribution and the "fast" qualifier are therefore asserted, not shown.
  - An aperiodic growth rate of 10²-10⁴ s⁻¹ is outside what a phasor model with algebraic network, no current limits and no delays can represent credibly.
- **Fix.**
  1. Identify the mechanism. Report participation factors (or state-group energy) of the dominant real eigenvalue, and track σ_min(g_z) along the replacement sequence. Test whether the eigenvalue passes through infinity, i.e. an SIB, as the replaced rating fraction ρ_i varies. `ReplacementPlan` already supports fractional ρ.
  2. Re-run the census for the new holdout with at least one of the following:
     - converter current limits;
     - a modulation/measurement delay;
     - voltage-dependent (ZIP) loads.

     Report whether the cliff and the FULL regret survive.
  3. Rewrite C2 as "a fixed list cannot anticipate an aperiodic high-penetration instability in this model". Replace "fast converter-control modes" everywhere unless (1) supports it. Cite the SIB/DAE-feasibility literature (e.g., Venkatasubramanian, Schättler and Zaborszky on local bifurcations and feasibility regions in DAE power systems).
- **New evidence needed:** yes. Item (1) needs only existing solves plus eigenvectors. Item (2) needs new runs.

### M2. The fixed-ranking regret metric and baseline overstate the failure, and "rankings do not transfer" is not supported

- **Quoted text.**
  - "This is the most favorable case for a node-only ranking."
  - "materially wrong in 0.38 of next-replacement decisions."
  - "Rankings do not transfer: a ranking fitted at one policy and used at another is materially wrong in 0.40 of decisions."
  - Discussion: "a replacement sequence evaluated with a fixed weak-bus list will choose a materially worse next unit in about a third of stable contexts."
- **What is wrong** (reviewer reanalysis of `H03_marginals.parquet` with the DP orders in `H04_ranking.csv`, new holdout, 19 policies):
  1. The DP maximizes pairwise agreement, not regret. The regret-optimal fixed order, found by enumerating all 9! orders, has median regret 0.345 rather than 0.376.
  2. The denominator includes contexts with no stable next replacement: a median 21.8 % of decisions. "Regret" there compares two unstable portfolios.
  3. Screening the list for feasibility reduces median material regret to 0.065 on contexts that have a stable option. The screened list takes the first ranked unit whose S∪{i} is STABLE.
  4. Transfer is nearly as good as fit. The off-diagonal regret is 0.402 against 0.376 on the diagonal, and off-diagonal accuracy is 0.725 against 0.763. The L6 rule is positive whenever L5 is, so it adds no independent evidence.
  5. "Weak-bus list" in the Discussion misdescribes the baseline: an in-sample oracle is not a weak-bus list.
- **Fix.**
  - Report regret for:
    - the DP order;
    - the regret-optimal order;
    - the stability-screened order;
    - the context-aware choice.
  - Restrict the denominator to contexts with at least one stable option, or report both.
  - Replace "do not transfer" with "transfer about as well as they fit; the failure is present in-sample".
  - Revise the Discussion sentence accordingly.
- **New evidence needed:** analysis only, on existing data.

### M3. The flagship reversal is not a same-mode reversal; the abstract misstates the direction; level B/C do not fully exclude mode switching

- **Quoted text.**
  - Sec. VI-A: "in S₂ = {31,34,37,38} ⊃ S₁ the same replacement moves **the same mode** right by 0.064 s⁻¹."
  - Fig. 1 panel titles: "(same tracked EM mode)".
  - Abstract: "In 17/19 base-stable new-holdout policies the same replacement moves the same tracked electromechanical mode **left in one stable portfolio and right in a portfolio that contains it**."
- **What is wrong.**
  1. In `H03_nested_pairs.parquet` the Fig. 1 pair (N24, unit 33, s2d) has hz₁ = 0.744 Hz and hz₂ = 1.310 Hz, and lvD = False. The critical modes of S₁ and S₂ are different modes, so the example is a level-C pair between two different EM-band modes. Overall, 17 % of level-C pairs are not level D. `strongest_new_C(1)` breaks a three-way tie at |Δ| = 0.042021 arbitrarily. Two of the tied pairs are level D.
  2. The abstract describes the s→d direction ("left ... then right in a superset"). Only 9/19 base-stable new policies have an s→d level-C reversal; 15/19 have d→s. Only 7/19 have both directions, which is what Prop. 1(c) needs to certify "neither submodular nor supermodular".
  3. Mode identity rests on three things:
     - MAC ≥ 0.8 of 78-component rectangular bus-voltage shapes;
     - a ±0.15 Hz window;
     - a candidate set restricted to band modes with Re ≥ −1 plus the critical mode (`_cdw.transverse_eval`).

     There is no participation check that the mode is electromechanical: "EM" is a frequency label, and AVR/PSS or outer-loop modes can sit in 0.1-2 Hz. There is also no continuation check. In 36 % of level-C pairs the real-part gap to the next eigenvalue (`gap2`) at one endpoint portfolio is below the reversal magnitude. Near such degeneracies MAC-based identity is unreliable.
- **Fix.**
  - Choose a level-D example for Fig. 1, or relabel the current one honestly. Break ties deterministically in favor of level D.
  - Rewrite the abstract and C1 in direction-neutral form, e.g. "moves the same tracked mode in opposite directions in two nested stable portfolios". State that both directions occur in 7/19.
  - Add a continuation test: sweep ρ_i from 0 to 1 for each level-C endpoint marginal and track the eigenvalue continuously.
  - Report the MAC margin (best minus second-best candidate).
  - Report rotor-speed/angle participation of the tracked mode.
  - Report coverage restricted to pairs with gap2 ≥ τ_mat at all four portfolios, as a sensitivity analysis analogous to the GOLD-B near-tie exclusion.
- **New evidence needed:** yes, for continuation and participation, on the existing policy set.

### M4. Model adequacy for the EM claims: damping, governors and loads

- **Quoted text.**
  - "Each SG uses a two-axis model on its own base with D = 0, constant mechanical power"; "with constant-power loads".
  - "Converter current limits, DC links and modulation delays are not modeled."
- **What is wrong.**
  - The level-C reversal magnitudes have quartiles 0.011 / 0.013 / 0.017 s⁻¹ (H03), i.e. at the materiality floor.
  - Replacing a machine that carries realistic damping (D of 1-2 pu on machine base) or a governor removes a stabilizing contribution of comparable or larger size from every marginal. This could erase the stabilizing side of most reversals.
  - The transverse-quotient construction exists only because D = 0 and there are no governors.
  - L7 already shows EM reversals disappear under converter-parameter draws (EC: 0/10). Their robustness to the omitted physics is unknown.
  - Matched dispatch (the GFL takes the SG's P and Q and rating) removes the static reactive-support consequence of retirement. That consequence is the dominant practical effect of SG retirement.
- **Fix.** Re-run H3 (levels C and D) on the 19 base-stable new policies with three additions:
  - D ∈ {0.5, 1, 2} pu;
  - a simple governor model;
  - ZIP loads.

  Report coverage and magnitudes. Discuss the matched-dispatch assumption and give one alternative, e.g. unity-PF GFL with the Q deficit moved to the remaining units.
- **New evidence needed:** yes.

### M5. GOLD-B: the baselines are condition-invariant by construction, the literal primary set is tiny, and T2 fails

- **Quoted text.**
  - "ranks finite reinforcements far better than nine static indicators".
  - "Dtot outranks every static index in every condition, including the best one chosen after the fact per condition."
  - "The sensitivity of the base portfolio has median ρ 0.034: the branch ranking belongs to the portfolio, not to the network."
- **What is wrong.**
  1. **The static indices never change.** `E34_sens.static_link_baselines` computes S1-S8 from the base network and base power flow, and S9 from the target set only. All nine static rankings are therefore identical in all 64 conditions. "Best static per condition" is a choice among nine fixed lists. They cannot, by construction, respond to the policy or draw that the paper says determines the ranking.
  2. **A fixed list learned from data is much closer.** The appropriate fixed-list competitor is a data-driven fixed branch ranking for V4. I used the mean truth rank over the other clusters, leave-one-cluster-out, from `H06_links_raw.parquet`:

     | set | fixed list ρ | fixed list top-5 | D_tot ρ | D_tot top-5 |
     |---|---|---|---|---|
     | all 64 conditions | 0.943 | 0.8 | 0.987 | 1.0 |
     | 24 new policies | 0.818 | 0.6 | 0.996 | 1.0 |
     | 40 draws at P4 | 0.959 | not computed | not computed | not computed |

     D_tot still wins in every condition. But the honest gap is 0.04-0.18 in ρ, not 0.49. On the new policies the fixed list (0.818) beats D_fro (0.712).
  3. **The "network" baseline is good where it matters.** The base-portfolio sensitivity D_conv has median top-5 = 0.60 and NDCG@5 = 0.91 (Table V). At the decision-relevant top of the list it is as good as effective resistance and better than |P_e| (0.20 / 0.61). The "belongs to the portfolio, not the network" conclusion rests on ρ over 46 branches, most of which have near-zero effect.
  4. **The literal primary set is small and mostly one policy.** It has 8 conditions: D01|EC|7, D01|EM-f|0, D01|EM-u|0, D01|EM-u|3, D01|EMC|6, D01|EMC|7, HARDENING_H09 and HARDENING_H18. Six of the eight are parameter draws at the discovery flagship, so evidence on held-out *policies* rests on 2 policies.
  5. **T2 fails and is not reported.** The Holm-family test T2 on this set gives p = 0.061 (Holm 0.061 > 0.05; `H19_evidence.json`). The paper reports only T3's p-value. Sec. V-E ("Three secondary tests form a Holm family") gives no hint that one of them failed.
- **Fix.**
  - Add the leave-one-cluster-out fixed branch list and a mode-shape-based branch proxy as baselines.
  - Make D_fro the headline competitor.
  - Report top-5 and NDCG@5 alongside ρ in the text.
  - Report T1, T2 and T3 with raw and Holm p.
  - Simplest remedy for eligibility: re-solve the 56 ineligible equilibria to a tighter Newton tolerance. The engine's own `solve` accepts tol = 1e-10. The preregistered R0 ≤ 1e-8 rule then applies literally to all 64 conditions, with no deviation or sensitivity analysis.
  - Mark n for the "port form, ρ elig." cell (2 conditions).
- **New evidence needed:** analysis plus a cheap re-solve.

### M6. Cross-model holdout: the metric hides the EM mode, test D uses an excluded policy, and "the method transfers" is marginal

- **Quoted text.**
  - "No stable-context reversal occurs: in every stable core context, every replacement leaves α_⊥ unchanged or increases it."
  - "The method transfers."
  - Table VI: D "corridor top-1 agreement 0.75, transfers".
- **What is wrong.**
  1. **The global abscissa is pinned.** In ALT, 62 of the 183 stable-context core marginals have Δα_⊥ exactly 0. In those contexts α_⊥ is set by a real pole at −0.100 s⁻¹ that no replacement moves (`alt/cases/*.json`). The global-α metric therefore cannot show a stabilizing marginal of the EM mode whenever that mode lies left of −0.1.
  2. **Stabilizing EM marginals do exist in ALT.** I re-evaluated on the MAC-tracked EM-band rightmost mode, the same device the paper uses for topology (H14). Stabilizing tracked-EM marginals occur at H02 (5), H07 (7) and HARDENING_H04 (2). HARDENING_H04 has 4 nested tracked-EM reversals, e.g. unit 33: ∅ (−) vs {30,35} (+). So the A/B verdict "model-specific (0)" depends on the metric. Under the paper's own EM-tracking philosophy it would be "partial". My implementation reuses the paper's thresholds (MAC 0.8, 0.15 Hz, τ_mat) but is not the preregistered one. The authors should redo it properly.
  3. **Test D's verdict depends on an excluded policy.** The prereg defines "tested policies" as those whose ALT base is STABLE. `H18_compare.py` applies this to A and B only. C, C2, D and E use all 8 policies, including HARDENING_H02, whose ALT base is UNSTABLE. On the 7 tested policies D = 5/7 = 0.71, which is PARTIAL, not TRANSFERS. C2 becomes 0.203, still barely above 0.20.
  4. **C2 is marginal.** It is descriptive in the prereg. Its median is 0.206 against 0.20, 3 of the 8 per-policy differences are below 0.20, and it uses 12 lines.
  5. **The holdout lacks the paper's key mechanism.** The frozen chain regulates constant Q (QFLAG = 0, plant Q-PI disabled), so it lacks the voltage control (g) that the paper's policy axis and interpretation revolve around.
  6. **The qualification rules were refined after one matrix datum was seen.** Refinements 1, 1b and 2 were defined after the H02/H4 matrix case had been seen, as the deviation log admits.
- **Fix.**
  - Recompute A and B on the tracked EM-band mode, labelled as a deviation.
  - Apply the "tested policies" rule uniformly.
  - Soften "the method transfers" to "a smooth α_⊥(γ) makes the first-order ranking accurate in ALT too (12 lines, descriptive)".
  - Add a voltage-controlling configuration of the same library chain (REECB1/REPCA1 voltage control enabled, gain matched to g) as a second alternative. Otherwise the reading "the stabilizing branch requires converter voltage control" remains untested.
- **New evidence needed:** yes (re-analysis plus one new configuration).

### M7. Audit trail: deviation-log timestamps contradict git and file times; post-hoc analyses are unlabelled; a preregistered phase is incomplete and silently omitted

- **Quoted text.** "All thresholds, holdout points and tests were frozen before the computations reported here; deviations are logged."
- **What is wrong.**
  1. **Entry times postdate the commits that contain them.** The entries timestamped "21:40" and "21:50" are already present in commit `4de412d2`, dated 21:32:13 −0500. The entries timestamped "22:05", "22:10" and "22:20" are in commit `eb57d907` (22:06:09). The log file itself was last modified at 22:00:41.
  2. **One "before" claim cannot be checked.** The 22:05 entry (refinement 1b) states it was "recorded BEFORE any H18 matrix case was run". Yet all 488 ALT case files carry file times 21:36:01-21:46:09 and already contain the `isolated_eigs` field that refinement 1b introduces. The implementation may well have preceded the runs, but the log timestamps cannot establish any "before" claim.
  3. **The eligibility rule was not in the committed analysis code.** The H07 code committed before the results computed the GOLD-B gate on all 64 conditions without the R0 filter (`GB_hardened_pass = GB_primary_pass`). It was changed to the literal eligibility rule in the results commit. The change is in the conservative direction, but its ordering relative to reading the 64-condition result is undocumented.
  4. **The 92 % anatomy result is post hoc and unlabelled.** `regret_anatomy_new` (the 2089/2265 = 92 % figure in the abstract logic and Sec. VI-C) was added to `H04_ranking.py` after the results commit. It is uncommitted in the working tree and absent from the deviation log. It is presented without an "exploratory" label.
  5. **H10 is incomplete and not disclosed.** `CDWH_RUN_STATUS.json` shows H10 (the corridor size-matched null) as RUNNING. `sec_corridors.tex` is a placeholder. `\corRobust{}` evaluates to "none" because the H10 dictionary is empty, not because the gate failed. The paper still uses "the top equal-budget corridor" in Sec. VI-G without defining corridors. The prereg requires negative results to be reported with the same prominence as positive ones.
- **Fix.**
  - Anchor every deviation to a commit hash, and record wall-clock times from the system clock.
  - Label the regret anatomy (and any other post-results analysis) as exploratory, and log it.
  - Complete H10 and report the corridor result, or state in the paper that the preregistered corridor phase is incomplete and why.
- **New evidence needed:** documentation, plus completion of H10.

### M8. Planner relevance and computational claims are not demonstrated

- **Quoted text.** "its cost is one factorization of the equilibrium Jacobian per portfolio and one back-substitution per candidate branch, instead of a re-solved equilibrium and eigenanalysis per candidate."
- **What is wrong.**
  - The implementation (`_sens.py`, `Engine.derivatives`) builds four finite-difference A matrices per candidate. Each is a full numerical DAE Jacobian plus a solve with g_z.
  - The equilibrium Jacobian R_w is itself built by finite differences.
  - No timings are given.
  - For unit sequencing, the recommended practice (evaluate in context) is exhaustive in-context eigenanalysis.
  - Only IEEE 39 is used, and there is no multi-step planning experiment.
- **Fix.**
  - Report wall-clock for finite re-solve and eigen per branch against D_tot per branch, on this system and on one larger system (e.g., a 2000-bus synthetic or the WECC 179-bus).
  - Use analytic or AD Hessian-vector products if the cost claim is kept.
  - Demonstrate one planning use, e.g. selecting reinforcements for a multi-step retirement path with and without first-order screening. Report decision quality against cost.
  - Add a paragraph answering Q5 explicitly.
- **New evidence needed:** yes.

### M9. Statistical reporting details that change the reading of the gates

- **Quoted text.** "(coverage 0.89; bootstrap interval 0.74-1.00) ... the gate (≥0.75) passes."
- **What is wrong.**
  1. **The gate is not robust to resampling.** The H3 gate passes on the point estimate, but the design-bootstrap interval crosses the gate. The exact binomial interval for 17/19 is [0.67, 0.99].
  2. **Pair counts mix populations and overlap.** "3263 new-holdout level-C pairs", the magnitudes and the direction counts (2172/1091) include 196 pairs from the base-unstable policies HARDENING_H07, H08 and H15, which the prereg excludes from reversal gates. Pairs overlap heavily (604 distinct S₁-marginals), so counts are not evidence weights.
  3. **p-values need exchangeability statements.** T1 (19/19 policies with p* < 0.90, p = 1.9e-6) and T3 are reported without the exchangeability statements the plan requires.
  4. **Coverage depends on the domain.** It depends on the domain and its measure (g = u²). A different but equally defensible parameterization (uniform g, or excluding AVR tunings that destabilize the all-SG base) could change coverage. No such sensitivity is shown.
- **Fix.**
  - State that the H3 gate passes on the point estimate only, and give the exact interval.
  - Restrict descriptive pair statistics to base-stable policies, and report distinct marginals and units.
  - Report T1-T3 with nulls.
  - Add a domain-measure sensitivity (coverage under uniform g; coverage excluding base-unstable-prone AVR regions).
- **New evidence needed:** analysis only.

---

## 5. MINOR comments

1. **Sec. VI-H.** "(50 of the V4 cases)": 50 is a count. Write "50 of 1296 V4 (policy, action) cases" (`H15_topology_gate.json`, `label_counts_H4`).
2. **Abstract.** It mixes condition sets: "0.985 on preregistration-eligible held-out conditions ... against 0.491 for [|P_e|] and 0.499 for the best of nine static indices chosen per condition". The 0.499 is from the 64-condition set. On the 8-condition set the oracle static is 0.491, identical to |P_e|. Use one set per sentence.
3. **Sec. VI-D.** "the threshold ... is below the solver's typical residual of 10⁻⁷". The H4 median R0 is about 5×10⁻⁸ and the maximum is 2.1×10⁻⁷. Say so.
4. **Sec. V-C.** "Reruns of selected tasks were byte-identical." The determinism check compares summary fields (16 H01 + 4 H06 tasks, `differing_keys: []`). Say "identical on every summary field".
5. **Limitations.** "Reversal magnitudes are a few hundredths of s⁻¹." The quartiles are 0.011-0.017 s⁻¹, i.e. about one hundredth, which is roughly 0.3 % damping ratio at 0.7 Hz.
6. **Fig. 4 (right).** It plots a quantity that satisfies the bound by Lemma 1 (the mean of m terms equals D/m), so every point lies above the line by construction. Replace it with the EM-clean witness share by direction (s→d 88 %, d→s 97.5 %) and the witness magnitude distribution.
7. **Sec. VI-I.** "attributes 19 % of the variance of the mixing change to the controller term and the rest to frequency movement". With cov(C, F) = −4.7×10⁻⁶ against a total variance of 3.1×10⁻⁵, the covariance term is −30 % and var(F)/var(C+F) ≈ 1.11. Report the three terms; the split is not additive.
8. **Sec. VI-H, "do not order actions within that distinction."** My re-computation from `H14_topology_em.parquet` supports this. The within-doublings median ρ is 0.00 (Fiedler), 0.03 (Kirchhoff) and −0.16 (gSCR). Within outages it is −0.19, 0.21 and −0.14. Losses reach |ρ| ≥ 0.6 in 19 % of policies within outages. Report these numbers.
9. **C4 provenance.** EM removal occurs only at D01 (discovery) and H06 (old holdout). Of the 45 EM creations, 42 are at discovery or old-holdout policies (D01: 20, D03: 14, H01: 4, H02: 4). Only 3 are on the new holdout. All creations are outages and all removals are doublings. Disclose this and temper C4's weight as a contribution.
10. **Sec. VI-A.** "Fig. 1 shows the strongest level-C nested reversal on the new holdout." The strongest (0.068 s⁻¹) is at base-unstable N08. Say "among base-stable new-holdout policies".
11. **Sec. VI-D.** "(the P4 incompatibility)" is undefined in this paper. Define it (V₄ unstable at P4 with a 0.62 Hz critical mode).
12. **Sec. V-D.** "validated against an EMT implementation in earlier work" applies to the TX3 composition. The ALT-WECC composition (library chain plus SG2AX transcription) is validated here only at the all-SG base (Q1) and by initialization residuals. Rephrase.
13. **Table I.** It omits the machine data needed to reproduce the benchmark: H, reactances, time constants, ratings, the bus-39 equivalent and the GFL ratings/base. "Harmonized benchmark" is undefined. Add an appendix or a data DOI.
14. **Data and code availability.** It names a local branch. Provide a public URL or archive (DOI), the environment lock files, and the TX4 `ibr_cycles` dependency version.
15. **"24 distinct optimal orders appear among the 52 base-stable policies."** With near-ties this is uninformative. Report the median Kendall distance between optimal orders.
16. **Figures.**
    - Fig. 1 text (panel titles, footer) is illegible at column width.
    - The Fig. 5 transfer heatmap has no policy tick labels or split separators.
    - In Fig. 6 (right), color points by whether the finite step caused a rightmost-mode switch (gap2). The off-diagonal cluster is a saturation/switch effect.
17. **Sec. VI-D.** "the ranking is unaffected" by saturation is too strong. Quantify top-5 changes between the first-order and the finite ranking for the saturating branches.
18. **Sec. II-C.** Specify that MAC uses the 78 rectangular bus-voltage components normalized to unit norm with phase fixed at the largest entry. Also specify that candidates are band modes with Re ≥ −1 plus the critical mode.
19. **Sec. V-E / Sec. VI.** Report T1 (sign test) with its null. The paper only states "19/19 policies have p* < 0.90". Also list the Holm family members and their adjusted p-values in one place.
20. **Table V.** "port form (policies only)" in the "ρ elig." column is computed on 2 conditions. Mark n.
21. **Discussion.** It could also cite work on the grid-synchronization stability of GFLs in weak grids and on the stability limits of low-inertia/high-penetration systems. This would help situate the penetration cliff reported here.
22. **"Held-out".** D_tot has no fitted parameters, so "held-out" has meaning only for the static-index selection (S1 chosen on discovery). Say so, to avoid implying out-of-sample validation of the derivative.
23. **Sec. VI-G.** "V₄ is unstable at all eight policies, with 0.3-0.8 Hz critical modes" is correct (0.31-0.81 Hz). Add that in ALT the V₄ target is unstable, so the ALT ranking ranks stabilizing reinforcements of an unstable portfolio.
24. **Sec. VI-F.** The fresh draws are all at P4, a discovery policy. Their "held-out" status is only with respect to parameter values. Make this explicit next to "held-out policies and parameter draws" in C3.

---

## 6. Numerical claim verification

Every macro-driven number I checked matches its source file unless flagged.

| claim (location) | file / field | file value | status |
|---|---|---|---|
| 17/19 level-C, 16/18 old, 14/15 disc. (Abstract, VI-A, Tab. II) | H03_gate.json `new_frac_C` | 0.895 (17/19) | match |
| coverage interval 0.74-1.00 | H05_stats.json `new_cov_C` | [0.737, 1.0] | match (crosses gate, M9) |
| level D 17/19 | H03 `new_frac_D` | 0.895 | match |
| both directions 0.37 (7/19) | H03 `new_frac_C_both_dirs` | 0.368 | match |
| 3263 pairs; d→s 2172, s→d 1091; quartiles 0.011/0.013/0.017 | H03 `new_C_*` | same | match, but includes 196 pairs from base-unstable policies (M9) |
| EM-clean witness 94 % | H03 `new_C_frac_em_clean_witness` | 0.944 | match |
| Fig. 1 example Δ −0.042 / +0.064 "same mode" | H03_nested_pairs (N24, 33, s2d) | d1 −0.0420, d2 +0.0644; hz 0.744 vs 1.310; lvD False | **numbers match, "same mode" wrong (M3)** |
| Abstract "left ... right in a superset" in 17/19 | H03_policy_summary `C_n_s2d>0` | 9/19 | **mismatch (M3)** |
| p* 0.75 (0.74-0.78), regret 0.38 | H04 / H05 | 0.749 [0.740, 0.783]; 0.376 | match |
| EM p* 0.97 / regret 0.03; SAME 0.98 / 0.02 | H04 | 0.973 / 0.034; 0.978 / 0.023 | match |
| shift −0.002 (−0.036 to 0.041), 0.007 | H05 `shift_*` | same | match |
| 19/19 p* < 0.90 | H05 `T1_sign_test` | 19/19, p = 1.9e-6 | match (p not reported) |
| transfer 0.40 off-diag; 24 orders; 52 policies | H04 | 0.402; 24; 52 | match (interpretation, M2) |
| 2089 of 2265 (92 %) | H04 `regret_anatomy_new` | same | match; **uncommitted post-hoc (M7)**; all 2089 are 0 Hz real eigenvalues (M1) |
| GOLD-B elig. 8 conditions: 0.985 / 0.913 / 0.491; diff 0.492 (0.477-0.608); top-5 1.000 | H07 `GB_prereg_*` | same | match |
| GOLD-B 64: 0.987 / 0.895 / 0.497; diff 0.489 (0.484-0.638) | H07 `GB_primary_*` | same | match |
| R0 max 2.1e-7 | H07 `R0_max_all` | 2.146e-7 | match |
| strata: EM 0.988 vs 0.498; SAME 0.988; γ 1.10: 0.999; 1.25: 0.995; V9 −0.183 vs 0.994; Dconv 0.034 | cdw_numbers / H07 | same | match |
| new policies D_fro 0.712, D_tot 0.996 | H07 `GB_new_policies_*` | same | match |
| Table V rows | recomputed from H06_metrics.csv | all cells reproduce to 2 d.p. | match |
| IV 100 % of 4048 | H07 `IV_*` | 1.0, 4048 | match |
| controller agreement 2.2e-7; loads 61 %; vset frozen 0 | H07 `node_*` | same | match |
| 28 % material; ~half ratio > 1; every condition has a flip | H07 `H8_*` | 0.281; 0.457; 1.0 | match |
| L7 1.00/0.90/0.80/1.00; EM 0.00/0.50/0.60/0.70 | H05 `L7_*` | same | match |
| fresh draws 0.986 vs 0.503 | cdw_numbers | same | match |
| ALT: 488/488 init ≤ 3e-13, no active limit, base within 2e-7 | H18 / deviation log | 2.49e-13; 0 limit; 1.7e-7 | match |
| ALT: A 0/7 (custom 0.71), Kendall 0.02, C2 0.99 vs 0.79, D 0.75, E 0.94 (35) | H18_gate.json | same | match; **D on tested set = 0.71 (M6)**; "unchanged or increases" is pinned-pole artifact (M6) |
| topology: 45 creations, 0 fast-only; P4 removals dbl0/1/13/43 = lines 1-2, 1-39, 8-9, 23-36 | H15 | same (branch map verified) | match |
| "50 of the V4 cases" | H15 `label_counts_H4` | 50 of 1296 | count, needs denominator |
| static topology ρ −0.50, 0.56, −0.54, −0.54; max frac 0.38 | H15 `H16` | same | match |
| mixing ρ 0.66 → 0.21; new 0.55, p 0.016, Holm 0.033; 19 % | H13 / H19 | same | match (variance split, minor 7) |
| 7 of 8 layers | H19 | 7 | match |
| T2 (GOLD-B) | H19 `F_conf_*` | p = 0.061, Holm 0.061 | **omitted (M5)** |
| corridor results | CDWH_RUN_STATUS H10 | RUNNING | **incomplete, omitted (M7)** |

---

## 7. Reviewer reanalyses (for the authors to replicate)

All of the following use only the result files listed and read-only access. No project file was modified.

1. **Out-of-band anatomy** (`H03_marginals.parquet`, new holdout, `lvA`). Take marginals with `hz_Si` outside [0.1, 2.0]. Tabulate `hz_Si`, `alpha_Si` and `status_Si`. Join with `gz_cond` from `raw/H_H01/*.json` records, and tabulate by portfolio size.
2. **Regret-optimal fixed order.** For each base-stable new policy, enumerate all 9! orders. For each context (mask) with ≥ 2 level-A units, the choice is the first available unit in the order, and material means Δ(choice) − min Δ ≥ 0.01. Minimize the material rate.
3. **Stability screen.** Take the DP order from `H04_ranking.csv` (FULL). The screened choice is the first unit with `status_Si == STABLE`, falling back to the DP choice. Report on contexts with at least one stable option.
4. **Leave-one-cluster-out fixed branch list** (`H06_links_raw.parquet`, family `link`, target H4). The truth is −(α(γ=1.5) − α₀). The fixed list is the mean truth rank over conditions of other clusters. Compare Spearman and top-5 with −d_total.
5. **ALT EM-tracked reversals** (`alt/cases/*.json`, net NOMINAL, core lattice, stable S). Take the rightmost mode with 0.1 ≤ f ≤ 2 Hz. Match it into S∪{i} by the best MAC among band modes within 0.15 Hz, and require MAC ≥ 0.8 and that the match is the EM-band rightmost of S∪{i}. Sign classes use τ_mat, and nestedness is checked on subsets.
6. **gap2 at level-C endpoints.** Take `gap2` of S₁, S₁∪{i}, S₂ and S₂∪{i} from the raw census, and compare the minimum with the pair magnitude.
7. **Audit trail.** Run `git show 4de412d2:<deviations>`, `git show eb57d907 --stat`, `git diff eb57d907` and the file modification times of `results/hardening/*` and `alt/cases/*`.
