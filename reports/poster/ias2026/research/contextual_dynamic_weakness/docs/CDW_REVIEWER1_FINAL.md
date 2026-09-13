# Reviewer 1 report — IEEE Transactions on Power Systems

**Manuscript:** "Portfolio-Conditioned Dynamic Weakness and Reinforcement Ranking in
Inverter-Rich Power Networks" (`reports/papers/cdw_contextual_dynamic_weakness/main.tex`,
compiled `main.pdf`, 10 pp.)

**Reviewer focus:** theory, novelty, logical validity of claims.

**Date:** 2026-09-12

**Recommendation:** Major revision. The paper is not acceptable in its current form.
Some of the required changes are rewording, but three need new evidence:
- a second network or a second converter model with voltage control (M5);
- a physical characterization of the fast aperiodic instabilities behind C2 (M4);
- tests of the reinforcement sensitivity on large, realistic actions (M8).

---

## 1. Summary of the submission

The paper studies α⊥(S), the transverse spectral abscissa of the IEEE 39-bus system when the
synchronous generators in S ⊆ V9 = {30,…,38} are replaced by a custom grid-following
converter at matched dispatch. It makes four claims.

- **C1.** Nested sign reversals of the marginal Δ_iα⊥(S) exist on a tracked
  electromechanical (EM) mode, at 17/19 base-stable new-holdout policies. The chain
  identity turns them into "certificates", so α⊥ is neither submodular nor supermodular.
- **C2.** An oracle fixed ranking of units is materially wrong in 0.38 of next-replacement
  decisions, mostly through transitions to fast converter-control modes. In the EM band it
  is right for 0.97 of pairs.
- **C3.** The classical re-equilibrated ("total") eigenvalue sensitivity to branch
  admittance ranks finite ×1.5 reinforcements with median Spearman ≈ 0.99. Static indices
  reach ≈ 0.5 and the frozen derivative ≈ 0.90.
- **C4.** Single branch actions remove and create EM incompatibilities, and seven static
  topology scores fail a preregistered prediction rule.

A cross-model check with a WECC library GFL chain reproduces the ranking method but not the
reversals.

Several aspects deserve credit:
- a frozen preregistration and a deviation log;
- stratification of every headline by modal band;
- negative results reported openly (mixing, cross-model);
- every number generated from the result files as a macro;
- deterministic reruns.

The lemma-level algebra is correct. My concerns are about novelty, the interpretation of
the algebra, a mislabelled flagship example, the physical validity of the modes that drive
C2, generalization, and several claims that the authors' own result files contradict.

---

## 2. Answers to the four key questions

**"Isn't this just eigenvalue sensitivity?"**
- For C3, yes. The paper says so itself (Eq. (4), l. 333–335). What remains is an evaluation
  of a known predictor on one benchmark.
- The predictor has no fitted parameters, so the "held-out" framing adds robustness across
  conditions but not out-of-sample validation.
- The comparison that matters in practice is total vs frozen: 0.99 vs 0.90, or 0.996 vs
  0.71 on the new policies. The static indices are not designed to rank modal damping and
  work as strawmen.
- The "port-form" row is not independent. `E34_sens.port_items` builds it from Jacobians of
  the re-solved SPR case at 1 ± 1e-3, so it is the same total derivative in a different
  algebraic form.

**"Isn't sign reversal just mode switching?"**
- The level B/C/D grading is the right tool, and in most level-C pairs the reversal stays
  on one mode.
- But the paper's own flagship example is a mode switch between the two contexts (M1). In
  the N24/unit-33 pair, the critical mode is 0.744 Hz at S1 and 1.310 Hz at S2. The pair is
  level C, not level D.
- 524 of the 3067 base-stable level-C pairs (17 %) are not level D. These include 334 of
  the 1077 s→d pairs, the direction that matters for non-submodularity.
- At level D, both directions occur in only 3/19 policies.
- No mechanism is offered for the genuine same-mode reversals, and the mixing mechanism
  failed its test. The paper therefore does not fully answer the question.

**"What is mathematically new?"**
- Nothing, and the paper says so (l. 243–245). Lemma 1 and Propositions 1–3 are one-line
  consequences of telescoping. Remark 1 is a textbook observation.
- The "certificate" is automatic. For any nested pair, Proposition 2 guarantees a
  certifying term, so its existence is not an empirical finding (M2).
- The conclusion "neither submodular nor supermodular" follows more directly from one
  one-step square of each sign. It has no algorithmic consequence here: α⊥ is not
  monotone, and greedy guarantees would not apply even under submodularity.

**"Why should IEEE-39 generalize?"**
- The paper gives no reason. There is one network, one custom GFL and one parameterization.
- The only out-of-model test (ALT-WECC) does not reproduce C1.
- The converter-parameter envelope EC removes all EM reversals at P4 (0/10 draws).
- The evidence therefore points to C1 being specific to the custom GFL's voltage control
  and parameters (M5).

---

## 3. Scores (1–10)

| criterion | score | justification |
|---|---|---|
| Novelty | 3 | No new mathematics (acknowledged). C3 is Smed 1993 / Nam 2000 applied to portfolio-specific operating points. C4 is classical N-1 small-signal behaviour. The gap is carved narrowly by a 2020–2026 window that excludes the most relevant prior work on location-dependent converter displacement (Gautam et al. 2009, which is in the authors' own gap matrix; Quintero et al. 2014). The genuinely new element is an empirical observation (same-mode nested reversals) without a mechanism, which fails the only cross-model test. |
| Mathematical correctness | 5 | Lemma 1, Propositions 1–3 and Remark 1 are correct. Checked: the counterexample's six second differences are 0, −3, −1, 0, −3, −1, and Eq. (4) is the standard IFT eigenvalue derivative. Four statements are wrong: (i) "a nested one is not [compatible with submodularity]" is false, and the paper's own counterexample contains nested d→s reversals; (ii) Remark 2 wrongly asserts ẇ = 0; (iii) the SPR description omits the slack and the gauge; (iv) the abstract's unqualified "neither submodular nor supermodular" does not follow from single-direction reversals. |
| Model adequacy | 3 | Phasor model with D = 0, no governors, constant-power loads and matched dispatch. The idealized GFL has no limits and no delays. The fast modes that drive C2 are real eigenvalues with median +447 s⁻¹ and up to 2.7×10⁴ s⁻¹ (my check), far outside phasor validity. There is one network and one converter model. |
| Statistics | 5 | Preregistration, cluster bootstrap and Holm correction are all good. Against that: T1 and T2 of the preregistered confirmatory family are not reported, and T2 (p = 0.061) fails. The primary GOLD-B set has 8 conditions in 6 clusters. Coverage intervals are percentile bootstraps on 17/19. There is no threshold-sensitivity analysis. Evidence layers L1–L4 are logically nested but counted as separate. |
| Evidence | 5 | Strong for "total beats frozen and static at small steps". Weak for C1's practical significance: median 0.013 s⁻¹, about 0.3 % damping ratio. C2's EM-band relevance is also weak, since the prereg A3 gate fails in the EM stratum. The claimed non-transfer of rankings is contradicted by the numbers. |
| Reproducibility | 7 | Frozen seeds, macros generated from JSON and determinism checks are excellent. But the repository is identified only by a branch name, Table I omits data needed to rebuild the model, and one macro (`\corRobust`) is a default value for a missing file. |
| Industrial relevance | 4 | "Recompute reinforcement rankings per portfolio" is useful, but it is what planners already do with modal sensitivity at the study case. The materiality threshold is an order of magnitude below planning damping criteria (3–5 %), and the unit-sequencing result is driven by aperiodic instabilities of an idealized converter. |
| Writing | 6 | Terse and carefully hedged in places. But internal jargon leaks throughout (H4, S1, TX3, A3 bar, "P4 incompatibility", HARDENING_H01, N24), notation is inconsistent (V4/H4), Fig. 3 is never referenced, and corridors appear without a definition. |
| Figures | 5 | Fig. 1 labels are illegible at print size, and the figure shows a pair that is not same-mode. Fig. 3 has the same problem. Fig. 4 (right) plots a theorem. Fig. 5 (right) has unlabeled policy axes and visually contradicts the "does not transfer" text. |

---

## 4. Verification log (numbers I recomputed from the result files)

All checks used the frozen outputs in `results/hardening/`:
- `H03_nested_pairs.parquet`, `H03_marginals.parquet`, `H03_policy_summary.csv`;
- `H04_ranking.csv`, `H07_goldb_gate.json`, `H18_crossmodel.csv`;
- `H13_mixing.json`, `H15_topology_gate.json`, `H19_evidence.json`.

The following numbers were checked.

1. **Flagship pair (Fig. 1; §VI-A, l. 414–418).** HARDENING_H24, i = 33,
   S1 = {31,34,38}, S2 = {31,34,37,38}.
   - d1 = −0.0420 and d2 = +0.0644: these match the paper.
   - hz1 = 0.744 Hz and hz2 = 1.310 Hz; `lvD = False`. **This is not a same-mode pair.**
   - Fig. 3 panel 3 (HARDENING_H09, i = 38, ∅ → {30,33,35}) is also `lvD = False`
     (0.643 → 0.453 Hz).
   - Panel 2 (H24, i = 35) is level D.
2. **Level-C pair population.**
   - The paper's 3263 pairs (2172 d→s, 1091 s→d, 94 % EM-clean) include the base-unstable
     policies H07, H08 and H15 (196 pairs).
   - Base-stable only: 3067 pairs (1990/1077), 95.7 % EM-clean.
   - These pairs come from only 87 distinct (policy, unit) combinations and 1513 distinct
     S2 marginals.
   - Chain length: m = 1 for 262 pairs. For these the "certificate" is identical to the
     reversal itself.
3. **Threshold sensitivity (base-stable new holdout).**

   | τ (s⁻¹) | level-C coverage | both directions (C) | level-D coverage | both directions (D) |
   |---|---|---|---|---|
   | 0.01 | 17/19 | 7/19 | 17/19 | 3/19 |
   | 0.015 | 17/19 | 4/19 | 14/19 | 2/19 |
   | 0.02 | 15/19 | 3/19 | 11/19 | 1/19 |
   | 0.03 | 10/19 | 1/19 | 5/19 | 0/19 |
   | 0.05 | 0/19 | 0/19 | 0/19 | 0/19 |

   - The largest level-C magnitude on base-stable new policies is 0.042 s⁻¹.
   - In damping-ratio terms (|Δσ|/ω), the quartiles are 0.27 %, 0.32 % and 0.39 %.
4. **Fast modes behind C2.**
   - On the base-stable new holdout, 9040 of 32663 stable-context replacements (27.7 %)
     give α⊥(S∪i) ≥ 1 s⁻¹.
   - 100 % of these have a purely real critical eigenvalue (hz = 0).
   - Their median α⊥ is 447 s⁻¹, the 95th percentile 5511 s⁻¹, and the overall maximum
     2.68×10⁴ s⁻¹.
   - They occur at post-replacement sizes 5–6.
   - `H04_gate.regret_anatomy_new`: 2032 of the 2265 material regrets involve an
     out-of-band choice with Δ ≥ 1 s⁻¹.
5. **Regret-optimal fixed ranking.** The paper's ranking maximizes pairwise agreement
   (`H04_ranking.py`, `optimal_linear_order`), not regret. An exact subset DP over the 2⁹
   placed prefixes gives the regret-minimizing fixed ranking:
   - FULL: median regret 0.345, vs 0.376 for the p*-optimal ranking;
   - EM: 0.030, vs 0.034.

   The qualitative conclusion stands, but "the best fixed ranking … 0.38" is not correct.
6. **Transfer (`H04_gate`).**
   - Median off-diagonal regret is 0.402 vs 0.376 on the diagonal; accuracy is 0.725 vs 0.763.
   - `H19_evidence.em_qualifiers.L5_EM_would_pass_A3_rule = false`: the L5 gate fails in
     the EM stratum.
7. **Cross-model (`H18_crossmodel.csv`).**
   - On the 7 ALT-tested policies, the *custom* GFL has a level-B/C core-lattice reversal
     at only 1/7 (HARDENING_H04). Test B therefore cannot discriminate between the models.
   - Test A is discriminative: custom 5/7 vs ALT 0/7.
   - Eight policies were preregistered, but 7 were tested: HARDENING_H02 has an unstable
     ALT base.
8. **Preregistered tests not reported.**
   - `H19_evidence.F_conf_raw_p`: T1 = 1.9e-6, T2 = 0.061, T3 = 0.016. Holm-adjusted:
     5.7e-6, 0.061, 0.033.
   - The paper reports only T3.
   - `H10_gate.json` does not exist, so `\corRobust{none}` is the default for a missing
     file (`CDWH_NUMBERS.py`, line 136).
9. **Topology creations (`H15_topology_gate.json`).**
   - All 45 EM creations are outages.
   - 34 of them are at the discovery policies D01 and D03; only 3 are at new-holdout
     policies.
   - EM removal occurs at P4 (4 doublings) and once at old-holdout H06.
10. **Items that check out.**
    - 17/19 level C; 7/19 both directions.
    - p* = 0.749 and regret = 0.376.
    - GOLD-B medians (0.985 / 0.913 / 0.491 on the 8 eligible conditions; 0.987 / 0.895 /
      0.497 on all 64).
    - IV 4048/4048; ALT init residual 2.5e-13; 32256 = 63 × 512.
    - Holm arithmetic; mixing ρ = 0.548 and p = 0.0163.

---

## 5. MAJOR comments

### M1. The flagship "same-mode" reversal is a between-context mode switch
- **Where:**
  - Abstract, l. 52–55: "the same replacement moves the same tracked electromechanical mode
    left in one stable portfolio and right in a portfolio that contains it";
  - §VI-A, l. 414–418: "in S2 … the same replacement moves the same mode right by
    0.064 s⁻¹";
  - Fig. 1 caption;
  - Fig. 3 (y-axis "on the tracked EM mode").
- **Problem:**
  - Level C only requires each endpoint marginal to track *its own* critical mode (l. 223–229).
  - In the flagship pair the critical mode is 0.744 Hz at S1 and 1.310 Hz at S2
    (`lvD = False`). Replacing unit 37 (S2∖S1) switches the critical mode, so unit 33 then
    acts on a different mode.
  - Fig. 3 panel 3 has the same defect.
  - The abstract and C1 use level-D language ("the same tracked … mode") for level-C counts.
  - l. 422–423 ("Every policy with a level-C reversal also has one at level D") is a
    policy-level statement. It does not make the displayed pairs level D.
  - 17 % of base-stable level-C pairs, and 31 % of the s→d pairs, are not level D.
- **Fix:**
  - Make level D the headline for any "same mode" claim.
  - Replace Fig. 1 and panels 1 and 3 of Fig. 3 with level-D pairs. For example, H18 unit 38:
    {35} → {30,33,35}, −0.031 → +0.045 s⁻¹, 0.66/0.68 Hz, MAC 0.91. H24 unit 35 is already
    level D.
  - Report the level-D direction split (1800 d→s, 743 s→d) and the per-policy both-direction
    count (3/19).
  - State in the abstract that both directions on one tracked EM mode occur in 3/19 policies.
- **Needs:** rewording and regenerated figures from existing data. No new simulations.

### M2. Misstatement of what nested reversals prove, and the "certificate" framing is vacuous
- **Where:**
  - l. 284–288: "A reversal between incomparable contexts is therefore compatible with
    submodularity; a nested one is not (Proposition 1)";
  - abstract l. 55–57; C1 l. 133–134;
  - §VI-B l. 462–465; Fig. 4 (right).
- **Problems:**
  1. The sentence is false for d→s reversals. The paper's own counterexample contains two
     nested d→s reversals in a submodular f:
     - ∅ ⊂ {a}: Δ_i f(∅) = +2, Δ_i f({a}) = −1;
     - {b} ⊂ {a,b}: Δ_i f({b}) = +1, Δ_i f({a,b}) = −2.

     Proposition 1(a) only excludes s→d reversals. Two thirds of the observed pairs are d→s
     (2172/3263), and these refute only supermodularity.
  2. The abstract derives "neither submodular nor supermodular" from "each such nested
     reversal". A single-direction reversal proves only one of the two failures. Both
     directions at a policy are needed. Per the results they occur in 17/19 policies at
     level A, 7/19 at level C and 3/19 at level D. The claim must be qualified per level and
     per policy.
  3. The "certificate" is automatic. By Lemma 1 and Proposition 2, *every* nested pair with
     D ≠ 0 has a term of the required sign and size, whatever f is. That "every one of the
     3263 … has a certifying curvature term" (l. 462–464) can only fail through a software
     bug. Fig. 4 (right) plots the inequality of Proposition 2, and for the 262 pairs with
     m = 1 the witness *is* the reversal. The only informative quantity is the EM-clean
     fraction.
  4. "Neither sub- nor supermodular" is shown more directly by one one-step square with
     d > 0 and one with d < 0 (the classical second-difference characterization). Such
     squares abound (Fig. 4 left). Non-submodularity of spectral and control-energy set
     functions is known (Olshevsky 2018; Summers & Kamgarpour 2019). Because α⊥ is not
     monotone, no greedy guarantee is lost either. The set-function section therefore has no
     algorithmic consequence in this paper.
- **Fix:**
  - Rewrite l. 286–288 as: "a nested s→d reversal is incompatible with submodularity, a
    nested d→s reversal with supermodularity".
  - Qualify the abstract and C1 claim with per-level, per-policy counts.
  - Remove "certify/certificate" as a contribution and describe Proposition 2 as an
    identity used to *locate* an interaction term.
  - Move Fig. 4 (right) to the supplement, or replace it with the EM-clean witness
    distribution.
  - If the set-function framing stays, give it content. For example, compute the
    submodularity ratio or weak-supermodularity constants of −α⊥ restricted to stable
    chains (Das & Kempe 2011; Bian et al. 2017; Summers & Kamgarpour 2019) and connect them
    to a sequencing guarantee.
- **Needs:** rewording. The submodularity-ratio analysis is optional new analysis on
  existing data.

### M3. Novelty and the gap statement
- **Where:** l. 106–127 (credited ingredients and gap), contributions C1–C4, Discussion
  l. 677–684, l. 704–707.
- **Problems:**
  1. The effect of displacing SGs by converters on specific oscillatory modes is known to
     depend on which units are displaced, and can be beneficial or detrimental:
     - Gautam, Vittal & Harbour, TPWRS 24(3), 2009. This paper is in the authors'
       `CDW_LITERATURE_GAP_MATRIX.csv` but is not cited.
     - Quintero, Vittal, Heydt & Zhang, TPWRS 29(5), 2014.

     Limiting the search to 2020–2026 (l. 121, l. 704) excludes exactly this line of work.
     Desai et al. 2026 (IEEE-39 with replaced SGs, dispatch-dependent stability) is also in
     the gap matrix but not cited.
  2. C3 is the Smed/Nam feasible sensitivity evaluated on one benchmark. The paper admits
     this (l. 333–335, l. 677–684). The "port form" in Table V is the same total derivative
     computed from re-solved Jacobians (`E34_sens.port_items` calls `spr_case` at 1 ± h).
     It is a code cross-check, not an independent method from [11].
  3. C4 amounts to "outages can create instability and reinforcements can remove it". All
     45 creations are outages, and 34/45 are at two discovery policies. This is N-1
     small-signal security, which the paper itself credits to [23], [24].
  4. C3's "the re-equilibration term … is what separates it from the frozen derivative"
     (l. 143–144) is a definition: D_tot − D_fro *is* that term. It is not a finding.
- **Fix:**
  - Cite and discuss the 2009–2016 displacement literature, and extend the search window.
  - Restate the contribution as an empirical case study.
  - Demote C4 to a supporting observation.
  - Delete the tautological clause in C3.
  - To earn a theory contribution, give a mechanism for the level-D reversals. For example,
    use a homotopy in the replaced unit's rating, so that Δ_iα⊥(S) becomes an integral of
    an eigenvalue sensitivity. Then show which change in the mode shape of the tracked mode
    between S1 and S2 flips the sign of the integrand.
- **Needs:** rewording and citations. The mechanism needs new analysis.

### M4. C2 is driven by aperiodic, extremely fast eigenvalues of an idealized phasor converter model, and fails in the EM band
- **Where:** abstract l. 57–61; C2 l. 135–139; §VI-C l. 496–507; Discussion l. 686–690;
  Conclusion l. 727–729.
- **Problems:**
  1. My check (§4, item 4):
     - 27.7 % of stable-context replacements on the base-stable new holdout drive α⊥ above
       1 s⁻¹ through a *purely real* eigenvalue;
     - median +447 s⁻¹ (time constant about 2 ms), 95th percentile 5.5×10³ s⁻¹, maximum
       2.7×10⁴ s⁻¹;
     - they occur when 5–6 of 9 SGs are converted;
     - 2032 of the 2265 material regrets are of this type.

     Real eigenvalues of this size in an index-one phasor DAE point to one of three causes:
     a near-singular algebraic Jacobian g_z (loss of the index-one property); an aperiodic
     PLL/weak-grid loss of synchronism; or an artefact of an algebraic network with
     current-controlled sources and constant-power loads. None of these is validated. The
     paper's own reference [16] (Markovic et al. 2021) shows that network dynamics must be
     modelled for converter-dominated fast modes. Calling them "fast converter-control
     modes" (l. 502) hides that they are non-oscillatory.
  2. The preregistration defines a result as surviving fast-mode removal "if its gate or
     effect direction holds in the EM stratum" (statistical plan §9). It fails:
     - `L5_EM_would_pass_A3_rule = false`;
     - EM regret 0.03, EM p* 0.97.

     The abstract still lists C2 as a contribution, and the conclusion (l. 733–734)
     generalizes from it.
  3. Aperiodic weak-grid GFL instabilities are the classical target of static strength
     indices. The paper never tests a *portfolio-conditioned* static screen for the
     next-replacement decision, for example the gSCR or minimum SCR of S∪{i}. If such a
     screen predicts the fast-mode jumps, the thesis that eigenanalysis is required for
     sequencing collapses.
- **Fix (new evidence):**
  - For these events, report σ_min(g_z), the participation factors of the real mode, and
    the sensitivity to PLL/inner-loop bandwidth and to load model (ZIP).
  - Spot-check them with network dynamics included or with an EMT model.
  - Test the gSCR(S∪{i}) screen on the sequencing decisions.
- **Fix (rewording):**
  - State in the abstract and C2 that the prereg A3 gate fails in the EM band.
  - Describe the fast events as real, aperiodic divergences of 10²–10⁴ s⁻¹.
  - Remove "Replacement planning should therefore rank interventions in the portfolio
    context" unless it survives item 3.

### M5. Generalization and cross-model evidence
- **Where:** l. 150–152; §V-D; §VI-G (l. 598–613); Table VI; Limitations l. 709–713.
- **Problems:**
  1. One network, one custom GFL, one parameter box. The new holdout interpolates within
     the same box in the same simulator, so it tests robustness to controller coordinates,
     not generality.
  2. The only model-level holdout (ALT-WECC) finds no stable-context reversal. In addition:
     - EC draws alone remove all EM reversals at P4 (0/10, l. 591–593);
     - ALT is run in constant-Q mode;
     - V4 is unstable at all 8 policies, and only 76/488 ALT cases are stable.

     The number of stable contexts per policy, which bounds the power of test A, is not
     reported.
  3. Tests B/C in Table VI are uninformative. At the same 7 policies the *custom* GFL has a
     same-mode/EM core-lattice reversal at only 1/7. "The reversals … do not [transfer]"
     (l. 152, l. 732) conflates the discriminative global test A (5/7 vs 0/7) with the
     non-discriminative tests B/C (1/7 vs 0/7).
  4. The sentence "consistent with, but does not establish, the reading that the
     stabilizing branch … requires converter voltage control" (l. 604–606) names a testable
     hypothesis. Enabling voltage control in REEC-B/REPC-A is a direct test.
  5. The H17 structural-zero refinements (Q4 rules 1 and 1b) and the limiter-flag
     exemptions were adopted after one ALT case (H02/H4) had been seen (deviation log,
     21:50 and 22:05). They determine "488 of 488 … no active limit" (l. 599–600) and
     should be disclosed in the paper.
- **Fix (new evidence):**
  - Rerun ALT-WECC with plant/electrical-control voltage regulation enabled.
  - Add a second network, for example the IEEE 68-bus NETS–NYPS system with sampled
    replacement chains instead of the full lattice.
  - Report stable-context counts per policy for both models.
- **Fix (rewording):** reword the transfer verdict so that only test A is called
  discriminative, and disclose the H17 refinements.

### M6. Practical significance and threshold dependence of C1
- **Where:** l. 212–215 (τ_mat), l. 427–432, abstract, C1.
- **Problem:**
  - The reversals are barely above τ_mat: quartiles 0.011/0.013/0.017 s⁻¹, or 0.27–0.39 %
    damping ratio.
  - The maximum on base-stable new policies is 0.042 s⁻¹ (below 1 % damping ratio).
  - Coverage depends steeply on τ (§4, item 3): at τ = 0.02, level C is 15/19 and level D
    11/19; at τ = 0.03, 10/19 and 5/19; at 0.05, zero.
  - Planning criteria (3–5 % minimum damping) are an order of magnitude larger.
  - The abstract and C1 do not convey this. The limitations paragraph says only "a few
    hundredths of s⁻¹".
- **Fix:**
  - Add a coverage-vs-τ curve for levels C and D, including both-direction counts.
  - State the maximum magnitude and the damping-ratio equivalent in the abstract.
  - Discuss whether any reversal changes a decision relative to a damping criterion.
- **Needs:** existing data only.

### M7. Preregistered analyses omitted from the paper
- **Where:** l. 403–407 ("Three secondary tests form a Holm family"); l. 643
  (`\input{sec_corridors.tex}`); l. 612–613 ("The top equal-budget corridor …").
- **Problems:**
  1. `sec_corridors.tex` is an empty placeholder, and `H10_gate.json` does not exist. The
     preregistered PRIMARY GATE H10 and analyses H9/H11 are missing without explanation.
     Yet §VI-G uses "the top equal-budget corridor" without ever defining corridors.
     `\corRobust{none}` is a default for a missing file, not a result.
  2. The statistical plan calls F_conf = {T1, T2, T3} the "Confirmatory family". The paper
     calls them "secondary" and reports only T3.
     - T2, the GOLD-B paired sign-flip test on the literal 8-condition set, has p = 0.061
       and fails at α = 0.05.
     - With 6 clusters the smallest attainable one-sided p is about 1/64, so the test was
       underpowered by design.
     - The prereg promises that negative results get "the same prominence as positive ones".
- **Fix:**
  - Add the corridor results, or a deviation entry declaring H10 incomplete.
  - Report T1–T3 with raw and Holm-adjusted p.
  - State the power limitation of T2.
  - Use the prereg's labels (confirmatory vs secondary).
- **Needs:** reporting of existing or missing results. H10 may need to be completed.

### M8. C3: over-reading of a small-perturbation result and missing tests on realistic actions
- **Where:** abstract l. 61–67; C3; §VI-D (l. 526–553); §VI-E; Discussion l. 691–694;
  Table V.
- **Problems:**
  1. At γ = 1.10 the median ρ is 0.999, which is Taylor's theorem. The planning-relevant
     actions are the ones evaluated in §VI-H: doublings (γ = 2), outages (γ = 0) and new
     lines. D_tot is never scored on them, although the data exist.
  2. "Held-out" (l. 125–126, l. 141–142) has no force for a predictor with no fitted
     parameters. Only the static baseline S1 was selected.
  3. "The ranking is unaffected" by saturation (l. 553) is contradicted by Fig. 6:
     - right panel: predictions of 0.25–0.38 have finite effects of 0.02–0.10;
     - left panel: several conditions have ρ_Dtot ≈ 0.65–0.8.

     Saturation probably reflects the max over modes (another mode becomes rightmost).
     Report top-1/top-3 errors per condition instead.
  4. Predictor and truth share the SPR semantics, which re-dispatch every reference after
     the change. The alternative "reference-preserving" semantics (theory note C4,
     old E34 RP files) is what happens if nobody resets setpoints. It is not reported, yet
     the re-equilibration term, the claimed distinguishing feature, depends on this choice.
  5. The cost statement (l. 692–694: "one factorization … and one back-substitution per
     candidate branch") does not describe the implementation. `_sens.Engine.derivatives`
     evaluates four central-difference Jacobians of the DAE per candidate. An analytic
     version needs Hessian-vector products of f and g, as well as the back-substitution.
  6. The abstract juxtaposes numbers from different sets: 0.491 is the eligible-set |P_e|,
     while 0.499 is the all-64 per-condition best static. On the eligible set, the best
     static is also 0.491.
- **Fix:**
  - New evidence: score D_tot on doublings and outages; report GOLD-B under RP.
  - Rewording: lead the abstract with the D_tot-vs-D_fro comparison; correct the cost
    statement; relabel the port-form row as a verification; replace "unaffected" with
    per-condition top-k statistics; use one condition set for every number in the abstract.

### M9. Model and SPR definitions are incomplete or incorrect in the paper
- **Where:** l. 196–199, l. 321–342 (SPR and Remark 2).
- **Problems:**
  1. "Every unit keeps its scheduled active power and terminal voltage" is impossible
     under a branch change, because losses change. In the implementation the bus-39 slack's
     P is free and its angle is fixed (theory note C4; `_sens.py` l. 148–149).
  2. With D = 0 and no governors, the equilibrium manifold has angle and frequency
     directions. R_w is nonsingular only once a gauge and a frequency are pinned. Neither is
     stated, and the invariance of D_tot to that choice is not shown.
  3. Remark 2 claims "the SPR equilibrium does not depend on a, so ẇ = 0". This is
     incorrect.
     - The freed references move: for an AVR gain, SPR sets v_ref = |V| + e_fd/K_A.
     - The correct statement, as in the theory note R1–R2: (x*, z*) do not move; the
       references do; they enter f and g affinely, so ∂_wA[ẇ] = 0 and D_tot = D_fro.
     - The same affinity argument explains why the frozen derivative for V_set is
       "identically zero" (l. 569–570). That is a structural fact, not a result.
  4. The transverse quotient, and the vector on which the MAC is computed, are not defined
     precisely enough to reproduce. The state dimension changes between S and S∪{i}, which
     makes the choice of vector for the MAC important.
- **Needs:** rewording and definitions.

### M10. Claims about fixed rankings that the data contradict or do not establish
- **Where:**
  - C2 l. 135–137 ("The best fixed ranking … is materially wrong in 0.38 … and does not
    transfer across policies");
  - §VI-C l. 490–494; Table IV (L6);
  - l. 521–523 ("Seven of eight are positive").
- **Problems:**
  1. The DP ranking maximizes pairwise agreement, not regret. The regret-minimizing fixed
     ranking has median regret 0.345 (FULL) and 0.030 (EM), not 0.376 and 0.034 (§4,
     item 5).
  2. Transferred rankings are barely worse than the in-policy oracle. Off-diagonal regret
     is 0.402 vs 0.376 on the diagonal (new holdout); accuracy is 0.725 vs 0.763 (all 52
     policies). The heatmap in Fig. 5
     is nearly uniform. The data say that fixed rankings are *uniformly* imperfect, not that
     they fail to transfer. The "24 distinct optimal orders" can be near-ties.
  3. Layers L1 ⇐ L2 ⇐ L3 ⇐ L4 are logically nested: a level-C nested reversal implies the
     A and B layers and an arbitrary-context reversal. L7 uses the arbitrary-context global
     rule, whose EM counterpart at EC is 0. "Seven of eight" therefore overstates how much
     independent evidence exists.
- **Fix:**
  - Report the regret-optimal ranking; the computation is trivial.
  - Replace "does not transfer" with the measured increment.
  - Present the layers as nested and give the EM version of L7.
- **Needs:** a cheap computation and rewording.

### M11. The framing of grid-strength indices is a strawman, and the sequencing comparison lacks a static baseline
- **Where:** abstract l. 45 ("Weak-bus and grid-strength indices rank locations once, for a
  given network"); Introduction l. 92–104; Discussion l. 661–675.
- **Problem:**
  - Multi-infeed indices (gSCR, gOSCR) are defined *for a given set of converter buses*,
    that is, per portfolio, and the paper's own ΔgSCR(H4) baseline is portfolio-conditioned.
  - In C2, no static index is evaluated at all. The only comparator is a fixed ranking
    fitted to the dynamics themselves.
- **Fix:**
  - Rewrite the framing to separate context-free weak-bus lists from portfolio-conditioned
    strength indices.
  - Add portfolio-conditioned gSCR as a sequencing baseline (see M4, item 3).
- **Needs:** rewording and a cheap computation on existing portfolios.

---

## 6. MINOR comments

1. **Population mismatch (l. 428–432, l. 462–465).** "3263 new-holdout level-C pairs",
   2172/1091 and 94 % include the base-unstable policies H07, H08 and H15. The coverage
   figures use base-stable denominators. Use one population: base-stable gives
   3067, 1990/1077 and 95.7 %. Pair counts are combinatorially inflated (87 distinct
   policy-unit combinations, 1513 distinct S2 marginals); report those counts too.
   *Rewording.*
2. **Abstract length.** The abstract is about 260 words; IEEE asks for at most 250.
   *Rewording.*
3. **Fig. 3 (`fig:examples`) is never referenced in the text.** *Rewording.*
4. **Missing denominator (l. 621–622).** "(50 of the V4 cases)" should read "50 of 1296
   (policy, action) cases". *Rewording.*
5. **Internal labels.** Replace or define H4 (vs V4), S1, TX3, "A3 bar", "P4
   incompatibility", HARDENING_H01, and N24 (the paper uses N24 while the data use
   HARDENING_H24). *Rewording.*
6. **Table VI.** The caption says "8 preregistered policies" while row A says "of 7".
   Explain that HARDENING_H02 has an unstable ALT base. Row C2 ("ALT total vs ALT finite,
   minus |P|") is cryptic, and the C2 margin (0.206 vs the 0.20 bar, with 3/8 policies
   at 0.196) should be stated. *Rewording.*
7. **Table IV, L7.** Add the EM same-mode coverages (EC 0.0, EM-f 0.5, EM-u 0.6, EMC 0.7)
   next to the global ones. *Rewording.*
8. **Coverage intervals.** A percentile bootstrap on 17/19 gives an upper limit of 1.00 and
   a lower limit (0.74) below the gate (0.75). Report exact binomial intervals
   (Clopper–Pearson [0.67, 0.99]) as descriptive, and state that the gate is on the point
   estimate. *Rewording.*
9. **Mixing (l. 650–656).**
   - "19 % of the variance … to the controller term and the rest to frequency movement"
     ignores Cov(C, F) ≠ 0; variance shares are not additive.
   - ρ = 0.548, p = 0.0163 is identical for four frequencies and for both n_rev and
     n_nested. Report the distribution and ties of n_rev.
   - With n = 19, the p ≤ 0.01 bar has low power at ρ ≈ 0.55, so "not confirmed" is more
     accurate than "negative".

   *Rewording.*
10. **Eq. (3).** The bracket notation [·]_{S0}^{Sm} is undefined, and d_{jj_r} clashes
    visually with d_{ij_r}. Define the bracket and use k_r for chain elements.
    *Rewording.*
11. **Proof of Proposition 1.** "(c) is the contrapositive" should read "(c) follows from
    the contrapositives of (a) and (b)". *Rewording.*
12. **Terminology.** "Discrete curvature" collides with the standard *total curvature* of
    submodular functions (Conforti & Cornuéjols 1984). Prefer "second difference" or
    "interaction term". *Rewording.*
13. **Remark 1.** Remark 1 shows that curvature heterogeneity is *necessary* for ranking
    failure, not sufficient. "What defeats a fixed ranking is the heterogeneity" should say
    so. The paper never measures heterogeneity; the distribution of d_ik − d_jk on EM-clean
    squares would link Remark 1 to the data. *Rewording or optional analysis.*
14. **"What the Laplacian does not capture" (l. 696–702).** This paragraph cites
    first-campaign exploratory numbers (Jaccard 0.2; 3 of 9 events) with no method or
    preregistration. Drop it or move it to the supplement. *Rewording.*
15. **Reproducibility of the model.** Table I lacks:
    - machine data (reactances, inertias, time constants) and their source;
    - the PSS block structure;
    - the bases of the PLL and PI gains;
    - the GFL rating convention;
    - the definition of the "harmonized benchmark".
    *Rewording or additional material.*
16. **Data availability (l. 736–739).** A branch name is not an archive. Provide a public
    URL and an archival DOI, and state the software versions (the run status lists
    Python 3.13, NumPy 2.5, SciPy 1.18). *Rewording.*
17. **"Byte-identical" (l. 389).** `CDWH_DETERMINISM.json` shows identity on summary fields
    for 20 tasks. Say "identical on all summary fields". *Rewording.*
18. **Figures.**
    - Fig. 1: node labels are illegible at column width.
    - Fig. 2: orange vs black crosses are unexplained.
    - Fig. 4 (left): the "49 % > +τ_res, 44 % < −τ_res" split is dominated by
      fast/unstable squares (|d| up to 10⁴). Show EM-clean squares separately.
    - Fig. 5 (right): policies are unlabeled and old/new are not separated.
    - Fig. 6 (right): mark the saturating branches.
19. **Topology section.**
    - All 45 creations are outages, and 34 of them come from the discovery policies D01 and
      D03; only 3 come from new-holdout policies.
    - EM removal is essentially a P4 phenomenon.

    State this, since the preregistration treats the new-holdout data as the confirmatory
    part. *Rewording.*
20. **l. 419–420.** "Every tested policy has a stable-context and a nested reversal of α⊥"
    should state that this is level A and holds for all 24 new policies, including
    base-unstable ones. *Rewording.*
21. **Table V.** The "ρ elig." entry for the port form is computed on 2 conditions (only two
    eligible conditions are policies); say so or omit it. *Rewording.*
22. **Missing references.**
    - Gautam et al. 2009 and Quintero et al. 2014 (TPWRS) on location-dependent
      displacement;
    - Desai et al. 2026 (arXiv 2601.05070) on IEEE-39 with replaced SGs;
    - Summers & Kamgarpour 2019 and Das & Kempe 2011 on non-submodular greedy selection.

    Gautam 2009, Desai 2026 and Summers & Kamgarpour 2019 are already in the authors' own
    gap matrix.
23. **Stable family and lattice structure.** The stable portfolios do not form a lattice or
    a down-closed family: unions and subsets of stable sets can be unstable. Say explicitly
    that sub/supermodularity is a property of α⊥ on all of 2^V, and that planners act only
    within the stable family, where matroid or greedy theory would not apply directly.
    *Rewording.*
24. **l. 568.** "agree to 2.2e-07 relative" refers to the controller coordinates. The
    V_set ratio behind Fig. 7 (right) (median ≈ 1.3e11 in `H07_goldb_gate.json`) divides
    by a frozen term that is numerically zero. Plot it as "frozen ≡ 0", not as a ratio.
    *Rewording.*

---

## 7. What would change my assessment

- **Towards acceptance:**
  - M1, M2, M6, M7, M9 and M10 corrected; these are mostly rewording on existing data;
  - a mechanism for level-D reversals (M3);
  - the fast aperiodic modes shown to be physical, or C2 restated in EM terms (M4);
  - one additional network or a voltage-controlling library converter reproducing
    same-mode nested reversals at material magnitude (M5);
  - D_tot evaluated on doublings and outages (M8).
- **Towards rejection:** if the fast aperiodic modes turn out to be numerical artefacts
  (near-singular g_z), and the same-mode reversals do not appear in a second system or in
  a voltage-controlling library model. In that case what remains is a known sensitivity
  evaluated on one benchmark.
