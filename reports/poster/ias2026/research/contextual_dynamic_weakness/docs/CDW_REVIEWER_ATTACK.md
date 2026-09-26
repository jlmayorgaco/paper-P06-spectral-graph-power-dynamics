# CDW reviewer attack (H22)

Each question is answered with the evidence file, or with the limitation it runs
into.

- Paths are relative to `contextual_dynamic_weakness/`.
- "New holdout" means the 24 preregistered maximin policies HARDENING_H01–H24 (19
  base-stable), frozen in `docs/CDW_HARDENING_PREREG_V1.md` (commit 05b507e3).

## Reviewer 1: theory and novelty

**Q1. Isn't this just eigenvalue sensitivity?**

Partly, and the paper says so.
- The total (re-equilibrated) derivative is classical: Smed 1993; Nam et al. 2000;
  Mendoza-Armenta & Dobson 2016; Li et al. 2019; the Joswig-Jones et al. 2026
  preprint for line admittances (`results/CDW_LITERATURE_GAP_MATRIX.csv`).
- What is not classical is its use as a portfolio-conditioned ranking of branch
  reinforcements, validated against finite re-solved outcomes:
  - on held-out policies and fresh envelope draws;
  - at three step sizes;
  - against a discovery-selected static index and the per-condition best of nine
    static indices.
- Results: `results/hardening/H07_goldb_gate.json`, `H06_metrics.csv`.
- The contextual sign-reversal results (C1, C2) are not sensitivity results. They
  come from the exact 512-portfolio census.

**Q2. Isn't sign reversal merely mode switching?**

No. Four levels are graded (`H03_gate.json`):
- A: global α.
- B: the critical mode of S is MAC-matched (≥ 0.8, abs(Δf) ≤ 0.15 Hz) to the critical
  mode of S∪i.
- C: B, with both modes in 0.1–2 Hz.
- D: C, with the two contexts' critical modes matched to each other.

At level C, 17/19 base-stable new-holdout policies have a nested reversal. Every
one of them also has a level-D reversal. A level-C marginal is by construction
the displacement of one tracked EM mode that remains rightmost.

Limitation: the intermediate squares of the certifying chain can still involve
mode switches. 94 % of the level-C reversals have an EM-clean witness.

**Q3. Isn't it caused by fast unstable modes?**

- **Reversals: no.** Level C restricts to the EM band.
- **Fixed-ranking failure: largely yes.**
  - In the EM stratum the oracle ranking orders 0.97 of pairs, with 0.03 material
    regret. Overall the figures are 0.75 and 0.38.
  - 92 % of the materially wrong decisions (2089/2265) are choices that push the
    rightmost eigenvalue out of the band (`H04_gate.json: regret_anatomy_new`).

  The paper states this as C2's anatomy, not as a caveat.

**Q4. Is the weak-corridor result just corridor size?**

This was tested with equal intervention budgets (Σ abs(Δγ) = 0.5) and size-matched
nulls: 500 arbitrary groups and all connected groups for k ≤ 4.
- Verdict: `results/hardening/H10_gate.json`.
- The unnormalized ×1.5 ranking is no longer used as evidence.

**Q5. Isn't the topology result already known from transmission switching?**

- **Known:** topology switching changes small-signal margins (Saric 2015; Li,
  Chiang & Du 2018), and Braess-type effects exist (Coletta & Jacquod 2016).
- **What is added:**
  - an EM-only reclassification (tracked pre-action EM mode; 50 fast-mode-dominated
    V4 cases excluded);
  - both removal and creation of EM incompatibility by single admissible actions;
  - seven static topology scores failing a preregistered prediction rule against the
    EM effect.
- **Files:** `H15_topology_gate.json`, `H16_static_prediction.csv`.
- **Scope:** no novelty is claimed for topology switching itself.

**Q6. What is actually mathematically new?**

Nothing in the algebra (chain identity; Propositions A–E). What the paper adds:
- It uses them as certificates: a nested reversal certifies curvature of the matching
  sign with a magnitude bound.
- It separates two things: sub/supermodularity, and fixed-ranking adequacy, which is
  governed by curvature heterogeneity (Remark 2 of the theory note).
- It shows a submodular function can have arbitrary-context reversals, which is why
  nestedness is tested.

See `theory/CDW_CONTEXTUAL_CURVATURE_THEOREM.md`.

**Q7. Why should results on IEEE-39 generalize?**

They are not claimed to generalize.
- The cross-model holdout shows the limits: the ranking *method* transfers to a
  frozen WECC library GFL chain (ρ 0.99 vs 0.79 for abs(P)), but reversals and the
  specific branch ranking do not (`H18_gate.json`).
- The claim is structural: rankings must be computed per portfolio. The phenomenon's
  existence depends on the converter control.

**Q8. Is the best fixed-ranking baseline fair?**

It is maximally favorable to fixed rankings.
- The DP ranking is fitted exactly on each policy's own stable contexts, which is
  in-sample.
- Out of sample (transfer across policies), material regret is 0.40 against 0.38
  in-sample (`H04_transfer_regret.csv`).

## Reviewer 2: numerics, statistics, modeling

**Q9. Are policies and envelope draws independent samples?**

No. They are designed points (maximin LHS; declared bounds).
- **Wording:** every fraction is "coverage across tested policies/draws".
- **Intervals:** bootstraps describe the stability of a summary across the design
  (`docs/CDW_HARDENING_STATISTICAL_PLAN.md`).
- **p-values:** secondary, with explicit nulls, in a Holm family.

**Q10. Doesn't the alternative converter invalidate the line ranking?**

It invalidates the specific ranking as a converter-independent fact: the median
Kendall τ between models is 0.02. It does not invalidate the method, since the ALT
total sensitivity predicts ALT outcomes. A planner must compute the ranking with a
validated model of the installed converters. This is consistent with the earlier
cross-tool finding of 7/12 signs.

**Q11. Why does a planner need contextual ranking instead of exhaustive eigenanalysis?**

Exhaustive re-solution is what the ranking is validated against.
- **Cost of the total derivative:** one Jacobian factorization per portfolio plus one
  back-substitution per candidate. Exhaustive evaluation costs one equilibrium solve
  and eigenanalysis per candidate and magnitude.
- **The finding:** the cheap *network-only* (static) or *base-portfolio*
  (conventional sensitivity) shortcuts are unreliable (ρ ≈ 0.5 and ≈ 0.03), while the
  portfolio-conditioned derivative is sufficient for ranking.

**Q12. Were thresholds changed after seeing results?**

No gate threshold was changed. The deviation log records every change.
- **Execution fixes:**
  - the relaunch;
  - the import order;
  - the bound-check tolerance;
  - the H06 worker-pool rerun.
- **H17 refinements, recorded before any H18 matrix case was computed:**
  - exact structural-zero removal;
  - limiter classification.
- **H6 eligibility:** the literal eligibility rule (R0 ≤ 1e-8) was *applied as written*.
  It reduced the primary GOLD-B set to 8 conditions, and the full 64-condition set is
  reported as a sensitivity analysis. Both pass.
- **Disclosed mistake:** one ALT datum (H02/H4) was seen during qualification.

**Q13. Model adequacy: no converter current limits, D = 0, constant-power loads, phasor only.**

These are limitations, stated in the paper.
- **Fast control-mode instabilities.** These drive the fixed-ranking failure. They
  are those of an idealized converter without current limits, so their practical
  weight is model-dependent.
- **D = 0 and no governors.** This creates exact structural modes, handled by the
  transverse quotient.
- **No EMT validation.** A previous EMT line on this benchmark was stopped by its
  own qualification rules.

**Q14. What is industrially useful?**

Two things:
- a warning that weak-bus lists are unreliable for replacement sequencing when
  control-mode instabilities are possible;
- a cheap and validated way to rank reinforcements for the portfolio actually
  planned.

The magnitudes of EM reversals are small (median 0.013 s⁻¹). The paper does not
claim that they change most EM decisions.

**Q15. Statistical strength of GOLD-B with 8 eligible conditions.**

- **Effect size:** the lower 95 % bootstrap bound of the median paired advantage is
  0.477, against a bar of 0.20.
- **Test:** the cluster sign-flip p is 0.061, because there are only 6 clusters. It is
  a secondary summary.
- **Full set:** over the 64 solved conditions the bound is 0.484. Dtotal beats every
  static index in every condition.
