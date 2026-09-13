# Author response to the two internal reviews (H31)

Manuscript: `reports/papers/cdw_contextual_dynamic_weakness/main.tex`.
Reviews: `docs/CDW_REVIEWER1_FINAL.md` (R1) and `docs/CDW_REVIEWER2_FINAL.md` (R2). Both recommend major revision.

Rule applied (H31): fix a comment only where the existing evidence permits, and add no unsupported claim. Every analysis added after the preregistered results is labelled post hoc in the paper and listed in the deviation log (entry of 2026-09-12T23:08:13-0500). No preregistered verdict was changed.

Each item below is marked:
- **FIXED**: changed in the paper, supported by existing data or by a post-hoc analysis on existing solves.
- **PARTLY FIXED**: reworded or narrowed; the full request needs new evidence.
- **OPEN**: needs new evidence that this campaign does not have. Listed as a limitation.

## 1. Changes that alter the headline claims

| # | comment | status | what changed |
|---|---|---|---|
| R1-M1, R2-M3 | Flagship pair is a between-context mode switch (0.744 vs 1.310 Hz, level D false) | FIXED | Fig. 1 is now the strongest level-D pair (N24, unit 33, {31} → {31,34,38}, +0.053 → −0.042 s⁻¹, 0.64/0.74 Hz, MAC 0.97). Fig. 3 shows three level-D pairs with a deterministic tie-break. The headline uses level D (17/19). |
| R2-M3 | Abstract describes one direction (s→d), which holds in 9/19 | FIXED | The wording is now direction-neutral ("opposite directions in two nested stable portfolios"). Direction counts are reported: level C 9/15/7 and level D 6/14/3 (s→d / d→s / both). |
| R1-M2 | "A nested reversal is not compatible with submodularity" is false for d→s; "certificate" is automatic | FIXED | Proposition 1 now reads: s→d contradicts submodularity and d→s contradicts supermodularity. The paper states that the bound holds for every nested pair by construction and is used only to *locate* an interaction term. The "certificate" contribution is removed. Old Fig. 4 (right), which plotted the bound, is replaced by the coverage-vs-τ curve. The left panel now shows EM-clean one-step squares only. |
| R1-M6 | Practical significance and τ dependence | FIXED | The abstract states the median magnitude (0.013 s⁻¹, about 0.3 % damping ratio) and that the effect vanishes above 0.05 s⁻¹. Fig. 5 (right) and supplement Table S-tau give coverage at τ = 0.01–0.05 for levels C and D and both directions. |
| R1-M4, R2-M1 | C2 is driven by real, aperiodic eigenvalues of 10²–10⁴ s⁻¹; mechanism not identified | PARTLY FIXED | Post-hoc characterization on existing solves: 99 % of these critical eigenvalues are real (5/50/95 % quantiles 20/176/3456 s⁻¹). They appear at 5–7 converted units. Participation is concentrated on the PLL angle (0.39) and q-axis current (0.36). In six fractional-replacement sweeps α⊥ jumps from about 0 to 3.8×10³–1.1×10⁵ s⁻¹ within one 0.025 step, exactly where σ_min(g_z) is 0.002–0.007. This is the signature of a singularity-induced bifurcation (Venkatasubramanian et al. 1995, now cited). C2 is rewritten as "a fixed list cannot anticipate this singularity crossing". The phrase "fast converter-control modes" is removed. **OPEN:** current limits, delays, ZIP loads and an EMT spot check. |
| R2-M2, R1-M10 | Regret baseline overstates the failure; "does not transfer" not supported | FIXED | Post hoc: regret-optimal fixed order 0.345 (DP over placed prefixes); stability-screened order 0.065 on contexts with a stable option (plain order 0.267); 22 % of contexts have no stable option. The paper says "the failure is present in sample and transfer adds little" (0.402 vs 0.376) and reports the median Kendall between optimal orders (0.72). |
| R1-M4 item 3, R1-M11 | No portfolio-conditioned static screen for sequencing | FIXED | Post hoc: choosing the unit that maximizes gSCR(S∪{i}) errs in 0.76 of decisions (EM 0.68); minimum SCR errs in 0.86. The introduction now separates context-free weak-bus lists from portfolio-conditioned strength indices. |
| R2-M5 | Static baselines cannot respond to conditions; a learned fixed list is the fair competitor | FIXED | Post hoc: a leave-one-cluster-out fixed branch list reaches 0.818 on the new policies and 0.959 on the draws, against D_tot 0.996 and 0.986. D_tot beats it in every condition, but the median paired margin is only 0.03 (draws) and 0.17 (new policies), and the paper says so. D_fro is now the lead comparison. The paper states that static indices do not change across conditions. |
| R1-M8 | D_tot never scored on large actions | FIXED | Post hoc, at six new-holdout policies: Spearman with the finite effect of doubling is 0.992 and of outage 0.885 (top-5 0.80). "The ranking is unaffected by saturation" is replaced by "the first-order prediction overstates the largest reinforcements". |
| R1-M8 item 5, R2-M8 | Cost claim misdescribes the implementation | FIXED | The cost claim is removed. The paper reports 1.46 s per branch for the FD-Jacobian total derivative against 0.61 s for a finite re-solve with eigenanalysis. The Discussion answers R2-Q5 directly: in-context unit selection costs nine eigenanalyses per step, and C2 is diagnostic. |
| R1-M7, R2-M5 item 5 | T1–T3 not reported | FIXED | Section V reports T1 (raw p 1.9e-6, Holm 5.7e-6), T2 (0.061, fails with 6 clusters) and T3 (0.016, Holm 0.033). |
| R1-M7, R2-M7 item 5 | H10 corridors incomplete and omitted | FIXED WHEN H10 COMPLETES | H10 was still running when the reviews were written. The corridor section and supplement tables are generated from `H10_gate.json` / `H11_corridor_table.csv` once the preregistered run finishes; see the final report for the verdict. |
| R2-M6, R1-M5 | Cross-model: pinned pole, non-uniform tested-policy rule, voltage control untested | FIXED (analysis); OPEN (second network) | The tested-policy rule is now applied to every test (D becomes PARTIAL, 0.71; C2 0.203 with 3/7 below the bar). The pinned −0.100 s⁻¹ REPCA1 s5_xi pole is disclosed. EM-tracked analysis: stabilizing marginals at 3/7, nested reversal at 1/7 (N04, the same policy where the custom GFL has a level-C reversal on that lattice). QFLAG=1 variant (library voltage control on, gains unchanged, outside the freeze): 0/7 global and 0/7 tracked, so the data do not support "voltage control enables reversals". **OPEN:** a second network (e.g. IEEE 68-bus). |
| R2-M7 | Hand-written deviation-log timestamps contradict git | FIXED | A correction entry with system-clock time reconstructs every earlier time from commits and file times. The original entries were not rewritten. The paper's limitations paragraph discloses it. |
| R2-M7 item 4 | Regret anatomy post hoc and unlabelled | FIXED | Labelled post hoc in the deviation log and in the paper ("post-hoc analyses requested by internal reviewers are labelled as such"). |

## 2. Definitions and correctness

| # | comment | status | change |
|---|---|---|---|
| R1-M9.1–2 | SPR ignores slack and gauge | FIXED | Section IV: the slack at bus 39 absorbs the loss change and fixes the angle reference; the frequency is nominal, which pins the symmetry directions. |
| R1-M9.3 | Remark 2 wrongly asserts ẇ = 0 | FIXED | The states do not move; the freed references do, and they enter f and g affinely, so ∂_wA[ẇ] = 0 and D_tot = D_fro. The zero frozen derivative of voltage setpoints is stated as a structural fact. |
| R1-M9.4, R2-minor 18 | MAC and transverse quotient under-specified | FIXED | Section II-C gives the 78 rectangular bus-voltage components, unit norm, phase fixed at the largest entry, the candidate set (0.1–2 Hz, Re ≥ −1 plus the critical one), and |Δf| ≤ 0.15 Hz. |
| R1-minor 13 | Remark 1: heterogeneity is necessary, not sufficient | FIXED | Worded as necessary. |
| R1-minor 23 | Stable family is not a lattice | FIXED | Stated in the Remark. |
| R1-minor 10–12 | Notation and "curvature" terminology | FIXED | "Second difference (interaction term)" replaces "discrete curvature". Chain elements use k_r. |
| R1-minor 15, R2-minor 13 | Table I incomplete; "harmonized benchmark" undefined | FIXED | Table I and Section II-A cite the ANDES IEEE 39-bus dynamic case as the source of machine, AVR and PSS data. The unit-specific ranges are given (K_A ∈ {10.1, 40}, T_E 0.25–0.95 s, M 4.9–8.4 s). The earlier single "T_E = 0.25 s" was wrong and is corrected. "Harmonized" is removed. |
| ALT description | "four decoupled states per converter" | FIXED | Three per converter (REPCA1 s2_xi, REECB1 PIQ_xi and PIV_xi), verified in the case files. |
| Intro citation | Gautam 2009 and Quintero 2014 characterized too strongly in a draft | FIXED | Gautam: some modes improve and others degrade. Quintero: penetration changes participation and adds converter-dominated modes. Neither is claimed to show dependence on which units are displaced. |

## 3. Reporting and presentation

- **Pair statistics:** base-stable only (3067 level-C pairs, 87 distinct policy–unit combinations). The strongest pair is chosen among base-stable policies. (R1-minor 1, R2-M9.2, R2-minor 10)
- **Coverage intervals:** the gate is stated to pass on the point estimate. The bootstrap interval reaching below it is reported, and the exact binomial interval [0.67, 0.99] is in the claim matrix. (R2-M9.1, R1-minor 8)
- **Robustness:** the gap-robust coverage (16/19) and the continuation check (60/60 endpoint marginals end on the MAC-matched mode) are added. (R2-M3)
- **Evidence layers (Table IV):** the table states that L1–L4 are nested, that L6 follows L5, and that L5 fails in the EM stratum. L7 now shows the EM values (EC 0.0, EM-f 0.5, EM-u 0.6, EMC 0.7). (R1-M10.3, R1-minor 7)
- **Cross-model table (Table VI):** 7 of 8 policies with the reason; C2 margin and 3/7 below the bar. (R1-minor 6)
- **Table V:** the port form is labelled a verification, with n = 2 for its eligible entry. The caption states the condition sets. (R1-minor 21, R2-minor 20)
- **Topology provenance:** 45 creations, all outages (34 discovery, 8 old, 3 new). The denominator "of 1296 V₄ policy–action cases" is given. Within-type correlations are about 0. C4 is demoted to a supporting observation. (R1-minor 4 and 19, R2-minor 1, 8, 9)
- **Mixing:** the three variance terms are reported (controller 24 %, frequency 106 %, 2cov −30 %). The result is described as "not confirmed" and the test's low power is noted. (R1-minor 9, R2-minor 7)
- **Laplacian paragraph:** dropped from the paper. (R1-minor 14)
- **"Byte-identical":** replaced by "identical in every summary field". (R1-minor 17, R2-minor 4)
- **"Held-out":** stated to qualify only the baselines, since D_tot has no fitted parameters. (R1-M8.2, R2-minor 22)
- **Abstract:** 232 words, one condition set per sentence. (R1-minor 2, R2-minor 2)
- **Jargon:** H4, S1, TX3, "A3 bar" and "P4 incompatibility" are removed or defined. (R1-minor 5)
- **Figures:** (R1-minor 18, R2-minor 16)
  - F1 is replaced by a level-D pair with larger labels, and its caption now explains the colors.
  - F2 uses white crosses on dark cells for contrast only, stated in the caption.
  - F4 heatmap has policy-group separators and labels.
  - F5 shows EM-clean squares only.
  - F7 no longer plots a ratio for V_set.

## 4. Open items that need new evidence (limitations in the paper)

1. **Second network** (R1-M5, R2-M8): e.g. IEEE 68-bus with sampled replacement chains. This is the single experiment most likely to change the assessment of C1.
2. **Omitted physics for the EM claims** (R2-M4): SG damping D ∈ {0.5, 1, 2}, governors, ZIP loads, unmatched reactive dispatch.
3. **Physical course of the singularity crossings behind C2** (R1-M4, R2-M1): converter current limits, measurement delay, network dynamics or EMT.
4. **Mechanism for the level-D reversals** (R1-M3): a homotopy in the replaced rating with the integrand's sign traced to the mode-shape change.
5. **Reference-preserving semantics for GOLD-B** (R1-M8.4): the RP files exist from the first campaign but were not re-evaluated on the new holdout.
6. **Analytic derivatives and a large system** (R2-M8): needed for any computational-cost claim.
7. **Domain-measure sensitivity** of coverage (R2-M9.4): uniform g, excluding AVR regions that destabilize the all-SG base.
8. **Data archive with DOI** (R1-minor 16, R2-minor 14): the repository branch is local; the review and submission bundles are the current archive.
9. **Submodularity ratio** of −α⊥ on stable chains (R1-M2 optional).

## 5. Expected effect on the scores

The rewording and post-hoc analyses address correctness, statistics, writing and figures. They do not address the two lowest scores, novelty (3/5) and model adequacy (3/3), which need items 1–4 of Section 4. The paper is accordingly framed as an empirical, preregistered case study on one benchmark.
