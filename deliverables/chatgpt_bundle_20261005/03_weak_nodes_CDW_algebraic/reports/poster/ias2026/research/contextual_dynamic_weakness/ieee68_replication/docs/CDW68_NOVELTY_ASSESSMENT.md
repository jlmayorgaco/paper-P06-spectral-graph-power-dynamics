# CDW68 novelty assessment (R20)

This update to `docs/CDW_NOVELTY_BOUNDARY.md` (hardening H20) takes the IEEE-68 replication into account.

**Evidence base.**
- The 57-entry hardening search in `results/CDW_LITERATURE_GAP_MATRIX.csv`.
- The 9 rows added in this campaign: `results/CDW68_LITERATURE_GAP_ADDENDUM.csv`. Each was checked through Crossref or arXiv on 2026-09-13.
- IEEE Xplore full-text search was again not available. Absence statements therefore read "to our knowledge, within this search". The word "first" is not used.

## 1. What the replication changed

| hardening item (H20 §2) | IEEE-39 status | IEEE-68 result (this campaign) | status after R19 |
|---|---|---|---|
| 1. Nested same-mode contextual sign reversal | 17/19 base-stable holdout policies (Model A); Model B: global 0/7, EM-tracked 1/7 | **0/16 level D and 0/16 level T for both models.** The EM-tracked marginals are almost all destabilizing (A: −0.0069 to +0.053 s⁻¹; B: −0.0008 to +0.057 s⁻¹). Nothing changes when governors are removed (NOGOV) or under the published simplification (SP33): 0/6 in each | **IEEE-39-specific evidence.** Not a general property of SG→GFL replacement |
| 2. Failure of the oracle fixed node ranking | material regret in FULL and EM strata | **Not replicated.** FULL has no material preference. TRACKED: regret 0 and p* = 1.0 for both models. The LEVEL-D and EM strata cannot be evaluated | **IEEE-39-specific.** On IEEE-68 a fixed order ranks the tracked EM effect without error |
| 3. Portfolio-conditioned total branch sensitivity as a ranking | ρ ≈ 0.996 on new policies (Model A) | **Replicates under both converters.** A: ρ 0.996, advantage 0.80, top-5 1.0. B: ρ 0.981, advantage 0.55, top-5 0.8. The preregistered target is an immaterial control pole (median maximum finite effect 4e-5 s⁻¹). The post-hoc EM-tracked ranking also holds (A: ρ 0.998, median maximum effect 0.020 s⁻¹, 5 material branches per condition in the median) | **The one finding carried across two benchmarks and two converter models**, with the materiality caveat |
| 4. EM-only topology effects | IEEE-39 | not in the R-scope; not tested | unchanged (IEEE-39 only) |
| 5. Cross-model check with a frozen library GFL chain | IEEE-39, 4 candidates | IEEE-68, qualified ANDES transcription (Q1–Q4 pass), 16 conditions | extended to a second network |

## 2. Prior work that the IEEE-68 result touches

- **Mode-dependent displacement effects.**
  - Gautam, Vittal and Harbour (2009) report that displacing SGs by DFIGs improves some electromechanical modes and degrades others.
  - Quintero et al. (2014) and Eftekharnejad et al. (2013) report system-dependent changes in damping and participation as converter penetration grows.
  - On IEEE-68, CDW68 finds replacement effects on the rightmost EM mode that are one-sided across contexts. This fits mode- and system-dependence and does not show context reversal.
  - None of these works tests nested portfolios, so they neither confirm nor contradict the IEEE-39 reversal.
- **Benchmark data.**
  - The published 68-bus benchmark has no governor model: Singh & Pal (2013) use constant mechanical torque, and Canizares et al. (2017) validate the same model.
  - The primary-frequency data come from PST (`data16m.m`, Chow & Cheung 1992; PSTess `tg.m`).
  - Reversal was absent both with and without these dynamics (REAL vs SP33). Governors and damping neither create nor remove it on this network.
- **Share-dependent stability boundaries.** Conte et al. (2026, preprint) study stability manifolds over controller parameters and inverter share. CDW68's conclusion that ranking quality depends on the portfolio sits alongside theirs. Neither study ranks branch reinforcements.
- **Total eigenvalue sensitivity.** The earlier conclusion stands: Smed 1993, Nam et al. 2000, Mendoza-Armenta & Dobson 2016 and Joswig-Jones et al. 2026 (preprint) already cover it. The IEEE-68 result adds no new formula.

## 3. What remains defensible as a contribution

1. **Transfer of the ranking.** A re-equilibrated (total) first-order sensitivity, taken at the planned converter portfolio (a known construction), ranks finite branch reinforcements of 1.10–1.50 with median Spearman 0.98–1.00 on held-out policies. This holds on two benchmarks and two converter models. It is compared head-to-head with nine static indices, the comparator being selected on discovery conditions.
2. **Materiality of that ranking.** The same statement on IEEE-68 carries its limit. The preregistered target is the rightmost transverse eigenvalue, a slow real control pole, and branch actions move it by at most 1.4e-3 s⁻¹ (A) and 3.4e-5 s⁻¹ (B). The EM-tracked ranking, which is material, was defined after the fact.
3. **A preregistered negative replication.** Nested contextual sign reversal and fixed-ranking failure appeared on IEEE-39 and did not replicate on IEEE-68. This was tested:
   - with governors and damping, with governors off, and without either;
   - under two converter models;
   - on 16 designed policies per model.

   The IEEE-39 findings are benchmark-specific evidence.

## 4. Wording after R19 (CASE B)

**Allowed**
- "On the IEEE 39-bus benchmark the same replacement moves the same electromechanical mode in opposite directions in nested stable portfolios. On the IEEE 68-bus benchmark with documented governors and damping this did not replicate (0/16 policies, both converter models)."
- "A re-equilibrated first-order sensitivity at the planned portfolio ranks finite branch reinforcements well on both benchmarks and under both converter models (median Spearman 0.98–1.00 on held-out policies)."
- "On IEEE-68 the preregistered target is a slow real control pole whose branch effects stay at or below 1.4e-3 s⁻¹. The ranking of the rightmost electromechanical mode, defined post hoc, reaches median Spearman 0.998."
- "Coverage over the designed conditions is X/Y." (deterministic fraction)

**Prohibited**
- "contextual weakness is a general property", "replicated reversal", "model-independent reversal".
- "fixed rankings fail" without the benchmark qualifier.
- "first", "novel sensitivity", "we introduce total sensitivity".
- Any probability of reversal, or any EMT claim.
- Any topology or Laplacian explanation of the IEEE-68 results. The Fiedler, effective-resistance and betweenness scores ranked finite effects with median ρ ≤ 0.29 on the discovery conditions of both models, and served only as baselines.
