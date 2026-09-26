# CDW hardening — experiment matrix (frozen with the preregistration)

Code lives in `experiments/cdw_hardening/`. The orchestrator is `CDWH_MASTER_RUN.py`.
Checkpoints are `raw/H_<phase>/<task>.json`, one file per task, resume-safe. Status
is in `results/hardening/CDWH_RUN_STATUS.json` and logs in `logs/hardening/`.

The model evaluations use `.venv/tx3-analysis`, except H17/H18 ALT-WECC, which uses
`.venv/xtool-andes-gfl` with a private home. BLAS is pinned to one thread.

| phase | purpose | inputs | tasks (approx.) | evaluations | outputs | depends on |
|---|---|---|---|---|---|---|
| H01 | new-holdout census | 24 HARDENING policies × 512 V9 portfolios, with mode shapes | 384 | 12 288 portfolio solves | `raw/H_H01`, `results/hardening/H01_*` | — |
| H01e | fresh-draw core lattice | P4 × 40 draws × 16 | 40 | 640 | `raw/H_H01e` | — |
| H03 | nested reversal A/B/C/D + curvature witnesses | E1 raw (39 policies) + H01 (24) | 0 (analysis) | — | `H03_nested_*.csv/json` | H01 |
| H04 | fixed ranking, regret, transfer matrix, stratified | E1 + H01 | 0 | — | `H04_*` | H01 |
| H05 | GOLD-A statistics | H03, H04 | 0 | — | `H05_*` | H03, H04 |
| H06 | link sensitivity + finite γ ∈ {1.10, 1.25, 1.50} | 24 policies × {H4, V9} + 40 draws × H4 (88 conditions); conv (64 pkeys); port (24) | ≈ 480 | ≈ 88 × 46 × (4 A-evals + 5 re-solves) | `H06_links.parquet`, `H06_metrics.csv` | — |
| H07 | GOLD-B paired test | H06 | 0 | — | `H07_gate.json` | H06 |
| H08 | frozen vs total decomposition; node coordinates | H06 + 24 policies × H4 node specs (SPR) | ≈ 96 | ≈ 24 × 43 specs | `H08_*` | H06 |
| H09 | equal-budget corridors (3 budgets × 11 names; identical sets solved once) | 101 conditions | 101 | ≈ 101 × 25 | `H09_corridors.csv` | — |
| H10 | size-matched nulls (families A, B; k ∈ {2,3,4,8,12}; k = 1 from links) | 96 holdout conditions × 4 108 groups (A: 2 500; B: 1 608; groups identical across families solved once) + 46 singles at the 32 old conditions | ≈ 8 000 (chunks of 50) | ≈ 395 000 solves (≈ 0.45 s each, ≈ 2.5 h on 20 workers) | `H10_null.parquet`, `H10_gate.json` | H06 (k = 1 new), H09 |
| H11 | corridor cross-validation table | H09, H10 | 0 | — | `H11_table.csv` | H09, H10 |
| H12 | mixing at fixed frequencies | 63 policies (H4), then θ_ref at every ω_c(θ) | 64 | 64 + 63 frequency evaluations | `H12_mixing.csv` | — |
| H13 | two-factor decomposition + correlations | H12, H03, H04, E1 | 0 | — | `H13_*` | H12, H03, H04 |
| H14 | topology EM tracking | 16 policies × 82 topologies × 16 portfolios, with modes | ≈ 1 312 | ≈ 21 000 | `H14_topology_em.parquet` | — |
| H15 | H_EM vs H under topology | H14 + old E7 census raw | 0 | — | `H15_*` | H14 |
| H16 | static topology scores vs Δα_EM | H14 | 0 (plus static PF per action) | 82 PF | `H16_*` | H14 |
| H17 | ALT-WECC qualification Q1–Q4 | 8 policies | ≈ 8 + checks | ≈ 30 ANDES cases | `H17_qualification.json` | — |
| H18 | cross-model matrix (both models) | 8 policies × (16 + 12×3 + 3 + 6 + 1) on each model | ≈ 16 | ≈ 8 × 62 per model | `H18_*` | H17 |
| H19 | master evidence table | all | 0 | — | `H19_evidence.csv` | H03–H18 |
| DET | determinism reruns | H01 at HARDENING_H01; H06 at HARDENING_H01/H4 | small | — | `H_DET.json` | H01, H06 |

## Figure definitions (frozen)

The main figures are:

- **Figure 1 (concept).** A schematic of one network with the same intervention i
  under S₁ (stabilizing) and S₂ (destabilizing), with Δ_iα(S₁) < 0 < Δ_iα(S₂). The
  illustrated pair is the strongest level-C nested reversal on the new holdout.
- **Figure 2 (contextual reversal).** A heatmap of policy × unit, showing the
  fraction of stable contexts in which unit i is destabilizing. Rows cover the
  discovery, old and new holdout policies. Level-C nested reversals are overlaid,
  and the nested coverage statistic is annotated.
- **Figure 3 (strongest same-mode examples).** Three panels, one per level-C nested
  reversal with the largest magnitude on the new holdout (distinct units). Each
  panel plots Δ_iα against context size and marks S₁, S₂ and the same tracked mode.
- **Figure 4 (fixed ranking).** The distributions of p* and material regret rate (old
  vs new), and the train × test transfer heatmap.
- **Figure 5 (curvature).** The distribution of d_ij(S) and the minimal
  positive/negative witnesses, with a schematic of the chain accumulation.
- **Figure 6 (GOLD-B).** Left: paired Spearman per new-holdout condition (S1_absP,
  oracle static, Dfrozen, Dtotal). Right: Dtotal-predicted first-order effect
  (0.5·dα/dγ) vs finite γ = 1.5 effect, pooled.
- **Figure 7 (why total matters).** Frozen, indirect and total for the three
  rule-selected branches and for the node coordinate classes.
- **Figure 8 (topology).** Only if question 13 is YES. The finite EM-mode effect for
  admissible actions, with the static scores.

The supplementary figures are:

- **S1:** equal-budget corridors.
- **S2:** size-null percentiles.
- **S3:** mixing decomposition.
- **S4:** cross-model transfer matrix.
- **S5:** uncertainty robustness.
- **S6:** fast vs EM stratification.

The format is vector PDF, SVG and PNG (300 dpi). No decorative plots; every main
figure is cited by a manuscript claim.
