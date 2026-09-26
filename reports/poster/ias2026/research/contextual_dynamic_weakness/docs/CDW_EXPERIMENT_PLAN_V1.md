# CDW experiment plan V1

This is the executable plan for the preregistered experiments. The rules are
in `docs/CDW_PREREG_V1.md`.
- **Orchestrator:** `experiments/cdw/CDW_MASTER_RUN.py`.
- **Status file:** `results/CDW_RUN_STATUS.json`.
- **Cost unit:** one portfolio evaluation, about 0.35–0.45 s on one thread
  (from the timing pilot). Workers are chosen at run time from the free RAM,
  with ≥ 2 CPU threads reserved and BLAS pinned.

| exp | inputs | outputs | frozen truth | metrics | baseline | cost (evals) | pass/fail |
|---|---|---|---|---|---|---|---|
| E1 | D01–D15, H01–H24; 512 V9 portfolios each | `results/CDW_E1_contextual_marginals.parquet`, `CDW_E1_summary.csv`, `CDW_E1_portfolios.parquet` | α_⊥ (SPR) | sign fractions, entropy, reversal indicators (global / stable / tracked), p*, regret | best node-only ranking (oracle-fitted) | 19 968 | A1, A2, A3 |
| E1b | E1 | `CDW_E1b_submodularity.csv` | E1 α | violations, magnitudes, minimal counterexamples | — | 0 | descriptive + HS |
| E2 | P4, G_S × envelope draws (TX4 4×100, CDW 4×40) × core; V9 × 10 CDW draws per envelope at P4 | `CDW_E2_robust_contextuality.parquet`, `CDW_E2_robust_summary.csv` | α_⊥ | sign and reversal persistence, Kendall, Spearman, top-3, witness change | nominal | 17 920 + 20 480 | A4 |
| E3 | 5 discovery + holdout conditions × {H4, V9}; NG, NP, NK, NV, NL | `CDW_E3_node_sensitivity.parquet` | finite re-solves (same semantics) | IV, frozen vs total sign and Kendall, prediction ρ and sign accuracy | frozen partial | ~(43 params × 5 solves) per condition | IV, H2-node |
| E4 | as E3; 46 branches | `CDW_E4_link_sensitivity.parquet`, `CDW_E4_baselines.csv` | finite ×1.5 (SPR) | Spearman, Kendall, top-5, sign accuracy | S1–S9, Dconv, Dfrozen, Dport | ~(46 × 5 solves) per condition | GOLD-B, H2-link |
| E5 | E1–E4 | `CDW_E5_static_vs_dynamic.csv` | — | categories, Kendall(static, dynamic) | static | 0 | descriptive |
| E6 | L_B partitions K2–K4, TX groups; E4 conditions | `CDW_E6_corridors.csv` | finite ×1.5 on C | rank robustness, NA index, Spearman | Σ static, Σ single-line | ~15 solves per condition | H3 |
| E7 | 46 outages + 46 doublings; core at 10 policies; V9 at P4 | `CDW_E7_topology.parquet` | α_⊥, H, κ after re-solve | removal / creation counts, static-score ρ | ΔFiedler, ΔKirchhoff, ΔgSCR, ΔminSCR | 14 720 + 47 104 | H4-topology |
| E8 | boundary θ* (and P4); 50 control–link pairs | `CDW_E8_exchange.csv` | paired finite re-solves | compensation ratio | — | ~300 | ≥ 80 % pairs ≤ 0.2 |
| E9 | T1–T4 × families a, b, c × single vs plan | `CDW_E9_design.csv` | exact subset re-evaluation | Φ_T, α(T), conflict witness | single-boundary tuning | ~6 iterations × Φ evaluations | GOLD-D |
| E10 | E1 | `CDW_E10_shapley.csv` | E1 α | Shapley, variance, hidden-reversal flag | — | 0 | descriptive |
| E11 | L_B, Kron, dynamic SVD | `CDW_E11_spectral.csv` | E1, E4 | predictive ρ of graph scores | — | small | "predicts" rule |
| E12 | core × D and H; FC18 events | `CDW_E12_mixing.csv` | — | μ_mix range, ρ with reversal count, support transitions, Jaccard | — | ~700 T(jω) evaluations | Q1–Q4 rules |
| E13 | GM, POD, TS families on holdout | `CDW_E13_reduction.csv` | full α_⊥, H, κ | accuracy, H match, speedup, state ratio | full model | ~1 000 | GOLD-C |
| E14 | E13 families on contour Γ | `CDW_E14_certificate.csv` | full verdict | coverage, false certifications | — | ~1 600 contour points per case | false = 0 |
| E16 | stable core portfolios, D and H | `CDW_E16_modal_energy.csv` | eigendata | fraction with argmax energy ≠ limiting mode | — | ~620 | "systematic" rule |
| E17 | repository search | `docs/CDW_AFRICANO_MISSING_INPUTS.md` or a provenance file | — | presence of every required input | — | 0 | gate |
| E23 | static-injection converter | `CDW_E23_crossmodel.csv` | α_⊥ (static model) | reversal presence, link Kendall, corridor overlap | GFL results | ~5 000 | transfer rule |

**Ordering and dependencies.**
- E1 → (E1b, E10) → E5.
- E2, E3, E4, E6, E7, E12, E13, E14, E16 and E23 are independent.
- E8 and E9 need the E3/E4 IV gate for their derivative parts.
- E11 predictive tests need E1 and E4.
- E22 is optional.

**First minimal set (if the campaign had to be cut).** E1 + E2 (GOLD-A), E4
(GOLD-B) and E9 T1 (GOLD-D).
