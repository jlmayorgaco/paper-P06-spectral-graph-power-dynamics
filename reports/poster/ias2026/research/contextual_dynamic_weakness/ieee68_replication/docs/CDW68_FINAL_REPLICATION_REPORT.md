# CDW68 final replication report: IEEE 68-bus confirmatory replication

- **Branch:** `research/cdw-ieee68-replication`, from `e03a5312`.
- **Preregistration:** commit `6f188531` (2026-09-13T09:19). Every section below cites the result file it reads.
- **Scope:** phasor-domain small-signal analysis. No EMT claim is made.

## 1. Executive decision

**CASE B**, by the exhaustive preregistered rule (prereg §17): the ranking method generalizes and the reversal does not. The gates came out as follows.

| component | gate | result |
|---|---|---|
| REV68A | level-D nested reversal in ≥ 50 % of base-stable policies, Model A REAL | **FAIL**: 0/16 (Clopper–Pearson 95 % upper bound 0.21) |
| REV68B | the same gate on Model B | **FAIL**: 0/16; level T also 0/16, so not PARTIAL |
| RANK68A | R12 gate, Model A | **PASS**: median ρ 0.996, advantage 0.80, top-5 1.0 |
| RANK68B | R12 gate, Model B | **PASS**: median ρ 0.981, advantage 0.55, top-5 0.8 |

**Consequence.** The paper pivots to portfolio-conditioned dynamic reinforcement ranking. The contextual reversal is reported as IEEE-39 benchmark-specific evidence.

**The caveat that travels with the pass.** On IEEE-68 the preregistered ranking target is the rightmost transverse eigenvalue.
- In every portfolio of both models it is a slow real control pole: the PSS washout at −0.067 s⁻¹.
- Branch reinforcements move it by at most 1.4e-3 s⁻¹ (A) and 3.4e-5 s⁻¹ (B), well below τ = 0.01.
- The electromechanical version of the ranking was defined after the fact. It passes as well (A: ρ 0.998, top-5 1.0) with material effects of up to 0.047 s⁻¹, but it cannot carry the preregistered verdict.

## 2. Why replication was required

The completed hardening review left one experiment blocking a claims paper: every CDW result came from the IEEE 39-bus system with D = 0, no governors and constant-power loads. The review asked two questions:
- whether the nested same-mode reversal (17/19 base-stable holdout policies on IEEE-39) and the portfolio-conditioned branch ranking survive a different network with documented primary-frequency dynamics;
- which of them survive a second converter model.

This campaign is that single test. No other network was searched.

## 3. IEEE-68 provenance

Full detail: `docs/CDW68_MODEL_AUDIT.md` and `results/CDW68_MODEL_MANIFEST.csv`.

- **Network** (Singh & Pal 2013, IEEE PES TF v3.3):
  - 68 buses, 83 branch rows, 16 of them transformers with taps; 100 MVA, 60 Hz.
  - Transcribed at Gate 3 (commit `d9fa097e`) and reproduced here: power flow within 5e-5 pu and 5e-5 deg; all 15 EM modes within 5e-4 Hz; Ybus rebuild error 0.0.
- **Machines:** 16 sub-transient machines.
  - Excitation: DC4B (G1–G8, G10–G12), ST1A (G9), manual (G13–G16).
  - PSS: speed washout plus three lead-lags on G1–G12.
  - Loads: constant impedance. Slack: G16.
- **Primary frequency** (not part of the published benchmark):
  - PST `data16m.m` (GridSTAGE copy, sha256 `c16e78f4…`) documents machine ratings (300–1900 MVA; equivalents 10–12 GVA), rotor damping d_o (zero on G1–G12, equal to H on G13–G16) and a PST tg model-1 governor block that is commented out in the file.
  - PSTess `tg.m` (sha256 `c1584021…`) gives the equations.
  - Frozen in `inputs/pst_primary_frequency_v1.json`.

## 4. Dynamic model realism

| variant | governors | rotor damping | base α⊥ (s⁻¹) | rightmost EM mode |
|---|---|---|---|---|
| REAL (confirmatory) | PST tg model 1 (1/R = 25 on rating, Ts 0.1, Tc 0.5, T3 0, T4 1.25, T5 5 s) on every surviving machine | PST d_o | −0.0669 (PSS washout, real) | −0.367 s⁻¹ at 0.555 Hz |
| NOGOV (ablation B) | none | PST d_o | −0.0461 (common-frequency mode, 0.002 Hz) | −0.371 at 0.518 Hz |
| SP33 (ablation C) | none | D = 0 (published) | −0.0556 | −0.118 at 0.520 Hz |

- **Common to all variants:** equilibrium residual 7.3e-13; identical network voltages (difference 1.1e-16); smallest documented-limit margin 3.00.
- **Damping effect.** The equivalents' documented damping moves the 0.52 Hz inter-area mode from −0.118 to −0.37 s⁻¹.
- **Consequence for the tests.** With it, the rightmost transverse eigenvalue is a slow real pole in every portfolio. This was established in R1 before any converter portfolio was evaluated, and it is the reason the EM-tracked level T was preregistered next to level D.

## 5. Converter models

- **Model A: the frozen TX4 GFL** (`69f200df`; device sha256 `bcf0c76b…`).
  - SRF-PLL 53/1400; power filter 0.03 s; outer PI 0.2/8; inner PI 0.25/6 behind x_f 0.15, r_f 0.01.
  - Leaky Q/V regulator with gain g (leak 0.05 rad/s).
  - Rating |S_gen|/0.8; matched dispatch.
  - No current limits, DC link or delays.
- **Model B: the WECC library chain TX3-GFL-0.1** (sha256 `f07a6a40…`).
  - Blocks: PLL2, BusFreq, REGCP1, REECB1, REPCA1; constant Q; forced-PQ equilibrium.
  - Runs in ANDES 2.0.0 on a new transcription SG68D/S/M of the Singh & Pal machine with the PST governor and damping (`code/andes68/cdw68_models.py`).
  - Not retuned.
- **Model B qualification** (`results/b68/B68_qualification*.json`), all criteria met before any B portfolio was run:

| check | result |
|---|---|
| Q1 | max \|Δα⊥\| against Model A's all-SG base: 9.9e-11 s⁻¹ over 16 k values, same status; EM top 3.1e-7 |
| Q2 | initialization residual ≤ 1.7e-13 |
| Q3 | descriptor QZ against reduced spectrum 7.7e-10 relative |
| Q4 | no active limiter |

- **Grid-forming model:** none is validated in the repository, so none was run.

## 6. Preregistration

Files: `docs/CDW68_PREREG_V1.md`, `docs/CDW68_STATISTICAL_PLAN.md` and `docs/CDW68_CLAIM_MATRIX_V1.csv`, committed at `6f188531` before any converter portfolio.

**Frozen:**
- the candidate set, the 16 policies, the draws and the TDS disturbance;
- τ = 0.01 s⁻¹;
- the level definitions A–D and T;
- the R11/R12 gates and the exhaustive CASE rule.

**Deviations** (`docs/CDW68_DEVIATIONS.md`, 15 entries):
- execution fixes (Model B check, JSON, ANDES shards, runtime registration, an R15 empty-store guard);
- two post-hoc analyses (EM-tracked branch ranking; exploratory TDS);
- a correction of an interpretation note (entry 11: NOGOV has three global marginals ≥ τ);
- layout-only figure fixes;
- paper and report production after the CASE decision (entries 13–14) and a runtime note (entry 15).

No gate, threshold, candidate, policy or model parameter was changed.

## 7. Candidate selection

**V68 = {G3, G4, G6, G9, G11, G12}.** These are the six physical plants with the largest scheduled active power (6.32–13.50 pu), with no tie. The rule extends the Gate 3 rule to both areas and uses no dynamic quantity.
- The area equivalents G13–G16 are excluded.
- Size: 2⁶ = 64 portfolios.
- Replacement: full, with matched dispatch.

The set was not changed after results.

## 8. Policy design

- **Model A: 16 maximin-LHS policies, P68_01–16** (seed 20260913, best of 2000 designs, minimum distance 0.1748).
  - Window: u ∈ [0, 1] with g = u², and k ∈ [0.5, 2.0].
  - P68_01–04 are discovery and are used only to select the static comparator. P68_05–16 are holdout.
- **Model B: 16 conditions B68_01–16**, carrying the k of the matching A policy.
- All 32 conditions have a stable all-SG base, so none was excluded. No policy search took place.

## 9. Portfolio census

Output: `results/CDW68_PORTFOLIOS.parquet` (A and B).

| set | portfolios | status | α⊥ range (s⁻¹) | α⊥ in the EM band |
|---|---|---|---|---|
| A REAL, 16 × 64 | 1024 | all STABLE | −0.0671 to −0.0581 | 0 % |
| A NOGOV, 6 × 64 | 384 | all STABLE | −0.0612 to −0.0461 | 0 % |
| A SP33, 6 × 64 | 384 | all STABLE | −0.0622 to −0.0539 | 0 % |
| B REAL, 16 × 64 | 1024 | all STABLE | −0.06695 to −0.06688 | 0 % |
| A draws, 40 × 64 | 2560 | all STABLE | max −0.0669 | 0 % |
| B draws, 20 × 64 | 1280 | all STABLE | max −0.0669 | 0 % |

- No portfolio is unstable, so H_RHP is empty and κ is undefined in every condition.
- **Rightmost EM-band mode:**
  - A REAL: 0.550–0.806 Hz (median 0.557 Hz), real part −0.373 to −0.263 s⁻¹;
  - B: 0.552–0.635 Hz, real part −0.369 to −0.251 s⁻¹.

## 10. Contextual reversal

**Primary result** (`results/CDW68_R07_A_gate.json`, `results/CDW68_R07_B_gate.json`).
- Level-D nested reversal coverage is **0/16 for Model A** and **0/16 for Model B** (95 % Clopper–Pearson upper bound 0.21 each).
- Levels A, B and C are also 0/16, at every τ from 0.01 to 0.05.
- The median-magnitude clause is undefined, because no level-D pair exists.

**Why levels A–D are empty.** The global marginals move a slow real pole.
- The largest level-A |Δα⊥| is 8.7e-3 s⁻¹ (A REAL, P68_16, replacing G12). For B it is 2.3e-5.
- No global marginal is material, so no sign class is non-zero.

**Preregistered secondary, level T (tracked rightmost EM mode):** 0/16 for both models.

| model | EM-tracked marginals | range (s⁻¹) | median | materially stabilizing (≤ −τ) | materially destabilizing (≥ τ) |
|---|---|---|---|---|---|
| A REAL | 2770 | −0.0069 to +0.0528 | +0.0020 | 0 | 437 |
| B REAL | 2804 | −0.0008 to +0.0570 | +0.0101 | 0 | 1429 |

Replacing a synchronous machine by either converter never materially stabilizes the tracked electromechanical mode on this network, in any context. A reversal needs a stabilizing side of at least τ, so none can exist. F2 and F3 show the one-sided pattern.

## 11. Discrete curvature

- **R8.** Chain-identity verification is vacuous: no nested reversal exists at any level.
- **Interaction terms still exist** (F4). On the rightmost EM mode, second differences have both signs beyond 1e-3 s⁻¹ in 13/16 A REAL policies and in 16/16 B policies. Their range is −0.0067 to +0.0457 s⁻¹ (A) and −0.0045 to +0.0562 s⁻¹ (B).
- **Reference policy P68_10:** 240 one-step squares with range −0.0052 to +0.0009 s⁻¹; 104 are below −1e-3 and none is above +1e-3. On the global α⊥ all |d| ≤ 2.9e-5.
- **Reading.** The set function is not modular, but its interactions do not change the sign of a replacement's effect: marginals shrink or grow with context and stay destabilizing.
- **Post-hoc TDS pairs** (§18). On their single-step chains the identity holds exactly (residual 0).

## 12. Fixed ranking

Source: `results/CDW68_R09_gate.json`.

| stratum | Model A | Model B |
|---|---|---|
| FULL | 57 decisions per policy (median), 0 material preferences: p* undefined, regret 0 | same |
| EM | 0 decisions: not evaluable | same |
| LEVEL-D | 0 decisions: not evaluable | same |
| TRACKED | 51.5 decisions, p* = 1.0, regret 0, transfer regret 0 | 46.5 decisions, p* = 1.0, regret 0, transfer regret 0 |

- The insufficiency bar (p* ≤ 0.90 and regret ≥ 0.10) is not met in any stratum.
- The T1 sign test cannot be computed.
- The headline "fixed ranking insufficient" is not allowed.
- On IEEE-68 the oracle fixed order of units reproduces every material preference on the tracked EM mode.

## 13. Total sensitivity validation (R11)

Source: `results/CDW68_R12_gate.json`.

| model | test | n | sign agreement | median relative error | Spearman | verdict |
|---|---|---|---|---|---|---|
| A | IFT total derivative (SPR-68) vs central re-equilibrated difference at h = 1e-4 | 80 | 1.000 | 1.3e-5 | 0.99995 | **PASS** |
| B | numerical total derivative, h = 1e-3 vs h = 1e-4 | 80 | 1.000 | 5.1e-7 | 0.99991 | **PASS** (numerical consistency only; B has no analytical engine) |

## 14. Reinforcement ranking (R12)

Holdout: 12 conditions per model, γ = 1.5.

| | Model A | Model B |
|---|---|---|
| median Spearman(D_tot, finite) | 0.9955 | 0.9810 |
| selected static (discovery) | S4 electrical distance, ρ 0.197 | S9 ΔgSCR, ρ 0.434 |
| median paired advantage (95 % bootstrap) | 0.799 (0.797–0.800) | 0.547 (0.546–0.547) |
| median top-5 precision | 1.0 | 0.8 |
| median Kendall / NDCG@5 | 0.958 / 1.000 | 0.896 / 0.998 |
| ρ at γ = 1.10 / 1.25 / 1.50 | 0.9993 / 0.9980 / 0.9955 | 0.9984 / 0.9939 / 0.9810 |
| T2 sign-flip p (Holm) | 0.0078 (0.0156) | 0.0898 (0.0898) |
| Wilcoxon p | 0.00024 | 0.00024 |
| gate | **PASS** | **PASS** |

- **T2B.** It does not reach 0.05. The preregistered statistic is the median under exact sign-flip enumeration. The twelve B advantages are nearly equal in size, and in that case the median statistic has little power. The gate itself is defined on medians and passes.
- **Materiality.** The median over conditions of the largest finite effect is 4.2e-5 s⁻¹ (A) and 3.3e-5 s⁻¹ (B). In the median condition no branch has a material effect.
- **Post-hoc EM-tracked ranking (A only).** Median ρ 0.998 and top-5 1.0. The median of the per-condition largest effect is 0.020 s⁻¹ (maximum 0.047). A median of 5 branches per condition exceed τ.
- **Top branches.**
  - EM-tracked, by median finite effect: lines 18–50, 51–50 and 18–49.
  - Preregistered target (both models): the G1 step-up transformer 1–54.

## 15. Static baselines

Discovery medians of Spearman with the finite ×1.5 effect (`results/CDW68_R10_static.csv`):

| index | A | B |
|---|---|---|
| S1 \|P\| | 0.043 | 0.202 |
| S2 \|S\| | 0.047 | 0.221 |
| S3 \|z\| | 0.021 | 0.230 |
| S4 electrical distance | **0.194** | 0.337 |
| S5 effective resistance | 0.116 | 0.294 |
| S6 Fiedler edge | 0.093 | 0.242 |
| S7 betweenness | 0.029 | 0.081 |
| S8 dV/dQ | 0.123 | −0.014 |
| S9 ΔgSCR | 0.067 | **0.434** |

- The static indices do not change across conditions.
- In the holdout, the best static index per condition has median ρ 0.197 (A) and 0.434 (B).
- The Laplacian scores (S5, S6) are baselines only; no topology explanation is offered.

## 16. Governor/damping ablation (R13)

Model A, P68_01–06, 64 portfolios each.

| variant | level D | level T | EM-tracked range (s⁻¹) | stabilizing ≥ τ | destabilizing ≥ τ | largest global \|Δα⊥\| |
|---|---|---|---|---|---|---|
| REAL | 0/6 | 0/6 | −0.0052 to +0.0393 | 0 | 267 | 1.5e-4 |
| NOGOV | 0/6 | 0/6 | −0.0102 to +0.0228 | 2 | 271 | 1.10e-2 (3 stabilizing, common-frequency mode, level B only) |
| SP33 | 0/6 | 0/6 | −0.0136 to +0.0274 | 6 | 277 | 4.7e-3 |

- Removing governors, or both governors and damping, makes a few EM-tracked marginals materially stabilizing. None of them forms a nested reversal.
- Reversal is absent in the realistic model and in both simplifications. On this network the documented primary-frequency dynamics neither create nor remove it.
- The simplification does shift the EM marginals slightly toward stabilizing values, which is the direction the IEEE-39 limitation paragraph anticipated.

## 17. Cross-converter transfer

- **What transfers between A and B:**
  - the absence of reversal (0/16 in both);
  - one-sided destabilizing EM effects (B's are larger: median +0.0101 against +0.0020);
  - an error-free fixed ranking on the tracked mode;
  - the ranking method (R12 PASS in both).
- **What does not transfer: the line list.** Across matched conditions, the Spearman between the two models' finite truths has median 0.50 (range 0.00–0.51), Kendall 0.37 and top-5 overlap 0.4. The 0.00 is P68_16, where Model A's critical pole is a different pole.
- **Which static index serves as comparator also depends on the converter:** electrical distance for A, ΔgSCR for B.

## 18. Nonlinear TDS

Source: `results/CDW68_R14_tds_posthoc.json`, figure F10.

**As preregistered.** No case exists: there is no level-D or level-T reversal and no level-D-eligible negative control. R14 is "not applicable".

**Post-hoc exploratory run.** The run took the two EM-tracked nested pairs with the largest same-sign marginals at the reference policy P68_10.
- Model A REAL; 0.5 pu reactor at bus 3 for 10 s; decay of the tracked modal coordinate fitted over 11–29 s.
- Integrator: BDF with a network solve at each step.

| pair | TDS decay change | linear decay change | sign agreement |
|---|---|---|---|
| G12, BASE → {11} | +0.011463 / +0.011467 | +0.011460 / +0.011466 | 2/2 |
| G4, {3} → {3,9} | +0.010275 / +0.010310 | +0.010273 / +0.010307 | 2/2 |

The traces confirm the linear signs and magnitudes to within 4e-6 s⁻¹. They say nothing about reversals, because none exists.

## 19. Uncertainty

Deterministic stress envelopes E05 (±5 %) and E10 (±10 %) at the reference policy (`results/CDW68_R15_uncertainty.json`).

| set | draws | level-D coverage | level-T coverage | draws with any materially stabilizing EM marginal | all 64 portfolios stable | smallest EM-tracked marginal |
|---|---|---|---|---|---|---|
| A E05 | 20 | 0/20 | 0/20 | 0/20 | 20/20 | −0.0055 |
| A E10 | 20 | 0/20 | 0/20 | 0/20 | 20/20 | −0.0056 |
| B E05 | 10 | 0/10 | 0/10 | 0/10 | 10/10 | −0.0009 |
| B E10 | 10 | 0/10 | 0/10 | 0/10 | 10/10 | −0.0014 |

- **Witness identity** is undefined: the nominal reference has no witness and no draw creates one.
- **Ranking under draws** (Model A; first 5 draws of each envelope, 83 branches, ×1.5): median ρ of D_tot 0.996 (10–90 % range over draws 0.9955–0.9960); post-hoc EM-tracked ρ 0.998 (0.9974–0.9981); median largest finite effect 4.3e-5 s⁻¹.

These are fractions of designed draws, not probabilities.

## 20. Cross-benchmark comparison

Source: `results/CDW_CROSS_BENCHMARK_MATRIX.csv`, figure F12.

| | IEEE-39 A | IEEE-39 B | IEEE-68 A | IEEE-68 B |
|---|---|---|---|---|
| governors / damping | none / D = 0 | none / D = 0 | PST tg / PST d_o | PST tg / PST d_o |
| candidates, policies | 9; 19 base-stable holdout | 4; 7 | 6; 16 | 6; 16 |
| level-D reversal | 17/19 (median 0.0127 s⁻¹) | global 0/7; EM-tracked 1/7 | 0/16 (T 0/16) | 0/16 (T 0/16) |
| fixed-ranking regret | FULL 0.38; EM 0.03 | not computed | FULL 0; TRACKED 0 | FULL 0; TRACKED 0 |
| total-sensitivity median ρ | 0.996 | 0.993 | 0.996 | 0.981 |
| static comparator median ρ | 0.445 | 0.790 | 0.197 | 0.434 |
| advantage | 0.49 | 0.20 | 0.80 | 0.55 |
| top-5 | 1.0 | n/a | 1.0 | 0.8 |
| TDS | not performed | not performed | post hoc, signs 4/4 | not applicable |

## 21. Negative results

All are preserved in the result files.
1. The primary reversal gate fails on both models: 0/16 and 0/16.
2. Level T fails on both models: 0/16 and 0/16.
3. No materially stabilizing EM-tracked marginal exists in REAL for either model, in the census or in any of the 60 draws.
4. Ablations: 0/6 reversals without governors and 0/6 without governors or damping.
5. The fixed-ranking failure of IEEE-39 does not replicate. Every IEEE-68 stratum has regret 0 or cannot be evaluated.
6. The preregistered branch-ranking target is immaterial: its effects are ≤ 1.4e-3 s⁻¹, and in the median condition no branch is material.
7. T2B is not significant (sign-flip p 0.090).
8. The line list does not transfer between converter models (median ρ 0.50).
9. The preregistered TDS holdout has no case.
10. Entry 7 of the deviation log was wrong for NOGOV and is corrected in entry 11.

## 22. Novelty implications

See `docs/CDW68_NOVELTY_ASSESSMENT.md`.
- Nested contextual reversal and fixed-ranking failure become **IEEE-39 benchmark-specific evidence**.
- The one finding carried across two benchmarks and two converter models is the portfolio-conditioned total branch sensitivity used as a ranking. It is a known construction (Smed 1993; Nam et al. 2000; Joswig-Jones et al. 2026, preprint), validated against finite reinforcements on held-out policies and compared with discovery-selected static indices.
- On IEEE-68 its preregistered target is immaterial. The material EM version is post hoc.
- The word "first" is not used.

## 23. Paper decision

**CASE B.** The CDW manuscript (`reports/papers/cdw_contextual_dynamic_weakness/`) is re-titled and its contributions reordered:
1. Portfolio-conditioned re-equilibrated reinforcement ranking on two benchmarks and two converter models.
2. When re-equilibration matters: on IEEE-68 the frozen derivative reaches ρ 0.55 against 0.996; total and frozen disagree in sign for a median 33 % of branches.
3. Contextual reversal on IEEE-39, with its failure to replicate on IEEE-68.
4. Boundaries: fixed rankings, cross-model transfer, and the immaterial IEEE-68 target.

**Readiness.** The manuscript is not yet ready for TPWRS as a claims paper. The IEEE-68 ranking pass rests on an immaterial pole. Its material EM version is post hoc and needs a preregistered confirmation, or at least a second EM-target test on a fresh set of conditions.

## 24. Reproducibility

- **Software:** Python 3.13; NumPy, SciPy and pandas from `.venv/tx3-analysis`; ANDES 2.0.0 in `.venv/xtool-andes-gfl`. Thread counts are fixed to 1.
- **Code** in `code/`:

| phase | scripts |
|---|---|
| inputs and audit | `R00_inputs.py`, `R01_audit.py`, `R01b_manifest.py` |
| Model B | `R02_export_b.py`, `andes68/run68.py` |
| census and reversal | `R06_census.py`, `R07_reversal.py` |
| fixed ranking | `R09_ranking.py` |
| branches and gates | `R10_branches.py`, `R10_static.py`, `R12_ranking_gate.py` |
| TDS, uncertainty, synthesis | `R14_tds.py`, `R15_uncertainty.py`, `R18_crossbench.py`, `R21_figures.py` |

- **Raw per-task JSON** in `raw/` is not committed. It is listed with sha256 in `CDW68_RAW_MANIFEST.csv` and packed in `CDW68_RAW_RESULTS_20260913.zip`.
- **Runtimes (wall clock):**

| phase | time |
|---|---|
| census, 1792 cases | 451 s |
| Model B census | 8 shards, 55 s |
| Model B branches | 6 shards, ≤ 919 s |
| Model B draws | 5 shards, ≤ 160 s |
| Model A branches, 128 tasks | 1862 s |
| Model A draws, 2560 cases | 703 s |
| ranking under draws, 40 tasks | 712 s |

- The campaign ran from the preregistration commit (09:19) through the last computation (R15L, finished 10:48:02) on 2026-09-13; the times of the document, paper and bundle commits are in git.
- **Figures:** F1–F12 in `figures/` (PDF, SVG, PNG).
