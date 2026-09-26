# Phase I — unit correction (MW vs MVA): data contract, audit and claim impact

Status: complete for the IAS2026 research tree, 2026-09-10. Branch
`ias2026/nonlinear-portfolio-v2`. Every correction is **append-only**. No
frozen table, figure, manifest or document (FREEZE_SHA256SUMS.txt,
JOURNAL_GATES_SHA256SUMS.txt) was rewritten. Corrected outputs are in
`results/UC/` (manifests in `outputs/unit_correction_v1/`). The superseding
vocabulary is `docs/POSTER_FINAL_SAFE_CLAIMS_v2.md`.

## 0. The error

`ReplacementCase.replaced_mw` (`src/ibr_cycles/models/ieee39_case.py:559`) returns

    sum over converter slots of  weight * 100 MVA,   weight = rho * Sn / 100,

that is, the **apparent-power rating** of the replacement converters. By the
frozen rating rule this equals `sum(rho_b * Sn_b)` of the retired machines. It
is **MVA, not MW**.

The value is read only after the case is solved. It is never an input to the
DAE, the Jacobian, the eigen-solver, the port operator or any certification
code (code-path audit, §2). So the spectral results cannot depend on it.
Everything that *used* it as megawatts can:

- predictors;
- matched-control selection;
- planning cost axes;
- the wording of claims.

## 1. The data contract (new; `src/ibr_cycles/units/`)

| canonical quantity | definition | unit | provenance |
|---|---|---|---|
| `replaced_sn_mva` | converter apparent-power rating as built (= legacy `replaced_mw`, to 1.8e-12) | MVA | slot weights |
| `retired_sg_sn_mva` | `sum rho_b Sn_b` of the retired machine share; NaN for IEEE-68, where `Sn` is a per-unit base | MVA | machine data |
| `replaced_pg_mw` | active power carried by the replacement converters at the solved equilibrium, **measured** from the device injections `Re(V conj(I))` | MW | equilibrium |
| `replaced_q_mvar` | reactive power carried by the converters | Mvar | equilibrium |
| `replaced_pmax_mw` | `sum rho_b Pmax_b` of the retired machines, **only where the source documents Pmax**; the synchronous machine's active limit, **not** a PV nameplate | MW | IEEE-39: ANDES `PV/Slack.pmax` of `ieee39_full.xlsx` (sha `9c2048dc…`); Kundur: placeholders (qmax 99 pu) → NaN; IEEE-68: none → NaN |
| `sg_sn_mva_left_at_replaced_buses` | machine share left at the replaced buses plus condenser ratings | MVA | plan |
| `condenser_sn_mva` | condenser rating added | MVA | plan |

Rules:

- No function converts MVA to MW.
- A **PV nameplate MW does not exist** in these benchmarks: the converter holds
  `P_ref` = the displaced dispatch, with an unlimited DC side. Therefore:
  - "PV capacity retained = X MW" is **never** written;
  - the permitted phrases are *converter active dispatch retained [MW]* and
    *replacement apparent-power rating [MVA]*.
- The legacy property is kept, unchanged, for the reproducibility of frozen
  results. New code must call `ibr_cycles.units.measure(case)`.

Tests (`tests/test_unit_contract.py`, 4/4):

- flagship: Sn 4270.7 MVA, Pg 2096.61 MW measured = base dispatch, Q 420.65 Mvar,
  Pmax 2983 MW;
- a 50 % replacement carries `rho Pg`;
- a condenser rating is a separate axis (the converter keeps all P, the
  condenser takes Q);
- Pmax is NaN where it is not documented.

Check over the 511 E12 census portfolios (UC01):

| check | result |
|---|---|
| `replaced_sn_mva` against the legacy column | 1.8e-12 MVA |
| measured `replaced_pg_mw` against the base-dispatch sum | 3.3e-11 MW |
| legacy/Pg ratio | 1.16 to 3.02 |
| Spearman correlation, legacy vs Pg | 0.964 |

## 2. Audit of every historical use

Search terms: `replaced_mw`, `MW`, `megawatt`, `4270.7`, `1040`, `970.2`,
`PV retained/lost/kept/forgone`, `removed/replacement MW`, `pv_mw`,
`MW-matched`.

Scope: the whole `reports/poster/ias2026` tree, including the poster sources.
`sections/*.tex` and `tikz/*.tex` contain no MW, MVA, AUC or PV-capacity claim,
so **the poster itself is unaffected**.

Legend:

- **Δnum**: the number changes.
- **W**: wording only.
- **Rerun**: needed a recomputation. "done" means it was redone in UC02/UC03.

| # | file | claim / table / figure | old quantity | old unit (as labelled → actual) | actual physical quantity | corrected value / unit | Δnum | W | Rerun |
|---|---|---|---|---|---|---|---|---|---|
| 1 | `src/ibr_cycles/models/ieee39_case.py:559` | property `replaced_mw` | Σ converter weight×100 | MW → **MVA** | converter rating | kept as legacy; superseded by `units.measure` | no | yes | no |
| 2 | `experiments/E11_replacement_atlas.py`, `results/tables/E11_replacement_atlas_{atlas,limits}.csv` | `replaced_mw`, `mw_replaceable = rho* × rating` | Sn-based | MW → MVA | rating replaceable | `UC03_E11_limits_reexpressed.csv`: `pg_mw_replaceable = rho* Pg` (e.g. bus 30: 1040 MVA / 436.1 MW; bus 38: 1684.1 MVA / 764.8 MW) | yes | yes | done |
| 3 | `experiments/E12_compatibility_census.py`, `E12_*.csv`, `results/manifests/E12_*.json` | census column `replaced_mw` | Sn sum | MW → MVA | rating | `UC01_census_quantities.csv` adds Pg, Q, Pmax per portfolio | column only | yes | done |
| 4 | `experiments/E13_baseline_challenge.py`, `E13_baseline_challenge_{auc,predictors,screening}.csv`, `results/manifests/E13_*.json`, `docs/PHASE_A_D_REPORT.md:127,141,163` | "replaced MW" predictor: AUC 0.813/0.773/0.861; screening at size 4: precision 0.300, 7 FN | Sn sum | MW → MVA | rating | Pg: AUC **0.572/0.636/0.727**; screening precision **0.100**, 9 FN. The Sn predictor keeps 0.813/0.773/0.861 as "retired rating [MVA]" | **yes** | yes | done (UC02) |
| 5 | `experiments/E14_negative_controls.py`, `E14_negative_controls_{all,N4_matched}.csv`, `docs/E14_NEGATIVE_CONTROLS.md:14-16`, CLAIMS O15, discovery_log row 23 | N4: 25 stable controls "matched on MW" 4665–5025 vs failing 4144–4985 | Sn sum | MW → MVA | rating | Pg matching: failing 2096.6–2683.8 MW; **107** stable in the ±5 % window; top 25 at 2533.9–2787.8 MW, **all stable**, worst abscissa **−0.1145**; min SCR 0.907–2.95 (no longer identical); 9/25 overlap with the frozen set | **yes** | yes | done (UC02) |
| 6 | `experiments/E14_negative_controls.py` N6 | MAC branch tracking into the 25 matched controls | set chosen on Sn | — | — | not recomputed for the Pg-matched set | — | — | **open** (non-blocking) |
| 7 | `experiments/E18_holonomy.py`, `E18_holonomy_matched_contrast.csv`, `docs/PHASE_E_REPORT.md:290-291`, CLAIMS O25/O26, discovery_log row 34 | 12 stable controls "matched on megawatts" (4176–4366); closure distance 0.142–0.741 vs flagship 1.58e-7; cycle ratio 0.16 | Sn gap | MW → MVA | rating | **recomputed** with Pg matching (2055.5–2162.5 MW): closure distance **0.156–0.863** vs 1.58e-7 (ratio 9.9e5); longest-cycle ratio **0.72** (< 1: the flagship cycle is still not larger than a control); 1/12 overlap with the frozen set | **yes** | yes | done (UC02) |
| 8 | `experiments/E21_minimal_repair.py`, `E21_minimal_repair_strategies.csv`, `results/manifests/E21_*.json`, `docs/PHASE_E_REPORT.md:316-337`, CLAIMS O27, `configs/ias2026/trackA_final_v1.yaml:61-80`, discovery_log row 35 | "PV kept / forgone" 4270.7 / 3230.7 / 3300.5; "costs 970–1040 MW of PV" | Sn | MW → MVA | converter rating kept or forgone | converter dispatch kept: flagship **2096.6 MW**; restore 30: 1660.5 (forgone 436.1); 33: 1444.6 (652.0); 35: 1409.6 (687.0); 37: 1775.1 (321.5); `UC03_E21_reexpressed.csv` | **yes** | yes | done (UC03) |
| 9 | `experiments/E23_closure_and_frontier.py`, `E23_closure_and_frontier_frontier.csv`, `results/manifests/E23_*.json`, CLAIMS O30 | `retune_pv_retained` 4270.7; `restore_pv_retained` 3300.5 | Sn | MW → MVA | rating | 2096.6 MW / 1775.1 MW; restore choice (bus 37) **unchanged** by the min-Pg rule | yes | yes | done (UC03) |
| 10 | `experiments/E34_mitigation_frontier.py`, E34 CSVs, `E34_MITIGATION_FRONTIER.md:15-45`, FINAL_TRACK_A C9, discovery_log row 80 | cost axis `pv_mw_forgone` (Sn × kept fraction); Pareto front | Sn | MW → MVA | rating | Pg axis measured on 23 plans: front membership **identical** (10 points); M5 costs **65.4 MW** (156 MVA) and **109.0 MW** (260 MVA); M4 still never efficient; zero-forgone strategies unchanged (M1–M3) | yes | yes | done (UC03) |
| 11 | `experiments/E38_bus30_collective_role.py`, `E38_*.csv`, `E38_bus30_collective_role.md:18,41-54`, CLAIMS N15, FINAL_TRACK_A C12, `FAILED_FINAL_VALIDATION.md:75-80` | "megawatts replaced alone 1040, 4th"; logistic adjusted OR 4.09 [0.14, 124.2], p 0.418 | Sn | MW → MVA | rating | bus 30 Pg **436.1 MW, 8th of 9**; adjusted for Pg: OR **0.053 [0.0007, 3.75], p 0.176**. Verdict (not supported) unchanged | **yes** | yes | done (UC02) |
| 12 | `experiments/E39_baseline_audit.py`, `E39_baseline_audit.csv`, `E39_figure_source.csv`, **`E39_ROC_PR.png`**, `E39_BASELINE_AUDIT.md:50,61-63`, CLAIMS N16, `FAILED_FINAL_VALIDATION.md:87`, `POSTER_FINAL_SAFE_CLAIMS.md:80,108-109` (S10), discovery_log row 89 | "replaced MW 0.81/0.77/0.86, good"; "megawatts are informative rankers" | Sn | MW → MVA | rating | Pg: 0.57/0.64/0.73, size-4 permutation **p = 0.46**; new figure `results/UC/UC02/UC02_ROC_PR_unit_corrected.png` | **yes** | yes | done (UC02) |
| 13 | `docs/PHASE_E_REPORT.md:214,327`, `docs/F7_POLICY_HYPERGRAPH.md:214`, CLAIMS O66, `trackA_mitigation_definitions.yaml:25`, `overnight_policies.yaml:39` | "keeps every megawatt of PV", "no change of PV megawatts", "converter carries every displaced megawatt" | dispatch statements | MW (correct) | active dispatch | unchanged: true in MW | no | no | no |
| 14 | `docs/PHASE_A_D_REPORT.md:47-56`, `E10_*`, `E20_*` | "rating MVA" columns | Sn | MVA (correct) | rating | unchanged | no | no | no |
| 15 | `POSTER_FINAL_SAFE_CLAIMS.md` S8, FINAL_TRACK_A C10, E34 | condenser "166–270 MVA" | condenser rating | MVA (correct) | rating | unchanged | no | no | no |
| 16 | BC00/BC Phase I docs, NL docs (this programme) | already use MVA/MW correctly | — | — | — | — | no | no | no |

## 3. Flagship physical numbers (UC01, measured)

| bus | Pg [MW] | Pmax [MW] (documented) | Qg [Mvar] | Sn [MVA] |
|---|---|---|---|---|
| 30 | 436.09 | 1040 | 92.70 | 1040.0 |
| 33 | 652.00 | 682 | 130.04 | 1174.8 |
| 35 | 687.00 | 697 | 225.53 | 1085.7 |
| 37 | 321.52 | 564 | −27.62 | 970.2 |
| **total** | **2096.61** | **2983** | **420.65** | **4270.7** |

**Correct statement.** The four replacements displace 2096.6 MW of
pre-replacement active dispatch (and 420.6 Mvar). This dispatch came from
synchronous machines rated 4270.7 MVA in apparent power, with a documented
active limit of 2983 MW. The replacement converters are rated 4270.7 MVA and
carry exactly that dispatch. No PV nameplate is modelled.

| case | converter Pg kept [MW] | converter Sn [MVA] | converter Q [Mvar] | synchronous MVA kept or added | statement |
|---|---|---|---|---|---|
| flagship unrepaired | 2096.6 | 4270.7 | 420.6 | 0 | as above |
| restore SG30 | 1660.5 | 3230.7 | 328.0 | 1040 (machine kept, Pmax 1040 MW, dispatching 436.1 MW) | restoring bus 30 keeps a 1040 MVA machine and gives up 436.1 MW of converter dispatch |
| restore SG33 | 1444.6 | 3095.9 | 290.6 | 1174.8 (652.0 MW) | — |
| restore SG35 | 1409.6 | 3185.0 | 195.1 | 1085.7 (687.0 MW) | — |
| restore SG37 | 1775.1 | 3300.5 | 448.3 | 970.2 (Pmax 564 MW, 321.5 MW) | the cheapest restoration in both MVA and MW |
| RB / RC converter retune | 2096.6 | 4270.7 | 420.6 | 0 | all four converters kept; controller gains change no power quantity (solved with the retuned controllers) |
| RD condenser 0.25 Sn at each retired bus | 2096.6 | 4270.7 | 0 (condensers take Q) | 1067.7 condenser | — |
| E34 condenser at bus 30 only, 0.16–0.26 Sn | 2096.6 | 4270.7 | 328.0 | 166.4–270.4 condenser | — |
| E34 M5 keep 0.15 / 0.25 of machine 30 | 2031.2 / 1987.6 | 4114.7 / 4010.7 | — | 156 / 260 (machine share, 65.4 / 109.0 MW) | — |

## 4. The E39 baseline audit with true active power (UC02)

Pipeline check: re-running the frozen E39 loop (same predictors, sizes, 5000
permutations, seed 20260917) reproduces the frozen CSV to **1.1e-16**. New
predictors used a fresh seed (20260918).

Raw-value ROC AUC; permutation p in brackets:

| predictor (expected: higher = more unstable) | size 4 (10/126 unstable) | size 5 (51/126) | size 6 (78/84) |
|---|---|---|---|
| retired rating Sn [MVA] (old "replaced MW") | 0.813 [0.0006] | 0.773 [<2e-4] | 0.861 [0.0014] |
| **removed dispatch Pg [MW]** | **0.572 [0.455]** | **0.636 [0.009]** | **0.727 [0.052]** |
| retired Pmax [MW], documented | 0.874 [<2e-4] | 0.723 [<2e-4] | 0.789 [0.018] |
| removed Q [Mvar] | 0.361 [0.151] | 0.430 [0.184] | 0.436 [0.612] |
| removed inertia fraction (unchanged) | 0.869 | 0.785 | 0.846 |
| min SCR, raw (low = risky) | 0.395 | 0.347 | 0.322 |
| additive / pairwise reconstruction (unchanged) | 0.144 / 0.153 | 0.146 / 0.147 | 0.090 / 0.105 |

Sign orientation: every expectation is "higher is more unstable". Pg and Pmax
keep that sign at every size (AUC > 0.5). Q is weakly inverted (0.36–0.44,
not significant).

PR AUC (base rates 0.08 / 0.40 / 0.93):

| predictor | size 4 | size 5 | size 6 |
|---|---|---|---|
| Sn | 0.328 | 0.747 | 0.988 |
| Pg | 0.104 | 0.567 | 0.972 |
| Pmax | 0.365 | 0.654 | 0.982 |

Precision / recall at the first `k` portfolios ranked:

| predictor | k = 10, size 4 | k = 10, size 5 | k = #unstable, size 5 | k = #unstable, size 6 |
|---|---|---|---|---|
| Sn | 0.30 / 0.30 | 1.00 | 0.61 | 0.96 |
| Pg | **0.10 / 0.10** | 0.70 | 0.51 | 0.95 |
| Pmax | 0.30 / 0.30 | 0.70 | 0.61 | 0.94 |

At size 4, k = 10 is also the number unstable, so the size-4 values double as
precision/recall at k = #unstable.

Why the ranking moves: the old column was a machine-**size** measure. Its
Spearman correlation with the removed inertia fraction is 0.84 at sizes 4–5
(Pmax 0.73), against **0.44** for dispatch. Removed inertia is `M = 2H Sn`.

Matched-portfolio analyses redone on Pg:

- **E14 N4.**
  - 107 stable four-replacement portfolios lie inside the failing Pg range
    ±5 %, and 106 of them also inside the failing min-SCR range.
  - The top 25 by Pg are all stable, worst −0.1145.
  - For **every** one of the 10 failing portfolios there is a distinct stable
    portfolio with the same removed dispatch within 5 %; the median gap is
    0.007 %.
- **E18.**
  - Controls matched on Pg: closure distance 0.156–0.863 against 1.58e-7 for
    the flagship.
  - Longest-cycle ratio 0.72: the flagship's cycle is not larger than the
    controls'.
- **E38.**
  - Adjusted for Pg instead of Sn, the bus-30 odds ratio moves from 4.09 to
    0.053 [0.0007, 3.75], p = 0.176.
  - That is still not significant, and the interval still spans more than three
    orders of magnitude.

**Answers.**

1. *Does true removed MW remain a useful conventional baseline?*
   - **No, not at the minimum failing order.** At size 4, removed active
     dispatch is indistinguishable from chance: AUC 0.57, permutation p = 0.46,
     precision@10 = 0.10. At sizes 5–6 it is moderate: 0.64 (p = 0.009) and
     0.73 (p = 0.052).
   - The informative "size" baselines are the retired **rating** (MVA), the
     documented **capability** Pmax (MW) and the removed **inertia**, all at
     0.72–0.87.
2. *Does the previous conclusion about nodal/simple metrics weaken, strengthen
   or stay qualitatively unchanged?*
   - **Qualitatively unchanged in its core, narrowed in one clause.** These
     statements still hold:
     - the short-circuit family is weak;
     - the lower-order reconstructions never fire (0 TP at every size);
     - no conventional measure identifies the irreducible fourth-order
       structure;
     - the matched contrasts (E14, E18) survive with Pg matching.
   - What changes: the clause "removed inertia **and megawatts** are
     informative rankers" is **invalid for megawatts of dispatch**. It becomes
     "removed inertia and retired machine rating/capability are informative
     rankers; removed active dispatch is weak at the minimum failing order".
   - C13's own numbers (inertia 0.78–0.87; reconstructions 0 TP) are
     unchanged. Its companion wording (N16, F15, S10) is corrected in
     `POSTER_FINAL_SAFE_CLAIMS_v2.md` and §6 below.

## 5. Planning and mitigation (UC03)

Every optimisation or table that preserved or minimised "PV MW" using `Sn` is
**INVALID as an active-power planning result**:

- E21 `pv_mw_*`;
- E23 `*_pv_retained`;
- E34 `pv_mw_*`, the "PV MW lost" axis;
- `trackA_final_v1.yaml` `pv_mw_*`.

Each was re-scored on four separate axes: converter Pg [MW], converter Sn [MVA],
synchronous MVA added, and controller change. The stability outcomes were not
re-optimised: a relabel of cost cannot change an abscissa. Results:

- **E34.**
  - The Pareto front recomputed with the frozen axes reproduces the frozen 10
    points exactly.
  - On the Pg axis the membership is **identical**. Only the M5 costs change:
    156 → 65.4 MW and 260 → 109.0 MW.
  - M4 (restore a whole machine) is still never efficient.
  - M1–M3 still forgo 0 MW.
  - S8 and C9/C10 stand as worded ("no photovoltaic energy forgone" is true in
    MW).
- **E23.** The restoration chosen at every target is bus 37, both by minimum
  MVA and by minimum MW forgone. The retained converter dispatch is 1775.1 MW,
  not "3300.5 MW".
- **E21.** Restoring 37 gives up the least dispatch (321.5 MW), restoring 30
  gives up 436.1 MW. The ordering is the same as on the Sn axis (970.2 < 1040
  MVA); the numbers are not.
- **E11.** `mw_replaceable` is a rating (MVA). Replaceable dispatch is
  `rho* Pg`; `rho* = 1` at every candidate.

No active-power optimisation in the tree changes its selected option. All
reported "MW of PV" figures change value.

## 6. Claim impact classification

| claim | class | corrected form |
|---|---|---|
| C1–C8, C10, C11; S1–S7, S9; spectral counts, `kappa`, `H` (F2C, F7, G1, BC01), policy regions (F2B, F7, F11), port closure (E16–E18 flagship values, E33, E41), F7 boundaries, symmetry (BC00, NL00), nonlinear DAE identities (NL01, NL02), O66 | **UNAFFECTED** | Verified by code path: `replaced_mw` is a read-only, post-solve property with no consumer in `src/` models, solvers, ports, certification, or any F- or G-experiment. Only E11–E14, E18, E21, E23, E34, E38, E39 and BC00 read it, all listed above |
| O27, O30, C9, E21/E23/E34 reports, `trackA_final_v1.yaml`, discovery_log rows 35 and 80, E11/E12 column names, PHASE_A_D line 163 header | **WORDING_CORRECTION** | "keeping all four converters: 4270.7 MVA of converter rating carrying 2096.6 MW of active dispatch" |
| E13 screening (PHASE_A_D lines 127, 141); O15 / E14 N4; O25 / O26 / E18; C12 / N15 / E38 (numbers); E21 / E23 / E34 "PV forgone" values; PHASE_E line 337 ("970–1040 MW of PV" → 321.5–436.1 MW of dispatch, 970.2–1040 MVA of rating) | **NUMERICALLY_CHANGED** (rerun done, conclusion unchanged) | values in §2–§5 |
| E14 N6 (branch MAC into the matched controls) | **REQUIRES_RERUN** (open, non-blocking: N4's stability outcome was re-established on Pg) | — |
| "Removed inertia **and megawatts** are informative rankers" (N16, F15, S10, E39 report); "megawatts replaced alone 1040, 4th" (E38); the "PV MW" cost axes as MW (E21, E23, E34); "matched on MW / megawatts" labels (E14, E18) | **INVALIDATED** as MW statements, superseded | §4–§5; `POSTER_FINAL_SAFE_CLAIMS_v2.md` |
| C13 | **UNAFFECTED** in its numbers, **companion wording corrected** | "Lower-order reconstructions: 0 true positives at every size. Removed inertia (AUC 0.78–0.87) and retired machine rating (0.77–0.86, MVA) rank instability; removed active dispatch does not at the minimum failing order (0.57, p = 0.46); the short-circuit family is weak" |

## 7. What this does not change

- The fourth-order phenomenon, its dispatch caveat (Q1) and its machine-data
  caveat (Q2).
- The port / closure mathematics.
- The certification programme's H, kappa, symmetry and zero ledger.

The unit error touched only the quantities used to **describe, match, rank or
price** portfolios, never the quantities used to **decide stability**.
