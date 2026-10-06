# Execution plan for Sonnet — poster v5 (diagram- and simulation-led), LaTeX

Written 2026-10-06. Self-contained. Goal: typeset the approved mock-up `output/pdf/MOCKUP_v5_36x48.png` (source: the mock-up script copied to `research_gold/mockup_v5.py`, see Task 0) as a clean LaTeX poster, 36 x 48 in. **No new science and no new numbers**: every number below already exists in `research_gold/checks/*.csv` or `research_gold/04_primary_recommendation.md`.

Repo root: `C:\Users\walla\Documents\Github\paper-P06-spectral-graph-power-dynamics` (Windows; Bash = Git Bash; PowerShell available).

## 0. Hard rules (read first)

1. **Write LaTeX only with the Write/Edit tools.** The Bash tool collapses `\\` and turns `\n`, `\b`, `\f` into control characters. In Python launched from Bash, build backslashes with `chr(92)`.
2. Do not modify earlier posters (`collective_interaction_*`, `one_mode_*`, `dense_v3_*`) or anything in `resources/`.
3. **No photo, no CV, no QR code** anywhere on the poster (user decision).
4. After every compile, open the PNG with the Read tool and crop panels with PIL to check for overflow. A clean log is not enough.
5. Keep all fonts ≥ 20 pt at final size (body ≥ 23 pt). Figures must be generated at sizes where their text stays ≥ 18 pt on the poster.
6. Do not copy anything from the user's reference photo except its visual language (hook band, block diagrams, small multiples, KPI strip, message banner). Its KPIs and references are not ours.
7. Commit only if the user asks, with explicit pathspecs.

## 1. Task 0 — Setup

```bash
cd reports/poster/ias2026
cp -r dense_v3_20261006 v5_20261006
rm -rf v5_20261006/build v5_20261006/panels/*
sed -i 's/IAS2026_Dense_v3_20261006/IAS2026_v5_20261006/' v5_20261006/compile_section.ps1
grep -c v5_20261006 v5_20261006/compile_section.ps1     # must print 2
cp "C:/Users/walla/AppData/Local/Temp/claude/C--Users-walla-Documents-Github-paper-P06-spectral-graph-power-dynamics/0f2abfe1-aeb7-4077-8780-c655a5d0ca0c/scratchpad/mockup_v5.py" ../../../research_gold/mockup_v5.py 2>/dev/null || echo "mock-up script not found; use the PNG as the reference"
```
The page size (91.44 × 121.92 cm) and header come from `poster_layout.tex`, `poster_common.tex`, `poster_design_system.tex` (the dense v3 design system already defines `DensePanel`, `DenseWords`, `DenseStep`, `\DenseBody`, `\DenseSmall`, `\DenseHead`, `\DenseEq`, `\DenseBig`). Reuse them. Header: `sections/header_white.tex` (logos) and `sections/header_green.tex`.

Compile: from `v5_20261006/`, `powershell -NoProfile -ExecutionPolicy Bypass -File compile_section.ps1 -Section full`. Output: `output/pdf/IAS2026_v5_20261006.pdf/.png`.

## 2. Figures (all already exist; regenerate only if a fix is needed)

Folder: `reports/poster/ias2026/one_mode_v2_20261006/generated/figures/` (PDF + PNG). Copy the ones used into `v5_20261006/generated/figures/`.

| Figure | Content | Script |
|---|---|---|
| `pll_loop` | block diagram: delay → PI → PLL → grid → phase error | `research_gold/checks/fig_v5_extra.py` |
| `splane4` | 4 s-planes, PLL cluster roots 40 ms (hollow) → 44 ms (filled), per strategy | same |
| `traces4` | 4 time traces, 5 Hz part of bus-30 frequency, bus 16 +100 MW, 44 ms | same |
| `price_prediction` | predicted vs exact worst root, 11 candidate modes (rank corr. 0.99) | `one_mode_v2_20261006/build_poster_figs.py` |
| `central` | rightmost root vs delay + event chips | `one_mode_v2_20261006/build_central_fig.py` |
| `share_delay_map` | share × delay map | `one_mode_v2_20261006/build_poster_figs.py` |
| `rank_two2` | δD_H eigenvalues per site | same |
| `decision` | pair 30–37 decision plot | same |

New figure to make (Task 2): **`h4_lattice`** — Hasse lattice of the 16 subsets of {30,33,35,37}: 15 blue nodes, the full set red, edges light grey, node size larger for the full set; title "15 stable / 1 unstable". Write it in `research_gold/checks/fig_v5_extra.py` (append) and save to the figures folder, figsize ≈ 6 × 4.5 in, no axes.

Strategy colours everywhere: no retune `#5b6b73`, low-frequency rule `#C99A20`, protect 4.91 Hz `#176494`, protect 1.04 Hz `#003B2D`; red `#B3363A` for unstable/fail.

## 3. Layout (mm, grid width 868.4, left margin 23, column gap 8, row gap 4)

Write a new `main.tex` (copy the header macro from `dense_v3_20261006/main.tex`). Below the header, stack these rows; each row is a `minipage` of the given height. Use `DensePanel` for numbered panels (it takes height, number, title). Navy title bars are optional; the existing green is fine.

| Row | Height | Content |
|---|---|---|
| A. Hook | 170 mm | full-width band titled **"A LOCAL RETUNE IS EXACT. THE WHOLE GRID PAYS ITS PRICE."** Three parts: (a) `pll_loop` + two lines "Each inverter sees the grid through a PLL whose measurement arrives τ late. Its only knobs: two gains."; (b) s-plane cartoon drawn in TikZ: a red vertical line (imaginary axis) with a pink "unstable" region on the right, one blue protected pole with a dashed ring labelled "protected λ: held exactly", four hollow grey poles with arrows (one arrow crossing into the red region), and the equation dμ/dh = −k_p(μ−λ)(μ−λ̄) ∂μ/∂k_I with the line "PLL: quadratic factor. Network: residue ∂μ/∂k_I"; (c) pale-green summary table "Same delay increase, 40 → 44 ms, ten PLLs" with columns rule / unstable roots / events: no retune 8 / aborted; low-frequency rule (usual) 12 / aborted; protect 4.91 Hz 0 / 2 of 5; protect 1.04 Hz 0 / 5 of 5; then in red bold: "The usual compensation is worse than none. Which mode you protect decides everything." |
| B. Law + price | 175 mm | two half-width panels. **2 THE EXACT LAW: TWO GAINS, ONE MODE**: transport equation (hero size), explicit K_p(h), K_I(h) formulas, the remainder Δκ(s), and a DenseWords box "No network model needed for the gains. Budget: h < φ_PI(ω)/ω (8.99 ms at 4.91 Hz)." **3 THE GRID SETS THE PRICE: CHOOSE THE MODE**: residue formula ∂μ/∂k_I = −ℓᴴZ_{k_I}(μ)r / ℓᴴZ_s(μ)r, `price_prediction`, the 4-step rule with `\DenseStep`, footnote "rank correlation 0.99: same order, not same values". |
| C. Simulation 1 | 180 mm | full width **4 SIMULATION 1 — WHERE THE ROOTS GO (IEEE-39, FULL DELAYED DAE, 204 STATES)**: `splane4` at full width; caption "open circles: 40 ms; filled: 44 ms; arrows: continuation of 8 catalogued PLL-cluster roots". |
| D. Simulation 2 | 140 mm | full width **5 SIMULATION 2 — WHAT THE GRID DOES IN TIME (BUS 16 +100 MW, 44 ms)**: `traces4` full width; two caption lines: "Unprotected, the 5 Hz oscillation grows until the delayed nonlinear integration diverges; protected, it dies out within 2 s." and "Full-spectrum count at 44 ms: 8 / 12 / 0 / 0 unstable roots. Events (5 load steps; limits on Δf, RoCoF, V, SG reserve): aborted / aborted / 2 of 5 / 5 of 5." |
| E. Four panels | 205 mm | four quarter-width panels: **6 PASS BOTH TESTS** — two cards (not the central figure): "Roots clean, events fail: protecting 4.91 Hz leaves 0 roots beyond the margin but only 2 of 5 events (0.505 Hz; SG reserve ≈ 0); k_I 246.7 → 137.8." and "Events clean, roots unstable: no retune at 42 ms has 4 unstable roots (+0.096 s⁻¹ at 4.8 Hz) yet meets the limits for 30 s; 0.5 s windows miss slow growth." **7 REACH: SHARE × DELAY** — `share_delay_map` + "Keeps the margin at 44/48 ms wherever 40 ms is feasible; no higher maximum share (90 % fails events even at 40 ms). Second design (85 %), rule unchanged: 5 of 5. Nine-bus bank: +0.886 → −0.131 s⁻¹ (re-run)." **8 BEYOND NODAL DAMPING** — D_H = ½(Z+Zᴴ), δD_H = ½(avᴴ+vaᴴ), `rank_two2`, "Never pure damping: one + and one − direction (60/60). Counting law: r bad directions need ≥ r retuned PLLs. Nodal +38.7 vs collective −737.5." **9 SAFE ALONE ≠ TOGETHER** — `h4_lattice`, det(I+Q_H) = det(I+Q_RR) det(I−R_{i|R}), "Pair retunes too: buses 30 & 37 hold −0.0657 alone, −0.0414 together; interaction +0.0243 ± 4×10⁻⁷ (interval-certified), repaired to −0.060 in 11 steps. Blind H4 prediction 16/16 (4 buses), 395/512 (9 buses)." |
| F. KPI strip | 62 mm | title bar "EVIDENCE AT A GLANCE"; six cells with a big number and two-line label: **12 → 0** (red) unstable roots at 44 ms, usual rule → our rule; **5 / 5** delayed nonlinear events passed; **0.99** rank correlation of the predicted price; **3.9×10⁻⁵** max error of the spillover identity, 219 mode pairs; **485** interval panels certifying one root cluster; **204** states in the full delayed IEEE-39 model. Thin vertical separators. |
| G. Banner | 26 mm | dark band, centred white bold: **CHOOSE THE MODE. COMPUTE THE GAINS. CHECK THE WHOLE GRID.** |
| H. Bottom | 105 mm | three panels: **10 THREE RULES FOR THE ENGINEER** (3 numbered rules + grey line "Not claimed: maximum share, optimal mode, superiority over other compensators, novelty of the algebra."); **11 REFERENCES** (the six references in `dense_v3_20261006/sections/footer.tex`); **12 CONTACT** (text only: name, jl.mayorga@uniandes.edu.co, Universidad de los Andes, Bogotá, Colombia, IEEE IAS Annual Meeting 2026, Vancouver). |

The heights sum to about 1063 mm plus header (140) and gaps; adjust row heights so the poster fills the page exactly, leaving no large empty band at the bottom. Remove the old green footer (`sections/footer.tex`) from `main.tex`: block H replaces it.

## 4. QA checklist

- [ ] Compile has no fatal error; inspect the full PNG and crops of every row.
- [ ] No text outside panels; nothing hidden under the next row.
- [ ] Strategy colours identical in the hook table, Simulation 1, Simulation 2, panel 6.
- [ ] Spot-check numbers against CSVs: `T3_full_spectrum_count.csv` (8/12/0/0 at 44 ms, as unstable roots; note T3 lists `protect_slow` = low-frequency rule), `T5_nonlinear_events.csv` (2 of 5), `T5c` / `T9_A_node1Hz_*` (5 of 5), `T4_predict_which_mode.csv` (Spearman 0.991), `T1_spillover_ieee39.csv` (max rel. error of well-scaled pairs ≈ 3.9e-5), `T9_M90_base_40.csv` (0 of 2), `T9_B_rule_node_44ms.csv` (5 of 5).
- [ ] No photo, no QR, no CV text.
- [ ] Never claims: novelty, a maximum share, an optimal mode, 48.9 ms as an operating limit.

## 5. Report to the user (Spanish, short)

1. Final PNG via SendUserFile.
2. What changed from the mock-up and why.
3. Any number that did not match its CSV.
4. Anything left unresolved.

## 6. Budget

Setup 5 min; H4 lattice figure 10 min; hook TikZ 30 min; rows B–H 60 min; QA 20 min. No Julia runs are needed.
