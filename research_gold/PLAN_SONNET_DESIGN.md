# Execution plan for Sonnet — poster redesign (phase 2)

Written 2026-10-05. Self-contained. Goal: rebuild the poster around a clear story with one new time-domain figure, claim-style titles and ~30 % less text. **Do not change any scientific claim or number**; every number must already exist in `research_gold/checks/*.csv` or be produced by Task 2 below.

Repo root: `C:\Users\walla\Documents\Github\paper-P06-spectral-graph-power-dynamics` (Windows; Bash = Git Bash; PowerShell available).

## 0. Starting point

- Current poster (keep it intact): `reports/poster/ias2026/one_mode_20261005/` → `output/pdf/IAS2026_Delay_Transport_20261005.pdf`.
- Evidence: `research_gold/04_primary_recommendation.md` (with Addenda 1–4 and the review note), scripts/CSV in `research_gold/checks/`.
- Previous Sonnet plan (done): `research_gold/PLAN_SONNET.md` — its sections 1 (hard rules) and 2 (tools) still apply. Read them first.

## 1. Hard rules (repeat)

1. **LaTeX only through the Write/Edit tools.** The Bash tool collapses `\\` and turns `\n`, `\b`, `\f` into control characters. In Python run from Bash, build backslashes with `chr(92)`.
2. Edit `bodies/NN.tex`, then `python splice_bodies.py`; never edit `sections/` directly.
3. Compile from the poster folder: `powershell -NoProfile -ExecutionPolicy Bypass -File compile_section.ps1 -Section all`. Then **look at the PNG** and crop panels with PIL to check for overflow. A clean log is not enough.
4. Keep the outer grid (`poster_layout.tex`) unchanged. Panel sizes: row 1 = 224 / 412.4 / 216 mm × 252 mm; strip 44 mm; panel 4 full width × 181 mm (`\PosterFourHeight`); panels 5–6 = 430.2 × 178 mm; panels 7–8 = 430.2 × 210 mm; panel 9 full × 84 mm.
5. No new claims. If something does not fit, cut words, not facts.
6. Commit only if the user asks, with explicit pathspecs.

## 2. Setup (Task 0)

```bash
cd reports/poster/ias2026
cp -r one_mode_20261005 one_mode_v2_20261006
rm -rf one_mode_v2_20261006/build
sed -i 's/IAS2026_Delay_Transport_20261005/IAS2026_Delay_Transport_v2_20261006/' one_mode_v2_20261006/compile_section.ps1
grep -c v2_20261006 one_mode_v2_20261006/compile_section.ps1   # must print 2
```
All work below happens in `one_mode_v2_20261006/`. Output: `output/pdf/IAS2026_Delay_Transport_v2_20261006.pdf/.png`.

## 3. Global design rules

- **Story order:** Problem → Surprise → Law → Why → In time → Spectrum ≠ events → Reach → Transfer and limits → Three rules.
- **Titles are claims** (exact titles given per block below).
- **One colour per strategy everywhere** (figures and text accents): no retune `#5b6b73` (grey), low-frequency rule `#C99A20` (gold), protect 4.91 Hz `#176494` (blue), protect 1.04 Hz `#003B2D` (green). Use the same legend names: "no retune", "low-frequency rule", "protect 4.91 Hz", "protect 1.04 Hz".
- **At most three big numbers on the whole poster:** "12 vs 0" (block 2), "0.99" (block 4), "5 of 5" (block 8).
- Body font: `\FigureSize` (24.7 pt) minimum; figure text ≥ 19 pt at final size (keep figsize ≈ 7.5–9.5 in wide so fonts stay large).
- Every block keeps at most one `InWords` box.

## 4. New material to produce

### Task 1 — Time-domain trajectories (the missing simulation)

Write `research_gold/checks/t15_trajectories.jl` (with the Write tool) based on `research_gold/checks/t9_events.jl`. Differences:
- After `met, r = DelayedEvents.adaptive(m; tau=tau, horizon=horizon, dtmax=.01, tol=1e-9)`, save the trajectory with
  `R.metrics(m, r; dt=.01, window=.5, monitor_buses=collect(1:39), savepath=joinpath(@__DIR__, "T15_" * label * "_trajectory.csv"))`
  (this call pattern is used in `experiments/interaction_decision_20261004/validate_julia.jl`).
- Arguments as in t9: `LABEL RULE NODE_RE NODE_IM H_MS RHO HORIZON EVENT_IDX`.

Run in background (from repo root), event 2 = bus 16 +100 MW, horizon 20 s, ρ = 0.875, delay 44 ms:
```
julia --project=. research_gold/checks/t15_trajectories.jl none44     none    0 1 4 0.875 20 2
julia --project=. research_gold/checks/t15_trajectories.jl lowfreq44  lowfreq 0 1 4 0.875 20 2
julia --project=. research_gold/checks/t15_trajectories.jl node104_44 exact  -16.34629153843205 6.543885486116494 4 0.875 20 2
julia --project=. research_gold/checks/t15_trajectories.jl node491_44 exact  -0.8015607744782098 30.850471693237324 4 0.875 20 2
julia --project=. research_gold/checks/t15_trajectories.jl none42     none    0 1 2 0.875 30 2
```
- `none44` and `lowfreq44` previously aborted with `DDE integration MaxIters`. If they abort again, the CSV may be missing: then plot whatever partial output exists. If none exists, rerun with horizon 8 s. As a last resort, show the abort time as a vertical "integration stopped" marker. Do not hide the failure.
- Inspect the trajectory CSV header first (`head -2`), then pick the frequency of one generator bus (prefer bus 30 or the bus with the largest deviation).

### Task 2 — Figures (add to `build_poster_figs.py`, same style)

1. **`time_traces.pdf/.png`** (block 5). Frequency deviation [Hz] vs time [s], 0–20 s, the four strategies at 44 ms in the strategy colours. Dashed ±0.5 Hz limits. If a run aborted, end its line and mark "integration stopped". Size ≈ 9.0 × 4.4 in.
2. **`slow_growth.pdf/.png`** (block 6, optional). `none42` trajectory over 30 s, showing the slow ≈ 4.8 Hz growth while the windowed limits are met. Size ≈ 6.5 × 3.6 in.
3. **`price_prediction.pdf/.png`** (block 4). Scatter from `research_gold/checks/T4_predict_which_mode.csv`: x = `predicted_worst`, y = `exact_worst`, one point per candidate protected mode, labelled with `protected_Hz`; diagonal y = x; text "rank correlation 0.99". Highlight 4.91 Hz (blue) and 1.04 Hz (green). Size ≈ 6.5 × 4.6 in.
4. Keep `central.pdf` (`build_central_fig.py`) and `share_delay_map.pdf` unchanged except for the colours if they differ from §3.

### Task 3 — Problem schematic (block 1, TikZ inside the body)

Horizontal chain drawn with TikZ (copy node styles from the loop diagram in `one_mode_20261005/bodies/03.tex` of the 20261004 version, or from `collective_interaction_20261004/bodies/03.tex`): box "SG" → box "inverter + PLL" with a small clock labelled "delay τ" on its input → cloud/box "IEEE-39 grid (10 sites)". A dashed return arrow from the grid to the PLL labelled "measured phase, τ late". Width ≤ 200 mm, height ≤ 70 mm.

## 5. Block contents (exact)

Write each body with the Write tool. Keep the outer `minipage`/panel wrapper exactly as in the current bodies (copy the first two and last two lines of the existing file).

**Header** (`sections/header_green.tex`, edit with Edit tool): eyebrow `BEYOND NODAL DAMPING`; title `PRESERVING ONE MODE DOES NOT PROTECT THE GRID`; subtitle `Exact PLL delay compensation and its network-wide price · IEEE-39`. If the title overflows at 90 pt, shorten it to `ONE MODE SAFE, GRID NOT`. Do not shrink the font.

**Block 1** — `bodies/01.tex`, title `A PLL THAT LISTENS LATE`
- Schematic (Task 3).
- Question, `\PosterParagraph` bold: "If a retune keeps one oscillation mode exactly where it was, is the grid safe?"
- Two lines: "IEEE-39, 87.5 % inverter share, ten PLLs. Their delay grows from 40 to 44 ms."
- No objectives list.

**Block 2** — `bodies/02.tex`, title `THE USUAL COMPENSATION MAKES IT WORSE`
- `central.pdf` as large as possible (≥ 270 mm wide).
- Side column: big number `12 vs 0` (`\PosterCondensed\bfseries 60pt`) with the caption "unstable roots at 44 ms: low-frequency rule vs protecting the right mode".
- InWords: "Which mode you protect decides everything. The low-frequency rule is worse than no retune; protecting the 1.04 Hz mode keeps the spectrum and passes 5 of 5 events."
- Remove the small table.

**Block 3** — `bodies/03.tex`, title `TWO GAINS PROTECT EXACTLY ONE MODE`
- Equation 1: transport `k_p(h)λ + k_I(h) = (k_pλ + k_I)e^{λh}`.
- Equation 2: remainder `Δκ(s) = −(s−λ)(s−λ̄)e^{−sτ}∫_0^h k_p(t)e^{−st}dt`.
- InWords: "Zero error at the protected mode, quadratic error everywhere else. No network model needed for the gains."
- Remove the derivation line.

**Strip** — keep as is.

**Block 4** — `bodies/04.tex`, title `WHY: THE NETWORK SETS THE PRICE`. Three columns:
- Left (~300 mm): spillover equation `dμ/dh = −k_p(μ−λ)(μ−λ̄) ∂μ/∂k_I`; one line "PLL: the quadratic factor. Network: the residue ∂μ/∂k_I (through Z(s))."
- Centre (~260 mm): `price_prediction.pdf` with big number `0.99`.
- Right (~260 mm): selection rule in 4 steps with `\StepBadge` (current text, shortened).
- Move "Three rules compared" and "Phase budget" out (to block 8).

**Block 5** — `bodies/05.tex`, title `IN TIME: WHAT THE GRID DOES`
- `time_traces.pdf`, full width of the panel.
- One line: "Bus 16 +100 MW at 44 ms. Delayed nonlinear DAE, method of steps."
- No table.

**Block 6** — `bodies/06.tex`, title `SPECTRUM AND EVENTS: PASS BOTH`. Two side-by-side cards:
- Card A, blue accent: "Protect 4.91 Hz: 0 roots beyond the margin, but 2 of 5 events (0.505 Hz; SG reserve ≈ 0). The retune spends integral gain." Source: `T5_nonlinear_events.csv`.
- Card B, grey accent: "No retune at 42 ms: 4 unstable roots, yet the event meets the limits for 30 s. 0.5 s windows do not see slow 5 Hz growth." Add `slow_growth.pdf` if produced. Source: `T5b_nonlinear_events_42ms.csv`, `T3_full_spectrum_count.csv`.

**Block 7** — `bodies/07.tex`, title `IT EXTENDS DELAY TOLERANCE, NOT THE MAXIMUM SHARE` — keep the current map and text; shorten to two sentences.

**Block 8** — `bodies/08.tex`, title `DOES IT TRANSFER? AND WHAT WE DO NOT CLAIM`. Three short items + one box:
- "Second design (85 %): rule applied unchanged → 0 roots, **5 of 5** events." (`T8_taylor_and_second_design.csv`, `T9_B_rule_node_44ms.csv`)
- "Nine-bus bank (assumed data): transport on three PLLs, −0.131 s⁻¹ (re-run); events reported."
- "First-order version of the formula gives the same counts: the choice of mode matters, not the truncation."
- InWords "Not claimed": optimal mode, maximum delay, maximum share, superiority over other compensators; root counts are floating point, not interval-certified.
- Add one line: "Phase budget h* = φ_PI(ω)/ω: 8.99 ms for the 4.91 Hz mode, a limit of the formula, not of operation."

**Block 9** — `bodies/09.tex`, title `THREE RULES FOR THE ENGINEER`. Three columns with `\StepBadge`:
1. "Do not compensate PLL delay at low frequency."
2. "Choose the protected mode with the spillover law: one modal analysis is enough."
3. "Check the spectrum and the events."

**Footer** — keep; ensure the takeaway reads "PRESERVING ONE MODE / DOES NOT PROTECT THE GRID." Reference [6]: confirm "F. Johansson, Arb: efficient arbitrary-precision midpoint-radius interval arithmetic, IEEE Trans. Computers 66(8), 2017" with one WebSearch; fix only if wrong.

## 6. QA checklist (Task 4)

- [ ] Compile has no fatal error; inspect the full PNG and crops of every panel.
- [ ] No text outside panels; no InWords box hidden under the next panel.
- [ ] Strategy colours identical in blocks 2, 5 and 6.
- [ ] Spot-check 6 numbers against CSVs: 12 vs 0 (`T3`/`T8`), 2 of 5 (`T5`), 5 of 5 at 85 % (`T9_B_rule_node_44ms`), 0.99 (`T4`), −0.131 (`T14` log), 0/2 at 90 % (`T9_M90_base_40`).
- [ ] The poster never says 48.9 ms is an operating limit, never claims novelty, never claims a maximum share.
- [ ] Read the poster top to bottom in 60 s following the story order; if a block needs more than two sentences to explain, cut text.

## 7. Report to the user (Spanish, short)

1. Whether the trajectories ran (and which aborted).
2. Final PNG via SendUserFile.
3. What was cut or could not fit.
4. Any number that did not match its CSV.

## 8. Budget

Task 1 (Julia, background) ~20–30 min wall time — start it first. Tasks 2–3 ~30 min. Block rewrite ~45 min. QA ~20 min.
