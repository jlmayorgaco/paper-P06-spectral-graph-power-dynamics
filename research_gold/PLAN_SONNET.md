# Execution plan for Sonnet — close the IAS 2026 poster

Written 2026-10-05. Self-contained: a new session can execute it without re-reading the conversation.
Repo root: `C:\Users\walla\Documents\Github\paper-P06-spectral-graph-power-dynamics` (Windows; Bash = Git Bash; PowerShell available).

## 0. Context in five lines

- **Poster to finish:** `reports/poster/ias2026/one_mode_20261005/` → `output/pdf/IAS2026_Delay_Transport_20261005.pdf` (+ `.png`). Earlier posters (`collective_interaction_20261003`, `collective_interaction_20261004`) must stay untouched.
- **Headline (decided, do not change):** exact delay transport of PLL gains, network-mediated modal spillover, choice of protected mode. Motto: "Preserving one mode does not protect the grid."
- **Evidence:** `research_gold/04_primary_recommendation.md` (and its three addenda), `research_gold/02_claims_evidence_matrix.md`, scripts and CSV in `research_gold/checks/`.
- **Rejected for the poster:** H4 / portfolios, maximum replacement, GSP / commutator, cooperative K_p-only compensation (negative on IEEE-39, Addendum 3), "stronger grid less stable".
- **Memory:** `gold-result-one-pll-one-mode` in the auto-memory directory.

## 1. Hard rules (read before touching anything)

1. **Never write LaTeX through a Bash heredoc or inline Python string**: the Bash tool collapses `\\` to `\` and `\n`/`\b`/`\f` become control characters. Write LaTeX with the **Write/Edit tools** only. In Python run from Bash, build backslashes with `chr(92)`.
2. Panel content lives in `one_mode_20261005/bodies/NN.tex` (`01`…`09`, `strip`, `footer`). After editing a body run `python splice_bodies.py` (it copies bodies into `sections/`). Never edit `sections/*.tex` directly.
3. Compile: from the poster folder, `powershell -NoProfile -ExecutionPolicy Bypass -File compile_section.ps1 -Section all`. Output goes to `output/pdf/IAS2026_Delay_Transport_20261005.pdf/.png`. Then **look at the PNG** (Read tool). A clean log does not prove the layout fits; overflow hides under the next panel.
4. Every number on the poster must exist in a CSV under `research_gold/checks/`. If a check fails, report it; do not tune until it passes.
5. Commit only with explicit pathspecs (never `git add -A`); the working tree has many unrelated untracked files. Do not push unless the user asks.
6. Long runs: start in background, poll with a bounded loop (`for i in $(seq 1 28); do …; sleep 10; done`), never a bare `sleep` > 60 s. Kill stray `julia`/`python` with `Stop-Process -Id`.

## 2. Tools already built (reuse, do not rewrite)

| Purpose | Command / file | Runtime |
|---|---|---|
| Full-spectrum root count (banded argument principle, validated vs eig at τ=0) | `from t3_full_spectrum_count import count` → `count(Model(p, tau), re0)` returns `(n, err)`; with `re0=-0.05` subtract 1 for the gauge root | ~20 s per call |
| Exact transport of gains | `from t1_spillover_ieee39 import transport, roots` | instant |
| Rule-based mode selection | `from t8_taylor_and_second_design import select_node, gains` | ~15 s |
| Share × delay spectral map | `research_gold/checks/t11_share_delay_map.py` → `T11_share_delay_map.csv` | ~6 min |
| Nonlinear delayed events (Julia, method of steps) | `julia --project=. research_gold/checks/t9_events.jl LABEL RULE NODE_RE NODE_IM H_MS RHO HORIZON [EVENT_IDX...]` with RULE ∈ `exact`/`lowfreq`/`none`; writes `research_gold/checks/T9_LABEL.csv` | ~10–12 min for 5 events × 60 s; run up to 3 in parallel |
| Central figure | `one_mode_20261005/build_central_fig.py` (event results typed in the `EVENTS` dict) | seconds |
| Other poster figures | `one_mode_20261005/build_poster_figs.py` | seconds |

Event indices for `t9_events.jl`: 1 = bus 8 −100 MW, 2 = bus 16 +100 MW, 3 = bus 16 −100 MW, 4 = bus 29 +100 MW, 5 = bus 29 −100 MW. Pass criteria: |Δf| ≤ 0.5 Hz, RoCoF ≤ 0.5 Hz/s, V ∈ [0.9, 1.1], SG reserve ≥ 0.002.
Rule-selected protected modes (from `T11_share_delay_map.csv`): ρ=0.80 → 0.697 Hz, 0.85 → 0.930 Hz, 0.875 → 1.041 Hz (`-16.34629153843205 + 6.543885486116494j`), 0.90 → 1.152 Hz. Get the exact complex node for ρ=0.90 with `select_node(p0)[0]` (p0 = P0 with `p0[:10]=0.90`).

## 3. Tasks, in order

### Task 1 — Nonlinear events for the share × delay map (required)
Goal: close the evidence for a "feasible region" panel.
1. Compute the ρ=0.90 node: small Python script in `research_gold/checks/` using `select_node`.
2. Run in background (three in parallel at most):
   - `t9_events.jl M90_rule_44 exact <re> <im> 4 0.90 60`
   - `t9_events.jl M90_rule_48 exact <re> <im> 8 0.90 60`
   - `t9_events.jl M90_none_44 none 0 1 4 0.90 30 2` (expected to abort)
   - `t9_events.jl M80_rule_44 exact <re80> <im80> 4 0.80 60`
3. Record pass counts (passed/5), worst |Δf|, worst RoCoF, min SG reserve.
**Acceptance:** CSVs exist; results are written into a new "ADDENDUM 4" of `04_primary_recommendation.md`, failures included.

### Task 2 — Replace panel 7 with the feasible-region map (required)
1. Add a function to `build_poster_figs.py` producing `generated/figures/share_delay_map.pdf/.png`. Grid: x = inverter share (80, 85, 87.5, 90, 92.5, 95 %), y = delay (40, 44, 48 ms). Each cell shows two marks: no retune vs rule (green = 0 roots beyond margin, red = >0, grey = base already fails at 40 ms). Overlay event results (passed/5) where Task 1 or earlier runs exist (design A ρ=0.875: 44 ms 5/5, 46 ms 5/5, 48 ms 4/5; ρ=0.85: 44 ms 5/5). Fonts ≥ 19 pt at figsize ≈ 7.5×4.5 in.
2. Rewrite `bodies/07.tex` (Write tool) with title `HOW FAR IT REACHES: INVERTER SHARE × DELAY`, the figure, one sentence: "The rule extends latency tolerance at every share where the 40 ms design is feasible; it does not raise the maximum share (92.5 % already fails at 40 ms)." Keep a one-line note that the damping-signature result (δD_H has one positive and one negative eigenvalue) is in the companion paper.
3. Splice, compile, inspect the PNG.
**Acceptance:** no overflow, labels readable, numbers match `T11_share_delay_map.csv` and T9 CSVs.

### Task 3 — Reproduce the nine-bus bank (recommended, cheap)
The code is at `research_gold/external_gpt_kp_only/BND_IEEE9_Python/` (`src/model.py`, `src/validate.py`, `data/config.json`).
1. Read `validate.py`; run the equilibrium, the three-PLL analytical transport 20→21 ms and the spectral check.
2. Expected (from the bank report): rightmost real part −0.130667 s⁻¹ for the analytical three-PLL design, +0.886150 s⁻¹ unchanged at 21 ms. Within 1e-4 → reproduced.
3. If reproduced, edit `bodies/06.tex`: change "Third model, reported not re-run" to "Third model, re-run here", keeping "assumed dynamic data".
**Acceptance:** numbers recorded in Addendum 4 with the command used.

### Task 4 — Verify the six footer references (required)
`bodies/footer.tex` lists: [1] Marković et al., IEEE TPWRS 2021; [2] Huang et al., arXiv:1903.05489; [3] Michiels & Gumussoy, arXiv:2003.05496; [4] Dörfler et al., IEEE TPWRS 2014 (sparsity-promoting wide-area control); [5] Bindel & Hood, SIAM J. Matrix Anal. Appl. 2013; [6] Johansson, IEEE Trans. Comput. 2017 (Arb).
Use WebSearch (standard mode) to confirm venue/year of each. Consider replacing [5] or [6] with a dominant-pole-placement-with-delay reference if one is confirmed (e.g. the PI–PR² dominant pole placement paper from ITU), since that is the closest prior art to the transport law. Fix only what is wrong.
**Acceptance:** a short table in Addendum 4: reference → confirmed / corrected (URL).

### Task 5 — Poster QA pass (required)
Open `output/pdf/IAS2026_Delay_Transport_20261005.png` and crop each panel (PIL) to check:
- no text outside a panel, no overlap of IN WORDS boxes with the next panel;
- every number traceable (spot-check five numbers against the CSVs);
- panel 8 still says 48.9 ms is spectral reach, not an operating limit;
- the central figure's `EVENTS` dict matches T5/T5b/T5c/T9 CSVs.
Fix with the Write/Edit tools, re-splice, re-compile.

### Task 6 — Optional, only if Tasks 1–5 are done and time remains
Cooperative compensation, one bounded attempt (Addendum 3 lists it as unexplored): protect the pole that actually crosses at 44 ms (find it from the no-retune 44 ms spectrum), allow three sites, keep all K_I fixed, minimise the worst catalogued root. Stop after one run of ≤ 15 min. Report the outcome in Addendum 4; **do not put it on the poster** unless it gives 0 roots beyond the margin AND 5/5 events.

### Task 7 — Commit (only if the user asks)
Pathspecs only:
```
git add research_gold/ reports/poster/ias2026/one_mode_20261005/ output/pdf/IAS2026_Delay_Transport_20261005.pdf output/pdf/IAS2026_Delay_Transport_20261005.png
```
Exclude `research_gold/checks/*_log.txt` if large. Message: "Add delay-transport poster and gold-result evidence". End with the attribution lines from the system reminder.

## 4. What to report back to the user (Spanish, short)

1. Task 1 table (share, delay, roots, events).
2. Whether the nine-bus bank reproduced.
3. Reference corrections.
4. Final poster PNG sent with SendUserFile.
5. Anything that failed, stated plainly.

## 5. Budget guide

Tasks 1 + 2 + 5 are the core (~1.5 h wall time, mostly waiting on Julia). Task 3 ~20 min, Task 4 ~15 min. Task 6 only with spare time.
