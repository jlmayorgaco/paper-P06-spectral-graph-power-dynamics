# 04 — Primary recommendation

**Decision: make the poster about the interpolation law of local PLL retuning — "One PLL, one mode".**

After reconstructing about 60 equations (ledger `01`) and discarding ten weaker narratives (list at the end), the strongest result is this:

> A delayed PLL has two real gains, so it can hide a change of delay from exactly one complex frequency. The gains that do it are closed-form and need no network model; the error left everywhere else is an exact quadratic factor times a network residue; and the standard choice of that frequency — the low-frequency one, which is what first-order delay compensation does — makes IEEE-39 unstable sooner than not compensating at all.

Everything in it was verified today on the exported 204-state IEEE-39 delay model (checks T1–T7 in `research_gold/checks/`). The missing test "R" (nonlinear events) was run; its outcome is reported below, including where it fails.

---

## THE RESULT

For a delay change `h` at PLL `i` and any complex node `λ = α + jω`:

1. **Transport (exact, local).** The unique real gains that keep the loop factor `κ(s) = (k_p s + k_I)e^{−sτ}` unchanged at `λ` are
   `k_p(h) λ + k_I(h) = (k_p λ + k_I) e^{λh}`.
2. **Remainder (exact, finite).** For every `s`:
   `κ(s;h) − κ(s;0) = −(s−λ)(s−λ̄) e^{−sτ} ∫_0^h k_p(t) e^{−st} dt`.
3. **Toll (exact, first order).** Every other mode `μ` of the grid moves by
   `dμ/dh = −k_p (μ−λ)(μ−λ̄) ∂μ/∂k_I`.
4. **Budget.** The transport stays inside positive gains while `h < h* = atan2(ω k_I, |λ|²k_p + α k_I)/ω`; for `α ≈ 0` this is the PI's own phase lag at `ω` divided by `ω`.

## THE EQUATION

```
κ(s;h) − κ(s;0) = −(s−λ)(s−λ̄) · e^{−sτ} ∫₀ʰ k_p(t) e^{−st} dt
```

One line says it all: zero at the node, quadratic everywhere else.

## THE THEOREM / PROPOSITION

**Proposition (interpolation law).** Consider the linearised delayed grid `Δ(s) = sE − A0 − Σ_i b_i(k_i) c_iᵀ e^{−sτ_i}` in which PLL `i` enters only through its loop factor `κ_i(s) = (k_{p,i}s + k_{I,i})e^{−sτ_i}` (equivalently `∂μ/∂k_p = μ ∂μ/∂k_I` at simple roots). Let `τ_i → τ_i + h_i`.

(a) *Complete parametrisation.* Every differentiable real retune is `k_p' = a k_p + k_I`, `k_I' = −b k_p` for some `(a,b) ∈ R²`, and then `∂_h κ_i = −(s² − a s + b) k_p e^{−s(τ_i+h)}`. The node is the root pair of `s² − a s + b`; `(a,b) = (0,0)` is first-order (moment) compensation, `(−k_I/k_p, 0)` is "do nothing".
(b) *Invariance.* If all retuned sites use the same node `λ` and `λ` is a simple root, then `λ, λ̄` stay roots for every set of sites and every heterogeneous `h_i`, with the same mode shape outside the PLL integrator states.
(c) *Toll.* For any other simple root `μ`: `dμ/dh_i = −k_{p,i}(μ² − aμ + b) ∂μ/∂k_{I,i}`.
(d) *Uniqueness.* If `∂λ/∂k_{I,i} ≠ 0`, no other real retune at site `i` preserves `λ`; in particular one PLL cannot hold two distinct modes.

*Proof.* (a) and the remainder: differentiate `κ_i` and integrate. (b): each `κ_i(λ)` is unchanged, so the reduced return difference at `λ` is unchanged. (c): chain rule with `μ_τ = −μ(k_pμ_{k_p} + k_Iμ_{k_I})` and `μ_{k_p} = μ μ_{k_I}`. (d): the real Jacobian of `(k_p,k_I) ↦ k_pλ + k_I` is non-singular for `ω ≠ 0`. ∎

**Proven:** (a)–(d) as identities of the linearised model. **Not proven:** that any node is optimal; that the first-order toll bounds the finite motion; anything about the nonlinear system. The half-plane/LP use of (c) for large steps **failed** (T7: an LP over 11 catalogued modes extrapolated to a design with 16 unstable roots) and must not be claimed.

## WHY THIS IS NEW (exact delta, stated conservatively)

Known: lead compensation; exact/dominant pole placement for PI with delay; PLL retuning guidance for multi-converter grids; Smith predictors for PLL delay. Not found in a targeted search (not a proof of priority): the finite remainder with its quadratic factor and the toll identity; the exact composition over sites and heterogeneous delays without a network model; and the measured reversal that moment-preserving compensation is worse than none in a converter-dominated grid. Label: **known components, possibly new synthesis — literature check still required.** Queries to run: `"PI" "time delay" "pole placement" remainder OR residual "(s-λ)"`; `delay compensation PLL "retuning" eigenvalue invariant multi-converter`; `"first-order" delay compensation destabilizing grid-following PLL`.

## WHY THIS IS IAS

It is a tuning rule a converter engineer can apply per device (own gains, own delay change, one target mode), with a stated limit and a stated price. It concerns latency in PLL-synchronised converters, the practical bottleneck of high inverter penetration, and it is validated on a nonlinear delayed IEEE-39 model at 87.5 % inverter share.

## WHY THIS CAN WIN

- One equation, one reversal, one figure, one rule.
- The reversal is concrete: the "obvious" compensation is the wrong one.
- It turns the project's earlier loose ends into consequences: decay-only retunes do not compose while full-pole retunes do (invariance (b)); moment cancellation went unstable at 41 ms (node at 0); an exactly assigned pole did not stabilise the 9-bus bank (uniqueness (d)); one delay creates a (+,−) damping pair (the same rank-one, two-gain structure seen in the damping operator).
- It is honest: the poster can show where the rule runs out (budget, frequency-security margin).

## EVIDENCE OBTAINED TODAY (IEEE-39, 204 states, ten PLLs, base delay 40 ms)

| Check | Result |
|---|---|
| T1 toll identity | 219 well-scaled mode pairs, max rel. error 3.9e-5, sign agreement 100 % |
| T2 invariance | protected pole drift 0.0 with all ten PLLs transported, h ≤ 10 ms |
| T3 full-spectrum count (validated vs `eig` at τ = 0) | no retune: 0 / 0 / **4** / **8** unstable roots at 40 / 41 / 42 / 44 ms; node at 0.05 Hz (moment-like): 0 / **4** / **12** at 41 / 42 / 44 ms; node at 4.91 Hz: **0 unstable and 0 beyond the −0.05 guard** at 42, 44, 46, 48, 48.9 ms |
| Budget | closed form 8.99 ms for the 4.91 Hz node; 50.2 ms for the 1.04 Hz node |
| T4 prediction | predicted vs exact worst root over 11 candidate nodes: Spearman 0.99 |
| T6 node choice | node at the heavily damped 1.04 Hz mode: 0 unstable, 0 beyond guard at 42–48 ms, and K_I falls only to 213 at 44 ms (138 for the 4.91 Hz node) |
| T7 retune plane | four special points agree with counts; LP extrapolation fails (negative result) |

**Nonlinear delayed events** (Julia, method of steps, frozen five events, limits 0.5 Hz, 0.5 Hz/s, V ∈ [0.9, 1.1], SG slack ≥ 0.002):

| Design | Delay | k_I | Outcome |
|---|---|---|---|
| Original gains | 44 ms | 246.7 | integration aborted (`MaxIters`) on bus 16 +100 MW; consistent with 8 unstable roots |
| Original gains | 42 ms | 246.7 | bus 16 +100 MW, 30 s: **meets the guards** (0.444 Hz, slack 0.0074) although the spectrum has 4 unstable roots (+0.096 s⁻¹ at 4.8 Hz). The 0.5 s-windowed metrics do not see a slow 4.8 Hz growth in 30 s |
| Node 4.91 Hz | 42 ms | 192.5 | **5/5 pass**: \|Δf\| ≤ 0.486 Hz, RoCoF ≤ 0.173 Hz/s, slack ≥ 0.0032 |
| Node 4.91 Hz | 44 ms | 137.8 | all five trajectories bounded; **only 2/5 meet every guard**: bus 16 −100 MW reaches 0.505 Hz; SG slack exhausted at bus 16 +100 MW (≈ 0) and bus 29 +100 MW (0.0015) |
| Node 1.04 Hz | 44 ms | 213.3 | **5/5 pass**: \|Δf\| ≤ 0.480 Hz, RoCoF ≤ 0.184 Hz/s, slack ≥ 0.0047 |

Reading: stability is recovered by any cluster-side node; the frequency-security reserve is spent in proportion to how much `k_I` the node consumes. That is the interpolation law again — one PI cannot hold both the mode and the low-frequency moments — and it gives the design rule below. Two cautions follow from the table: the event guards alone do not detect a slowly growing 5 Hz instability (always pair them with the root count), and the 4.91 Hz node at 44 ms is a genuine failure of the frequency-security contract, not of stability.

## THE STORY

**For a fourth-year student.** Each inverter listens to the grid through a PLL. If the PLL hears the grid a few milliseconds late, the grid can start to oscillate. You can retune the PLL to cancel the lateness, but only for one oscillation at a time. Cancel it for the wrong one and you make things worse; cancel it for the right one and the grid tolerates almost 9 ms more delay.

**For a controls/power PhD.** A PI with delay has two real degrees of freedom, so its loop factor can interpolate the delay-free one at one conjugate pair. The interpolation remainder is `(s−λ)(s−λ̄)R_h(s)`; because each PLL enters the grid characteristic matrix as a rank-one term, the motion of every other root is that factor times the residue `∂μ/∂k_I`. Moment matching is the node at the origin; in a grid whose fragile modes sit at 4–5 Hz that node carries the largest remainder exactly where the residues are largest.

**For a skeptical IAS professor.** Nothing here needs a new theorem of control; the claim is a complete and exact account of what a local PI can do about latency in a networked converter system, checked on a full-order delayed model with root counts and nonlinear events, including the cases where it fails. The measured reversal and the closed-form limit are the contributions.

## THE ONE FIGURE

`research_gold/FIG_GOLD_one_pll_one_mode.png` (draft built from T2/T3): rightmost root versus PLL delay from 40 to 50 ms for three strategies — no retune (unstable from 42 ms), first-order/moment compensation (unstable sooner and faster), closed-form transport at the 4.91 Hz node (flat below the guard up to the dotted budget line at 48.99 ms) — with the full-spectrum root counts written on the curves. Add the 1.04 Hz node as a fourth curve and an inset of the s-plane showing the node and the quadratic toll.

## THE ONE EXPERIMENT

Done: T3 (root counts) + T5/T5b/T5c (nonlinear events). To finish before printing: (1) nonlinear events for the 1.04 Hz node at 46 and 48 ms, to find where its frequency-security reserve ends; (2) the same three-strategy sweep from a second base design (e.g. the 88.455 % design) to show it is not a single-case effect; (3) a longer or unwindowed nonlinear run of the untuned 42 ms case to exhibit the 4.8 Hz growth directly.

## THE DESIGN RULE

1. Never compensate PLL latency at low frequency (`k_p += h k_I`).
2. Compensate at a node on the PLL-cluster side; among nodes that keep the spectrum safe, pick the one that spends the least `k_I` (a strongly damped node acts mostly as the gain reduction `e^{αh}`).
3. Check the phase budget `h*` before promising a latency.
4. Each device applies the rule with its own delay; no communication and no network model are needed for the gains — the network is needed only to choose the node (one eigen-analysis).

## THE ONE-LINE TAKEAWAY

**"A PLL can hide its delay from one mode only — choose the 5 Hz one, not the slow one, and the grid gains 9 ms."**

Suggested title: *One PLL, One Mode: the exact price of hiding latency in inverter-dominated grids* (keep "Beyond Nodal Damping" as the series tag: the price is set by a network residue, not by a nodal coefficient).

---

## DO NOT PUT THESE AT THE CENTRE OF THE POSTER

1. **Maximum SG→GFL replacement / delay frontier.** No upper bound exists; M9–M17 blocked.
2. **"Safe alone can fail together" (buses 30–37).** One pair, a stable pole, 2 violations in 88 resolved cases. Use as a witness for invariance (b).
3. **H4 and other subset-stable witnesses.** Strong pictures, no predictive theorem; different model (TX4). Background only.
4. **Möbius / Shapley / cycle attributions.** Representation-dependent; retire as physics.
5. **"Stronger grid, less stable".** Known phenomenon; negligible on the full model; surrogate tuned.
6. **Commutator / GSP mode-mixing as a stability criterion.** A linear-algebra identity with no stability consequence.
7. **Network inertia floor / area conservation.** Elegant, but reduced model only and not a grid-code quantity. One supporting sentence at most.
8. **Counting law for damping (m ≥ r).** Rigorous and worth a side panel as the damping-domain face of the same structure, but positivity at one frequency is not what stabilises the grid.
9. **Non-normal growth, hidden margins, cycle-space geometry, inertia-as-separator.** Textbook bounds or unreplicated scalings.
10. **`K_p* = 2rL⁺`, exact action-space rank ≤ 30, pair determinant.** Correct algebra that serves the result; not a headline.

---

## ADDENDUM (5 Oct, after external review of this recommendation)

Corrections adopted:
- Motto changed from "One PLL, one mode" to **"Preserving one mode does not protect the whole grid."** The result is about *preserving* one prescribed conjugate pair with two real gains; a retune moves many other poles and some improve.
- The interpolation law, the rank signature of the damping change and the support-invariant response are three complementary consequences of the local architecture, not one theorem.
- The transport formula needs no network model; **choosing** the protected mode does.
- 48.9 ms is the spectral reach of the 4.91 Hz choice, not an operating limit.

New checks (scripts and CSV in `checks/`):
- **T8 comparator.** The first-order modal rule (`k_p + h(2αk_p + k_I)`, `k_I − h|λ|²k_p`) gives the same root counts as the exact transport at 42–48 ms for both protected modes (0 unstable / 0 beyond margin). The low-frequency rule gives 14 / 14 at 46 ms. So the failure is the choice of interpolation point, not the truncation.
- **T9 events, design A.** Protect 1.04 Hz: 5/5 at 46 ms; 4/5 at 48 ms (bus 16 +100 MW, SG reserve 0.0018 < 0.002). Low-frequency rule at 44 ms: integration aborted on both events tried.
- **T8/T9 second base design (85 % inverter share at every site).** Rule fixed on design A and applied unchanged selects the 0.93 Hz mode. At 44 ms: no retune 10 unstable / 12 beyond margin and integration aborted; low-frequency rule 12 / 12; rule-selected transport 0 / 0 and **5/5 events** (|Δf| ≤ 0.404 Hz, RoCoF ≤ 0.185 Hz/s, SG reserve ≥ 0.027).

Still open: the nine-bus bank is cited from its report, not re-run here (its code package is not in this repository); no interval-certified counts; no comparison with other compensators (e.g. Smith predictor); one network.

Poster built with this hierarchy: `reports/poster/ias2026/one_mode_20261005/` → `output/pdf/IAS2026_Delay_Transport_20261005.pdf`.

## ADDENDUM 2 (5 Oct): comparison with the portfolio and co-design lines; T11

- **Portfolio line (TX3/TX4, PD39).** H4 = {30,33,35,37} is a clean witness (15 safe / 1 unsafe; blind prediction 16/16 on four buses), but its own reports limit it: the nine-bus transfer census is 395/512 with 116 false-safe verdicts; the PD39 255+1 and blocker-atlas campaigns both closed as CASE C; H4 is a minimal blocker in only about 11.5 % of sampled conditions; the model has no governor and no delay. It explains *that* sets fail, not which or how to prevent it.
- **Co-design / maximum replacement line.** No upper bound and no validated frontier (mega experiment M9–M17 blocked; the 94.2 % candidate retracted; 88.455 % is one point on one search path).
- **T11 (spectral screen only), `checks/T11_share_delay_map.csv`.** Uniform inverter share 80, 85, 87.5, 90 %: baseline gains are inside the −0.05 margin at 40 ms and have 8–12 roots beyond it at 44 ms and 16 at 48 ms; the rule-selected transport has 0 at both delays for all four shares (protected mode 0.70–1.15 Hz, k_I ≈ 213 at 44 ms and ≈ 183 at 48 ms). At 92.5 % and 95 % the 40 ms base itself is outside the margin (1 and 3 roots): the rule extends latency tolerance across shares; it does not raise the maximum share. No nonlinear events were run for T11.

## ADDENDUM 3 (5 Oct): cooperative K_p-only compensation with all K_I fixed (external proposal)

Math checked: (i) one PLL cannot keep a complex pole by changing only its own K_p (the required rate `λk_p + k_I` is not real); (ii) with two K_p the determinant is exactly `c0 + c1 x + c2 y + c12 xy`, so real solutions come from a real quadratic; (iii) with all K_I fixed the low-frequency moments H0, H1 are unchanged under the moment-law hypotheses. All three are correct.

Reproduced here: the nine-bus case of the proposal (`external_gpt_kp_only/check_remote_repair.py`): delay of PLL 1 from 20 to 21 ms gives +0.275 s⁻¹ and 2 unstable roots; changing only K_p at PLLs 2 and 3 gives −0.128 s⁻¹ and 0. Same numbers as reported.

IEEE-39 tests (`checks/T12_*.csv`, `T12b_single_site_delay.csv`), target pole 4.91 Hz, spectral counts:
- **All ten delays 40 → 44 ms.** Of 45 pairs only (37,38) has admissible real solutions; they keep the target pole but leave 8 and 10 unstable roots. The minimum-norm solution over all ten K_p leaves 10. (Exact per-site transport: 0.)
- **One delay only, +8 ms at bus 30:** no admissible two-K_p solution among the other nine sites. **+8 ms at bus 37:** pair (30,39) keeps the pole (gain change 1.3 %) but 2 roots stay beyond the margin, the same as with no retune. Local transport at the delayed site: 0 beyond the margin.

Reading: on IEEE-39 the cooperative K_p-only construction preserved the chosen pole and did not restore the margin in any case tested. It is one more instance of the poster's statement, not a repair. Not explored: other target poles (e.g. the pole that actually crosses), three or more sites with an objective on the cluster, or mixed K_p/K_I allocations. Status: works in the nine-bus bank; **negative on IEEE-39 so far**; keep out of the poster.


## ADDENDUM 4 (5 Oct): share × delay events, nine-bus re-run, reference check

**Nonlinear events, rule-selected protected mode (5 events × 60 s; limits |Δf| ≤ 0.5 Hz, RoCoF ≤ 0.5 Hz/s, V ∈ [0.9, 1.1], SG reserve ≥ 0.002).** Nodes from `checks/t13_nodes.py`; CSV `checks/T9_M*.csv`.

| share | delay | node (Hz) | spectral roots (T11) | events | worst \|Δf\| | worst RoCoF | min SG reserve |
|---|---|---|---|---|---|---|---|
| 0.80 | 44 ms | 0.697 | 0 | 5/5 | 0.310 | 0.174 | 0.0551 |
| 0.90 | 44 ms | 1.152 | 0 | **0/5** | 0.596 | 0.176 | ≈ 0 |
| 0.90 | 48 ms | 1.152 | 0 | **0/5** | 0.610 | 0.156 | ≈ 0 |
| 0.90 | 44 ms, no retune, event 2 | – | 8 unstable | aborted (DDE MaxIters) | – | – | – |

Reading: at 90 % the full-spectrum count is clean but the nonlinear events fail the frequency and SG-reserve limits in all five events, so spectral margin is necessary, not sufficient. This is shown on the poster (gold cells). Earlier results (ρ = 0.875: 44 ms 5/5, 46 ms 5/5, 48 ms 4/5; ρ = 0.85: 44 ms 5/5) unchanged.

**Nine-bus bank re-run** (`checks/t14_ieee9_bank.py`, log `T14_ieee9_bank_log.txt`; dynamic data assumed). Gains from `data/config.json` (ρ = 0.75), 20 → 21 ms, spectral only:
- unchanged gains at 21 ms: rightmost root +0.886150 s⁻¹ (expected +0.886150);
- three-PLL exact transport toward the root −0.5279 + 49.58j: rightmost −0.130667 s⁻¹ (expected −0.130667). Reproduced within 1e-4.
- The 12/12 nonlinear events of the bank report were **not** re-run; the poster says so. Transport toward the other catalogued roots did not restore the margin (+0.87 to +0.99 s⁻¹ when protecting the 12.2, 15.1 or 2.1 rad/s mode), i.e. the choice of protected mode matters here too.

**Reference check (WebSearch, 5 Oct).**
| ref | result |
|---|---|
| [1] Marković et al., TPWRS 2021 | confirmed: Understanding Small-Signal Stability of Low-Inertia Systems, vol. 36 no. 5, 3997–4017, doi 10.1109/TPWRS.2021.3061434 |
| [2] Huang et al., arXiv:1903.05489 | confirmed (Impacts of Grid Structure on PLL-Synchronization Stability); footer topic corrected from "multi-converter PLL retuning" to "PLL stability and grid structure" |
| [3] Michiels & Gumussoy, arXiv:2003.05496 | confirmed |
| [4] Dörfler et al., TPWRS 2014 | confirmed (Sparsity-Promoting Optimal Wide-Area Control, vol. 29, 2281–2291, arXiv:1307.4342) |
| [5] Bindel & Hood, SIMAX 2013 | confirmed (arXiv:1303.4668) |
| [6] Johansson, IEEE Trans. Comput. 2017 | Arb paper (arXiv:1611.02831) confirmed; search showed ARITH-24 2017 only, the IEEE Trans. Comput. venue was **not** verified |

No replacement with a dominant-pole-placement-with-delay reference: none was confirmed.

Task 6 (single bounded cooperative-compensation attempt) was not run.

**Review of Addendum 4 (5 Oct):** the 90 % event failure is caused by the share, not by the retune. Base design at 90 % with original gains and 40 ms (`checks/T9_M90_base_40.csv`, events bus 16 +100 MW and bus 8 −100 MW): 0 of 2 pass, |Δf| 0.541 and 0.577 Hz, SG reserve ≈ 0 on bus 16 +100. Panel 7 text and map were corrected accordingly ("necessary, not sufficient" removed; the 90 %/40 ms cell now shows 0/2).

## ADDENDUM 5 (6 Oct): poster v2 (design phase) — `reports/poster/ias2026/one_mode_v2_20261006/`
Output `output/pdf/IAS2026_Delay_Transport_v2_20261006.pdf/.png`. New: time-domain figure from `checks/T15_*_trajectory.csv` (amplitude of the 4.2–5.6 Hz band of the 39 bus frequencies, event bus 16 +100 MW, 44 ms, ρ = 0.875). Findings: no retune and low-frequency rule grow (4e-5 → 5e-3 Hz in ~6 s, and 1e-4 → 3e-3 in 3 s) until the integrator diverges; protecting 4.91 or 1.04 Hz decays below 1e-6 Hz within ~7 s. The 20 s runs of the two unstable cases diverge and save nothing, so shorter horizons were used (8 s no retune, 5 s low-frequency); the figure marks the divergence. The 42 ms no-retune 30 s trajectory shows slow growth (6e-5 → 6e-4 Hz) while event limits are met. Also new: "price prediction" figure (T4). Honest note on that figure: the zoom shows prediction gets the ORDER of the nine safe modes right but not their exact values; the headline number is a rank correlation (0.99, dominated by separating the two slow modes from the nine safe ones).
