# FINAL REPORT: feedback-cycle falsification in the full IEEE-39 SG/GFL model (2026-10-05)

## 1. Executive verdict
**PARTIAL_PASS.** The first half of the tested statement survives; the second half is refuted.
- Survives: SG→GFL replacement creates a distinct PLL/controller family (4.5–5.0 Hz, PLL-state participation 0.90–0.99) that does not exist in the all-SG model; its first delay crossing is a collective closure (eig(Q) = −1 with all ten local factors regular); the dominant pair cycles at the leading critical roots are 35–36, 30–37 and 33–34, as predicted.
- Refuted: targeted sparse PLL retuning does not break that core. The preregistered S0 = {30,33,36,37} fails in Ki-only, Kp-only and joint form (worst root stays +0.65 to +0.69 s⁻¹ at 44 ms). Every four-site support fails. A feedback-core-guided ordering needs all ten sites; a physical-coupling ordering needs six (with large gain cuts and 4/5 events); random orderings need 8–10 (mean 8.8).
- What does repair it: dense changes. One common Kp reduction of 14 % (28.27 → 24.19) restores the full margin at 44 ms and passes 5/5 nonlinear events.
Decision gate: PARTIAL_PASS (not STRONG_PASS). Exact frequency (7.879 Hz) and the reduced-model numbers were not matched; the full-model family is at 4.65 Hz.

## 2. Frozen model provenance
HEAD `8f21ec0b`, branch `codex/collective-interaction-bounds-20261003` (six tracked files dirty, none used). Julia 1.11.9, PowerDynamics 5.0.0, NetworkDynamics 1.3.0, OrdinaryDiffEqRosenbrock 2.7.1; Python 3.11.10. Model: `ReducedDAE.jl` (`:physical_supply`, 204 states, algebra eliminated exactly), exact-exponential reader `interaction_decision_20261004/model.py`, nonlinear `DelayedEvents.jl`. Baseline rho = 0.875, Kp = 28.2743, Ki = 246.740, 40 ms; stress case 44 ms. Hashes: `ENVIRONMENT.md`, `PREREGISTRATION.sha256` (e327c2f1…), amendments 01–05 with their own hashes. Baseline parity (TABLE02): equilibrium residual 2.8e-9, Julia Jacobian equals the Python reader to 1e-16 (rho 0.125–0.875), baseline roots reproduced to 1e-10, 0 roots beyond margin at 40 ms and 8 unstable / 10 beyond margin at 44 ms (as in earlier work), baseline five events 5/5. Label: POWERDYNAMICS_VALIDATED (spectral parity), NONLINEAR_TDS_VALIDATED (events).
Scope note: linear analysis uses the exact linearisation of the existing model, not a new model. The all-SG model is the same code with rho = 0 (114 states, no PLL states).

## 3. Preregistered Python predictions (outcomes)
P1 all-SG has no PLL family: PASS. P2 distinct converter family after replacement: PASS (at 4.5–5 Hz, not 7–8 Hz; frequency not required). P3 localisation 35–36, 30–37, 33–34: PASS at the leading root (top-2 contains frozen pairs); 5 of 6 critical roots; one root (4.637 Hz) closes on 36–37. P4 Ki more authority than Kp: REFUTED (Kp is 1.7–2.0 times larger per unit log-gain, threshold was ≤ 0.5). P5 S0 first: done; S0 FAILED. P6 forced response: PARTIAL (below). P7 distinguishable from all-SG resonance: PASS.

## 4. All-SG reference
19 oscillatory physical modes, 0.05–1.5 Hz, all participation 1.0 on SG states, rightmost real part −0.098 s⁻¹ (TABLE03, FIG01). No mode above 1.5 Hz. Load-to-bus-frequency response (exploratory, TABLE03b) is small (≤ 4 mHz/MW).

## 5. SG→GFL modal emergence
At rho = 0.875 the 40 ms baseline has nine PLL-family roots at 4.5–5.0 Hz (PLL participation 0.90–0.99, damping 0.03–0.17) plus SG modes at 0.8–1.9 Hz (PLL participation < 0.07) (TABLE04a). Tracked by Newton continuation and MAC (TABLE_C01, FIG02) the family exists for every rho > 0 (it is ten local PLL loops at small share, split by network coupling into 4.5–5.0 Hz as rho grows). All 19 all-SG branches were matched at rho = 0.125 (MAC ≥ 0.7). Two PLL branches (PLL2, PLL5) have min MAC 0.76/0.64 so are reported as a family. Label: NUMERICALLY_VERIFIED.

## 6. Delay dependence
Exact exponential delay, no Padé in any reported root (Padé only seeds). Uniform delay 40 → 50 ms (TABLE05, FIG03): the four rightmost branches (4.91, 4.95, 4.99, 4.81 Hz at 40 ms) cross the margin at 41.57, 41.68, 42.27, 43.26 ms and the axis at 41.68, 41.79, 42.39, 43.38 ms; two more cross near 43.9 and 44.0 ms. Counts (argument principle, floating point, not certified): N_unstable/N_margin = 0/0 (40, 41 ms), 4/4 (42), 6/6 (43), 8/10 (44), 14/14 (46), 16/16 (48, 50). Hidden-block condition numbers recorded (TABLE05). A first-pass tracking bug that merged two branches was found and fixed (NEGATIVE_RESULTS).

## 7. Feedback-return construction
G(s) = C (sI − A0)⁻¹ B (ten PLL ports, grid closed), K = I − G E, local factors L_i = 1 − G_ii e^{−sτ_i}, Q_ij = −G_ij e^{−sτ_j}/L_i. Determinant identity det Δ = det(sI − A0) ∏L_i det(I+Q) verified at 24 random complex points, max error 1.1e-12 (TABLE_E00). EXACT_IDENTITY (numerically verified). At each critical root min|1 + eig(Q)| ≤ 4e-10 while |L_i| ∈ [0.011, 0.28]: local loops are regular, the closure is collective. Caveat: |L_i| is small (0.01–0.03) at its minimum, never near zero; "regular" is relative. On the imaginary axis σ_min(I+Q) at 40 and 44 ms is in FIG04.

## 8. Dominant SCC / cycle
Leading 44 ms root (+0.969, 4.649 Hz): top pairs 35–36 (|1−p| = 0.18), 30–37 (0.21), 30–35, 30–36, 33–34. Next root (+0.929, 4.681 Hz): 35–36 with |1−p| = 0.008; root 4.566 Hz: 33–34 (0.018); root 4.637 Hz: 36–37 (0.005); 4.713 Hz: 35–36 (0.27). 40 ms baseline root: 35–36 (0.19), 30–37 (0.22). SCCs (edge threshold 0.1 max|Q|) are almost whole (7–10 buses), so the localisation lives in the cycles, not in the strongly connected components (TABLE_F01–F04, FIG05/06). Label: NUMERICALLY_VERIFIED.

## 9. Physical graph versus dynamic graph
The three frozen pairs are exactly the top-3 physical Kron-coupling pairs (33–34, 35–36, 30–37). Spearman(physical coupling, |p_ij|) = 0.39–0.61. So at the leading choices the dynamic ranking adds no information beyond the physical ranking. Physical ordering also reached a successful sparse support earliest (k = 6). Consistent with "physical graph = candidate locality"; does not support "dynamic cycle ranking = better guide" (TABLE_G01/G02).

## 10. Kp/Ki authority
Exact simple-root sensitivities (finite-difference check 6e-8) at the five critical roots (TABLE08, FIG07): the leading root is controlled almost entirely by site 30 (dRe/dlogKp = +8.5, dRe/dlogKi = +4.6 s⁻¹), then 37, 39; sites 35 and 36 (the cycle pair) have authority ≤ 0.05. Kp raw authority per unit log-gain is 1.7–2.0× Ki, so P4 is REFUTED. Kp-only fails when sparse and succeeds when dense (uniform Kp 28.27 → 24.19).

## 11. Sparse repair (tau = 44 ms, bounds 0.25–4× nominal, margin −0.05, exact spectrum after each step)
| design | result |
|---|---|
| S0 Ki-only / Kp-only / joint | fail: worst +0.688 / +0.654 / +0.658 s⁻¹ |
| sens4, phys4, core4 (four sites, joint or Ki) | fail |
| size scan, joint: physical order | success at 6 sites ({30,33,35,37,39,32}) |
| sensitivity order | success at 9 |
| feedback-core order | success at 10 |
| 5 random orders | 9, 8, 9, 10, 8 |
| uniform Kp / Ki / joint (1–2 scalars) | success (Kp −14.4 %, Ki −21.5 %, joint −9 %/−5.6 %) |
| all-ten nodal Kp / Ki / joint | success in 2–4 iterations |
Engine v1 (no step acceptance) gave the same qualitative picture; kept in `derived/engine_v1/`. Support search is rule-ordered, not exhaustive: a smaller successful support may exist (UNRESOLVED). Final margins are −0.052 to −0.053 (just inside the declared margin).

## 12. Forced-response experiment (5 MW-equivalent sinusoid, bus 29 by the frozen residue rule, 30 s; TABLE11, FIG09)
- all-SG: no PLL; bus-frequency deviation ≤ 0.7 mHz at all three frequencies.
- 40 ms baseline (stable): PLL error peak 0.17°, 0.38°, 0.28° at 0.7 f*, f*, 1.3 f* (resonant ratio 2.2 and 1.4: below the preregistered ≥ 3).
- 44 ms stress (unstable): integration aborts (retcode MaxIters) at 9.9–11.9 s with PLL error ≈ 40°; time to 5° is 4.31 s at f* versus 6.73 and 6.41 s off resonance. Peak criterion not meaningful once growth saturates.
- sparse-repaired (physical 6-site, stable): peak 2.67° at f* versus 0.22° and 0.28°, i.e. 12× and 9.6× resonant ratio; the design is bounded, but its f* response is 7× larger than the 40 ms baseline.
P6 is therefore PARTIAL: resonance at f* is real in the stable cases, the literal criterion on the stress case cannot be evaluated. Label: NONLINEAR_TDS_VALIDATED (small amplitude, exploratory design).

## 13. Nonlinear validation (frozen five events, TABLE12, FIG10)
Baseline 40 ms 5/5 (max |Δf| 0.471 Hz, min SG slack 0.0074). Uniform Kp, Ki, joint and nodal all-ten Kp/Ki/joint: 5/5 each (max |Δf| 0.471–0.486 Hz, min slack 0.0032–0.0073). Physical 6-site sparse design: 4/5 (min SG slack 0.0014 < 0.002). The two sparse failures-by-spectrum (S0, core4) are not run nonlinearly because they fail the spectral margin. A finite event set is a declared security campaign.

## 14. RHP-zero check (exploratory, TABLE_L01)
No right-half-plane zeros in the ten PLL return channels (SISO or 10×10 MIMO): NEGATIVE_RESULT. Load-to-generator-bus-phase channels have RHP zeros already in the all-SG model (10–12 real zeros, 7–20 rad/s) and more (10–30, including 1.1 Hz complex pairs for loads at 8 and 16) at rho = 0.875. The confirmation column in the table reuses the same pencil, so this is NOT independently confirmed and is not a claim.

## 15. Negative results (all kept)
S0 fails in all three variants; four-site supports fail; core-guided order is the worst of the three rules; P4 refuted; first-pass tracking bugs (Phase C/D) and engine v1 documented; 7–8 Hz family not present (4.65 Hz); local factors never near zero; SCC not localised; frozen-pair dominance is explained by physical coupling; stress-case forced response saturates; RHP-zero search negative for PLL channels and not independently confirmed elsewhere. See `NEGATIVE_RESULTS.md`.

## 16. What is actually new
Candidate contribution, as tested: (a) exact-delay determinant factorisation det Δ = det(sI − A0) ∏L_i det(I+Q) on the full 204-state model with a numerically verified sign convention; (b) its use to show the leading delay crossings are collective closures with regular local factors, whose closest pair cycles sit on the electrically strongest neighbourhoods; (c) a fair comparison showing that this localisation does NOT translate into a smaller retuning support in this model. (c) is a negative result and is the main outcome.

## 17. What remains classical
Nyquist/return-ratio and loop-dressing ideas, modal sensitivity, graph/GSP descriptions, sparse control, small-signal analysis, and delay-compensating retuning (a uniform 14 % Kp cut or the exact-transport rule from the poster work). Uniform retuning being sufficient makes the sparse-repair framing unnecessary here.

## 18. Poster recommendation
Do not print the claim that targeted sparse retuning breaks the core. A defensible panel: "the leading delay crossings are collective: local PLL loops stay regular while network coupling closes a cycle (35–36, 30–37, 33–34), the same neighbourhoods as the strongest electrical coupling; yet the repair is global, not local". Keep the one-mode delay-transport poster headline.

## 19. Exact claims safe to print
- In the IEEE-39 model with 87.5 % inverter share, the delay-induced instability at 41.6–44 ms comes from the 4.5–5 Hz PLL family, absent from the all-SG model (exact delay, floating-point counts: 4/6/8 unstable roots at 42/43/44 ms).
- At these roots all ten local PLL factors are regular while eig(Q) = −1; closest pair cycles are 35–36 (|1−p| = 0.008), 33–34, 36–37.
- Retuning four sites (including {30,33,36,37}) does not restore the margin at 44 ms; a common Kp reduction of 14 % does, and passes the five frozen events.
- Frozen five events are a finite security campaign.

## 20. Claims that must NOT be printed
"Certified" root counts; the 7–8 Hz family; 35–36 as a universal cause; sparse core-guided retuning as a repair; Ki dominating Kp; the RHP-zero findings; any universal delay limit or optimal support; EMT validity; that the dynamic cycle graph improves on the physical graph.

## Artifact index
Preregistration and amendments: `PREREGISTRATION.md`, `AMENDMENT_01…05_*.md` (+ `.sha256`). Code: `src/`. Tables: `derived/TABLE01…TABLE12*.csv`, `TABLE_C01/D01/E00/E01/F01-F04/G01/G02/L01*.csv`, `derived/engine_v1/`. Figures: `figures/FIG01…FIG10` (png + pdf). Raw: `raw/jac`, `raw/forced`, `raw/events`, `raw/cases`. Baseline event copy: `research_gold/checks/T9_FBK_base40.csv`. Logs: `logs/`.
