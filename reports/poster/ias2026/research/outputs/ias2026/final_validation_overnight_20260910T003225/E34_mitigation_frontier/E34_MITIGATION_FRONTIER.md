# E34 — mitigation Pareto frontier at four declared margins

Four margins declared before optimizing — **−0.02, −0.05, −0.10 and the base case
itself at −0.126478** — and five strategy families costed on three axes reported
**separately**, never collapsed into one scalar. 198 candidates, **198 meet their
target**, **10 Pareto-efficient points**.

The constrained quantity is the **spectral abscissa**, which is what
`trackA_mitigation_definitions.yaml` froze and what E21 used. The tracked
inter-area family envelope is reported beside it in every row, per the same
frozen rule.

## The Pareto front

| target | strategy | detail | PV MW lost | sync MVA | controller change | abscissa | zeta_min | m4 |
|---|---|---|---|---|---|---|---|---|
| −0.02 | **M1 converter only** | PLL wn ×1.11, zeta ×0.81, P bw ×0.49, Q bw ×1.16, filter ×1.26 | **0** | **0** | 0.803 | −0.0200 | 0.0057 | 0.0190 |
| −0.02 | M2 condenser | rating 0.16 at bus 30 | 0 | 166 | 0 | −0.0225 | 0.0063 | 0.0385 |
| −0.02 | M5 fractional machine | keep 0.15 of machine 30 | 156 | 156 | 0 | −0.0326 | 0.0091 | 0.0360 |
| −0.05 | **M1 converter only** | PLL wn ×1.05, P bw ×0.43, Q bw ×1.05, filter ×1.05 | **0** | **0** | 0.849 | −0.0500 | 0.0144 | 0.0377 |
| −0.05 | M2 condenser | rating 0.20 at bus 30 | 0 | 208 | 0 | −0.0678 | 0.0190 | 0.0501 |
| −0.10 | **M1 converter only** | PLL wn ×1.09, zeta ×1.04, P bw ×0.36, Q bw ×1.08, filter ×1.11 | **0** | **0** | 1.049 | −0.1000 | 0.0289 | 0.0704 |
| −0.10 | M2 condenser | rating 0.24 at bus 30 | 0 | 250 | 0 | −0.1110 | 0.0241 | 0.0682 |
| **base** | **M1 converter only** | PLL wn ×0.95, zeta ×1.06, **P bw ×0.32**, filter ×1.08 | **0** | **0** | 1.150 | **−0.1265** | 0.0340 | 0.0854 |
| **base** | M2 condenser | rating 0.26 at bus 30 | 0 | 270 | 0 | −0.1317 | 0.0241 | 0.0781 |
| base | M5 fractional machine | keep 0.25 of machine 30 | 260 | 260 | 0 | −0.1374 | 0.0265 | 0.0813 |

## What the frontier says

**A converter retune alone reaches every declared margin, including the full base
margin, with zero megawatts of PV forgone and zero synchronous capacity added.**
The price is a controller change norm rising from 0.80 to 1.15, and the dominant
move is always the same: **weaken the outer active-power loop**, ×0.49 at the
easiest target and ×0.32 at the hardest. This improves on E21, which reached only
−0.050 with a converter retune; the difference is the physical parameterization —
PLL natural frequency and damping ratio as coordinates rather than raw gains.

**The cheapest synchronous mitigation is a condenser at bus 30 alone**, 166 MVA
for −0.02 rising to 270 MVA for the base margin, with no PV lost and no
controller change. That the single cheapest location is bus 30 is an independent
corroboration of the descriptive finding in E38, reached by a different route.

**Restoring a whole machine is never Pareto-efficient.** M4 meets the targets and
is dominated everywhere — by M5 at a fraction of the capacity, and by M2 at less
capacity still.

## Does repair move the interaction spectrum away from minus one?

**In proportion to the margin achieved, and not always away.** Median `m4` after
repair, by target: 0.065 (−0.02), 0.075 (−0.05), 0.085 (−0.10), **0.111 (base)**,
against **0.0898 unrepaired**.

A repair that only just clears −0.02 leaves `m4` at 0.019 to 0.065, **below** the
unrepaired 0.0898. That is not a failure of the invariant; it is what the
invariant means. E33 established that `m4` tracks `abs(alpha)`, the distance to
the boundary in either direction. The unrepaired flagship sits 0.145 past the
boundary; a repair to −0.02 sits 0.02 before it, and is therefore *closer to
closure*. Only repairs reaching the base margin push `m4` above the unrepaired
value.

**Required reading of `m4`:** a boundary-proximity invariant, not a stability
score. Its value carries no sign information; the system eigenvalue does.

Sanity check: the base case against itself gives `m4 = 1.0000` exactly, because
with no replacement the interaction operator is identically zero.

## A correction made during the run

The first pass reported M1 as failing every target. The optimizer converges *onto*
the constraint and lands within its own tolerance of it; judging that at `1e-9`
reported a successful repair as a failure. The test is now `1e-6` and the
shortfall is recorded per row. No target, cost axis, model or coordinate was
changed. The verdict moved from "no converter-only solution exists" to "a
converter-only solution exists at every target", which is why the tolerance
mattered enough to fix and re-run.

## Files

`E34_mitigation_frontier.csv` (+ `.parquet`), `E34_pareto_front.csv`,
`E34_Pareto_margin_vs_controller_change.png`,
`E34_Pareto_margin_vs_syncMVA.png`, `E34_repair_closure_before_after.png`,
`E34_figure_source.csv`, `manifest.json`.
