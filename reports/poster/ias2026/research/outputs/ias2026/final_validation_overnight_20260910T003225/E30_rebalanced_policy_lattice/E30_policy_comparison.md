# E30 — the order-4 mechanism under four reactive policies

**Gate G1: QUALIFIED, not failed.** Re-equilibrated engineering policies do not
destroy the Track-A mechanism entirely — under unity power factor the flagship is
still unstable, and worse. But the **irreducibly fourth-order character is
specific to the reactive-matched dispatch**, and under voltage regulation the
instability does not occur at all.

This is the single most consequential result of the overnight run and it changes
what may be claimed.

Capability constants were frozen in `configs/ias2026/overnight_policies.yaml`
before any policy case was solved. Observable: the frozen v2C inter-area family
envelope, unchanged, with the comparison basis pinned to the 12 flagship
survivors so the overlap threshold means the same thing at every node.

## Result

| policy | class | effective order | `α_base` | `α_exact` | `≤1` | `≤2` | `≤3` | `μ₄` | proper subsets stable | flagship band RHP |
|---|---|---|---|---|---|---|---|---|---|---|
| **P1 Q-matched** | **A** | **4** | −0.1265 | **+0.1447** | −0.2078 | −0.3352 | −0.2884 | **+0.4331** | **yes** | 1 |
| **P2 unity PF** | **B** | **3** | −0.1265 | **+0.4436** | −0.2245 | −0.4514 | **+0.2309** | +0.2127 | **no** (30+33+35 = +0.1407) | 1 |
| **P3 V-reg, 0.95 pf capability** | **D** | — | −0.1265 | **−0.0953** | −0.0473 | −0.1222 | −0.0901 | −0.0053 | yes | **0** |
| **P4 V-reg, S-capability** | **D** | — | −0.1265 | **−0.0953** | −0.0473 | −0.1222 | −0.0901 | −0.0053 | yes | **0** |

Reconstruction residual `0.0` for every policy; all 16 subsets tracked in all
four lattices, zero tracking failures.

## What each policy says

**P1, Q-matched — category A.** The frozen result reproduces exactly: all fifteen
proper subsets stable with zero right-half-plane modes in the band, the portfolio
unstable at `+0.1447` with exactly one, orders 1–3 stable, `μ₄ = +0.4331`,
`κ = 4`. The operating point is identical to base by construction, drift `0.0`.

**P2, unity power factor — category B.** The portfolio is unstable and *more*
so, `+0.4436` against `+0.1447`. The same family carries it, 0.581 Hz. But the
triple **30+33+35 is already unstable at `+0.1407`**, so "every proper subset is
safe" is false here, and the order-3 reconstruction already predicts instability
at `+0.2309`. **The effect is third order under this dispatch, not fourth.**
The operating point genuinely moves: worst voltage drift 0.090 pu, minimum
voltage 0.954.

**P3 and P4, voltage regulation — category D.** The flagship is **stable**,
`α_IA = −0.0953`, zero right-half-plane modes in the band, every subset stable.
`μ₄ = −0.0053`, negative and negligible. Replacement still costs damping
(−0.1265 → −0.0953) but never crosses. The family frequency rises to 0.725 Hz,
the same signature the repair studies report.

P3 and P4 are numerically identical here because **no plant saturates under
either capability rule**: worst reactive utilisation is 0.632 of the 0.95-power-
factor capability and 0.268 of the apparent-power capability. The two rules
differ in principle and did not differ at these operating points. That is a
result about this benchmark, not a general equivalence.

Under voltage regulation the family often has a **single** member rather than two
(family sizes {1, 2}, minimum overlap 0.873): regulating the bus voltage pulls
the second descendant out of the band or below the overlap floor.

## Consequences for what may be claimed

1. **"Irreducible fourth-order interaction" must be stated with its dispatch.**
   It is a property of the Q-matched mechanism-isolation configuration. Under
   unity power factor the same physical family goes unstable at third order;
   under voltage regulation it does not go unstable at all.
2. **"All proper subsets are stable" is policy-dependent.** True under P1, P3 and
   P4. **False under P2**, where the triple 30+33+35 is already unstable.
3. **The engineering headline strengthens.** Deploying these plants in
   voltage-regulating mode with realistic reactive capability removes the
   instability entirely at this operating point — the mitigation does not have to
   be retrofitted, it can be a commissioning requirement.
4. **The worst engineering case is unity power factor**, which is also the
   cheapest and the most common default for battery-free PV. That is the
   industrially relevant warning and it is a *lower-order*, more easily predicted
   failure.

## Recorded corrections

The classifier as first written decided category D on whether damping worsened
relative to base, which labelled the voltage-regulating policies "B" although
their flagship is stable with no right-half-plane mode in the band. That is a
definitional error in the classifier, not a threshold change: category D is
"mechanism disappears", and a portfolio with no instability is exactly that. The
classifier now decides D on the absence of an instability. No threshold, capability
constant, observable or model was altered.

## Files

- `E30_rebalanced_policy_lattice.csv` — all 64 cases with every requested metric
- `E30_policy_summary.csv` — per-policy reconstruction and classification
- `E30_family_boundary_by_Qpolicy.png` — the figure
- `E30_figure_source_panel_a.csv`, `E30_figure_source_panel_b.csv` — exact figure data
- `manifest.json` — config, commit, environment, verdict
