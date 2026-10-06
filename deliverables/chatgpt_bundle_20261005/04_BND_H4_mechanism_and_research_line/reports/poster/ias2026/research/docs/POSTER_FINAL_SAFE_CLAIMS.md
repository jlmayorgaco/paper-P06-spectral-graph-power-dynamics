# Poster-safe claims after the final validation run

Every statement below is backed by a manifest in
`outputs/ias2026/final_validation_overnight_20260910T003225/`. The poster itself
has **not** been edited; this is the permitted vocabulary, not a draft.

## Safe to state

**S1 — the phenomenon, with its dispatch named.**
> With the four replaced plants on fixed reactive setpoints, replacing the
> synchronous machines at buses 30, 33, 35 and 37 with grid-following PV leaves
> every one of the fifteen proper subsets stable and the full portfolio unstable.
> A model retaining interactions only to third order predicts stability; the
> exact fourth-order interaction predicts the crossing.

E30 category A, `μ₄ = +0.433`, orders 1–3 stable at −0.208, −0.335, −0.288 and
the exact portfolio at +0.145. Band right-half-plane count 0 for all fifteen
proper subsets and 1 for the portfolio.

**S2 — the branch-independent version, which is stronger.**
> No proper subset has any unstable mode in the 0.3–1.5 Hz band. The full
> portfolio has exactly one.

An integer count, dependent on no tracking decision.

**S3 — it is not caused by the converter tuning.**
> Over 1500 draws spanning ±20 % in PLL natural frequency, 0.5 to 1.0 in PLL
> damping and ±25 % in three outer bandwidths, the pattern held in every single
> draw.

E36, 1500 of 1500, `α` from +0.080 to +0.209.

**S4 — the operating condition decides whether it happens.**
> Across a held-out envelope of 1004 dispatchable operating points the portfolio
> is unstable at about a quarter of them; when it is unstable, it is an
> irreducible fourth-order failure 87 % of the time.

E35, 0.239 [0.213, 0.267] and 0.871 [0.822, 0.911], fresh seed 20260913.

**S5 — the port-space closure predicts the boundary.**
> The gauge-invariant closure distance reaches its minimum at the operating point
> where the inter-area family crosses the axis, to within a quarter of one
> percent of load, and at the frequency of the crossing mode, to a median of
> 0.005 Hz.

E33, 290 continuation points, Spearman +0.861, p = 1.5e-86.

**S6 — the reduced representation is exact enough to trust.**
> An operator on eight action coordinates reproduces the eigenvalues of the
> 110-state model to 1e-9 over the whole subset lattice and reaches the same
> stability verdict in every case.

E41, 18 cases.

**S7 — the nonlinear check.**
> Nonlinear phasor-domain DAE simulation reproduces the predicted frequency to
> 0.003 Hz and the growth or decay direction in all 25 runs, over four
> disturbances and seven configurations.

E32. Call it **nonlinear phasor-domain DAE time-domain simulation**, never EMT.

**S8 — mitigation, stated as a frontier.**
> At this operating point a converter retune reaches the full base-case margin
> with no photovoltaic energy forgone and no synchronous capacity added, and a
> 166–270 MVA synchronous condenser at a single bus does the same without
> touching the converters.

E34, 10 Pareto-efficient points at four declared margins.

**S9 — condenser robustness.**
> The synchronous condenser restored stability at every one of the 240 unstable
> held-out operating points.

E35, 240 of 240, exact interval [0.985, 1.000].

**S10 — the conventional-baseline comparison, in its corrected form.**
> A model keeping interactions only to third order does not merely rank these
> portfolios badly — under the declared threshold it never predicts instability
> at all, declaring every portfolio stable including the 139 that are not.
> Removed inertia and replaced megawatts rank reasonably well, AUC 0.78–0.87;
> the short-circuit family is weak, 0.59–0.72.

E39, sizes 4, 5 and 6, permutation null 0.50002.

## Must be stated alongside, not omitted

**Q1 — the dispatch caveat.** The fourth-order character belongs to the
reactive-matched configuration. At unity power factor the effect is third order
and more severe; under voltage regulation with realistic reactive capability
there is no instability at this operating point. (E30)

**Q2 — the machine-data caveat.** With ±20 % inertia, ±10 % transient reactance,
±15 % AVR and ±20 % PSS, the full pattern holds in 42 % of draws. The result is
conditional on this machine data. (E37)

**Q3 — the independence caveat.** An independent implementation reproduces the
base inter-area mode but **not** the machine-removal ordering. (E31)

**Q4 — the repair caveat.** The *fixed* converter retune restores stability in
117 of 240 held-out unstable cases: all mild ones, none of the severe. (E35)

## Forbidden wording

- the fourth-order claim **without** its reactive dispatch — Q1
- "robust to machine parameters" — Q2 measures 42 %
- "independently validated" for anything but the base inter-area mode — Q3
- "a converter retune repairs the failure" as a general statement — Q4
- "conventional baselines are anti-predictive" — E39 refutes this for inertia and
  megawatts; only the short-circuit family is weak
- "bus 30 is the collective enabler" — E38's adjusted test does not support it.
  The permitted form is descriptive: *bus 30 appears in every one of the twelve
  genuine inter-area cores while ranking fifth to seventh on individual measures*
- EMT, electromagnetic transient, or any claim of industrial validation
- "mode splitting" — the supported wording is that the base inter-area family
  develops two closely related descendants
- any power law or fitted exponent for `m₄`
- "H5B was refuted" — it remains NON-EXECUTABLE
- any v2B `alpha_Q` boundary number

## Figures cleared for use

| figure | experiment |
|---|---|
| `E30_family_boundary_by_Qpolicy.png` | the dispatch dependence, both panels |
| `E33_closure_tracks_boundary.png` | the main mathematical figure |
| `E33_mu_approaches_minus_one.png` | the interaction eigenvalue reaching −1 |
| `E34_Pareto_margin_vs_controller_change.png` | the industrial frontier |
| `E34_Pareto_margin_vs_syncMVA.png` | condenser sizing |
| `E32_flagship_vs_repair_TDS.png` | nonlinear confirmation |
| `E41_port_vs_full.png` | theory against brute force |
| `E39_ROC_PR.png` | baseline comparison |

`E31_ANDES_vs_internal.png` may be used **only** to show the base-mode agreement,
never the ordering panel, unless the disagreement is stated in the caption.
