# TRANSACTIONS FINAL RESULT AUDIT

Scope: IAS 2026 Track A research line, after the F7 freeze
(`IAS2026_TRACKA_F7_POLICY_HYPERGRAPH_FREEZE`, commit `f12ae3a0`) and the phases
run after it: F8/F8B/F8C (service attribution), F10 (baselines), F11 (scaling),
F12 (second benchmark, preregistered in commit `34a64bdc`). Sources:
`docs/F7_POLICY_HYPERGRAPH.md`, `docs/F8_SERVICE_ATTRIBUTION.md`,
`docs/F10_F11_BASELINES_AND_SCALING.md`, `docs/F12_SECOND_BENCHMARK.md`,
`docs/CLAIMS.md`. The frozen F7 labels were re-solved after all later library
changes: 120 of 120 identical (`results/F7_regression_after_freeze.json`).

Every IEEE-39 statement below is conditional on: the harmonized first-order
excitation model, stabilizer on, matched dispatch (common operating point), core
30/33/35/37, protected band 0.3–1.5 Hz, and the leaky Q/V regulator coordinate
(boundary locations move with its leak, the phenomena do not).

---

## 1. What physical services reshape `H_Gamma`?

- **Electromagnetic presence at the retired buses** — a voltage behind transient
  reactance. Restoring it with a condenser stripped of inertia (1 %), EMF
  dynamics, AVR, PSS, damping and reactive output empties `H` at every unsafe
  point (four of five at 25 % rating, the fifth at full rating); along the
  condenser rating the hyperedges vanish one by one, reversing the F7
  contraction, at 0.25–19.5 % of the retired rating.
- **Excitation dynamics of the surviving fleet.** Slowing every surviving AVR
  moves the same fleet through `kappa = 1 -> 2 -> 3 -> 4 -> inf` before the base
  loses viability; manual excitation is not a viable configuration of this grid
  (aperiodic base instability). The excitation mean dominates the margin;
  heterogeneity at matched mean also changes `H` (causal, secondary, sign
  depending on the parameter and the policy).
- **Converter reactive policy** (F7): a partial, non-monotone substitute for the
  lost voltage stiffness.
- **Surviving-fleet inertia** changes `H` and is **destabilizing** here (x0.5
  empties `H` at four of five unsafe points; x2 worsens or destabilizes the base).
- **Margin-only in almost all cases:** condenser PSS (never changes `H`), and at
  most points condenser inertia, damping, reactive share and EMF dynamics once
  electromagnetic presence is there. Services that change the witness coalition
  itself: condenser inertia/EMF/damping/AVR/reactive share at P1, surviving-fleet
  inertia, damping and PSS at P1–P3 and P_fold.
- **Not a restoring service:** synthetic inertia without electromagnetic
  presence (fixes only the marginal fold point).

## 2. Causal or merely correlated?

**Causal, in the interventional sense, for the services above.** Each was changed
by a physical model intervention at fixed operating point, with everything else
held, and the effect was reproduced in both directions where physically
meaningful: presence restored at unsafe points empties `H`, removed from safe
condenser configurations returns the unsafe `H`, and shrinking the rating
reverses the hyperedge sequence at the same thresholds; surviving inertia acts
the same way whether reduced or increased; AVR speed orders the regions along a
continuous path. The mechanism reading — **loss of local voltage stiffness at the
retired buses interacting with the excitation dynamics of the remaining
machines** — is consistent with every intervention. It is causal within the
model; no field or hardware evidence exists, and the analysis is small-signal
only.

## 3. What does generalized Nyquist already explain?

**All of the per-subset stability information.** `det T_S / det T_0 =
det(I + M[S,S])` is the return difference of the portfolio loop; its winding on
`dGamma` with the open-loop pole correction is a generalized Nyquist count, and
the `-1` crossing is the classical Nyquist/impedance condition. Applied to every
principal sub-loop it recovers `H` exactly (52/52 points, 6/6 same-`kappa`
regions, 36/36 policy-line points) — because it is the same computation. The
`-1` crossing is **not** new and is not claimed. Novelty is explicitly
downgraded to what follows in 4 and 5.

## 4. What new information do `H_Gamma` / `kappa_Gamma` provide?

Not a new stability test — a new **object and its behaviour**:

- the minimal incompatible coalitions themselves, their maximal guaranteed-safe
  complements (planning content without assuming hereditary stability), and the
  partition of physical parameter space they induce, with its local-constancy,
  boundary, fold and topology theory (F2B–F2D);
- the empirical finding that policy and excitation reorganize the coalitions
  (36 hypergraphs over three slices; up to 15 behind one `kappa`; contraction
  `{30,33,35,37} -> {30,33,35} -> {30,33}` robust along every excitation axis and
  reversible by policy), and a mechanism for it (F8);
- what practitioners' methods short of exhaustive per-subset Nyquist miss:
  flagship-only Nyquist (never the coalition), single-port impedance screening
  (singletons only), first-order modal sensitivity (`H` right at 22/52 points,
  3/36 on the policy line), extrapolation from exhaustive single and pair studies
  (38/52; misses both order-3 regions), SCR/gSCR (policy-blind), structured mu
  (never certifies, never identifies).

## 5. What computational advantage does the port formulation provide?

Organizing all `2^n` Nyquist tests as principal minors of **one** operator: two
operator builds per parameter point, no equilibrium solves, all portfolios at
once. At a common operating point, with zero count or `H` errors against the full
model:

| lattice | full model (M1) | localized assembly (M2) | port core (M3) |
|---|---|---|---|
| 16 subsets | 2.9 s / point | 0.018 s (165x) after 7.8 s init | 0.13 s (22x), no init |
| 512 portfolios | 118 s / point | 0.68 s (173x) after 338 s init | 0.99 s (119x), no init |

M3 wins for few points, M2 for dense maps (break-even about 70 and 1 100
points). **Neither advantage exists under re-equilibration** (unity power factor:
naive reuse gives the wrong `H` at 4 of 6 points); no speed-up is claimed there.
Port closure also marks every boundary whatever its physical cause: 29 851 of
29 851 policy/excitation crossings (F7) and 1 511 of 1 511 service-induced
crossings (F8C) port-visible and detected; no port-invisible case occurred.

## 6. Is policy-dependent composability reproduced on a second benchmark?

**Yes, under the preregistered rule, with a strong qualifier.** On Kundur
two-area (documented data, preregistered protocol) reactive policy changes `H` on
every pure-policy line of both slices through port-detected imaginary-axis
crossings, and is non-monotone (`inf -> 2 -> 1 -> 2 -> inf`). But every pair of
Kundur replacements diverges aperiodically outside the band, so no policy makes
the Kundur fleet fully composable, and the band-limited `H` agrees with the
full-RHP `H` at only 0.3–8 % of nodes (post-hoc sensitivity; on IEEE-39 it agrees
at 96.5–100 %). Policy still moves the full-RHP hypergraph on every line.

## 7. Which results are IEEE-39-specific?

- `kappa = 4`, the flagship coalition, the specific witness chain
  `{30,33,35,37} -> {30,33,35} -> {30,33}`, the 36 hypergraphs and their areas;
- the ~0.68 Hz instability tongue and its two folds;
- a policy-reachable fully safe region;
- surviving-fleet inertia being destabilizing; the dose thresholds of
  electromagnetic presence; the sign of the heterogeneity effect;
- agreement of band-limited and full-RHP hypergraphs.

## 8. Which results appear framework-general?

- **Theorems** (F2B, F2C, F2D): local constancy of `kappa` and of the whole `H`;
  boundaries only at spectral boundaries; elementary moves; `H`-free portfolios
  safe without heredity; non-monotone and disconnected regions possible; tongue
  tips are structurally stable folds.
- **Reproduced on two independently parameterized grids:** policy reorganizes
  `H`; voltage regulation is non-monotone for composability; every `H` change is
  an imaginary-axis crossing where the port closure reaches `-1`, apart from
  band-edge re-classifications, which are separable.
- **Structural:** the principal-minor port core is exact wherever the operating
  point is common to all portfolios, and invalid otherwise.

## 9. What is safe for the IAS poster?

Wording that stays inside the evidence:

- "On IEEE-39, with a documented first-order excitation model, the set of
  replacement coalitions that destabilize the inter-area band is not fixed: the
  plant voltage-regulator gain alone moves the same four-unit fleet from a
  four-unit coalition to a single pair and back to safe." (O66)
- "Electromagnetic presence at the retired buses, not inertia, is what restores
  composability: a small condenser with almost no inertia removes the incompatible
  coalitions one by one." (F8)
- "Voltage regulation is not monotonically stabilizing." (N22)
- "Every boundary is a return-difference eigenvalue reaching -1 — a generalized
  Nyquist condition applied to every sub-portfolio; the contribution is the map,
  not the test." (F10)
- Hero figure: `figures/F7A_hypergraph_map.png`, with its model conditions in the
  caption.

Not safe for the poster: any unconditional IEEE-39 statement; "first", "new
stability criterion", "~100x" without "at a common operating point"; any claim
that the Kundur fleet can be made safe; any claim from the band-limited analysis
where outside-band modes exist.

## 10. What is safe for TPWRS?

The above, plus, as a full paper:

- the F2B–F2D theory as a formal framework with the generalized-Nyquist lineage
  stated in the introduction and the novelty placed on the object, its parameter
  dependence and its mechanism;
- F7 with its validation (Safeguards A–D, random direct-path checks, leak
  sensitivity), F8 as the causal mechanism study, F10 as an honest comparison in
  which exhaustive sub-loop Nyquist is shown to be equivalent, F11 as a scoped
  computational result, F12 as a replication with its limitation stated in the
  abstract.

Still required before submission: (i) an IEEE 68-bus/NETS-NYPS replication with
documented dynamic data, because Kundur is a small benchmark that reproduces the
policy dependence but not a safe region; (ii) a full-RHP (or real-axis-inclusive)
protected region reported beside the band on every benchmark; (iii) an EMT or at
least nonlinear time-domain check of the tongue and of the electromagnetic-presence
threshold; (iv) the ANDES cross-check extended to the service interventions. No
TPWRS claim may omit TX3's negative result on the matched three-GFL benchmark
(`CLAIMS.md`, relationship to the frozen TX3 envelope).
