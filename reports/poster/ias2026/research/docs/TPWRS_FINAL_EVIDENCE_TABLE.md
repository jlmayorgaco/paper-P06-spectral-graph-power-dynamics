# TPWRS final evidence table

Sources:

- the frozen F7 results (`IAS2026_TRACKA_F7_POLICY_HYPERGRAPH_FREEZE`);
- the frozen F8–F12 results (`IAS2026_TRACKA_F8_F12_POST_F7_FREEZE`);
- the journal gates: G1 `docs/G1_WHOLE_RHP_COMPOSABILITY.md`, G2
  `docs/G2_TDS_VALIDATION.md`, G3 `docs/G3_IEEE68.md` (preregistered, commit
  `d9fa097e`), G4 `theory/F2E_whole_rhp_and_mechanism_projections.md`, G5
  `docs/G5_NOVELTY_STATEMENT.md`;
- the ledger `docs/CLAIMS.md` (rows cited).

Conventions:

- "RHP" means `Gamma_RHP = {Re s > 0, |s| > 1e-3}`, the planning object.
- "Band" means the 0.3–1.5 Hz inter-area projection.
- "TDS" means the nonlinear phasor-domain simulation, which is **not EMT**.
- Every IEEE-39 result is conditional on the harmonized first-order AVR,
  matched dispatch, core 30/33/35/37, and the leaky Q/V coordinate.

| # | claim | theorem / empirical | IEEE-39 | Kundur | IEEE-68 | full-RHP validated? | TDS validated? | conventional baseline? | limitation | final wording |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | `H` is an antichain of minimal unsafe portfolios; a portfolio containing no hyperedge is stable; `H` is locally constant and changes only at spectral boundaries | **theorem** (F2C Props 1–3; F2E §3) | all maps consistent; 0/1 153 spot errors (F7) | 0 spot errors (F12) | 60/60 spot checks (G3) | yes (F2E: the proofs hold for any `Gamma`) | n/a | clutters and argument principle are standard | heredity fails, so the stable-set family is not a simplicial complex; the converse "contains a hyperedge implies unsafe" is false | "Portfolios free of every hyperedge of `H_RHP` are small-signal stable; `H_RHP` changes only where an eigenvalue crosses the imaginary axis or leaves the origin." |
| 2 | Equal `kappa` can carry different `H`; regions can be non-monotone and disconnected | **theorem** (F2D, one-mode constructions) + empirical | up to 15 `H` behind one `kappa`; the tongue; 30 / 36 / 16 distinct `H_RHP` (F7A/B/C) | 4–5 `H_RHP` over `kappa` in {1, 2} | not observed (map empty; 3 family points) | yes (F7A/C identical, F7B 96.5 %) | the tongue: U-S-U-S along `k = 1.30` | `kappa` alone (a scalar order) does not show it | IEEE-39 only for the rich structure | "The scalar order hides which coalitions fail; the hypergraph does not." |
| 3 | `H_RHP = min(∪ H_m)`; every band hyperedge contains an RHP hyperedge; **`H_Gamma = EMPTY` does not imply stability unless `Gamma` covers every mechanism** | **theorem** (F2E Props 1–3) + counterexamples | undamped condenser: band empty, RHP not | band empty at 556/374 nodes, RHP never empty | 12-plant, `g = 0.25`: only hyperedge is a 2.34 Hz PLL mode | this *is* the RHP result | condenser swing (8.65 Hz) and Kundur aperiodic pair diverge in TDS | — | mechanism projections still useful as explanation | "The planning object is `H_RHP`; mechanism-resolved `H_Gamma` explain transitions and never certify safety." (O82) |
| 4 | `H_RHP` depends on the converter reactive policy at fixed network and dispatch | empirical | yes: every slice | yes: 61/61, 36/36 lines change; never empty | preregistered 4-candidate map: **NOT REPRODUCED** (empty everywhere). 12-plant: `kappa` 6 → 9 → 11 at 3 points | yes | yes: tongue sides, Kundur `g = 0.08` S vs `g = 0.11` U | F10: flagship-only Nyquist, SCR/gSCR and single-port impedance miss the coalition changes | system-dependent; on 68-bus only at ≥ 6/12 plants and at three points, not a map | "Policy-dependent on IEEE-39 and Kundur; on the documented 68-bus system only at high penetration." (O80, O86, N27) |
| 5 | The witness coalition contracts `{30,33,35,37} → {30,33,35} → {30,33}` as excitation strengthens | empirical | yes (O71) | n/a (3 candidates) | n/a | yes (F7A/C identical) | `kappa = 4` (flagship U, triple S) and `kappa = 2` (pair U, singles S) confirmed | per-subset GN reproduces it (it is the same computation, O78) | IEEE-39 only | "The minimal incompatible coalition shrinks as excitation gain rises." |
| 6 | Voltage support is not monotonically stabilizing: an instability tongue with two folds | empirical (+ F2D construction) | yes (O70, N22) | yes (the Kundur `2` boundary: stable at g = 0.08, unstable at g = 0.11) | not observed | yes | yes: U-S-U-S at `g = 0.020/0.042/0.102/0.173`; Kundur sides | non-monotone damping is known and is not claimed as new (G5) | leak-dependent boundary positions (F14) | "Intermediate reactive gains can destabilize a composable fleet; full gain restores it." |
| 7 | Every oscillatory boundary is a zero of the reduced return difference (port closure) | empirical, on the classical GN identity | 29 851/29 851 (F7), 1 511/1 511 (F8C), 128/128 RHP | 473/473 RHP | not tested | yes: 601/601 RHP | n/a | **it is generalized Nyquist on each sub-loop (O78, N23)** | aperiodic crossings excluded (row 8); speed only at a common operating point | "The classical return difference localizes the hypergraph's oscillatory boundaries at reduced port order." |
| 8 | Aperiodic boundaries through the origin are not port-certifiable | theorem-level remark (F2E §3) + empirical | none occur | 194 located | none on the map | yes | Kundur pair diverges along the predicted real mode (+0.77 against +0.685 s^-1) | — | must be certified on the full model | "Aperiodic transitions require the full model; the port test does not apply at s = 0." |
| 9 | Electromagnetic presence at the retired buses restores composability, dose-dependently, **if the condenser damps its own swing mode** | empirical (causal, both directions) | RHP thresholds 0.24 % (P_fold) to 45 % (P1); needs `D = 2` or flux dynamics | not tested | not tested | yes (G1; N25 restates O72) | threshold at P4: TDS 2.47 % against 2.48 % | inertia-based screening does not identify it (row 11) | IEEE-39 only; condenser per-unit services held fixed along the rating | "A condenser with a damped swing mode empties `H_RHP` above a point-dependent rating." (O84) |
| 10 | An undamped low-inertia classical condenser is itself unstable | empirical | 7.7–14.3 Hz swing mode, `Re` up to +0.50 | — | — | only visible in RHP (band blind) | yes: 8.65 Hz, +0.034 s^-1 | — | idealized condenser (D = 0, EMFs frozen) | "Minimal condensers must be damped." (F19) |
| 11 | Synthetic inertia without electromagnetic presence does not restore; surviving-fleet inertia is destabilizing for this mechanism | empirical | yes (O73, O74) | not tested | not tested | yes (those configurations identical) | no | inertia-based screening would rank these wrongly | IEEE-39 only; not a general inertia statement | "Inertia is not the missing service here." |
| 12 | Surviving-AVR speed orders the regions | empirical | yes (O75) | not tested | not tested | **partly**: order holds, but P3 is never RHP-composable (N26) | no | — | slow-AVR end leaves the band (0.21–0.27 Hz) | "AVR speed moves the fleet between regions; it cannot be removed and does not always reach composability." |
| 13 | Excitation heterogeneity is causal at matched mean, but the mean dominates | empirical | yes (O76) | not tested | not tested | not recomputed (F7C slice: `H_RHP = H_IA` 100 %) | no | — | definition-dependent shares | "Heterogeneity is a genuine but secondary cause." |
| 14 | Conventional shortcuts miss the coalition structure; per-subset GN recovers it exactly | empirical | yes (O78) | not tested | not tested | band study | n/a | **this row is the baseline comparison**: GN per sub-loop = ours; flagship Nyquist, modal sensitivity 22/52, additive 30/52, pairwise 38/52, SCR/gSCR, μ-bound fail | GN equivalence means no new stability test | "The contribution is the hypergraph and its behaviour, not the test." (N23) |
| 15 | Localized / port-core evaluation is 22–173x faster | empirical | yes (O79) | — | exact affine assembly used (relative error `3e-16`) | n/a | n/a | full re-solve is the baseline | only at a common operating point (N24) | "Faster only at a common operating point." |
| 16 | The Kundur fleet cannot be made composable by policy | empirical | — | yes: `H_RHP` never empty; pairs diverge aperiodically | — | yes | yes (aperiodic pair; oscillatory boundary sides) | — | Kundur alternative benchmark | "On Kundur no reactive policy makes the fleet composable." (O80) |
| 17 | The documented 68-bus model is reproduced | empirical | — | — | power flow to `5e-5`; 15/15 modes to `5e-4` Hz | n/a | no | the report's own tables | converter rating rule declared (`|S|/0.8`) | "Our 68-bus model reproduces the IEEE PES benchmark report." (O86) |
| 18 | The nonlinear model does what the eigenvalues say | empirical | 26 runs | 6 runs | not run | n/a | **32/32 verdicts**; oscillatory `|df| <= 2.5e-4` Hz, `|d alpha| <= 6e-4` s^-1 | E25/E32 earlier checks | not EMT; small disturbances only; no limits | "Small-signal predictions are confirmed in nonlinear phasor-domain simulation." (O85) |

## Bottom line for the manuscript

- **Supported:**
  - the object and its theory (rows 1–3, 8);
  - policy dependence, witness changes and non-monotone regions on IEEE-39
    and Kundur (rows 2, 4–6);
  - oscillatory port localization (row 7);
  - the IEEE-39 service mechanism with the damped-condenser restatement
    (rows 9–11);
  - time-domain confirmation of every declared case (row 18);
  - a faithful 68-bus reproduction (row 17).
- **Restricted:**
  - policy dependence on the 68-bus system: absent for the preregistered
    candidates, present only at ≥ 6/12 plants (rows 4, 17);
  - the AVR window (row 12);
  - the speed-up (row 15).
- **Not claimed:** generalized Nyquist, the `-1` crossing, hypergraphs, or
  non-monotone damping as new; any EMT, transient-stability or field
  validation.
