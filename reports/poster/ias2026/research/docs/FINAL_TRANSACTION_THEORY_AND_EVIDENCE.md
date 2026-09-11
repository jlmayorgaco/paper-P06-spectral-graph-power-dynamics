# Final theory and evidence (closure campaign)

**Run.** `outputs/ias2026/final_math_nonlinear_validation_20260910T231539/`
(stamp `outputs/ias2026/FINAL_CLOSURE_CURRENT_RUN`).
Branch `ias2026/final-closure-campaign`. The preregistrations are committed
before the runs they govern: 9f581244, a489e4fe, 7e4e792b, d79fb295. The
per-claim table is `results/FINAL_EVIDENCE_TABLE.csv`.

**Status vocabulary.**
- **PROVED**: a mathematical proof exists (it may be classical).
- **INDEPENDENTLY REPRODUCED**: another implementation or code path gives the
  same result.
- **NUMERICALLY VALIDATED**: tested beyond the case it was designed on.
- **BENCHMARK-SPECIFIC**: observed on one frozen model, under its declared
  conditions.
- **FAILED**: the tested hypothesis does not hold.
- **OUT OF MODEL SCOPE**: not decidable in the implemented model.

**Conventions.**
- "Stable" means **transversely** stable: `alpha(A_perp) < 0`, the quotient by
  the exact center subspace (`theory/TRANSVERSE_STABILITY_QUOTIENT.md`).
- "TDS" is the nonlinear phasor-domain DAE, **not EMT**.
- MW is active dispatch (`replaced_pg_mw`) and MVA is rating (`replaced_sn_mva`);
  they are never merged.

**Mandatory model-scope sentence.**

> The frozen IEEE-39 model has no primary frequency restoration and is analysed
> in relative / transverse coordinates. Absolute common-frequency restoration is
> outside this benchmark.

---

## A. PROVED MATHEMATICAL RESULTS

None of these is new mathematics. Each row states what is classical.

| # | statement | proof | classical source |
|---|---|---|---|
| A.1 | **Transverse quotient.** `C = span{R_x, w}` is exactly A-invariant with Jordan block `[[0, omega_B], [0, 0]]` whenever `D = 0` and there is no governor. `det(sI − A) = s^2 det(sI − A_perp)`. The rotation is a gauge; the drift is a physical missing-restoration property. The reduced field satisfies `r(x + phi R_x) = r(x)` and `r(x + c w) = r(x) + c omega_B R_x`, so `y' = Z^T r(x* + Z y)` is the exact nonlinear transverse system. With governors or `D != 0`, `C = span{R_x}`. | `TRANSVERSE_STABILITY_QUOTIENT.md` §2–3 | quotient by an invariant subspace; symmetry reduction |
| A.2 | **Minimal incompatibility hypergraph.** `H` is an antichain. A portfolio containing no hyperedge is stable. `H` is locally constant in the policy and changes only where a transverse eigenvalue crosses the imaginary axis. `H_RHP = min(∪_m H_m)`. A band-limited `H_Gamma = EMPTY` does not imply stability unless `Gamma` covers every mechanism. | F2C Props 1–3; F2E Props 1–3 | clutters; the argument principle |
| A.3 | **Order and structure are distinct.** One-mode constructions give equal `kappa` with different `H`, and non-monotone, disconnected composability regions. | F2D | elementary |
| A.4 | **Planning hierarchy.** C (any order safe) ⇒ B (one-at-a-time safe path) ⇒ A (final state stable); neither converse holds. With the base stable, `S` is C iff it contains no hyperedge. | `FINAL_COMBINATORIAL_THEOREMS.md` §3 | elementary |
| A.5 | **Monotone class.** A common order-preserving Metzler similarity implies that `alpha` is monotone on the lattice, stability is hereditary, and A = B = C. One robust re-stabilizing pair excludes every such realization. | ibid. §2 | Perron–Frobenius (positive systems) |
| A.6 | **Complexity.** Minimum-cardinality Boolean spectral destabilization is NP-complete, already for `D_delta W D_delta − (k−1) I`. | ibid. §1 | a corollary of the clique → sparse-PCA reduction (Magdon-Ismail 2017); robust stability NP-hard (Poljak–Rohn 1993; Nemirovskii 1993) |
| A.7 | **Resilience complex.** `K(rho)` is a simplicial complex. Its minimal non-faces are `H_NL(rho)`. Each `|H| = k ≥ 2` hyperedge induces the boundary of `Delta^{k−1}`, a local `S^{k−2}`. `K` is a decreasing filtration in `rho`, `R_k` is non-increasing, and `H_NL(0+) = H_RHP_perp`. | ibid. §4 | Stanley–Reisner; sublevel filtrations |
| A.8 | **Principal-minor structure.** `det(I + M_SS) = Σ det M_UU`. The Möbius coefficients are the block-touching minor sums. These sums, `det(I + Q_SS)` and the connected (log) coefficients are invariant under port-basis change; individual block-splitting minors are not. The interaction degree `d*` and `kappa` bound neither each other. | `PRINCIPAL_MINOR_PORTFOLIO_STRUCTURE.md` §1–5 | multilinearity; Möbius inversion |
| A.9 | **Zero-frequency port closure.** Two rank-one relocations move exactly the two center eigenvalues, so `sigma(A##) = sigma(A_perp) ∪ {−beta, −beta2}`. `det T##(0) = det g_z det A## / det f_x##`. Its sign changes iff the parity of the positive real transverse eigenvalues changes. | `SYMMETRY_DEFLECTED_PORT_CLOSURE.md` §2 | Brauer (1952); Schur determinant |

---

## B. NUMERICALLY VALIDATED GENERAL MECHANISMS

| # | result | evidence | status |
|---|---|---|---|
| B.1 | The transverse structure holds in every audited case: `dim C = 2`, Jordan `J_12 = omega_B`, invariance ≤ 1.1e-11, spectral identity ≤ 3.1e-8; the old count equals `N_RHP_perp` in 106/106. | FC01 (IEEE-39, Kundur, IEEE-68) | NUMERICALLY VALIDATED |
| B.2 | Transverse re-audit of all **345 229** F7 points: **0** label changes at resolved points. 5 128 remain unresolved (1.5 %): an oscillatory pair within the assembly error, never the old disc. After direct-path recomputation, 16 unresolved points remain at `g ≤ 0.005`. | FC01, `results/TSQ_ieee39_reaudit.csv` | NUMERICALLY VALIDATED |
| B.3 | Zero-frequency port closure, frozen before the holdout. Kundur holdout **28/29**; the miss is a real-pair coalescence (3 → 1 real, 0 → 2 complex), so 28/28 zero crossings. **0/445** false positives (109 oscillatory, 336 no-crossing), 0 beta-inconsistent, identity ≤ 8.8e-12. Synthetic 200/200 and 0/200. The BC01 194/194 is method development only. | FC02, `results/zero_frequency_port_validation.csv`, FIG6 | NUMERICALLY VALIDATED (Kundur + synthetic) |
| B.4 | Principal-minor identities hold to 1.3e-13. The full-order Möbius coefficient is non-zero on all three benchmarks: no exact low-degree representation exists. | BC02b | NUMERICALLY VALIDATED |
| B.5 | The oscillatory boundaries are zeros of the reduced return difference: 29 851/29 851 (F7), 601/601 RHP. This is generalized Nyquist per sub-loop. | F7/F8C/G1 | NUMERICALLY VALIDATED; classical |
| B.6 | **Second-order DAE curvature** (M0 / M1 / M2): the full second-order model is 170–3400× more accurate than linear; without the KCL curvature the frequency prediction is about 10³× worse than linear. Details in the subsection below. | FC07, FIG7 | NUMERICALLY VALIDATED |
| B.7 | **Hopf type** on the exact transverse system, with transversality checked and two step sizes: flagship boundary on the P4 line (`g* = 0.2077`, 0.706 Hz) `l1 = +0.0095`; tongue boundary `k = 1.30` (`g* = 0.0509`, 0.639 Hz) `l1 = +0.0250`; damped-condenser boundary (rating 2.483 %, 0.621 Hz) `l1 = +0.026` (spread 0.2 %; nearest other mode `abs(Re)` = 1.4e-3, so the normal form holds only in a small neighbourhood). **All subcritical.** | FC13 | NUMERICALLY VALIDATED (diagnostic) |
| B.8 | The resilience-complex statements hold on the measured thresholds in 6/6 (point, family) cases: downward closed, non-faces = minimal non-tolerating sets, local spheres, filtration, `R_k` monotone, `H_NL(0+) = H_RHP_perp`. | FC06 | NUMERICALLY VALIDATED (A.7 is a tautology of the definition) |
| B.9 | Small-disturbance agreement between the TDS and the eigenvalues: 32/32 verdicts; `abs(df) <= 2.5e-4` Hz; `abs(d alpha) <= 6e-4` s^-1. | G2 | NUMERICALLY VALIDATED |
| B.10 | **Certificates.** Generic sampled-curvature Lyapunov bound: `r_cert / r_TDS` = 1e-20 to 1e-11 (`kappa(P)` about 1e9). Energy certificate: NOT_APPLICABLE (no energy function for the GFL controls on a lossy network). BC03: sampled screening only. Small gain `rho(Rbar)` = 3.76 (P4), 12.3 (P2), 2.95 (P_inf) > 1 fails. Top-sum screening gives `kappa ≥ 2` at all three points (exact 4, 2, inf). | FC08, FC09 | **FAILED** as useful certificates |

### B.6 Second-order KCL-curvature audit (FC07)

**Setup.**
- Five IEEE-39 cases: BASE, single 30, triple 30+33+35, P4 flagship
  (transversely unstable) and P_inf flagship.
- Seven D1 load-pulse amplitudes from 2 to 128 MW; reference is the full DAE
  (BDF, 1e-10).
- The models are compared on the COI frequency, `abs(V_16)` and the bus-30
  current.
- **Run history.** The first run (v1, kept in `v1_foh_bug/`) had a
  first-order-hold coefficient error (1/dt too large) and one mismatched
  pulse-end sample. Both were found before any interpretation, fixed, and the
  run repeated (v2). Only v2 is reported.

**Error relative to the linear model M0, at the smallest amplitude (2 MW).**

| model | frequency | `abs(V_16)` | current |
|---|---|---|---|
| M1: full second order | **170–450× smaller** | **750–3400× smaller** | 60–490× smaller |
| M2: KCL curvature `F_z D^2 psi` dropped | **570–1700× larger** | 0.1–0.5× (below M0 at small eps, worse above 8–16 MW) | 1–60× larger |
| M2t: M2 with the center-subspace component of x1 removed | 2–3× larger | ≈ 1× | 0.9–1.4× |

**Slopes on 2–16 MW.**
- M0: `p` = 2.00 in all observables and all cases.
- M1, voltage: `p` = 3.0–3.2.
- M1, current: `p` = 2.5–3.0.
- M1, frequency: 1.9–2.8. A discretization floor (about 0.4 % of the
  second-order term at `dt = 1e-3`) dominates below 8 MW.
- The flagship points at 128 MW lie outside the valid region: every model
  fails there alike, so they are not used.

**Symmetry check.** Evaluating the full forcing on the transverse part of x1
leaves x2 unchanged to 2–7e-4 relative, which is the FD accuracy. This
confirms `D^2 r[zeta, c] = 0` for `c in C`. Without the KCL term, x2 is 48–228×
too large. That is the broken rotation / drift symmetry acting on the secular
common-angle drift of the `D = 0` benchmark. With the drift removed (M2t), x2
is still wrong (0.29–0.69 of the true size) and the model is **no better than
linear**. `‖F_z D^2 psi‖ / ‖D^2 F‖` has a median of 1.00–1.09: the two terms
nearly cancel.

**Status: NUMERICALLY VALIDATED** (IEEE-39, five cases). This is a consistency
requirement. It is not new (Taylor), and it is not a headline.

**Allowed wording:** "A second-order reduced SG–IBR response is consistent only
if the algebraic-manifold (KCL) curvature is retained. Dropping it breaks the
exact rotation / frequency-drift symmetry at second order and makes the
second-order frequency prediction about 10³ times worse than the linear one.
With that artifact removed, the device curvature alone is still no better than
the linear model."

---

## C. IEEE-39 BENCHMARK RESULTS

These results are conditional on:
- the harmonized first-order AVR (declared MODEL_CHANGED);
- matched dispatch and core buses 30/33/35/37;
- the leaky Q/V converter coordinate;
- `D = 0`, no governor, transverse stability.

Policy coordinates: `g` Q/V gain; `k` excitation-gain scale; `t` excitation
time-constant scale; `h` heterogeneity.

| # | result | numbers | status |
|---|---|---|---|
| C.1 | **Policy-dependent minimal incompatibility.** At fixed network and dispatch, the converter reactive policy changes *which* coalitions fail, not only whether. | 30 / 36 / 16 distinct `H` on F7A/B/C (345 229 points, 0 label changes under A.1); up to 15 `H` behind one `kappa`; witness contraction `{30,33,35,37} → {30,33,35} → {30,33}` as excitation strengthens; the tongue U-S-U-S on `k = 1.30`, also seen in TDS | BENCHMARK-SPECIFIC, NUMERICALLY VALIDATED |
| C.2 | **Flagship.** 2096.6 MW of active dispatch on 4270.7 MVA of converter rating. | at P4 (`g = 0.03625`, `k = 1.425`): `alpha_perp` = +0.127, `H_RHP_perp` = {30+33+35+37}, `kappa` = 4; at `P_inf` (`g = 1`, `k = 0.5`): EMPTY | BENCHMARK-SPECIFIC |
| C.3 | **The failing mode is the inter-area family, not dispatch.** E14 N6 redone with Pg-matched controls and the frozen modal-family definition. | 25/25 controls tracked, family `alpha <= −0.114` (damping ≥ 3.9 %); 10/10 failing portfolios tracked, 8/10 with that family unstable | BENCHMARK-SPECIFIC |
| C.4 | **Service mechanism.** A condenser with a damped swing mode empties `H_RHP_perp` above a point-dependent rating. An undamped condenser is itself unstable. Synthetic inertia alone does not restore composability. Track A is not controller-caused but is controller-repairable. A fixed retune has bounded authority. | P4 threshold 2.48 % (106.0 MVA) spectral, 2.47 % TDS; undamped swing 7.7–14.3 Hz; E35 117/240 | BENCHMARK-SPECIFIC, NUMERICALLY VALIDATED (TDS) |
| C.5 | **Robustness to documented primary frequency control** (`ieee39_governed_documented_v1`, TGOV1N from the same workbook, untuned). Policy-dependent incompatibility **survives**; the P4 `kappa = 4` coalition **does not**. | See the list below this table. | BENCHMARK-SPECIFIC; conditional statement |
| C.6 | **Monotone class excluded.** Adding bus 34 re-stabilizes an unstable portfolio. | 3 robust counterexamples, e.g. `alpha(30+33+35+37)` = +0.145 → `alpha(+34)` = −0.237 (E12 census policy) | NUMERICALLY VALIDATED |
| C.7 | **Nonlinear composability within model scope** equals the spectral structure. | See the list below this table. | NUMERICALLY VALIDATED **negative** (criteria A/B not met) |

**C.5 governed replication (FC03).**
- F8 points: `H` survives at P3 unchanged (`30+33+35|30+33+37`). It becomes
  EMPTY at P1, P2, P4 and P_fold.
- P4 flagship: `alpha` = +0.127 → **−0.0745**. The flagship stays the
  least-damped portfolio but no longer crosses.
- P4 line: `kappa = 4` only for `g` in [0.016, 0.023] (frozen [0.016, 0.20]).
  The tongue on `k = 1.30` disappears.
- Plane (399 points, `t = 1.5`, `h = 1`): 10 distinct `H` governed (frozen 26).
  - 74 non-empty points: `kappa` = 4 (24), 3 (32), 2 (18), and none of order 1.
  - The governed base is never unstable (frozen: 21 base-unstable points).
- Condenser: governed EMPTY at every rating. Frozen threshold between 2.48 %
  and 3 %, consistent with G1.

**C.7 nonlinear composability within model scope (FC05/06).** Three frozen
disturbance families are used: D1 (bus-16 load pulse, MW), D2 (SG36 `Pm` dip,
MW) and D3 (PV dip, fraction).
- On `[0, rho_scope)`, `kappa_NL = kappa_RHP_perp` and
  `H_NL = H_RHP_perp` for all six (point, family) cases.
- `rho_scope` at P4: D1 130.5 MW, D2 97.7 MW, D3 0.091.
- `rho_scope` at P_inf: D1 88.3 MW, D2 360.8 MW, D3 0.308.
- Of 96 thresholds:
  - 84 are OUTSIDE_MODEL_SCOPE;
  - 8 are NONE_IN_RANGE;
  - 2 are VACUOUS (the BASE in D3);
  - 2 are FAILS_DECLARED_SECURITY (the transversely unstable P4 flagship).
- **No transversely stable portfolio fails inside the model's validity.**
- The binding guards are the GFL terminal-voltage band [0.9, 1.1] pu and the
  documented PSS output limit.

**C.8 — first-event (limiter-activation) complex.** This is a scope
diagnostic, not composability. The order in which portfolios first leave the
declared envelope as `rho` grows is itself a structured hypergraph:
- at P4, D1: `{30+33+35+37}` until 130 MW, then `30+33+35`, then
  `30+33+35|33+35+37`, … `33|35`, and BASE at 681 MW;
- at P_inf, D1: EMPTY until 88 MW, then the flagship first.

It says which coalitions first need the omitted limiters. It says nothing about
their stability beyond that point.

---

## D. KUNDUR / IEEE-68 CORROBORATION

**Kundur (alternative benchmark, same code; not an independent
implementation).**
- Policy dependence: 61/61 and 36/36 lines change `H_RHP`, which is never
  empty.
- Aperiodic Q/V-integrator crossings: 194 located; `H_phys = H_G1` at
  5 917/5 917 points.
- The TDS confirms the aperiodic pair: +0.77 against +0.685 s^-1 predicted.
- The port-closure holdout (B.3) lives here.

**IEEE-68 (documented PES benchmark, reproduced).**
- Power flow reproduced to 5e-5; 15/15 modes to 5e-4 Hz.
- The preregistered 4-candidate policy map is **EMPTY everywhere: NOT
  REPRODUCED.**
- At ≥ 6/12 plants, `kappa` goes 6 → 9 → 11 at three policy points.
- No real transverse crossing, so no port-closure case.
- No Hessian or second-order claims (NL02 accuracy about 1e-3).

**Independent implementation (ANDES, FC14 and F1/E31).**
- Equation-equivalent pieces only:
  - IEEE-39 power flow;
  - the base inter-area mode (within 4.1 %);
  - the direction of the governor effect: TGOV1N moves the rightmost ANDES
    inter-area pair from −0.222 to −0.770 s^-1 at 0.62–0.64 Hz, matching FC03
    (frozen base −0.144 → governed −0.239).
- ANDES does **not** implement the L0 GFL. No portfolio, policy or nonlinear
  GFL result is ANDES-validated.
- The documented IEEEX1 is itself unstable in ANDES (six real modes near
  +1.03), which is why the frozen model uses the harmonized AVR.
- The E31 machine-removal ordering was not reproduced (C4 FAILED).

---

## E. ENGINEERING / PLANNING RESULTS

Objectives are active dispatch [MW] and synchronous support [MVA], on separate
axes.

| # | result | numbers | status |
|---|---|---|---|
| E.1 | P1 spectral any-order-safe plan on the 16-subset core | P4: `30+33+35`, 1775.1 MW (flagship excluded). P_inf: the flagship, 2096.6 MW | BENCHMARK-SPECIFIC |
| E.2 | P1 on the 512-portfolio census (E12 policy) | optimum `31+32+34+35+37+38`, 3652.3 MW / 6499.9 MVA; it is A, B and C, so the order constraint does not bind the optimum | BENCHMARK-SPECIFIC |
| E.3 | **Sequencing matters for 3 of 327 stable targets.** They are final-stable, have a safe one-at-a-time path, but are not any-order safe. There is no final-stable target without a safe path. | e.g. `30+32+33+34+35` (3008.1 MW / 5224.4 MVA). Its only unstable subset is `30+32+33+35`, so exactly the orders that convert 34 last are unsafe. The same pattern holds in all three cases: the unstable subset is the four-set without 34 | BENCHMARK-SPECIFIC (A.4 realized) |
| E.4 | P2 nonlinear-safe planning (`rho*` preregistered: D1 200 MW, D2 200 MW, D3 0.35) | P4: all three families OUT_OF_MODEL_SCOPE (`rho* >= rho_scope`). P_inf D2: evaluated, **identical** to P1 (flagship). The decision is unchanged wherever decidable | criterion C **not met** |
| E.5 | P3 minimum damped-condenser support for the P4 flagship | transverse stability: 2.48 % = **106.0 MVA**. Staying inside the declared envelope at `rho*`: D1 18.6 % = **792.4 MVA**; D2 16.0 % = **684.0 MVA**; D3 not reached at 25 % (1067.7 MVA). Every sub-threshold failure is an OUTSIDE_MODEL_SCOPE guard event | BENCHMARK-SPECIFIC; **downgraded (rule 24)**: an envelope / ride-through requirement, not a stability requirement |
| E.6 | Documented governors remove the P4 support need | governed P4 flagship stable with no condenser; governed EMPTY at every rating | conditional on the governor model |
| E.7 | The reactive policy that repairs spectral composability narrows one envelope | P_inf vs P4 flagship `rho_scope`: D1 88.3 vs 130.5 MW (narrower); D2 360.8 vs 97.7 MW and D3 0.308 vs 0.091 (wider) | BENCHMARK-SPECIFIC, envelope-conditioned |

---

## F. FAILED HYPOTHESES

| hypothesis | result | source |
|---|---|---|
| simple-cycle magnitude or holonomy explains the failure | falsified; the Pg rerun gives ratio 0.72 | E18, UC02 |
| bus 30 is the collective enabler | not supported; p = 0.18 after adjustment | E38, UC02 |
| replaced MW (or MVA, inertia) identifies the minimum failing coalition | no: true-dispatch AUC 0.57 (p = 0.46) | UC02 |
| a band-limited certificate suffices | no: `H_Gamma` empty with `H_RHP` non-empty (condenser, Kundur, 12-plant) | F2E, G1 |
| the plain port detects zero-frequency crossings | 24/194 | BC01 |
| a low-degree Boolean representation is exact | no: `d* = m` on all three benchmarks | BC02b |
| the reduced `A(delta)` is affine | no: rational, full degree; T4 NOT_APPLICABLE | BC02 |
| **nonlinear composability is stricter than spectral inside model scope** (`kappa_NL < kappa_RHP` or `H_NL != H_RHP`) | **not observed** in 6/6 cases; every finite threshold of a stable portfolio is a scope event | FC05/06 |
| a nonlinear / energy certificate predicts a useful fraction of the threshold | no: ratio ≤ 1e-11; energy NOT_APPLICABLE | FC08 |
| low-order principal-minor information certifies all subsets | no: SG `rho` = 3.76–12.3; top-sum screening only `kappa ≥ 2`; sampled only | FC09 |
| IEEE-39 is in the monotone (Metzler) tractable class | no: 3 robust counterexamples | FC10 |
| the IEEE-39 `kappa = 4` coalition is robust to primary frequency control | no: governed P4 flagship `alpha` = −0.0745 | FC03 |
| policy dependence on IEEE-68 (preregistered candidates) | not reproduced | G3 |
| the structure is universal | no | G3, FC03 |
| the undamped condenser restores composability | no: its own swing mode is unstable | G1 |
| a fixed converter retune repairs the failure in general | bounded authority, 117/240 | E35 |
| ANDES reproduces the machine-removal ordering | no (C4 FAILED) | E31 |

---

## G. MODEL-SCOPE LIMITATIONS

1. **Phasor-domain DAE, not EMT.** L0 has no current limiter, PLL or
   ride-through logic, anti-windup, DC link, modulation delay or ZIP loads.
   The guards only detect when an omitted limiter would activate. Beyond them
   the model is not valid, and thresholds there are OUTSIDE_MODEL_SCOPE.
2. **No primary frequency restoration** (`D = 0`, no governor). All stability
   is transverse. The documented-governor variant changes the P4 conclusion
   (C.5).
3. **Harmonized first-order AVR** (MODEL_CHANGED). The documented IEEEX1 is
   unstable in ANDES.
4. **No PV nameplate.** `P_ref` is dispatch with an unlimited DC source. Only
   dispatch or documented SG Pmax can be MW objectives. Kundur and IEEE-68 have
   no usable Pmax.
5. **Numerical, not interval, guarantees.** `d_axis` is a numerical criterion.
   Curvature is sampled (no AD or interval derivatives). 1.5 % of F7 points
   are unresolved.
6. **Three disturbance families**, each at one location. `rho*` is
   preregistered, not optimized.
7. **IEEE-68 derivative accuracy** is about 1e-3, so there are no second-order
   claims there.
8. **ANDES does not implement the L0 GFL.** The independent validation covers
   the network, the SG/exciter modal reference and the governor direction only.
9. **Kundur and IEEE-68 run through the same code.** They corroborate
   benchmark generality, not implementation correctness.

---

## H. WHAT IS ACTUALLY NOVEL (novelty audit)

**Comparison against the required method families.**

| family | what it already gives | relation to this work |
|---|---|---|
| generalized Nyquist / return ratio | per-configuration stability test; the −1 crossing | our labels are GN per subset (B.5); **no new test** |
| impedance stability | single-port criteria | single-port / flagship-only criteria miss the coalition changes (F10); not superseded, just a different question |
| passivity / small gain | sufficient all-configuration certificates | our SG/Perron bound *is* small gain, and it fails on the benchmark (B.10) |
| structured singular value / robust stability | worst case over structured sets; NP-hard in general | the `mu` bound fails (F10). Boolean robust stability is NP-hard (Poljak–Rohn), so A.6 adds nothing |
| sparse / structured stability radius | minimum structured perturbation to instability (Katewa–Pasqualetti 2020) | continuous perturbations of given structure; our object is Boolean replacement over a physical policy family. A.6 is the sparse-PCA corollary: **no priority** |
| minimal cut sets / reliability, N-k | minimal failing sets on hypergraphs (N-k defective k-sets) | the *static* analogue of `H`. Ours is dynamic, transverse-spectral and policy-dependent. The hypergraph itself is not new |
| hypergraph methods | clutters, minimal transversals | not new |
| continuation / bifurcation | boundaries, folds, Hopf type | our boundaries are Hopf (subcritical) and the tongue is a fold pair; standard tools, no new theory |
| transient energy / Lyapunov | energy-function margins | no energy function for the GFL model; the generic bound is useless (B.10) |
| nonlinear ROA methods | basin estimates (SOS, sampling) | our `r_S` is a TDS bisection radius on declared families, not an ROA; scope-censored |

**Verdict on the six candidate novelties of the brief.**

1. **Policy-dependent minimal dynamic incompatibility structure over physical
   binary SG→IBR portfolios.** *Retained: the primary contribution,
   application-level.*
   - The test is classical.
   - The object `H_RHP_perp(theta)` (the minimal transverse-unstable
     replacement coalitions as a function of converter policy) and its
     measured behaviour are not found in the audited literature:
     policy-dependence, witness contraction, equal-`kappa`-different-`H`,
     non-monotone regions, the damped-service mechanism.
   - The behaviour is established on IEEE-39 and Kundur. It is conditional on
     `D = 0` and no governor (weakened, not removed, with governors). It is
     absent on the preregistered IEEE-68 map.
2. **Disturbance-dependent nonlinear composability** (`H_NL`, `kappa_NL`,
   `R_k`). *Definitions retained as a theoretical extension; no benchmark
   headline.* Inside model scope `H_NL = H_RHP_perp` everywhere (hard criteria
   A and B fail), and D fails.
3. **Coalitions ↔ resilience complex.** *Not new mathematics.* It is
   Stanley–Reisner and a filtration. It is a correct organizing statement,
   verified in 6/6 cases.
4. **Cardinality-constrained binary destabilization and a tractable positive
   subclass.** *Not novel.* It is the sparse-PCA clique reduction plus
   Perron–Frobenius. It is useful only to explain why exhaustive witness
   search is needed and why IEEE-39 is outside the easy class.
5. **Symmetry-deflated port closure with zero-frequency crossings.** *Modest
   method-level novelty.* The ingredients are classical. The exact application
   to the structural center subspace of the phasor DAE, validated on a
   preregistered holdout (28/29, 0/445), is the new part. It is relevant to
   aperiodic (Kundur-type) crossings only.
6. **Conversion into safe-sequencing and planning constraints.** *Elementary
   (A.4), with benchmark-specific content.* The content is 3/327 order-dependent
   targets, and the P3 envelope gap (106 vs 684–792 MVA), downgraded to an
   envelope requirement.

**Headline gates (§23).**
- Nonlinear composability as a main contribution: **not met** (A–D all fail
  inside scope).
- Complexity as a headline: **not met** (literature equivalent).
- Certificate as a headline: **not met** (no all-subset guarantee, no pruning,
  screening only).
