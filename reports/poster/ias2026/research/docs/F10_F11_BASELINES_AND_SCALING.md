# F10 — Existing-method baselines · F11 — Computational scaling

Post-freeze. Code: `experiments/F10_baselines.py`, `experiments/F11_scaling.py`,
`src/ibr_cycles/models/port_core.py` (tested in `tests/test_port_core.py`). Data:
`results/F10/`, `results/F11/`. Same model scope as F7/F8.

## F10 — what conventional methods recover

**The central point first.** For a portfolio `S` the port operator gives
`det T_S / det T_0 = det(I + M(s)[S, S])`, the return difference of the
interconnection between the replaced devices and the rest of the grid. Its
winding on the boundary of the band region, corrected for open-loop device
poles, **is a generalized Nyquist count**, and the `-1` crossing of an eigenvalue
of the return ratio is the classical Nyquist/impedance condition. Nothing here
claims otherwise. Applied to every principal sub-loop, generalized Nyquist
therefore recovers **exactly** the same combinatorial information as the
incompatibility hypergraph — because it *is* the same computation.

Baselines (52 points: the six F8 points, the two largest same-`kappa` regions for
`kappa = 1, 2, 3`, and 40 random base-stable F7A points; plus a 36-point
pure-policy line, F7B `t = 0.852`, `g = 0 ... 0.35`):

| method | what it needs | A first witness | B all minimal subsets | C same `kappa`, different `H` | D policy-driven `kappa` (line) | E boundary frequency |
|---|---|---|---|---|---|---|
| B1 generalized Nyquist, full portfolio only | 2 operators per point | no (says "flagship unstable", correct at 52/52) | no | no: "unstable" in every same-`kappa` region | only flagship yes/no | yes, where the flagship crosses |
| **B2 generalized Nyquist, every sub-loop** (= this project's port closure) | 2 operators per point, all `2^n` minors | **yes** | **yes, 52/52** | **yes, 3/3 pairs, 6/6 regions** | **yes, 36/36** | **yes**, to `1e-7` Hz at every located boundary (F7, F8C) |
| B3 single-port impedance screening | per-bus minor loops | only singleton witnesses (singletons exact 52/52) | 31/52 (only where `H` has no multi-bus hyperedge) | no | — | per bus only |
| B4 modal participation (base inter-area mode) | 1 eigen-analysis | ranking only; not a subset | no | no | no | mode frequency only |
| B5 first-order modal sensitivity (1 % replacement per bus, additive) | `1 + n` full models | `H` exact 22/52, `kappa` 25/52 | 22/52 | 1 of 3 pairs; `H` exact in 1 of 6 regions (predicts `EMPTY` in 5) | `kappa` exact 3/36 | approximate: median error 0.03 Hz, worst 0.21 Hz |
| B5b additive extrapolation from exhaustive single studies | `1 + n` full models | `H` 30/52, `kappa` 34/52 | 30/52 | 2 of 3 pairs; `H` exact in 1 of 6 regions | 5/36 | — |
| B5b pairwise extrapolation from exhaustive single + pair studies | `1 + n + n(n-1)/2` = 11 of 16 lattice members | `H` 38/52, `kappa` 41/52 | 38/52 | 2 of 3 pairs; `H` exact in 3 of 6 regions; predicts `EMPTY` in both `kappa = 3` regions | 27/36; lags the true sequence and declares safety early | — |
| B6 SCR / gSCR | network only | ranking, static | no | **impossible**: identical at every policy point | **impossible** | none |
| B7 structured robustness (D-scaled mu upper bound over `diag(delta_a I_2)`) | 2 operators + optimisation | no | no | no | no | — |

B7 never certified a safe point (0 of the safe points) and never falsely
certified an unsafe one: over the fractional-replacement box the bound is far
above 1 (e.g. 6.3 at P2). It relaxes exactly the combinatorial vertex structure
that `H` keeps, so it cannot recover which subsets fail.

**F. Evaluation count.** Exhaustive direct evaluation needs `2^n` full-model
solves per parameter point (16 here, 512 for the nine-candidate fleet).
Generalized Nyquist on every sub-loop, done through one base operator and the
principal minors, needs **two** operator builds per point and no equilibrium
solves (matched dispatch). Pairwise extrapolation needs 11 of the 16 full solves
and is still wrong on a quarter of the points and on the order-3 structure.

**Downgrade, stated plainly.** Generalized Nyquist applied to every principal
sub-loop provides the same combinatorial information with the same effort as the
port-closure computation, because it is the same computation. The `-1` crossing
is not new, and neither is the per-subset test. What remains specific to this
work is (i) the object — the incompatibility hypergraph and its partition of
physical parameter space, with the local-constancy, boundary and topology
results of F2B–F2D; (ii) the observation that physical policy reorganizes it
(F7), with a mechanism (F8); (iii) organising all `2^n` Nyquist tests as
principal minors of one operator, which is what makes the map cheap (F11). The
methods that practitioners actually use short of that — flagship-only Nyquist,
single-port impedance screening, modal sensitivity, extrapolation from single and
pair studies, SCR/gSCR, mu — do not recover `H`, and none of them tracks the
policy dependence.

## F11 — computational scaling

Three methods, CPU seconds per task (so that parallel execution flatters none),
accuracy always against M1:

- **M1** full DAE build, equilibrium, central-difference Jacobian and eigensolve
  per portfolio and per point;
- **M2** localized `A_red` assembly (three M1 solves per portfolio at
  initialisation; then only the converter rows and each machine's efd row move);
- **M3** reduced port core (one base and one all-replaced operator per point;
  every portfolio a principal minor; counts by the argument principle).

| lattice | points | M1 per point | M2 per point | M2 initialisation | M3 per point | M3 initialisation | count errors M2 / M3 | `H` errors M2 / M3 |
|---|---|---|---|---|---|---|---|---|
| 16 core subsets | 24 | 2.90 s | **0.018 s** (165x) | 7.8 s | 0.13 s (22x) | none | 0 / 0 | 0 / 0 |
| 512 portfolios, nine candidates | 6 | 118 s | **0.68 s** (173x) | 338 s | 0.99 s (119x) | none | 0 / 0 | 0 / 0 |

Other measured quantities: M1 state dimension 72–110; M2 reference storage 1.6 MB
(16) and 64 MB (512); M3 operator dimension 78, action dimension 8 (core) and
18 (fleet), 240–400 contour points, peak memory about 23 MB; M2 band-eigenvalue
error against M1 at most `4.3e-7`; M3 winding residuals zero (integer to
round-off).

**Break-even.** M2 is faster per point but pays `3 x 2^n` full solves up front;
M3 pays nothing up front. M2 overtakes M3 after about 70 points for the 16-subset
core and about 1 100 points for the 512-portfolio fleet. For a handful of
operating points or a single new policy M3 is the cheaper method; for dense maps
on a fixed lattice M2 is.

**Regime.** Both speed-ups exist **only** where the operating point is common to
every portfolio (matched dispatch, frozen equilibrium). Under a re-equilibrating
policy (unity power factor) each portfolio has its own operating point; reusing
the base operator anyway gives the wrong hypergraph at 4 of 6 points
(`results/F11/F11_reequilibrated.csv`). There, only M1 is valid and no speed-up is
claimed. The ~100x figures must never be quoted without "at a common operating
point".
