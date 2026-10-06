# F7 — The policy-dependent spectral incompatibility hypergraph

**Verdict.** Physical policy does not merely move the scalar order of
non-composability. On all three declared slices it reorganizes the **minimal
spectral incompatibility hypergraph** `H_Gamma` of the same installed fleet.
36 distinct hypergraphs occur, many sharing a `kappa`, and every change is carried
by an identifiable witness coalition crossing the imaginary axis, except 128
band-edge re-classifications that are reported separately. The contraction
`{30,33,35,37} -> {30,33,35} -> {30,33}` is robust along every excitation axis,
and the plant voltage-regulator gain alone drives the same fleet through
`kappa = 4 -> 3 -> 2 -> 3 -> 4 -> infinity`. Voltage regulation is **not**
monotonically stabilizing.

Theory: `theory/F2C_incompatibility_hypergraph.md`,
`theory/F2D_topology_of_composability_regions.md`. Code:
`experiments/_f7_common.py`, `F7_policy_hypergraph.py`, `F7_tongue.py`,
`F7_report.py`, `F7_figures.py`, `F7_safeguard_A_audit.py`,
`F7_safeguard_A_maps.py`, `F7_leak_sensitivity.py`. Data: `results/F7/`.
Figures: `figures/F7A_hypergraph_map.*`, `figures/F7_three_slices.*`.

Model and scope, stated once for every number below: IEEE-39 with the harmonized
first-order excitation model, stabilizer on, matched dispatch (the operating point
is identical for every subset and every parameter value), the four-bus core
`A = {30, 33, 35, 37}`, tested family its power set,
`Gamma = { Re s > 0, 0.3 <= f <= 1.5 Hz }`. `H = EMPTY` means "no coalition of the
core", never "no coalition".

## 1. Coordinates and slices

| symbol | meaning | range |
|---|---|---|
| `g` | reactive-policy gain: leaky Q/V regulator `q = q_ref + g kp_v e + x_v`, `x_v' = g ki_v e - w (x_v - q_ref)`, `kp_v = 2`, `ki_v = 20`, `w = 0.05 rad/s`. `g = 0` is exact fixed-Q (matched), `g = 1` the E30 voltage regulator up to the leak | 0 to 1 |
| `k` | excitation-gain scale on every machine's own `KA` | 0.5 to 2.3 |
| `t` | excitation time-constant scale on every machine's own `TE` | 0.5 to 3.0 |
| `h` | heterogeneity amplitude: `TE_i(h) = G (TE_i/G)^h`, `G = 0.565 s` the geometric mean; `h = 1` native, `h = 0` uniform | 0 to 2 |

| slice | plane | fixed |
|---|---|---|
| F7A | `g x k` | `t = 1.5`, `h = 1` (extends the F2 line, which is its `g = 0` edge) |
| F7B | `g x t` | `k = 1`, `h = 1` |
| F7C | `g x h` | `k = 1`, `t = 1` |

No 4-D Cartesian grid was formed.

## 2. Safeguards

**A — regular policy path.** The naive coordinate (PI gains scaled to zero) is
singular at `g = 0`: four exact zero eigenvalues, one decoupled marginal
integrator per converter, so `g = 0` would have to be excluded. The leaky
coordinate is regular on `[0, infinity)`:

| check | answer | evidence |
|---|---|---|
| state dimension constant? | yes | per subset constant over all 1200 direct spot-check points; flagship 86 states at every `g` |
| `g_z` regular? | yes, and parameter-independent | relative spread of `cond(g_z)` below `8e-9` on every map; `2.79e3` for the flagship |
| zero/marginal artificial state absent? | yes | at `g = 0` the spectrum is matched-plus-`(-0.05)x4` to `2.4e-6`; at all 1 200 direct map points every subset has exactly two eigenvalues below `1e-3`, the structural angle-reference pair (`results/F7/F7_safeguard_A_maps.json`) |
| `A_red` continuous? | yes, affine in `g` | residual `5.4e-11` |
| baseline count unchanged for structural reasons? | yes | the base case has no converter; its matrix is bit-identical across `g` |

**Leak sensitivity — boundary locations move, phenomena do not.** The leak `w`
sets the regulator's finite DC droop and is a modelling constant. 360 map points
(half adjacent to a boundary, half interior, `g > 0`) re-labelled on the direct
path: agreement with `w = 0.05` is 89 % at `w = 0.02` (interior 99.4 %, boundary
79 %) and 72 % at `w = 0.2` (interior 88 %, boundary 56 %)
(`results/F7/F7_leak_sensitivity.json`). So boundary positions and area
fractions are **conditional on `w`**. Every qualitative F7 result survives a
tenfold range, `w = 0.02, 0.05, 0.2` (`F7_leak_qualitative.py`): the F7A
tongue line `k = 1.30` gives `kappa 3 -> 4 -> inf -> 4 -> inf` at all three; the
pure-policy `4 -> 3 -> 2 -> 3 -> ... -> inf` paths of F7B (`t = 0.852`) and F7C
(`h = 0.5`) persist at all three; F7B `t = 0.68` gives `2 -> 1 -> 2 -> inf` at all
three; the coarse F7A grid carries 23, 25 and 24 distinct hypergraphs.

**B — which part of dGamma.** Every crossing is matched across its bisection
bracket and labelled by the constraint of `Gamma` its safe-side partner violates;
the direct right-half-plane count on both sides is recorded as an independent
check. Only `IMAGINARY_AXIS_*` events are called dynamic-stability composability
transitions.

**C — tongues and islands.** Every cell is sampled at corners, edge midpoints and
centre; it is homogeneous only if all nine share `H`, `kappa` and the full vector
of 16 protected counts, no subset has a band mode within `0.01 s^-1` of the axis,
and no unstable mode is within `0.02 Hz` of a band edge; otherwise it is split to
depth 4 (finest cell `1/16` of a coarse cell, sampled at `1/32`). Hanging nodes
reopen contradicted leaves. The map is evaluated on an **exact parametric
assembly** of `A_red` (equilibrium, `g_z` and state dimension are fixed along
every coordinate, so only the converter rows, affine in `g`, and each machine's
efd row, `(K_i/T_i) u_i + (1/T_i) w_i`, move). Against the direct path it agrees
to `1.6e-9` in the matrix and `3e-8` on band eigenvalues.

**D — the narrow region near 0.68 Hz.** Section 6.

### Random-map validation (pass rule frozen before the run: 95 % Clopper–Pearson upper bound of the label error in homogeneous leaves below 2 %)

400 random points per map (half uniform in `g`, half uniform in `sqrt g`),
evaluated on the **direct** path, independent of both the fast assembly and the
refinement.

| slice | map points | spot points in homogeneous leaves | errors | 95 % upper bound | in boundary cells | nearest-node errors (all 400) | fast vs direct label disagreements | ISLAND_OR_TONGUE cells | verdict |
|---|---|---|---|---|---|---|---|---|---|
| F7A | 128 005 | 376 | **0** | 0.79 % | 24 | 1 | 0 | 0 | **PASS** |
| F7B | 88 163 | 386 | **0** | 0.77 % | 14 | 0 | 0 | 1 | **PASS** |
| F7C | 145 811 | 391 | **0** | 0.76 % | 9 | 1 | 0 | 0 | **PASS** |

The one ISLAND_OR_TONGUE cell (F7B, `g = 0.0094`, `t = 2.75`) is a single node on
the band-edge sliver of §5, not a separate island. The two nearest-node errors
fall in finest boundary cells, i.e. within `1/32` of a coarse cell of a located
boundary. The hero figure is released on this basis.

## 3. Boundary report

Every located crossing is stored in `results/F7/{F7A,F7B,F7C}_events.csv` with:
parameter coordinates, witness subset, `kappa` and `H` before/after, the crossing
eigenvalue and frequency, boundary type, the direct right-half-plane count on
both sides, port visibility, closure distance `min |mu + 1|` over `sigma(Q)`, the
refined closure minimum and port frequency, the factorization residual of
identity (8), the `K(s)` solve residual, and the direct-path real part at the
located crossing.

| | F7A | F7B | F7C |
|---|---|---|---|
| subset crossings located | 11 586 | 10 928 | 7 465 |
| of which change `H` | 9 931 | 9 823 | 7 246 |
| `IMAGINARY_AXIS_TRANSVERSAL` | 11 586 | 10 800 | 7 465 |
| `IMAGINARY_AXIS_TANGENCY` on the mesh | 0 | 0 | 0 |
| `LOWER_BAND_EDGE` | 0 | **128** | 0 |
| `UPPER_BAND_EDGE` | 0 | 0 | 0 |
| `CORNER_OR_MULTIPLE` | 0 | 0 | 0 |
| `DAE_OR_EQUILIBRIUM_FAILURE` | 0 | 0 | 0 |
| unresolved / replay mismatch | 0 / 0 | 0 / 0 | 0 / 0 |
| direct RHP change at axis crossings | always 1 | always 1 | always 1 |
| direct RHP change at band-edge crossings | — | **always 0** | — |
| crossing frequencies (axis) | 0.316–0.714 Hz | 0.300–0.719 Hz | 0.494–0.727 Hz |
| base-stability (admissibility) boundary events | 609 | 609 | 0 |

No mesh edge was classified as a tangency: at the lattice resolution every edge
crossing near the tongue tip is transversal, because the two crossings of a
line just past the tip fall in different edges. Tangencies are
codimension one in the family of paths; they were located directly (§6).

## 4. Answers

**1. kappa regions observed.** `{1, 2, 3, 4, infinity}` on F7A and F7B,
`{2, 3, 4, infinity}` on F7C. Physical-area fractions of the base-stable domain
(linear in `g`, so dominated by the large-`g` safe region):

| `kappa` | F7A | F7B | F7C |
|---|---|---|---|
| infinity | 76.3 % | 89.7 % | 79.3 % |
| 4 | 2.9 % | 2.2 % | 6.2 % |
| 3 | 3.4 % | 1.8 % | 5.4 % |
| 2 | 6.8 % | 2.4 % | 9.1 % |
| 1 | 10.6 % | 4.0 % | — |

**2. Distinct hypergraphs.** 30 (F7A), 33 (F7B), 16 (F7C); **36 distinct over
the three slices**. Per `kappa`: F7A `{1: 10, 2: 13, 3: 5, 4: 1}`, F7B
`{1: 10, 2: 15, 3: 6, 4: 1}`, F7C `{2: 10, 3: 4, 4: 1}`. `kappa = 4` always has the
single hypergraph `{30,33,35,37}`; below it, one `kappa` value covers up to 15
different witness structures. Full tables: `results/F7/*_labels.csv`, with the
maximal guaranteed-safe portfolios (F2C Proposition 4) of every region.

**3. Connected or fragmented.** Counted on the refined lattice with
8-connectivity; a component is substantial with at least 25 nodes.

- F7A and F7C: every `kappa` level is connected.
- **F7B: `P_3` has two substantial components** (12 988 and 11 310 nodes): a
  fast-exciter component (`t = 0.74-1.00`, up to `g = 0.37`) and a slow-exciter
  component near fixed-Q (`t = 1.64-3.0`, `g < 0.016`), separated by the `kappa = 4`
  band. So along `t` at fixed Q, `kappa` goes `... 3 -> 4 -> 3`.
- Fragmented hypergraph regions: F7A `30|33` (2) and `30|33|37` (2); F7B
  `30+33+35` (3), `30+33+35|30+33+37` (3), `30|33` (2), `30|33|37` (2); F7C none.
- Heredity fails widely: closure defect (a safe superset of a hyperedge) at
  7 928 (F7A), 3 099 (F7B) and 247 (F7C) points, all in `kappa <= 2` regions.
  The stable-set family is **not** a simplicial complex there.

**4. Genuine tongues or islands.** One genuine **instability tongue** per slice:
a lobe of the flagship's inter-area family, characterized on F7A (§6); on F7B and
F7C its crossings lie in the same band (0.63-0.70 Hz) at the same policy range
(`g = 0.03-0.18`), consistent with the same lobe but not separately tracked. It is a
**prong** of the connected `kappa = 4` region, not an island; the safe region has
the matching notch. Its trace on pure-policy lines: F7A `k = 1.250-1.306`
(`kappa: 3 -> 4 -> inf -> 4 -> inf`), F7B `t = 1.090-1.184`
(`4 -> inf -> 4 -> inf`), F7C `h = 1.26-1.52` (`4 -> inf -> 4`). No isolated
island was found by refinement or by the 1 200 random checks.

**5. Band artifacts.** 128 of 29 979 crossings (all on F7B): triples `30+33+37`
and `30+35+37`, at `t = 2.64-3.0` and `g < 0.0095`, carry an **already unstable**
mode (`Re = +0.001` to `+0.057`) that slides across the 0.30 Hz lower edge. Each
changes `H` while the direct RHP count stays 1 -> 1. These are
band-classification transitions and the corresponding labels
(`30+33+35|30+33+37`, `30+33+35|30+35+37`, `30+33+35|30+33+37|30+35+37` in that
corner) depend on the band convention. Every other change is an
imaginary-axis crossing.

**6. Witnesses on the major regions** (area at least 0.5 %; complete tables in
`results/F7/*_labels.csv`):

| region `H` | `kappa` | maximal guaranteed-safe portfolios | F7A | F7B | F7C |
|---|---|---|---|---|---|
| `30+33+35+37` | 4 | the four triples | 2.9 % | 2.2 % | 6.2 % |
| `30+33+35` | 3 | `30+33+37`, `30+35+37`, `33+35+37` | 0.7 % | 0.5 % | 0.7 % |
| `30+33+35 | 30+33+37` | 3 | `30+35+37`, `33+35+37`, `30+33` | 1.8 % | 1.0 % | 3.6 % |
| `30+33+35 | 30+33+37 | 33+35+37` | 3 | `30+35+37`, `30+33`, `33+35`, `33+37` | 0.8 % | — | 1.1 % |
| `30+33` | 2 | `30+35+37`, `33+35+37` | 2.7 % | 0.8 % | 2.3 % |
| `30+33 | 33+35+37` | 2 | `30+35+37`, `33+35`, `33+37` | 0.5 % | — | 1.3 % |
| `30+33 | 33+35 | 33+37 | 30+35+37` | 2 | `30+35`, `30+37`, `35+37`, `33` | 1.0 % | — | 1.3 % |
| `30` | 1 | `33+35+37` | 1.7 % | 1.2 % | — |
| `30|33|35|37` | 1 | none | 2.6 % | 0.9 % | — |

Every `kappa <= 2` witness set on every slice contains bus 30 or bus 33, and every
`kappa = 2` region (38 region labels over the three slices) has `30+33` among its
hyperedges.

**7. Can policy move the same installed fleet?** Yes, in every direction, with
no change of topology, of the replaced buses or of PV megawatts; only `g`.

| move | where (pure `g`, other coordinates fixed) |
|---|---|
| unsafe -> safe | every base-stable line of all three slices (292, 396, 163 lines); 559 / 642 / 724 edge events |
| safe -> unsafe (increasing `g`) | the tongue: 19 / 25 / 83 edge events |
| `kappa 4 -> 3 -> 2` | F7B `t = 0.832-0.871` and F7C `h = 0.41-0.59` (increasing `g`); F7A `k = 1.653-1.672` (decreasing `g`) |
| `kappa` falls then rises | F7B `t = 0.66-0.70`: `2 -> 1 -> 2 -> inf` |

The strongest single path, F7B at `t = 0.852`, as `g` rises from 0:

    kappa  4 -> 3 -> 2 -> 3 -> 4 -> inf
    H      {30,33,35,37} -> {30,33,35} -> {30,33,35 | 30,33,37} -> {30,33,35 | 30,33,37 | 33,35,37}
           -> {30,33 | 33,35,37} -> {30,33} -> {30,33,35 | 30,33,37} -> {30,33,37} -> {30,33,35,37} -> {}

The reactive policy first contracts the incompatibility hypergraph to the pair
`30+33` and then expands it back and empties it.

**8. Heterogeneity alone changes H.** Yes. Along `h` at fixed Q (`g = 0`):

    h: 0 -> 2    H: {30,33 | 33,35,37} -> {30,33,35 | 30,33,37 | 33,35,37} -> {30,33,35 | 30,33,37}
                    -> {30,33,35} -> {30,33,35,37}                (kappa 2 -> 3 -> 4)

and at every `g` up to 0.3 the same contraction appears as the fleet is made
uniform. Along `h` every `kappa` change is an increase (2 462 `H` changes at
constant `kappa`, none from safe to unsafe): **in this slice heterogeneity is
protective**, and a uniform fleet at the geometric-mean time constant is *more*
fragile (`kappa = 2`), not less.

**9. Does closure detect every port-visible boundary?** Yes.
**29 851 of 29 851** imaginary-axis crossings are port-visible and detected; none
is port-invisible, so F6's case did not arise here. Multi-port witnesses
(24 887): closure distance at most `2.6e-6` (median `8e-8`), port frequency within
`1e-7` Hz of the crossing. Single-port witnesses (4 964): the individual factor
`det(I + M_aa)` at most `1e-6`. Factorization residual at most `6.9e-7`, solve
residual at most `6e-17`, direct-path `|Re lambda|` at the located point at most
`3.5e-9`. The 128 band-edge events are not on the imaginary axis and are not
claimed as closure detections.

**10. Random-map validation.** Section 2: 0 errors in 1 153 homogeneous-leaf
points (each map's 95 % bound below 0.8 %), 0 fast/direct disagreements in 1 200,
one boundary-sliver ISLAND_OR_TONGUE cell. **PASS.**

**11. Is the 4 -> 3 -> 2 witness contraction robust?** As a **kappa-witness**
chain, yes: `{30,33,35,37} -> {30,33,35} -> {30,33}` appears in order along the
excitation axis at every policy value up to `g = 0.25-0.30` on all three slices
(44 of 60 lines), towards higher gain (F7A), faster exciters (F7B) and a more
uniform fleet (F7C). It is absent only where the flagship region itself has
vanished (`g >= 0.3-0.4`). As a statement about the **whole** hypergraph it is
not a simple contraction: on the `g = 0` edge of F7A the order falls in 14
elementary moves,

    {F} ->C {30,33,35} ->I +{30,33,37} ->I +{33,35,37} ->I +{30,35,37}
        ->C {30,33 | 30,35,37 | 33,35,37} ->C ... -> {33 | 30,35 | 30,37} ->C {30 | 33} -> ... -> {30 | 33 | 35 | 37}

(C contraction, I insertion), reproducing F2's two transitions exactly
(`k* = 1.1488`, bracket 1.125-1.150; `k* = 1.8168`, bracket 1.800-1.825). And along
the policy axis the contraction reverses (item 7). **Not monotone.**

**12. Candidate IAS hero figure.** `figures/F7A_hypergraph_map.png`: the F7A plane
with `kappa` as an ordinal fill, exact `kappa` boundaries in black, hypergraph
boundaries at constant `kappa` in white, the κ-witness coalitions written on the
major regions, both folds marked, and an inset of the inter-area real part along
the pure-policy line `k = 1.30`. It shows in one frame that the order is
policy-dependent, that `H` carries structure `kappa` cannot, and that the policy
effect is not monotone. `figures/F7_three_slices.png` is the journal version.

## 5. Where the band convention matters

Only on F7B's slow-exciter, near-fixed-Q corner (§4 item 5). Any claim about the
labels there must say "inside the 0.3–1.5 Hz window". Everywhere else the
hypergraph boundaries are imaginary-axis crossings and do not depend on the
window except through which modes are counted.

## 6. Safeguard D — the narrow instability region near 0.68 Hz

`experiments/F7_tongue.py`, `results/F7/F7_tongue.json`; slice F7A, flagship.

**Which case.** **A — a genuine right-half-plane instability region.**

- **B (band artifact): no.** Crossings at 0.605, 0.639 and 0.697 Hz, far from both
  edges; the direct RHP count moves with `N_Gamma`.
- **C (tracking artifact): no.** `N_Gamma` is label-free; in addition the mode was
  tracked continuously from `g = 0` along 41 lines of constant `k` with minimum
  consecutive eigenvector MAC **0.99995**.
- **D (numerical near-tangency): no** away from the tip: the largest interior
  maximum of `Re lambda` is `+0.048 s^-1` against a fast/direct band error of
  `3e-8`.

**Which mode.** The flagship's own **inter-area family**: the track starts at the
`g = 0` inter-area mode (0.53 Hz at `k = 1.30`) and climbs continuously to
0.68-0.72 Hz as `g` rises; MAC 0.970 to the frozen v2C inter-area reference and
0.957 to the `g = 0` mode on machine coordinates; participation led by the rotor
states of machines 38, 34, 39 and 36. It is not a new mode.

**Along the pure-policy line `k = 1.30`:**

| segment | `g` | `H` | `N_Gamma(flagship)` | direct RHP (flagship) |
|---|---|---|---|---|
| near fixed Q | 0.017 | `{30,33,35,37}` | 1 | 1 |
| stable gap | 0.042 | `EMPTY` | 0 | 0 |
| inside the region | 0.102 | `{30,33,35,37}` | 1 | 1 |
| after | 0.173 | `EMPTY` | 0 | 0 |

Boundaries of the region at `g = 0.0509` (0.639 Hz, enters) and `g = 0.1532`
(0.697 Hz, leaves). Port closure at both: closure distance `3.1e-8` and
`1.1e-8`, port frequency matching to `2e-9` Hz, factorization residual `5e-16`,
individual factors 0.51 and 0.40: **collective, port-visible** on both sides.

**The two folds** (nondegenerate, Proposition 6 of F2D):

| | `g*` | `k*` | frequency | `R` | `dR/dg` | `dR/dk` | `d2R/dg2` |
|---|---|---|---|---|---|---|---|
| tip of the instability tongue | 0.10257 | 1.24954 | 0.681 Hz | `-4e-15` | `-2.1e-7` | `+0.359` | `-16.3` |
| closing point of the stable gap | 0.04092 | 1.30710 | 0.622 Hz | `+3e-14` | `+7.8e-10` | `+1.186` | `+230.6` |

At both, `Re lambda = 0` and the derivative along the policy direction vanishes
while the in-plane gradient does not. A pure-policy path at exactly `k*`
therefore meets the axis **tangentially** — an `IMAGINARY_AXIS_TANGENCY`, a
**boundary fold**, not an F2B Proposition-3 simple crossing for that path. In the
`(g, k)` plane the boundary remains a regular curve: one S-shaped curve with two
turning points relative to `g`, so that lines between `k = 1.2495` and `1.3071`
cross it three times. Proposition 6 of F2D makes this structurally stable.

Wording: "narrow instability region" or "instability tongue", descriptively.
Nothing here supports "resonance tongue" or "Arnold tongue".

## 7. What changes in the claim ledger

New observations O62–O70 and the corrections N21–N22 are in `docs/CLAIMS.md`.
The most consequential correction: **O59 is restated.** Making the fleet uniform
at its geometric-mean exciter time constant makes it *more* fragile, so the
phenomenon does not "require" heterogeneity; O59's stable uniform case was a fast
common exciter (`T = 0.25 s`), which is an exciter-speed effect.

## 8. Limits

- One network, one excitation model family, one operating point, one core of
  four candidates. Proposition 5 of F2C makes the core hypergraph the exact trace
  of the fleet hypergraph, but nothing outside the core was tested.
- Area fractions are in the declared coordinates and change with them.
- The leak `w` is a modelling constant: boundary locations and area fractions
  move with it (72-89 % label agreement over a tenfold range), the qualitative
  results do not (§2).
- Features thinner than `1/32` of a coarse cell can escape any finite lattice;
  the near-axis monitor and the random checks bound, but do not exclude, them.
- Not EMT, not time-domain; small-signal only.
