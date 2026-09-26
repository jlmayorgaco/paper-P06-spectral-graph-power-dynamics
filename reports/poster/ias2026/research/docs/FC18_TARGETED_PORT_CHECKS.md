# Targeted port-structure checks on the frozen IEEE-39 (FC18)

**Scope and sources.**
- Script `experiments/final_closure/FC18_targeted_port_checks.py`.
- Output `outputs/.../FC18_targeted_port_checks/`: `FC18_summary.json`,
  `FC18_events.csv`, `FC18_walk_rho_band.json`.
- Frozen model only: core 30/33/35/37, matched dispatch, no GFM, no DC link,
  no limits, no governor.
- No parameter search. The ten boundary points are the nine witness changes of
  the frozen F7B line (`t = 0.8515625`, `k = h = 1`, κ: 4 → 3 → 2 → 3 → 4 → ∅)
  plus the FC13 flagship boundary on the P4 line.
- Each boundary was re-located by a 40-step bisection on the direct path.

**Objects.**
- Common descriptor realization C2 (BC02).
- Port `T_S(s)`.
- `M = D K`, with `D = blkdiag(ΔY_i)` and `K = Eᵀ T_0⁻¹ E`.
- The factorization `I + M = (I + D K_d)(I + Q)` with
  `Q = (I + D K_d)⁻¹ D K_o`, where `K_o` holds the off-diagonal network-closure
  blocks between candidate buses.
- Hence `det(I + M_SS) = Π_{i∈S} det(I + M_ii) · det(I + Q_SS)`, and `Q` has
  zero diagonal blocks.

## 1. Descriptor interaction lifting

**Affinity.** In all four descriptor blocks `f_x, f_z, g_x, g_z`, every Möbius
coefficient of order ≥ 2 over the 16 vertices is ≤ 4.5e-16 relative, at all ten
points. Each replacement block is local to its bus:

| block | nonzero entries of a singleton increment | what they are |
|---|---|---|
| `g_z` | 4 | the 2×2 bus block |
| `g_x` | 12 | injection rows × device states |
| `f_z` | 28–29 | device rows × bus voltage |
| `f_x` | 59–61 | device-internal |

**Port identity.** `det T_S / det T_0 = det(I + M_SS)` holds to ≤ 1.4e-13 for
all 16 subsets at all ten boundaries (160 checks). The critical eigenvalue of
`M_S` sits at −1 to ≤ 8.5e-14. The identity is Sylvester's determinant lemma;
it is not claimed as new.

**Expansion at the boundary.** At every boundary the zero of `det(I + M_SS)`
comes entirely from the interaction factor. The local factor stays at 0.08–0.31.

| boundary | S | local factor | interaction factor | truncation of det(I+Q) at order ≤ 2 / ≤ 3 / ≤ 4 |
|---|---|---|---|---|
| E1 (κ 4→3), 0.628 Hz | 30+33+35 | 0.315 | 1e-14 | 0.336 / **0** |
| E2 | 30+33+37 | 0.198 | 2e-13 | 0.396 / **0** |
| E3 | 33+35+37 | 0.080 | 3e-13 | 0.385 / **0** |
| E4 (κ 3→2), 0.656 Hz | 30+33 | 0.081 | 1e-13 | **0** |
| E5 | 33+35+37 | 0.135 | 9e-14 | 0.296 / **0** |
| E6 (κ 2→3) | 30+33 | 0.104 | 1e-13 | **0** |
| E7 | 30+33+35 | 0.130 | 7e-15 | 0.372 / **0** |
| E8 (κ 3→4) | 30+33+37 | 0.132 | 4e-13 | 0.293 / **0** |
| E9 (→ ∅), 0.718 Hz | 30+33+35+37 | 0.134 | 8e-14 | 0.590 / 0.132 / **0** |
| flagship, P4 line, 0.706 Hz | 30+33+35+37 | 0.089 | 8e-14 | 0.605 / 0.161 / **0** |

At every one of the ten boundaries, the characteristic zero is reached only by
the full-order network-closure term of the coalition whose label changes. Any
truncation below order `|S|` leaves 0.13–0.60. If network closure is removed
(`K_o = 0`), no boundary exists at all.

## 2. Rank / effective dimension of the interaction operator

Along the path, κ goes 4 → 3 → 2 → 3 → 4 → ∅, but the effective rank of `Q`
stays full:
- 7–8 at relative threshold 1e-2 (8 at 1e-3);
- 6–8 at threshold 1e-1.

The smallest singular values are 0.04–0.06, and no dimension is lost when the
witness contracts. The only rank bound (order ≤ rank = 8) is vacuous for four
candidates.

What does change is the **magnitude** of `Q`: σ₁ is 6.1–8.9 at the κ-changing
contractions E1 and E4, and 1.7–1.8 at the four-bus boundaries. This is
descriptive only (ten points, not tested).

**Verdict:** effective rank is unrelated to the observed contraction. Interaction
strength, not dimension, varies with the witness order.

## 3. Boundary-gradient validation

`ds*/dθ = −(p* F_θ q)/(p* F_s q)` with `F = I + M_SS` was compared with direct
full-DAE finite differences of the tracked transverse eigenvalue, over all four
policy coordinates (g, k, t, h) and in box-normalized coordinates, at all ten
boundaries:

| quantity | error |
|---|---|
| normal-vector angle | ≤ 4.3e-4 ° |
| normal magnitude | ≤ 4.3e-6 relative |
| imaginary part (frequency drift) | ≤ 9e-6 relative |
| `dRe s*/dg` | equal to 4 digits at every point, e.g. −0.4610 vs −0.4610 at the flagship |

**Retuning.** Target: the minimum Q/V-gain change that makes the P4 flagship
stable (a pure `g` move; the other parameters fixed). Reference: direct
continuation, `g* = 0.2076814045`, i.e. `Δg = +0.1714`. The existing retune
optimizers (E21/E23/E34) act on a different parameter vector at a different
operating point, so they are not a like-for-like comparison.
- **One-step prediction from P4:** `Δg = 0.0518`, only 30 % of the required
  change. `dRe s/dg` drifts from −2.45 at P4 to −0.46 at the boundary, so a
  single linear step from a deep unsafe point is not accurate.
- **Port-derivative Newton iteration:** converges in 5 steps to
  `g = 0.2076814043`, an error of 2e-10 against direct continuation.

## 4. Dynamic network-walk decomposition

The splitting keeps the network variables: `T_0 = D_g − N`, with `D_g` the 2×2
bus-diagonal blocks, and `T_0⁻¹ = Σ (D_g⁻¹ N)^k D_g⁻¹`.

| case | `rho(D_g⁻¹ N)` at the evaluation point | over the 0.3–1.5 Hz band |
|---|---|---|
| flagship boundary (P4 line, 0.706 Hz) | 1.44 | 1.003–1.556 |
| Pg-matched stable control 31+32+33+38 (E12 policy, rightmost band mode 1.03 Hz) | 1.10 | 1.003–1.554 |

**NON-EXECUTABLE.** The Neumann condition fails at every frequency of the band
for both cases, so no walk or path interpretation is given.

What remains well defined is the share of the critical coupling that comes from
network closure: `|p* D K_o q| / |p* M q|` = 0.68 at the flagship boundary. For
the control the ratio exceeds 1, because the local and network parts cancel; it
is not interpretable.

## 5. Second-order curvature spot check (FC07 v2)

Cases: baseline, P4 flagship, and the policy-repaired P_inf flagship.

**Complete second order vs linear, at 2 MW.**

| observable | error vs linear | slope on 2–16 MW |
|---|---|---|
| bus voltage | 750–3400× smaller | cubic remainder: `p` = 3.03–3.15 |
| current | 60–490× smaller | 2.0–3.0 |
| frequency | 170–450× smaller | about 2 at the smallest amplitudes |

In current and frequency a discretization floor (`dt = 1e-3`) keeps `p` near 2
at the smallest amplitudes.

**Removing `F_z D²ψ`.** The frequency error becomes 570–1700× the *linear* one.
With the rotation / drift artifact removed, the device curvature alone is still
no better than linear.

**Status:** methodological support only.

## 6. Frozen-benchmark rule

Complied with: nothing was added to IEEE-39. A mixed SG/GFL/GFM benchmark would
be a separately versioned case.

## Summary table

The mixed-device column records what was reported for the six-bus SG–GFM–GFL
case; it was not re-run here.

| claim | independent mixed-device validation | IEEE-39 result | pass/fail | paper role |
|---|---|---|---|---|
| Descriptor blocks affine in the binary replacement variables | reported (six-bus) | order ≥ 2 Möbius ≤ 4.5e-16 in `f_x, f_z, g_x, g_z`; increments bus-local | PASS | model lemma (support) |
| `det T_S / det T_0 = det(I + M_SS)` | reported | ≤ 1.4e-13, 160 checks at 10 boundaries | PASS | classical tool (Sylvester); not novel |
| Higher-order portfolio terms arise from network closure of local blocks | reported (principal-minor lifting) | local factor 0.08–0.31 ≠ 0; interaction factor 1e-13; truncation below order \|S\| leaves 0.13–0.60 at all 10 boundaries; no boundary without `K_o` | PASS | main mechanistic statement (exact bookkeeping of the witness structure) |
| −1 boundary closure | reported | critical eigenvalue of `M_S` at −1 to ≤ 8.5e-14 | PASS | tool (classical GN form) |
| Rank / order limitation | reported | effective rank 7–8 throughout κ 4→3→2→3→4; no dimension loss | FAIL (unrelated / vacuous on IEEE-39) | negative remark |
| Port boundary sensitivity | reported | normal angle ≤ 4.3e-4 °, magnitude ≤ 4.3e-6, 10 points | PASS | tool: exact boundary normals |
| Port-guided Q/V retuning | — | one step: 30 % of `Δg`; Newton: `g*` to 2e-10 in 5 steps | PASS (iterated) / FAIL (single step) | engineering remark: local, iterative only |
| Dynamic network-walk expansion | reported | `rho` = 1.003–1.556 over the band, both cases | NON-EXECUTABLE | none (future case) |
| KCL curvature at second order | reported | cubic in V; `F_z D²ψ` removal 570–1700× worse than linear in frequency | PASS | methodological support |
| Frozen benchmark untouched | — | nothing added | complied | rule |

## Answers

1. **Does descriptor network closure explain the observed high-order portfolio
   structure?** Yes, as exact bookkeeping.
   - The replacement blocks are exactly affine and bus-local.
   - All non-additivity sits in the network-closure operator `Q`.
   - At every one of the ten frozen boundaries, the zero is produced only by
     the full-order `|S|` closure term of the coalition that changes label.
   - It explains *where* the higher-order terms come from. It does not predict
     κ without computing `Q`.
2. **Does interaction-space dimension constrain the observed order usefully?**
   No. The effective rank is 7–8 throughout the 4 → 3 → 2 contraction, and the
   bound is vacuous. What varies is interaction magnitude (descriptive only).
3. **Does the port derivative predict F7 boundary motion accurately enough to
   guide retuning?**
   - Locally, yes: boundary normals agree to 4e-4 ° and 4e-6.
   - For retuning, only iteratively: Newton reaches the direct boundary to
     2e-10 in 5 steps. A single linear step from P4 recovers only 30 % of the
     required Δg.
4. **Can critical port coupling be tied to interpretable physical network
   paths?** No. The walk expansion is NON-EXECUTABLE: `rho ≥ 1.003` at every
   band frequency for both the flagship boundary and the matched control.
