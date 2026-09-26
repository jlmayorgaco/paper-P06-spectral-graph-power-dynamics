# F2E — Whole-RHP incompatibility and mechanism-resolved projections

This is the Gate 4 reframing. It is written after the post-F7 freeze. The frozen
notes `F2B`, `F2C` and `F2D` are not edited. Their statements hold for any
protected region `Gamma` with a piecewise-C1 boundary, so they hold for both
regions defined here. This note fixes which region carries the planning claim
and how the two are related.

Evidence cited: `docs/G1_WHOLE_RHP_COMPOSABILITY.md` (IEEE-39, Kundur, F8 under
RHP), `docs/G3_IEEE68.md` (68-bus), `docs/G2_TDS_VALIDATION.md` (nonlinear
time-domain).

## 1. The primary object

    Gamma_RHP = { s : Re s > 0, |s| > eps }       eps = 1e-3 (the reference double zero)
    N_RHP(S; theta) = # eigenvalues of A_red(S; theta) in Gamma_RHP
    U_RHP = { S : N_RHP(S) > 0 },    H_RHP = min U_RHP,    kappa_RHP = min_{E in H_RHP} |E|

`H_RHP` is the **primary** object. Take a tested portfolio that contains no
hyperedge of `H_RHP`. It is small-signal stable, since an unsafe portfolio
contains a minimal unsafe subset (F2C). The converse is not guaranteed. A
portfolio that contains a hyperedge can still be stable, because `U_RHP` need
not be upward closed; the closure defect `up(H) \ U` of F2C measures this. So
`H_RHP = EMPTY` means that every tested portfolio is small-signal stable. Every
planning statement in this project uses `H_RHP`, and it uses it as a
sufficient condition for safety ("no hyperedge contained").

The `eps` disc is the only departure from the open half-plane. It removes the
two zero eigenvalues of the angle reference, which are structural and present
in every portfolio. A real eigenvalue that crosses through the origin
registers when it leaves the disc, at `+eps` rather than at `0`. That shifts
the parameter value of an aperiodic boundary by `O(eps / |d lambda / d theta|)`.
The shift is declared rather than hidden.

## 2. Mechanism-resolved projections

A **declared mechanism** `m` is a subregion `Gamma_m` of `Gamma_RHP`, chosen
because one physical mode family lives there. An example is the inter-area band

    Gamma_IA = Gamma_RHP ∩ { 0.3 <= |Im s| / 2pi <= 1.5 Hz }.

    N_m(S) <= N_RHP(S),    U_m ⊆ U_RHP,    H_m = min U_m.

`H_m` is a **projection** of the incompatibility structure onto mechanism `m`.
It answers "which coalitions destabilize *this* mode family". It does not
answer "which coalitions are unsafe".

**Proposition 1 (domination).** Every hyperedge of `H_m` contains a hyperedge
of `H_RHP`. Hence `kappa_RHP <= kappa_m`, and `H_RHP = EMPTY` implies
`H_m = EMPTY` for every `m`.

*Proof.* Take `E` in `H_m`. It lies in `U_m`, which is contained in `U_RHP`.
Every member of a finite family contains a minimal member, so `E` contains
some `F` in `min U_RHP = H_RHP`. The inequality follows by taking sizes. QED

**Proposition 2 (the converse fails).** `H_m = EMPTY` does not imply
`H_RHP = EMPTY`. So

> **`H_Gamma = EMPTY` does NOT imply global (small-signal) stability unless
> `Gamma` covers every instability mechanism of interest.**

*Counterexamples, both from documented benchmarks and both on the full
nonlinear-model path:*

- **Kundur two-area.** `H_IA = EMPTY` at 556 of the K12A nodes and 374 of the
  K12B nodes. `H_RHP` is non-empty at every node of both slices. Pairs of
  grid-following replacements diverge aperiodically, with a real eigenvalue,
  which the band excludes by construction. G2 shows the nonlinear trajectory
  of such a pair leaving the operating point along the predicted real mode.
- **68-bus, 12-plant family.** At `g = 0.25` the only hyperedge (every NETS
  plant) destabilizes a converter PLL mode at 2.34 Hz, outside the band, so
  `H_IA` cannot contain it.
- **IEEE-39 with an undamped condenser.** A condenser with 1 % inertia, frozen
  EMFs and `D = 0` gives `H_IA = EMPTY` and `H_RHP != EMPTY` at every F8 point.
  The cause is the condenser's own swing mode at 7.7–14.3 Hz. G2 reproduces
  it in the time domain at 8.65 Hz.

**Proposition 3 (recomposition).** Let the declared mechanisms partition
`Gamma_RHP` into disjoint `Gamma_1, ..., Gamma_r`, for example "inter-area
band", "aperiodic, near the real axis", and "other oscillatory". Then
`N_RHP = sum_m N_m`, `U_RHP` is the union of the `U_m`, and

    H_RHP = min ( H_1 ∪ H_2 ∪ ... ∪ H_r ).

*Proof.* `N_RHP(S) > 0` exactly when some `N_m(S) > 0`, which gives the union
of the families.

- Let `E` be minimal in the union of the `U_m`, and say `E` lies in `U_m`. Any
  proper subset of `E` lying in `U_m` would also lie in the union, so `E` is
  minimal in `U_m`, that is `E` is in `H_m`. `E` is also minimal within the
  union of the `H_m`, because that union is contained in the union of the
  `U_m`.
- Conversely, let `E` be minimal in the union of the `H_m`, and suppose some
  `F` in the union of the `U_m` is a proper subset of `E`. Then `F` contains a
  minimal member `F'` of its own family, and `F'` is in the union of the
  `H_m`. Since `F'` is a proper subset of `E`, this contradicts the minimality
  of `E`.

QED

The planning object is therefore the minimal envelope of the mechanism
projections, not their union and not any single one. The projections remain
the right tool for **explaining** a transition. For example, G1 classifies
every discrepancy `H_RHP \ H_IA` by the mechanism whose modes cause it.

## 3. What carries over from F2B–F2D

These hold on `Gamma_RHP` without change of proof:

- `H_RHP` is an antichain (F2C Prop. 1).
- `H_RHP` is locally constant on the parameter set where no eigenvalue of any
  tested portfolio lies on `∂Gamma_RHP` (F2C Prop. 3, by the argument
  principle).
- `H_RHP` can change at fixed `kappa_RHP` (F2C, F2D). On IEEE-39 F7A and F7C,
  `H_RHP` equals the band hypergraph, which has up to 15 hypergraphs behind
  one `kappa`.
- The regions `{theta : H_RHP(theta) = H}` can be non-monotone and
  disconnected (F2D; the construction there uses a single mode and applies
  verbatim).
- The stable-set family is **not** called a simplicial complex. Heredity fails
  in general (F2C §5), and the 12-plant 68-bus family was not tested for it.

**Boundaries simplify.** `∂Gamma_RHP` consists of the imaginary axis and the
`eps` circle, so Safeguard B's band-edge types disappear. Every change of
`H_RHP` is exactly one of:

| boundary | where | port closure |
|---|---|---|
| oscillatory imaginary-axis crossing | `lambda = i omega`, `omega > 0` | evaluable. The return difference `det(I + M_S(i omega))` vanishes when `i omega` is not a device pole, as for the band. |
| aperiodic crossing through the origin | a real eigenvalue leaving the `eps` disc | **not evaluable**: the reduced port operator `T_0(s)` carries the reference zero at `s = 0`. These boundaries are certified on the full model only. |
| DAE or equilibrium failure | loss of the operating point or of index 1 | not a spectral boundary. It is reported separately. |

Evidence: 601 of 601 oscillatory RHP boundaries on IEEE-39 and Kundur are
port-visible. The largest multi-port closure distance is `1.0e-6`, and the
largest single-port individual factor is `7.8e-7`. The 194 aperiodic Kundur
boundaries are located and classified but not port-certified.

## 4. When the band object is the planning object

Only when it coincides with `H_RHP`, and that is an empirical statement for
each case and each region:

| case | `H_RHP = H_IA` |
|---|---|
| IEEE-39 F7A, F7C | 100 % |
| IEEE-39 F7B | 96.5 % (band-edge exits of the same branch) |
| IEEE-39 F8 services | 85.4 % (undamped condensers) |
| Kundur | 0.3–8 % |
| 68-bus, 4-candidate map | 100 % (both empty everywhere) |
| 68-bus, 12-plant family | the band misses the `g = 0.25` and `g = 1` hyperedges (converter PLL modes, 2.3–3.5 Hz) |

The band study of the declared inter-area mechanism is legitimate as a
**mechanism** study: F7's contraction of the witness coalition, the
instability tongue and its folds, and F8's service attribution of the
inter-area modes. It is not a safety study.

## 5. Wording fixed by this note

- "Composable" means `H_RHP = EMPTY`.
- "Composable with respect to the inter-area mechanism" means `H_IA = EMPTY`,
  and must be written that way.
- A `kappa` without a subscript refers to `kappa_RHP`.
- Never write "`H = EMPTY` implies stability" without naming `Gamma_RHP`.
