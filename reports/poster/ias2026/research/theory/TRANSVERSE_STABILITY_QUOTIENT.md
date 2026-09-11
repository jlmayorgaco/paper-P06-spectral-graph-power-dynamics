# Transverse (relative) stability of the frozen IEEE-39 benchmark

Status: definition plus proof (algebra); numerical verification in
`outputs/ias2026/final_math_nonlinear_validation_*/FC01_transverse_quotient/`
(summary in `results/TSQ_ieee39_reaudit.csv`). Code:
`src/ibr_cycles/certification/transverse.py`.

## 1. Model-scope statement (to be stated in every paper)

> The frozen IEEE-39 model has no primary frequency restoration (every machine
> has D = 0 and no governor) and is analysed in relative / transverse
> coordinates. Absolute common-frequency restoration is outside this benchmark.

"Stable", "Hurwitz" and "whole-RHP" always refer to the **transverse**
dynamics defined below. For IEEE-39 the planning object is written
`H_RHP^perp`, and the region is the **transverse whole-RHP**.

## 2. The exact structural center subspace

Let the phasor DAE be `x' = f(x, z)`, `0 = g(x, z)` (`M = I` in L0) with
reduced linearization `A = f_x − f_z g_z^{-1} g_x` at an equilibrium.

**Rotation (gauge symmetry).** Every device angle state (rotor angle `delta`,
PLL angle `theta_pll`) enters the equations only through `v e^{-j angle}`.
Every other device state is written in its own rotating frame, and the network
and the PQ / Z loads are rotation-covariant.

- Therefore `f(x + phi R_x, e^{j phi} z) = f(x, z)` and the KCL rotates
  covariantly, for every `phi`.
- Differentiating at `phi = 0` gives `f_x R_x + f_z R_z = 0` and
  `g_x R_x + g_z R_z = 0`, hence `A R_x = 0`.
- Here `R_x = 1` on every `delta` and `theta_pll`, and `R_z = jV`.
- This is an invariance of the equations: **no physics**.

**Common-frequency drift (physical marginal mode).** Let `w` be the uniform
frequency shift: machine speed `+1`, PLL integrator `+omega_B`,
synthetic-inertia filter `+1`, and on the 68-bus system the speed-input PSS
washout `+K`. The one-parameter family

    x(t) = x0 + c w + c omega_B t R_x,   z(t) = e^{j c omega_B t} z0

solves the nonlinear DAE exactly whenever the following all hold:

- every machine has `D = 0`;
- the mechanical power is constant (no governor);
- the network is frequency-independent (quasi-static phasors).

The reason is that the electrical power is rotation-invariant and the power-form
swing equation then has `omega' = 0` at speed `1 + c`. Differentiating in `c`
gives

    A w = omega_B R_x.

The drift is a **physical** property of the modelled plant: there is no
primary frequency restoration. It is not a gauge.

**Structure.** `C = span{R_x, w}` is A-invariant with

    A [R_x  w] = [R_x  w] J,   J = [[0, omega_B], [0, 0]].

Because `omega_B != 0`, this is a **Jordan chain of length 2** (one eigenvector
`R_x`, one generalized eigenvector `w`), not two independent zero modes:

- `dim C = 2`;
- the zero eigenvalue has geometric multiplicity 1 and algebraic multiplicity 2
  inside `C`.

**Left basis.** `W_c` spans the left generalized null space `null((A^2)^T)`,
normalized so that `W_c^T [R_x w] = I`. Then `W_c^T A = J W_c^T`.

**Invariance for every portfolio and policy point.** The derivation uses only:

- the device frames;
- `D = 0`;
- constant mechanical power;
- the frequency independence of the network.

It does not use any parameter value. `C` is therefore structurally invariant for
every frozen portfolio, reactive policy, AVR speed and condenser configuration
with `D = 0`. It fails, and must be recomputed, exactly when one of those
assumptions fails. For example:

- the G1 condensers with `D = 2` leave `C = span{R_x}`;
- the governed model of §6 also leaves `C = span{R_x}`.

## 3. Transverse dynamics

Let `U` be an orthonormal basis of `C` and `Z` an orthonormal basis of its
complement. Invariance gives `Z^T A U = 0`, so

    [U Z]^T A [U Z] = [[U^T A U, U^T A Z], [0, A_perp]],   A_perp := Z^T A Z,
    det(sI − A) = s^2 det(sI − A_perp).

`A_perp` is the matrix of the induced map on the quotient `X / C`. It is the
dynamics of every relative coordinate: angles relative to the rigid rotation,
and speeds relative to the common drift. It does not depend on the choice of
`Z`. Any other orthonormal complement `Z' = Z O` gives `O^T A_perp O`, the
same spectrum.

**Definitions.**

    alpha_perp(S, theta) = max Re sigma(A_perp(S, theta)),
    N_RHP^perp(S, theta) = #{lambda in sigma(A_perp) : Re lambda > 0},
    transverse stable  iff  alpha_perp < 0,

- `H_RHP^perp(theta)`: the minimal portfolios with `N_RHP^perp > 0`.
- `kappa_RHP^perp`: their minimum cardinality.

No eigenvalue is removed by magnitude. A physical eigenvalue of `A_perp`
arbitrarily close to 0 is counted with its sign, provided that sign is resolved
by the frozen four-state classifier (`d_axis > 10 eps`). Otherwise it is
`BOUNDARY_OR_UNRESOLVED`.

**Equivalence with Phase I.** The BC01 construction (quotient by `R_x`, then
deflation of the exact eigenvector `Z_1^T w`) produces the same spectrum as
`A_perp`. Both are quotients by `C`, and the spectrum of a quotient map is
basis-independent.

**Relation to the old disc rule.** `|lambda| <= 1e-3 → ignore` coincides with
the transverse count exactly when `A_perp` has no eigenvalue in that disc. The
disc is not a valid definition: it would hide a physical real crossing (the
Kundur Q/V integrator crosses at `g` of about 3e-4 with eigenvalue about 1e-2
beside it).

## 4. Numerical verification (FC01 part 1, direct Jacobians)

On 106 cases (six F8 points × 16 flagship subsets, a Kundur node × 8, and
IEEE-68 base plus 4 candidates):

| check | result |
|---|---|
| `dim C` | 2 in every case |
| Jordan block | `J_12 = omega_B = 376.99` in every case |
| invariance `‖A Z_c − Z_c J‖` | ≤ 1.1e-11 relative |
| coupling `‖Z^T A U‖` | ≤ 4.8e-11 relative |
| spectral identity `sigma(A) = {0, 0} ∪ sigma(A_perp)` | ≤ 3.1e-8 |
| left basis condition | 3–26 (IEEE-39, Kundur); about 2200 (IEEE-68, stiff scaling) |
| transverse modes with `abs(lambda) < 1e-2` | none at the audited points |
| old-cutoff RHP count vs `N_RHP^perp` | equal in 106 of 106 |

The label re-audit of all 345 229 F7 points is in `FC01_summary.json`.

## 5. Kundur and IEEE-68

Both frozen models also have `D = 0` and no governor, so the same `C` (dimension
2, Jordan chain) applies. BC01 verified the partner in 100 % of 49 272 cases. The
transverse definitions are therefore the correct ones for them as well. This is
a property of those models as frozen, not an IEEE-39 convention.

## 6. The governed variant (`ieee39_governed_documented_v1`)

With the documented TGOV1N governors (same source workbook), the rotation
survives (`A R_x = 0` to 1e-10). The partner does not: `‖A w − omega_B R_x‖` is
0.10–0.72 relative. The center subspace is `C = span{R_x}` (dimension 1), and
transverse stability is the plain T1 quotient.
