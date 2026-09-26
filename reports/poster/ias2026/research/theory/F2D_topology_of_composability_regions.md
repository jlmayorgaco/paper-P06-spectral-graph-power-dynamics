# F2D — Composability regions need not be monotone, nested or connected

Companion to F2, F2B and F2C. Every statement here is mathematics about analytic
linear time-invariant families; none is drawn from IEEE-39. The constructions
are `tests/test_composability_topology.py`.

**Positive result, unchanged.** Under (H0)–(H3) of F2B, `N_Gamma(S; ·)` for every
tested `S`, hence `kappa_Gamma` and the whole hypergraph `H_Gamma`, are locally
constant on `Theta_reg` and constant on each of its connected components
(F2B Theorem 1, F2C Theorem 6).

**Negative result.** Nothing more global follows without further assumptions:
along a one-parameter path `kappa` can fall and rise again, a level set `P_k` can
have several components, and `H` can move while `kappa` stands still. So there is
no general monotone ordering of physical parameter space by composability order.

## 1. The model class

All constructions use one oscillatory mode per replacement set,

    A(S; theta) = sigma(S; theta) I_2 + omega_0 J,     J = [[0, 1], [-1, 0]],
    omega_0 = 2 pi 0.6 rad/s,

whose eigenvalues are `sigma(S; theta) +- i omega_0`. With the project region
`Gamma = { Re s > 0, 0.3 <= |Im s| / 2pi <= 1.5 Hz }`, `N_Gamma(S; theta) = 1` iff
`sigma(S; theta) > 0`, and a crossing of `dGamma` is always an imaginary-axis
crossing at 0.6 Hz. The effect of the actions is **additive**,

    sigma(S; theta) = -1 + sum over a in S of c_a(theta),

so `A(S; theta) = A_0 + sum_{a in S} c_a(theta) I_2`: a common baseline with one
rank-two update per action, and no interaction term at all. Every obstruction
below therefore arises without any genuine joint effect. (H0)–(H2) hold
trivially (no algebraic variables, fixed dimension two), and (H3) fails only on
the zero sets of the `sigma(S; ·)`, which are finite (one parameter) or curves
(two parameters).

Two states is the minimum: a real one-state system has no eigenvalue in the band.

## 2. Counterexample 1 — kappa falls, then rises

Two actions, `theta` in `[0, 1]`,

    c_1(theta) = 0.6 + 2.4 theta (1 - theta),        c_2(theta) = 0.6.

Then `sigma({2}) = -0.4` always, `sigma({1,2}) = 0.2 + 2.4 theta(1 - theta) > 0`
always, and `sigma({1}) = -0.4 + 2.4 theta(1 - theta)` is positive exactly on
`(theta_-, theta_+)`, `theta_+- = (1 +- 1/sqrt 3) / 2 = 0.2113, 0.7887`.

| `theta` | unsafe sets | `H_Gamma` | `kappa_Gamma` |
|---|---|---|---|
| `[0, theta_-)` | `{1,2}` | `{ {1,2} }` | 2 |
| `(theta_-, theta_+)` | `{1}`, `{1,2}` | `{ {1} }` | **1** |
| `(theta_+, 1]` | `{1,2}` | `{ {1,2} }` | **2** |

Both transitions are transversal simple crossings
(`d sigma({1}) / d theta = 2.4 (1 - 2 theta) = +-1.386` at `theta_+-`), so they
are exactly the transitions F2B Proposition 3 describes, and still `kappa` is not
monotone along the path.

Minimality. With one action `kappa` takes values in `{1, infinity}`, so a
decrease followed by an increase between **finite** orders needs two actions.
With one action `kappa` can still go `infinity -> 1 -> infinity`
(`c_1(theta) = 6 theta (1 - theta)`, unsafe on `(0.211, 0.789)`).

## 3. Counterexample 2 — a level set with several components

**One parameter.** In Counterexample 1, `P_2 = [0, theta_-) union (theta_+, 1]`:
two components, separated by `P_1`.

**Two parameters, islands.** Two actions, `theta = (x, y)` in `[-2, 2]^2`,

    c_1(x, y) = 1.3 - 4 [ (x^2 - 1)^2 + y^2 ],         c_2 = 0.6.

Write `q = (x^2 - 1)^2 + y^2`. Then `sigma({1}) = 0.3 - 4q` and
`sigma({1,2}) = 0.9 - 4q`, so

- `P_1 = { q < 0.075 }`: two disjoint disks around `(+-1, 0)` (at `x = 0`, `q = 1`),
  i.e. **two islands** of order one;
- `P_2 = { 0.075 < q < 0.225 }`: **two components**, each an annulus around one
  island, so each is **not simply connected**;
- `P_infinity = { q > 0.225 }`: one component, surrounding both.

A `kappa` phase diagram is therefore a genuinely topological object: its regions
can be fragmented and can have holes, and the number of components is part of
what a map must report.

## 4. Counterexample 3 — H moves while kappa stands still

F2C Proposition 8, restated in this model class: three actions,

    c_1 = 0.6,    c_2(theta) = 0.6 - 0.4 theta,    c_3(theta) = 0.3 + 0.4 theta,

`H` goes `{ {1,2} } -> { {1,2}, {1,3} } -> { {1,3} }` at `theta = 1/4, 1/2` while
`kappa = 2` throughout.

## 5. Counterexample 4 — composability is not monotone in the actions either

`c_1 = c_2 = 0.6`, `c_3 = -0.5`: `{1,2}` is unsafe (`+0.2`) and `{1,2,3}` is safe
(`-0.3`). The unsafe family is not upward closed, the closure defect of F2C is
nonempty, and the safe family is not a simplicial complex. One stabilizing action
suffices; no interaction is needed.

## 6. What does force monotonicity

**Proposition 5 (monotone destabilization).** Let `theta(s)`, `s` in `[0, 1]`, be
a path in `Theta` along which, for every tested `S`, the indicator
`1[N_Gamma(S; theta(s)) > 0]` is non-decreasing. Then `U_Gamma(theta(s))` is
non-decreasing under inclusion and `kappa_Gamma(theta(s))` is non-increasing.

*Proof.* `U(s_1)` is contained in `U(s_2)` for `s_1 <= s_2`, and `kappa` is the
minimum size over `U`, which cannot increase when `U` grows. QED

In the additive one-mode class a sufficient condition is that every `c_a` is
non-decreasing along the path. The counterexamples each violate the hypothesis
in the smallest possible way: one subset's indicator switches on and then off.
`H` is not monotone even then — insertions can contract or reorganize it — but
`kappa` is.

**Remark.** A monotone map would require the physical parameter to act
monotonically on every subset's band abscissa. Nothing in the physics of a
voltage regulator, an exciter or a heterogeneity amplitude guarantees that, and
F7 records where it fails.

## 7. The generic geometry of a tongue

The mechanism behind Counterexample 1 is a real part with an interior extremum
along the path. In two parameters that is a structurally stable feature, not an
accident.

**Proposition 6 (fold of the boundary with respect to one coordinate).** Let
`R(g, y) = Re lambda(g, y)` be the real part of a simple eigenvalue branch with
`omega` in the open band, `C^2` near `(g*, y*)`, and suppose

    R = 0,    R_g = 0,    R_y != 0,    R_gg != 0      at (g*, y*).

Then

1. the boundary `{ R = 0 }` is, near `(g*, y*)`, a `C^1` curve (the in-plane
   gradient `(0, R_y)` is nonzero), so the F2B boundary is regular in the plane;
2. it is tangent to the `g` direction at `(g*, y*)` and lies on one side of the
   line `y = y*`, locally
   `y - y* = -(R_gg / 2 R_y) (g - g*)^2 + o((g - g*)^2)`;
3. the path `y = y*` meets `dGamma` **non-transversally**: along it
   `R(g, y*) = (R_gg / 2)(g - g*)^2 + o(...)`, a double root. F2B Proposition 3
   does not apply to that path;
4. lines `y` slightly beyond `y*` on one side cross the boundary twice, on the
   other side not at all: the unsafe set of `R` is a **tongue** with its tip at
   `(g*, y*)` when `R_gg < 0`, and the safe set has a notch there when `R_gg > 0`;
5. the configuration persists under `C^2`-small perturbations of `R`: the tip
   solves `R = R_g = 0`, whose Jacobian in `(g, y)` is
   `[[R_g, R_y], [R_gg, R_gy]] = [[0, R_y], [R_gg, R_gy]]` with determinant
   `-R_y R_gg != 0`, so the implicit function theorem moves it continuously.

*Proof.* (1) is the regular value theorem. For (2) solve `R(g, y) = 0` for `y`
by the implicit function theorem, `y = phi(g)`, with `phi'(g*) = -R_g / R_y = 0`
and `phi''(g*) = -R_gg / R_y` from differentiating twice. (3) is Taylor's theorem
at fixed `y = y*` with `R = R_g = 0`. (4) follows from (2): `R` changes sign across
the parabola, so a line on the concave side meets it twice. (5) is stated. QED

So a narrow instability region with a tip is a **nondegenerate fold of the
boundary relative to the policy coordinate**: codimension one in the family of
policy paths (only the path through the tip is tangent) and codimension zero as a
feature of the plane (it survives perturbation). Along a pure-policy path passing
near the tip, a small change of the second coordinate turns zero crossings into
two, and `kappa` or `H` changes twice where a naive monotone picture expects at
most once. Proposition 6 is what the F7 "tongue" is tested against
(`results/F7/F7_tongue.json`); nothing in it needs, or justifies, the words
"resonance" or "Arnold tongue".

## 8. What this does and does not establish

It establishes that the local theorems of F2B/F2C are the strongest general
statements available: without a monotonicity hypothesis such as Proposition 5,
no ordering of parameter space by `kappa` is implied, and a composability map
must report components, not just labels. It does **not** say that non-monotone
behaviour is common in power systems, and it does not predict where F7 will see
it. That is empirical.
