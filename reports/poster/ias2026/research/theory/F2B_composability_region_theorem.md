# F2B — The composability-region theorem

Companion to `theory/spectral_composability_order.md`. Numerically corroborated
by `experiments/F2_composability_continuation.py`.

## 1. Regular parameter set

Fix the protected region `Gamma` and the tested family `Script_S`, both finite
and fixed. Let `Theta_reg` be a subset of parameter space such that for **every**
`S` in `Script_S` and every `theta` in `Theta_reg`:

- **(H0) fixed structure.** The state dimension `n(S)` does not depend on
  `theta`, and `f`, `g` are C1 in `(x, z, theta)`.
- **(H1) equilibrium branch.** An equilibrium `(x*(S; theta), z*(S; theta))`
  exists and depends continuously on `theta`.
- **(H2) index-1 regularity.** `g_z(S; theta)` is invertible at that equilibrium.
- **(H3) boundary avoidance.** `sigma(S; theta)` has no point on `dGamma`.

(H2) excludes impasse surfaces; (H1) excludes equilibrium folds and the
non-existence region; (H3) is the condition that actually fails at a transition.

## 2. Theorem

**Theorem 1 (composability regions).** Under (H0)–(H3), for each fixed `S` the
count `N_Gamma(S; ·)` is locally constant on `Theta_reg`, hence constant on every
connected component of `Theta_reg`. Consequently `kappa_Gamma(·)` is constant on
every connected component of `Theta_reg`, and parameter space decomposes as

    Theta_reg = union over k of P_k,      P_k = { theta : kappa_Gamma(theta) = k },

each `P_k` a union of connected components, with `P_infinity` the set where no
tested subset disturbs `Gamma`.

*Proof.* Fix `S` and `theta_0` in `Theta_reg`.

**(a) Continuity of the reduced matrix.** By (H1) the equilibrium is continuous
in `theta`; by (H0) the Jacobian blocks `f_x, f_z, g_x, g_z` are continuous in
`(x, z, theta)` and therefore continuous in `theta` along the branch. By (H2)
`g_z` lies in the general linear group, on which inversion is continuous. Hence
`theta -> A_red(S; theta)` given by (1) is continuous at `theta_0`, and its size
is fixed by (H0).

**(b) A contour that stays valid.** Eigenvalues of a matrix are bounded by its
norm, so choose `R` larger than `sup` of `||A_red(S; theta)||` over a compact
neighbourhood `V` of `theta_0`, and set `Gamma_R = Gamma intersect { |s| < R }`.
No eigenvalue has modulus `R` or more on `V`, so `N_Gamma = N_{Gamma_R}` there,
and `dGamma_R` is a compact piecewise-C1 curve.

**(c) The count is a contour integral.** By the argument principle applied to the
characteristic polynomial `p(s; theta) = det(s I - A_red(S; theta))`, whose zeros
are the eigenvalues with algebraic multiplicity,

    N_{Gamma_R}(S; theta) = (1 / 2 pi i) * contour integral over dGamma_R
                            of  p'(s; theta) / p(s; theta)  ds.                 (6)

**(d) The integrand is continuous.** `p` and `p_s` are polynomials in `s` whose
coefficients are continuous in `theta` by (a). By (H3) at `theta_0` no eigenvalue
lies on `dGamma`, so `|p|` has a positive minimum on the compact set `dGamma_R`
at `theta_0`; by continuity and compactness there is a neighbourhood `U` of
`theta_0` inside `V` on which that minimum stays positive. On `U` the integrand
of (6) is jointly continuous and the contour is fixed, so the integral is
continuous in `theta`.

**(e) Integer plus continuous equals locally constant.** The left side of (6) is
an integer for every `theta` in `U`. A continuous integer-valued function on a
connected neighbourhood is constant, so `N_Gamma(S; ·)` is constant on `U`.

**(f) From local to global, and from counts to kappa.** Local constancy on a
connected component gives constancy there by a standard connectedness argument
(the level set is open and closed and non-empty). Finally `kappa_Gamma` is
determined by (4) from the finitely many integers `{N_Gamma(S; theta)}`, each
constant on the component; a function of finitely many constants is constant. QED

**Remark.** Nothing in the proof is deep. It is the argument principle plus
"integer and continuous implies constant". The content of F2B is not the
technique but the **object**: that a replacement-portfolio order defined this way
is forced to be piecewise constant, so a composability *phase diagram* is a
well-posed thing to draw rather than an artefact of gridding.

## 3. Where the order can change

**Corollary 2 (boundary characterisation).** If `theta_1` and `theta_2` lie in
`Theta_reg` and `kappa_Gamma(theta_1) != kappa_Gamma(theta_2)`, then every
continuous path from `theta_1` to `theta_2` leaves `Theta_reg`. If (H0)–(H2) hold
along the path, the failure is (H3): some `S` in `Script_S` has an eigenvalue on
`dGamma` at some point of the path.

So the order changes **only across spectral boundaries**, and the boundary is
carried by an identifiable witness subset.

**Proposition 3 (generic codimension one).** Suppose at `(omega*, theta*)` a
subset `S` has a **simple** eigenvalue `lambda(theta*) = i omega*` on the
imaginary-axis part of `dGamma`, with `omega*` in the open band. Then near
`theta*` there is a C1 branch `lambda(theta)` with `lambda(theta*) = i omega*`,
and if the **transversality** condition

    grad_theta ( Re lambda ) (theta*) != 0                                      (7)

holds, the set `{ theta : Re lambda(theta) = 0 }` is locally a C1 hypersurface of
codimension one in `Theta`.

*Proof.* Simplicity gives `p_s(i omega*; theta*) != 0`, so the implicit function
theorem applied to `p(lambda; theta) = 0` yields a C1 branch `lambda(theta)`. A
complex-conjugate pair moves together, so `Re lambda = 0` is a **single** real
equation. Under (7) the regular value theorem makes its zero set a codimension-one
C1 manifold. QED

The same argument applies at the band edges, where the single real equation is
`|Im lambda| = omega_lo` or `omega_hi`.

## 4. The port form of the same boundary

Let `T_0(s; theta)` be the baseline port operator and `T_S` the operator after
replacing the ports in `S`. Where the factorisation is regular — that is, where

  (P1) `T_0(s)` is invertible,
  (P2) each individual block `I + M_aa(s)` is invertible,
  (P3) no pole-zero cancellation occurs between `det T_S` and `det T_0`,
  (P4) the port update has the assumed rank at every replaced port —

the frozen factorisation gives

    det T_S(s) / det T_0(s)
      = [ product over a of det( I + M_aa(s) ) ] * det( I + Q_S(s) ).           (8)

Hence a zero of `det T_S` that is **not** a zero of `det T_0` and not attributable
to an individual block satisfies

    det( I + Q_S(i omega; theta) ) = 0,  equivalently  -1 in sigma( Q_S(i omega; theta) ).  (9)

**Corollary 4.** Under (H0)–(H3) and (P1)–(P4), and when the crossing mode is
visible at the action ports, the boundary of Proposition 3 coincides with the
locus (9): the collective composability boundary is where an interaction
eigenvalue reaches `-1`.

`sigma(Q)`, `det(I + Q)` and `dist(-1, sigma(Q))` are invariant under the per-port
basis changes that leave the physical operator unchanged, because those act on
`Q` by similarity. `||Q||` and `sigma_min(I + Q)` are **not** invariant and may
not be used. That is the frozen E17/N8c result and F4 will restate it as a
theorem.

## 5. Exceptions, stated rather than buried

| exception | what breaks |
|---|---|
| **multiple simultaneous crossings** | two or more witnesses reach `dGamma` together; `kappa` may jump by more than one and the boundary set need not be a manifold at the intersection |
| **defective eigenvalues** | at a non-semisimple eigenvalue the branch is only Puiseux-continuous, not C1; Proposition 3 fails, though Theorem 1 survives because it needs only continuity of the count |
| **DAE singularity** | `g_z` singular, (H2) fails, `A_red` undefined |
| **equilibrium fold or non-existence** | (H1) fails — in the power-system instance this is the non-dispatchable region, which is why it is excluded explicitly rather than counted as stable |
| **structure change along the path** | (H0) fails if a controller state appears or disappears; the reduced matrix is then not even the same size, and continuity is meaningless |
| **baseline pole on the contour** | `det T_0(i omega) = 0`, (P1) fails and (8) is undefined — this is exactly why the closure is evaluated on the imaginary axis and never at an eigenvalue of either system |
| **normalisation singularity** | some `det(I + M_aa) = 0`, (P2) fails and the individual/collective split is undefined |
| **port-invisible modes** | a crossing mode not observable at the action ports changes `N_Gamma` without any eigenvalue of `Q` reaching `-1`; (9) then detects nothing. F6 constructs this case |
| **finite tested family** | `kappa = infinity` means "no tested subset", never "no subset" |

## 6. Relation to generalised Nyquist — what is and is not new

**Not new.** The condition `det(I + L(i omega)) = 0`, equivalently
`-1 in sigma(L(i omega))`, is the classical generalised Nyquist / return-difference
criterion, and the eigenloci of a return ratio are standard multivariable control.
Impedance-based stability analysis applies the same object to a source/load
partition. Nothing in Corollary 4 claims otherwise, and the `-1` crossing itself
is explicitly **not** claimed as a contribution.

**Not new either.** Theorem 1 is the argument principle plus a connectedness
argument. Any competent reader will regard the proof as routine.

**What is new** is the composite object and its use:

1. the integer `kappa_Gamma` defined on a **combinatorial family of
   interventions** against a common baseline, rather than on a single loop;
2. the resulting **partition of physical policy space** into composability
   regions `P_k`, together with the statement that the partition is well posed —
   piecewise constant with codimension-one boundaries — so that a composability
   phase diagram is a mathematical object and not a gridding artefact;
3. the identification of the **witness subset** carrying each boundary, which is
   what a planner acts on;
4. the exact **individual/collective factorisation** (8), which separates "this
   one plant is marginal" from "these four are jointly marginal" — generalised
   Nyquist applied to the whole system gives a yes/no answer and does not, by
   itself, produce the subset structure;
5. the empirical result that the order is **movable by physical policy**, which is
   a statement about `P_k` and has no counterpart in a single-portfolio analysis.

An honest reviewer response is that (1)–(4) are a reformulation and organisation
of known ingredients, and that the substance is (5) plus the fact that the
reformulation makes (5) expressible. F10 will test claim (4) against a
generalised-Nyquist baseline implemented on the same flagship; if that baseline
recovers the subset structure as cheaply, the novelty claim reduces to (5).

## 7. Numerical corroboration

`experiments/F2_composability_continuation.py`, 61 points, excitation-gain scale
from 0.50 to 2.00 in steps of 0.025 at a fixed time-scale of 1.5, every point
base-stable, all sixteen subsets evaluated at every point.

| | |
|---|---|
| values of `kappa_Gamma` seen | **4, 3, 2** |
| maximal constant runs | **3** |
| transitions | **2** |
| transitions explained by a witness reaching `dGamma` | **2 of 2** |

| transition | gain scale | witness subset |
|---|---|---|
| `kappa` 4 → 3 | 1.125 → 1.150 | **30+33+35**, a triple crosses |
| `kappa` 3 → 2 | 1.800 → 1.825 | **30+33**, a pair crosses |

Exactly one witness at each transition, so both are simple crossings and
Proposition 3 applies. The order falls as the excitation gain rises, because
progressively smaller subsets become sufficient — which is the mechanism the
partition is meant to expose.

Median closure margin by region: `0.376` at `kappa = 2`, `0.273` at `kappa = 3`,
`0.0713` at `kappa = 4`, reproducing on this independent one-dimensional sweep
the ordering F1B found on the two-dimensional map. This remains an **empirical
observation**; no theorem here predicts it, and section 7 of the F7 brief is the
place it will be tested rather than asserted.
