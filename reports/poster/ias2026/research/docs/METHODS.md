# METHODS

Equations and numerical procedures, in the order the pipeline executes them.

## 1. Model and reduction

Nonlinear semi-explicit DAE, `xdot = f(x,z,theta)`, `0 = g(x,z,theta)`. At an
equilibrium both residuals vanish; `ibr_cycles.dynamics.equilibrium` certifies
`||f||_inf` and `||g||_inf` separately and returns a terminal status rather than
a silently converged point.

Linearization gives `dx_dot = fx dx + fz dz`, `0 = gx dx + gz dz`. For index-1
cases with nonsingular `gz`,

    Ared = fx - fz gz^-1 gx

`reduce_index_one` refuses the reduction when `cond(gz)` exceeds `1e12`, because
past that point the index-1 assumption itself is what fails, not the arithmetic.
Jacobians are central differences with an adaptive step `eps^(1/3) max(|v|, scale)`;
`jacobian_agreement` compares two independent estimates.

## 2. Controller self-energy

Partition `Ared` into retained `r` and hidden `h` coordinates:

    Sigma(s) = Arh (sI - Ahh)^-1 Ahr
    Teff(s)  = sI - Arr - Sigma(s)

**EXACT IDENTITY** (`reduction/schur.py`):

    det(sI - A) = det(sI - Ahh) det(Teff(s))

`Sigma` is never replaced by a constant nodal damping matrix; its frequency
dependence is the mechanism under study. `reduction/residues.py` expands
`Sigma(s) = sum_j R_j/(s - p_j)` when `Ahh` is diagonalizable and reports the
reconstruction error so a caller can reject the expansion instead of trusting it.

## 3. Action Green operator

Every action is written on the operator, not the state matrix:

    T_S(s) = T_0(s) + sum_{a in S} U_a C_a V_a^H,   T_0(s) = sI - A_0

so a state-matrix update contributes `dT_a = -dA_a`. `factorize_delta` measures
the numerical rank by SVD against a relative tolerance and stores the singular
values and the reconstruction error. **Rank one is a measurement, never a
declaration.**

By the matrix determinant lemma,

    det(T_S)/det(T_0) = det(I + C K(s)),   K(s) = V^H T_0(s)^-1 U

`K_ab(s)` reads: action `b` injects into the controller-dressed grid and is
observed by action `a`.

## 4. Isolated versus collective effects

With `M = C K` and `D_self` the block diagonal of `I + M`,

    det(I + M) = [prod_a det(I + M_aa)] * det(I + Q),   Q = D_self^-1 (I+M) - I

The individual factor is *exactly* the product of the single-action determinant
ratios, so the collective factor carries everything no single action explains.

> **Caveat (blocks of rank above one).** The block-diagonal split is a choice.
> `det(I+Q)` depends on it. The gauge-invariant objects are the cycle
> holonomies, not `Q`. For rank-one actions the split is unambiguous.

## 5. Closed intervention cycles

For a directed cycle `a1 -> a2 -> ... -> ap -> a1`,

    H_gamma = K_a1a2 K_a2a3 ... K_apa1

Under an internal basis change `U_a -> U_a S_a`, `V_a^H -> S_a^-1 V_a^H`, the
blocks transform as `K_ab -> S_a^-1 K_ab S_b`, so `H_gamma` transforms by
similarity. Only `trace`, `det` and `eig` may be reported; the entries may not.
`cycles/graph.py` builds the directed action graph from the block norms, finds
strongly connected components and enumerates simple cycles.

## 6. Connected spectral cumulants

With `Phi(theta,s) = log det(I + Theta K(s))` and one scalar activation per
action,

    d2 Phi / dtheta_a dtheta_b  = -tr(K_ab K_ba)      = -tr(H_ab)
    d3 Phi / dtheta_a dtheta_b dtheta_c
        = tr(K_ab K_bc K_ca) + tr(K_ac K_cb K_ba)     = tr(H_abc) + tr(H_acb)

Verified against central mixed differences over a step ladder; the study reports
the whole ladder, not just the best point, so truncation and roundoff regimes
stay visible. These are **infinitesimal-amplitude** quantities and are a
different object from Section 7.

## 7. Finite-amplitude Moebius interactions

For each subset, `v(S,s) = log det(T_S(s)/T_0(s))`, then

    mu(S) = sum_{R subset S} (-1)^(|S|-|R|) v(R)

Independent principal-branch logarithms are wrong exactly when the ratio winds
around the origin, so `cycles/mobius.py` follows the branch by homotopy
continuation from amplitude 0 to 1 and flags any step turning by more than
`pi/2`.

> A nonzero `mu_ABC` does not by itself prove that the pure third-order term
> caused instability. `MobiusDecomposition.truncated(members, k)` exists for
> exactly this test: compare the order-1, order-2 and exact reconstructions and
> report which one crosses the axis.

## 8. Mode provenance by winding

For a contour `Gamma` bounding a region `Omega`,

    Delta N_Omega = wind_Gamma det(T_A/T_0)

computed independently for the full ratio, the individual factor and the
collective factor. Each result carries the contour, sample count, minimum
magnitude on the contour, raw value, nearest integer, integer residual, and the
value at twice the resolution. A winding is reported as converged only when it
is integer *and* stable under refinement.

Two conventions are fixed and must not drift:

1. **The count is per eigenvalue, not per conjugate pair.** With `Omega` equal
   to the right half plane, one unstable oscillatory mode gives `+2`.
2. **`Omega` must be admissible.** `check_admissibility` requires that `Omega`
   contain at least one portfolio mode and *no* base mode and *no* mode of any
   lower-order portfolio. This is not a numerical nicety: on the toy case a
   circle of radius 0.30 around the created mode swallows the single-action mode
   of A and moves the count from the collective factor to the individual one.
   Both answers are arithmetically correct; only the admissible one is
   interpretable. `Provenance.verdict` returns `INADMISSIBLE_REGION` rather than
   a provenance label.

## 9. Minimum destabilizing order and repair

`kappa_Omega` is the smallest portfolio size producing a new unstable mode. It
is a property of the **frozen action amplitudes**, not of the system: scaling
every action changes it, so it is always reported with the amplitude vector.

Repair minimizes `||d theta||_W` subject to a spectral-abscissa or damping
target. The comparison that matters is not whether a repair exists but whether a
repair confined to the identified feedback core beats removing a useful asset or
retuning globally at equal or smaller change norm.

## 10. Open modelling decision: frozen versus re-equilibrated operating point

The low-rank action model assumes the intervention changes only the linearized
operator at a **fixed** equilibrium. For the toy this is exact. For a real
network it is not: changing a susceptance, a PLL gain or machine damping moves
the AC power flow and therefore moves `A` through `(x*, z*)`. Re-equilibrating
makes `dA` generically full rank and the determinant lemma loses its advantage.

The two readings answer different questions and may not be mixed in one claim:

- **frozen operating point** — "controller-mediated interaction between actions",
  `K(s)` exact, low-rank structure available;
- **re-equilibrated** — "engineering actions in a re-equilibrated grid", stronger
  physically, no low-rank structure. This is what the TX3 campaign did.

IEEE-39 will be built both ways as a cross-control (decision of 2026-09-09), and
every reported number carries which reading produced it.
