# F2 — Spectral composability order

A branch-independent replacement-order measure. Implemented in
`src/ibr_cycles/diagnosis/composability.py`, tested in
`tests/test_spectral_composability.py`.

## 1. Setting

A semi-explicit index-1 differential-algebraic system

    xdot = f(x, z; theta),        0 = g(x, z; theta)

with parameter `theta` in a parameter set `Theta`. `theta` may carry the
operating point, the reactive or voltage-control policy, the excitation
parameters, the machine parameters, or anything else that is held fixed while a
portfolio is evaluated.

Let `A` be a finite set of **candidate replacement actions** — in the power-system
instance, the buses at which a synchronous machine may be replaced. For a
**replacement set** `S` in a tested family `Script_S` of subsets of `A`, the
system has an equilibrium `(x*(S;theta), z*(S;theta))` and, where `g_z` is
invertible, the index-1 reduced matrix

    A_red(S; theta) = f_x - f_z g_z^{-1} g_x                                  (1)

evaluated at that equilibrium. Write `sigma(S; theta)` for its spectrum as a
multiset in the complex plane, with algebraic multiplicity.

The state dimension `n(S)` generally **depends on `S`**: replacing a synchronous
machine by a converter changes how many states the device carries. Two spectra
belonging to different `S` therefore live in different spaces and cannot be
matched eigenvalue by eigenvalue. This is not a technical nuisance; it is the
reason a label-based order is the wrong object.

## 2. The protected region and the count

Fix an open set `Gamma` in the complex plane with piecewise-C1 boundary
`dGamma` — the **protected region**, the part of the plane that must stay empty.

**Definition 1 (spectral count).**

    N_Gamma(S; theta) = # { lambda in sigma(S; theta) : lambda in Gamma }        (2)

counted with algebraic multiplicity.

**Definition 2 (spectral composability order).** With `Script_S` a finite family
of tested subsets containing the empty set,

    Delta N_Gamma(S; theta) = N_Gamma(S; theta) - N_Gamma(empty; theta)         (3)

    kappa_Gamma(theta) = min { |S| : S in Script_S, Delta N_Gamma(S; theta) != 0 }  (4)

with `min` of the empty collection written `infinity`. The sets attaining the
minimum are the **witnesses**.

In the power-system instance used throughout this project,

    Gamma = { s : Re s > 0,  0.3 Hz <= |Im s| / 2 pi <= 1.5 Hz }

with one representative kept per conjugate pair. That convention halves every
count and changes no verdict, but it is fixed once and stated.

`kappa_Gamma = infinity` is a statement about the **tested family**. With
`Script_S` the full power set of `A` it says no subset of `A` disturbs the
region; it never says anything about actions outside `A`.

## 3. Label invariance

**Proposition 1.** `N_Gamma(S; theta)`, and therefore `kappa_Gamma(theta)`,
is unchanged by

  (i) any relabelling or reordering of eigenvalues;
  (ii) any invertible change of differential-state coordinates,
       `A_red -> T^{-1} A_red T`;
  (iii) any diffeomorphic change of algebraic-variable coordinates.

*Proof.* (i) `sigma` is a multiset and (2) counts membership of a set; no
enumeration is involved. (ii) Similar matrices have the same characteristic
polynomial, hence the same spectrum with multiplicity. (iii) Let `ztilde = phi(z)`
with `J = phi'(z)` invertible, and `gtilde(x, ztilde) = g(x, phi^{-1}(ztilde))`.
Then `gtilde_ztilde = g_z J^{-1}`, `gtilde_x = g_x`, `ftilde_ztilde = f_z J^{-1}`,
so

    ftilde_ztilde gtilde_ztilde^{-1} gtilde_x
      = f_z J^{-1} (g_z J^{-1})^{-1} g_x
      = f_z J^{-1} J g_z^{-1} g_x
      = f_z g_z^{-1} g_x,

and `A_red` is not merely similar but **identical**. QED

The contrast with the observable this project used previously is the point.
`alpha_IA` required selecting one eigenvalue of the replaced system by shape
similarity to a reference eigenvector. That selection depends on the inner
product used, on which coordinates are retained in the comparison, and — as
`docs/V2B_RESULTS.md` records — it is not even single-valued when the reference
mode has two descendants of nearly equal assurance. Definition 2 contains no
eigenvector, no matching and no ordering, so none of those failure modes exist.

## 4. Composability order is NOT interaction order

Let `F : Script_S -> R` be a set function and let its Moebius transform be

    mu(S) = sum over R subset of S of (-1)^{|S| - |R|} F(R).                    (5)

The **irreducible connected interaction order** is the largest `|S|` carrying a
non-vanishing `mu(S)`. This is the quantity the earlier Track-A work reported as
"orders 1 to 3 stable, exact order 4 unstable".

**These are logically independent.** Both implications fail.

**Counterexample A — `kappa = 3` with no interaction at any order.** Three
actions, and a set function that is exactly additive:

| set | value | Moebius term |
|---|---|---|
| empty | −3.0 | −3.0 |
| each single | −1.8 | **+1.2** |
| each pair | −0.6 | **0** |
| the triple | **+0.6** | **0** |

Every term of order two and three vanishes, so there is no joint effect
whatsoever. Yet nothing of size one or two crosses zero and the triple does, so
under the threshold count `kappa_Gamma = 3` while the interaction order is 1.
`kappa` is a **threshold** notion, not an interaction notion.

**Counterexample B — `kappa = 1` with a large third-order term.** Take
`F(empty) = -1`, `F({1}) = +1`, `F({2}) = F({3}) = -1`, `F({1,2}) = F({1,3}) = +1`,
`F({2,3}) = -1`, `F({1,2,3}) = +5`. Then `mu({1,2,3}) = +4`, a large genuine
triple interaction, while a single action already changes the count, so
`kappa_Gamma = 1`.

Both are in `tests/test_spectral_composability.py`.

**Why the conflation went unnoticed.** On the Track-A flagship the two coincide:
all fifteen proper subsets stay out of `Gamma`, the portfolio enters it, *and*
the fourth-order Moebius term is large (`mu_4 = +0.433`). That is a genuine
coincidence of two independent properties at one parameter point, not an
identity. The third test in the file pins exactly that.

**They also answer different questions.** `kappa_Gamma` answers *how many
simultaneous replacements a planner may make before the protected region is
disturbed*. The interaction order answers *how much of the change is genuinely
joint rather than accumulated*. A planner needs the first; a mechanism claim
needs the second. Reporting one and calling it the other is the error this
section exists to prevent.

## 5. Conventions that must be fixed once

- **The region.** `Gamma` is part of the definition. A different band or a
  different half-plane offset defines a different order.
- **The tested family.** `kappa_Gamma = infinity` is relative to `Script_S`.
- **Conjugate pairs.** Counted once here.
- **Zero modes.** The two reference zeros of this DAE are excluded by a magnitude
  floor; they lie outside `Gamma` in any case.
- **Fixed structure along a parameter path.** The continuity results of F2B need
  `n(S)` constant as `theta` varies. A policy coordinate that switches a
  controller state in or out of existence violates this. The practical
  consequence for F7 is stated there: parameterise the reactive policy *inside* a
  fixed converter structure, with a gain going to zero, rather than by adding or
  removing an integrator.
