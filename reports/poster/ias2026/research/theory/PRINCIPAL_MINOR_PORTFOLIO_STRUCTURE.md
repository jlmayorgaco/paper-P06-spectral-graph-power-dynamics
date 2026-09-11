# Principal-minor structure of the portfolio family

Status: statements 1–5 are proved below (linear algebra). The numerical
verification is BC02b (`results/BC/BC02b/BC02b_summary.json`). The certificate
question in §6 is tested in FC12 (BC03 gate). Nothing here is claimed as a new
theorem of linear algebra. The identities are classical; what is new is their
use as the organizing object of the replacement family.

## 0. Setting and the correct Boolean object

Candidates `V = {1, …, m}`, each acting on its bus through `p_i = 2` real
voltage channels. The hierarchy established in Phase I (BC02, BC02b) is:

    descriptor pencil (affine in delta)
      -> port / action factorization
      -> principal submatrices M_SS
      -> det(I + M_SS)
      -> transverse spectral count, H, kappa

At a common operating point (matched dispatch; the C2 common realization of
BC02):

    T_S(s) = T_0(s) + sum_{i in S} E_i dY_i(s) E_i^T      (exact, BC02: 1e-14)
    M(s)   = D(s) K(s),  D = blkdiag(dY_i),  K = E^T T_0(s)^{-1} E
    det T_S / det T_0 = det(I + M_SS)                      (Sylvester)
    det P_S = h_S det T_S,   h_S = h_0 prod_{i in S} kappa_i(s)   (multiplicative)

The reduced matrix `A(delta)` is **not** used as a Boolean object. It is
rational in `delta` and has a full-degree Möbius spectrum (BC02b).

## 1. Principal-minor expansion (classical)

For any square `X` of order `n`,

    det(I + X) = sum_{U subset {1..n}} det X_UU     (det of the empty minor = 1).

*Proof.* Multilinearity of the determinant in the columns of `I + X`. Each
column is `e_k + X_k`; choosing `X_k` for `k in U` and `e_k` otherwise, the
determinant is `det X_UU`. ∎

## 2. Möbius coefficients are block-touching minor sums

Let `r(S) = det(I + M_SS)` for `S subset V`, and let `blocks(U)` be the set of
actions whose channels meet `U`. Then

    r(S) = sum_{T subset S} mu_T,    mu_T = sum_{U : blocks(U) = T} det M_UU,

and `mu_T` is the Möbius transform of `r`:
`mu_T = sum_{R subset T} (−1)^{|T|−|R|} r(R)`.

*Proof.*

1. Apply §1 to `X = M_SS`, whose index set is the channels of `S`.
2. Group the principal minors by the set of blocks they touch. This gives the
   first identity.
3. The Möbius transform on the Boolean lattice is the unique inverse of the zeta
   (sum-over-subsets) transform, so the coefficients of that expansion are the
   Möbius coefficients. ∎

**Numerical check.** Vertex values against `det(I + M_SS)`: ≤ 1.3e-13.
Möbius coefficients against the minor sums: ≤ 1.3e-13 (BC02b, three families,
four values of `s`).

## 3. Representation invariance

Let `S_b = blkdiag(S_1, …, S_m)` be any admissible port-basis change (one
invertible 2×2 per action), under which `M -> S_b^{-1} M S_b`.

1. `r(S)` is invariant for every `S`, because `det(I + (S_b^{-1} M S_b)_SS) =
   det(S_S^{-1} (I + M_SS) S_S)`.
2. Hence every Möbius coefficient `mu_T` is invariant: it is a fixed linear
   combination of invariant numbers.
3. The *individual* scalar minors `det M_UU` with `U` splitting a block are
   **not** invariant. Only their block-touching sums are.
4. The normalized collective factor `det(I + Q_SS)` with
   `Q = blkdiag(I + M_ii)^{-1}(I + M) − I` is invariant. Its coefficients satisfy
   `mu~_{i} = 0` for every singleton, because `Q_ii = 0`: the individual effects
   are exactly factored out (eq. 20 of the v1 note).
5. The "connected" coefficients `c_T`, defined as the Möbius transform of
   `log det(I + M_SS) = tr log(I + M_SS)`, are invariant. They collect the
   closed walks of `M` whose block support is exactly `T`: the cycle
   (holonomy) terms of E18.

## 4. Exact polynomial degree

`r` is a multilinear polynomial in the indicators `delta` of exact degree
`d* = max{|T| : mu_T != 0}`, with `d* <= m`. Generically `d* = m`, because the
fully-touching minor sum has no structural reason to vanish. BC02b measured
full-order coefficients at every tested `s`:

| benchmark | full-order coefficient relative to `r_0 = 1` |
|---|---|
| IEEE-39 P4, order 4 | 0.012–0.062 |
| Kundur, order 3 | 0.18–1.31 |
| IEEE-68 candidates, order 4 | 0.028–0.21 |

**No low-degree representation is exact on these benchmarks.** The IEEE-39
coefficients decay with order: 0.53, 0.20, 0.05, 0.012 at `s = 0.5 + 2j`. That
decay is structure a certificate may exploit. It is not degree reduction.

## 5. Algebraic interaction order is not kappa

`kappa` (the minimum destabilizing cardinality) and `d*` (the algebraic
interaction degree) are different quantities. Neither bounds the other in
general:

- **Additive but high kappa.** A scalar family that is affine in `delta` has
  `d* = 1`, yet `kappa = m` (proposition T2 of the v1 note).
- **Full degree, low kappa.** A family can have `d* = m` and `kappa = 1`: a
  destabilizing singleton together with arbitrary higher-order terms.

The fourth-order Möbius statement of Track A (`mu_4 = +0.433` of the abscissa)
is a statement about one functional. It is not implied by `kappa = 4` and does
not imply it.

## 6. The certificate question

> Can low-order principal-minor information certify the exclusion of
> higher-order closure, under additional structure?

The candidates tested in FC12 use only first- and second-order blocks of `Q`
(`Q_ij`, `i != j`) to bound every `det(I + Q_SS)`.

| candidate | condition on the contour | conclusion |
|---|---|---|
| small gain / Perron–Frobenius (T3) | `rho(Rbar) < 1`, with `Rbar_ij >= sup ‖Q_ij‖` | all subsets safe |
| top-sum comparison (T3a) | `c_r(d) < 1` | every portfolio of size ≤ r safe, so `kappa >= r + 1` |
| H-matrix / block comparison matrix | comparison matrix of `I + Q_SS` is an M-matrix | nonsingular for all S; the same bound as Perron-weighted T3a |
| block Gershgorin (row sums) | `max_i sum_{j != i} ‖Q_ij‖ < 1` | the unweighted special case |

All of them need two things:

- **(a)** individual actions transversely stable, with `N_h = 0` (BC02 ledger);
- **(b)** bounds valid on the **whole** contour.

A frequency grid gives lower bounds of the suprema, so a gridded result is a
**SCREENING**, never a **CERTIFICATE**.

**Structural obstacle identified before running.** `K = E^T T_0^{-1} E` has a
pole at `s = 0`: `T_0` is singular there because of the center subspace.
Whether that pole cancels in `Q` is an empirical question for each benchmark.
If it does not, every `H_inf` bound is infinite and the gate is vacuous, which
is a valid negative result. Relocating the center subspace (the relocated port)
removes the singularity but couples all devices through `R_x`, which destroys
the locality of §0. The two cannot be combined.
