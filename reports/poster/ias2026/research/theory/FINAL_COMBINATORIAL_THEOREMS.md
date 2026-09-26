# Combinatorial theorems of the final campaign (proofs and literature status)

Notation:

- `V = {1, …, m}` is the candidate set, `S ⊆ V` a portfolio, `delta in {0,1}^m`
  its indicator.
- `alpha(S)` is the transverse spectral abscissa (`theory/TRANSVERSE_STABILITY_QUOTIENT.md`).
- "Unsafe" means `alpha(S) >= 0`, i.e. not transversely stable.
- `H` is the family of minimal unsafe portfolios and `kappa = min |H|`.

## 1. Complexity of minimum-cardinality Boolean spectral destabilization

**Problem MCBSD(k).** Given rational matrices defining a family `A(delta)` and an
integer `k`, is there `S` with `|S| <= k` and `alpha(A(delta_S)) >= 0`?

**Theorem 1.** MCBSD is NP-complete already for the symmetric, degree-2
(multilinear) family

    A_G(delta) = D_delta W D_delta − (k − 1) I,

where `W` is the adjacency matrix of a simple graph `G` and `D_delta = diag(delta)`.

*Proof.*

1. **Structure.** For `delta = delta_S`, `A_G` is block diagonal:
   `W[S] − (k−1) I` on `S` and `−(k−1) I` on the complement. Hence
   `alpha = max(lambda_max(W[S]), 0) − (k−1)`, and for `k >= 2` the condition
   `alpha >= 0` is equivalent to `lambda_max(W[S]) >= k − 1`.
2. **A clique gives `alpha >= 0`.** A graph on `s` vertices has
   `lambda_max <= Delta_max <= s − 1`. If `G[S]` contains `K_k` with `|S| = k`,
   then `G[S] = K_k` and `lambda_max = k − 1`. (If `S` contains a `k`-clique and
   `|S| <= k`, then `S` is exactly that clique.)
3. **`alpha >= 0` gives a clique.** Suppose `|S| = s <= k` and
   `lambda_max(W[S]) >= k − 1`. Then `s − 1 >= lambda_max >= k − 1`, so `s = k`
   and `lambda_max = s − 1`. The Perron value `s − 1` equals the maximum degree
   only for a connected `(s−1)`-regular graph, i.e. `K_s`. So `G[S] = K_k`.
4. **Reduction.** `G` has a `k`-clique iff MCBSD(k) holds for `A_G`. The map
   `G -> (W, k)` is polynomial, and CLIQUE is NP-complete (Karp 1972), so
   MCBSD is NP-hard.
5. **Membership in NP.** Given `S`, `alpha >= 0` holds iff `(k−1) I − W[S]` is
   not positive definite. This is decidable in polynomial time by exact
   rational `LDL^T` (Bareiss-type elimination with polynomially bounded
   entries). ∎

The minimization version (`kappa` itself) is therefore NP-hard.

**Literature status — NOT new.** The reduction is the classical clique
reduction for cardinality-constrained principal-submatrix eigenvalue
maximization, i.e. the NP-hardness of sparse PCA (Magdon-Ismail, *Inf. Process.
Lett.* 2017; also Moghaddam–Weiss–Avidan). Robust stability and robust
nonsingularity of affinely parameterized and interval families are also
NP-hard (Poljak–Rohn 1993; Nemirovskii 1993). Structured and sparse stability
radii are studied in Katewa–Pasqualetti (*Automatica* 2020) and related work
under continuous structured perturbations.

The restricted binary / symmetric / degree-2 statement is a direct corollary of
the sparse-PCA reduction, so **no priority is claimed**. Its role in the paper is
explanatory only: exact minimal-witness search over replacement portfolios is,
in general, combinatorially hard. It does **not** say that IEEE-39 itself is
hard: `m = 4` there, and exhaustive search takes 16 evaluations.

## 2. A monotone tractable class

**Theorem 2.** Suppose there is one invertible `T` (common to all `S`) such that
every `B(S) = T^{-1} A(S) T` is Metzler, and `S ⊆ S'` implies
`B(S) <= B(S')` entrywise. Then `S ⊆ S'` implies `alpha(A(S)) <= alpha(A(S'))`.

*Proof.* `alpha` is similarity invariant. For a Metzler `B` and
`c >= max_i |B_ii|`, the matrix `B + cI` is nonnegative and, by
Perron–Frobenius, `alpha(B) = rho(B + cI) − c`. From `0 <= B(S) + cI <=
B(S') + cI` and the monotonicity of the spectral radius on nonnegative matrices,
`rho(B(S) + cI) <= rho(B(S') + cI)`. ∎

**Corollary 2.1 (heredity).** Stability of the full portfolio implies stability
of every subset. The safe family is downward closed, and `H` is its set of
minimal non-faces.

**Corollary 2.2 (no re-stabilization).** If `alpha(S) >= 0` and `S ⊆ S'`, then
`alpha(S') >= 0`. Adding a replacement cannot re-stabilize an unstable
portfolio. Superset pruning is then valid **for stability labels**, and
"final-state stable", "a safe path exists" and "any order is safe" coincide.

**Status.** Classical: Perron–Frobenius monotonicity of the spectral abscissa of
Metzler matrices. Not new; it is the positive-systems literature applied to
portfolios.

**Diagnostic.** One robust pair `S ⊂ S'` with `alpha(S) > 0 > alpha(S')` at the
same policy excludes **every** common order-preserving Metzler realization,
because `alpha` is similarity invariant (FC10). Such an example does **not**
show that "graph frustration causes instability". It only shows that the
benchmark is outside this class.

## 3. Planning hierarchy (final state, safe path, any order)

For a target `T`:

- **A**: `alpha(T) < 0`.
- **B**: there is a chain `∅ = S_0 ⊂ S_1 ⊂ … ⊂ S_{|T|} = T`, adding one element
  at a time, with every `S_j` stable.
- **C**: every `R ⊆ T` is stable.

**Proposition 3.** C ⇒ B ⇒ A. Neither converse holds in general.

*Proof.* C ⇒ B: any chain works. B ⇒ A: the last element of the chain is `T`.
Counterexamples to the converses:

- the diagonal family of v1 §12.2: `A_0 = −I`, `A_{1} = diag(2, −3)`,
  `A_{2} = diag(−3, 2)`, `A_{12} = −2I`. It is A but not B;
- any target containing an unstable subset that is avoidable along some other
  chain is B but not C. ∎

**Theorem 3.1 (v1 T6).** With every portfolio feasible and the base stable,
`S` is in `C` iff no hyperedge of `H` is contained in `S`. `C` is the largest
hereditary subfamily of the safe family, and it is exactly the set of targets
whose every one-at-a-time implementation order passes only through stable
stationary states.

The model results for IEEE-39 are in FC10 (census) and FC11/FC12
(16-subset lattice, spectral and nonlinear).

## 4. Resilience complex (nonlinear composability)

For a disturbance family `d` with nested sets `D_rho`, and `r_S` the critical
radius:

    K(rho) = { S : r_R > rho for every R ⊆ S }.

**Proposition 4.**

1. `K(rho)` is downward closed, hence an abstract simplicial complex on `V`.
2. Its minimal non-faces are exactly the minimal non-tolerating portfolios,
   `H_NL(rho)`.
3. For `H in H_NL(rho)` with `|H| = k >= 2`, the induced subcomplex `K(rho)[H]`
   is the boundary of the simplex `Delta^{k−1}`, a triangulated `(k−2)`-sphere.
   Its reduced homology is the field in degree `k − 2` and zero otherwise.
4. `rho_2 > rho_1` implies `K(rho_2) ⊆ K(rho_1)` (a decreasing filtration), and
   `R_k = min_{|S| <= k} r_S` is non-increasing in `k`.
5. As `rho -> 0+`, `H_NL(0+) = H_RHP_perp`, provided every transversely stable
   portfolio tolerates small disturbances.

*Proof.*

1. Downward closure holds by definition: the condition is imposed on all
   subsets.
2. `S` fails to be a face iff some `R ⊆ S` fails to tolerate `rho`. A
   non-face all of whose proper subsets are faces is therefore a minimal
   non-tolerating set, and conversely.
3. Every proper subset of `H` is a face and `H` is not, which is the boundary of
   the simplex. Its homology is standard.
4. Follows from `r_S > rho_2 > rho_1`. For `R_k`, the minimum is taken over a
   larger set.
5. By the definition of `r_S = 0` for transversely unstable `S`, and by local
   asymptotic stability of the transversely stable equilibria, which gives
   `r_S > 0`. ∎

**Status.** Standard combinatorial topology (Stanley–Reisner correspondence
between a complex and its minimal non-faces; sublevel filtrations as in
persistent homology). It is an exact organizing consequence of the safety
definition, not new mathematics. The **global** Betti numbers of `K(rho)` are
not implied by the local spheres and are reported separately (FC06).
