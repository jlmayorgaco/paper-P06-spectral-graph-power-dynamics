# Finite-amplitude connected cumulants of the replacement characteristic function

Status: definitions and proofs (Phase 1 of the 2026-09-11 extension). All
statements are checked numerically in `tests/test_connected_cumulants.py` (74
tests) and `experiments/connected_cumulants/CC02_toy_falsification.py`
(`results/CC/CC02/`).

The code is in `src/ibr_cycles/cycles/connected.py`. This object is **not** the
infinitesimal log-det cumulant of `src/ibr_cycles/cycles/cumulants.py` (§6).

**Classical ingredients, credited:**

- the Möbius function of the set-partition lattice (Rota 1964);
- the moment–cumulant formula (Leonov–Shiryaev 1959; Speed 1983);
- the exponential formula for set partitions (Stanley, *EC2*, §5.1);
- the Leibniz expansion and the cycle decomposition of permutations.

What is specific here is the *application* to the replacement set function, plus
four elementary consequences:

- the connectivity theorem for minimal coalitions (§2.3);
- the boundary identity (§2.4);
- the block-port expansion with its refutation of the naive holonomy formula (§3.2);
- the separation of three orders (§4).

None of these is claimed as deep new mathematics.

## 0. Setting

- `V = {1,…,m}` is the set of candidate replacements. Each candidate acts
  through a port with `p` channels (`p = 1` scalar, `p = 2` for the phasor
  ports of the benchmark).
- `Q(s)` is the normalized network-closure operator of FC18:
  `Q = (I + D K_d)^{-1} D K_o`. It has zero `p×p` diagonal blocks.
- The collective characteristic set function is

      F_s(S) = det(I + Q_SS(s)),   F_s(∅) = 1.

- The full characteristic function factors as

      det T_S / det T_0 = Π_{i∈S} det(I + M_ii) · F_s(S)     (FC18, Sylvester).

- `Π(S)` is the set of set partitions of `S`, and `|π|` is the number of blocks.

**Definition (connected cumulant).** For non-empty `S ⊆ V`,

    chi_S(s) = Σ_{π∈Π(S)} (−1)^{|π|−1} (|π|−1)! Π_{B∈π} F_s(B).

The coefficient is the Möbius function `μ_Π(π, 1̂)` of the partition lattice.
For a singleton, `chi_{i} = F(i) = det(I + Q_ii) = 1`. Every statement in §§1–2
holds for an arbitrary set function `F` with `F(∅) = 1` over a commutative
`Q`-algebra, and in particular pointwise in `s`.

## 1. Moment–cumulant inversion (C1)

Let `Λ_V = C[x_v : v ∈ V] / (x_v^2 : v ∈ V)`, with `x^A = Π_{v∈A} x_v`. Then
`x^A x^B = x^{A∪B}` if `A ∩ B = ∅`, and `0` otherwise.

**Lemma 1.1 (exponential formula; classical).** Put
`Φ(x) = Σ_{A⊆V} F(A) x^A`, so that `Φ − 1` is nilpotent. Then

    log Φ = Σ_{∅≠S⊆V} chi_S x^S.

*Proof.*

1. `log Φ = Σ_{k≥1} (−1)^{k−1}/k (Φ − 1)^k`. The sum is finite because
   `(Φ−1)^{m+1} = 0`.
2. The coefficient of `x^S` in `(Φ − 1)^k` is a sum over ordered `k`-tuples of
   non-empty, pairwise disjoint sets with union `S`. Overlapping products
   vanish, and disjoint ones multiply to `x^S`.
3. Each unordered partition `π` of `S` with `|π| = k` arises from exactly `k!`
   ordered tuples. Hence `[x^S](Φ−1)^k = k! Σ_{π∈Π(S),|π|=k} Π_{B∈π} F(B)`.
4. Multiplying by `(−1)^{k−1}/k` and summing over `k` gives `chi_S`. ∎

**Theorem C1 (inversion).** For every non-empty `S`,

    F(S) = Σ_{π∈Π(S)} Π_{B∈π} chi_B.

The family `(chi_B)` is the **unique** family with this property.

*Proof.*

1. Put `Ψ = log Φ = Σ chi_S x^S`. In a commutative `Q`-algebra,
   `exp(log(1+u)) = 1 + u` for nilpotent `u` (a formal power-series identity
   with finitely many non-zero terms). So `Φ = exp Ψ = Σ_k Ψ^k / k!`.
2. As in Lemma 1.1, `[x^S] Ψ^k = k! Σ_{|π|=k} Π chi_B`. Summing over `k` gives
   the stated identity.
3. Uniqueness: the identity reads `F(S) = chi_S + (terms in chi_B with |B| < |S|)`.
   This system is triangular with unit diagonal, so it determines `chi` from
   `F` recursively. ∎

**Lemma 1.2 (recursion used in the code).** Fix `i_0 ∈ S`. Then

    F(S) = Σ_{B : i_0 ∈ B ⊆ S} chi_B F(S ∖ B).

Hence `chi_S = F(S) − Σ_{i_0∈B⊊S} chi_B F(S∖B)`.

*Proof.* Group the partitions of `S` by the block `B` containing `i_0`. The
remaining blocks form an arbitrary partition of `S ∖ B`, and by C1 their sum is
`F(S ∖ B)` (with `F(∅) = 1`). ∎

**Numerical check (CC02).**

- Inversion residual ≤ 1.9e-13 (relative) over 120 random set functions,
  `n ≤ 6`.
- Recursion versus the partition definition ≤ 1.6e-15.

## 2. Factorization annihilation (C2) and minimal coalitions

**Theorem C2.** Let `S = S_1 ⊔ S_2` with both parts non-empty. Suppose
`F(A) = F(A ∩ S_1) F(A ∩ S_2)` for every `A ⊆ S`. Then `chi_T = 0` for every
`T ⊆ S` that meets both `S_1` and `S_2`; in particular `chi_S = 0`.

*Proof.*

1. Let `Φ_R = Σ_{A⊆R} F(A) x^A`. The factorization gives
   `Φ_S = Φ_{S_1} Φ_{S_2}`: expand the product over pairs `(A_1 ⊆ S_1, A_2 ⊆ S_2)`
   and use `x^{A_1} x^{A_2} = x^{A_1 ∪ A_2}`.
2. Both factors are `1 + nilpotent` and commute, so
   `log Φ_S = log Φ_{S_1} + log Φ_{S_2}`.
3. `log Φ_{S_i}` contains only monomials `x^T` with `T ⊆ S_i`.
4. By Lemma 1.1 (applied to `F` restricted to subsets of `S`), the coefficient of
   `x^T` in `log Φ_S` is `chi_T`. Mixed `T` therefore have coefficient 0. ∎

**Proposition 2.2 (converse).** If `chi_T = 0` for every `T ⊆ S` meeting both
`S_1` and `S_2`, then `F(A) = F(A ∩ S_1) F(A ∩ S_2)` for all `A ⊆ S`.

*Proof.* Then `Ψ_S = Ψ_{S_1} + Ψ_{S_2}`, so
`Φ_S = exp Ψ_{S_1} · exp Ψ_{S_2} = Φ_{S_1} Φ_{S_2}`. Compare coefficients. ∎

So exact factorization across `{S_1, S_2}` is equivalent to the vanishing of
**all** mixed cumulants.

The vanishing of the single top cumulant `chi_S` is strictly weaker. For the path
coupling `1–2–3–4` (§4, example b), `chi_1234 ≡ 0`, yet `F` does not factorize
across any bipartition (test `test_C2_top_cumulant_zero_does_not_imply_factorization`).

**Theorem 2.3 (minimal coalitions are cumulant-connected).** Suppose

    F(H) = 0   and   F(B) ≠ 0 for every ∅ ≠ B ⊊ H.

Then the hypergraph `𝒞(H) = {B ⊆ H : |B| ≥ 2, chi_B ≠ 0}` is connected and spans
`H`.

*Proof.* Suppose not. Then `H = H_1 ⊔ H_2`, with both parts non-empty, such that no
`B ∈ 𝒞(H)` meets both. (Take `H_1` to be one connected component; an isolated
vertex is a component.) Every mixed cumulant is then zero. By Proposition 2.2,
`0 = F(H) = F(H_1) F(H_2)`, which contradicts `F(H_i) ≠ 0`. ∎

The theorem is purely algebraic and holds pointwise in `s`. To read it
dynamically we need one hypothesis.

**Admissibility (A).** On the closed right half-plane `Ω̄`, less the structural
centre (the transverse quotient removes it), the following hold:

- the base and every singleton are stable;
- `s ↦ F_s(S)` is analytic;
- the zeros of `F_s(S)` in `Ω̄` are exactly the transverse eigenvalues of
  portfolio `S` there.

On the benchmark:

- the ratio identity is verified at 160 checks (FC18);
- the local factors are non-zero at the boundaries (0.08–0.31);
- κ ≥ 2 throughout;
- in the lag-network toy, (A) holds exactly.

**Corollary 2.3′.** Under (A), let `H` be a minimal unstable coalition, and let
`z ∈ Ω̄` be any unstable eigenvalue of `H`. This includes the critical `s* = jω*`
at a boundary. Then:

- `F_z(H) = 0`;
- `F_z(B) ≠ 0` for every proper `B`, because proper subsets are stable;
- hence `𝒞_z(H)` is connected and spans `H`.

A minimal coalition is never a disjoint union of non-interacting clusters at its
own characteristic zero.

The corollary does **not** say that `chi_H ≠ 0`. The path and star examples have
κ = 4 with `chi_H ≡ 0`: they are connected through their pair cumulants.

**Corollary 2.4 (boundary identity).** At any zero of `F(H)`,

    chi_H(s*) = − C_H(s*),   C_H := Σ_{π∈Π(H), |π|≥2} Π_{B∈π} chi_B,

so `|chi_H(s*)| = |C_H(s*)|` identically. The modulus of `chi_H` at a boundary is
therefore not an independent diagnostic. The informative question is whether the
composite part nearly cancels on its own.

**Definition 2.5 (connected share and deletion displacement).** At a zero `s*` of
`F(H)`, with `t_π = Π_{B∈π} chi_B(s*)`:

    nu_H = |chi_H(s*)| / Σ_{π∈Π(H)} |t_π| ∈ [0, 1/2].

- The bound `1/2` follows from `|chi_H| = |Σ_{π≠1̂} t_π| ≤ Σ_{π≠1̂} |t_π|`.
- `nu_H = 1/2` when the composite terms are aligned and `chi_H` alone cancels
  them.
- `nu_H → 0` when the composite terms cancel among themselves.
- **For `|H| = 2`, `nu_H = 1/2` is forced.** `F(H) = 1 + chi_H`, so at a zero
  `chi_H = −1`. A pair boundary is connected by necessity and uninformative.

For `B ⊆ H`, the partitions of `H` containing the block `B` sum to
`chi_B F(H∖B)`. Deleting `chi_B` therefore changes `F(H)` by
`−chi_B F(H∖B)`. At a simple zero, the first-order displacement of the
characteristic zero is

    delta_B = |chi_B(s*) F(H∖B)(s*)| / |∂_s F(H)(s*)|    [s^-1].

This puts all clusters on one physical scale (eigenvalue units).

## 3. Cycle expansions (C3)

### 3.1 Scalar ports: the conjectured formula is TRUE for |S| ≥ 2

**Theorem 3.1.** Let `Q ∈ C^{n×n}` have arbitrary diagonal, and `F(S) = det(I + Q_SS)`.
Then:

- for every `S` with `|S| ≥ 2`,

      chi_S = (−1)^{|S|−1} Σ_{γ} Π_{(i,j)∈γ} q_ij,

  where `γ` runs over the `(|S|−1)!` spanning directed cycles of `S` (the cyclic
  permutations of `S`; the two orientations are different cycles when
  `|S| ≥ 3`);
- for singletons, `chi_{i} = 1 + q_ii`. With `Q_ii = 0` this is `1`, not the
  value `0` the cycle formula would give. So the formula holds exactly for
  `|S| ≥ 2` and fails for `|S| = 1`.

*Proof.*

1. By Leibniz, `det(I + Q_SS) = Σ_{σ∈Sym(S)} sgn σ Π_{i∈S} (δ_{iσ(i)} + q_{iσ(i)})`.
2. Write `σ` as a product of disjoint cycles, with `sgn σ = Π_c (−1)^{|c|−1}`.
   The product over `i` factorizes over the cycles:
   - a fixed point `{i}` contributes `1 + q_ii`;
   - a cycle `c` with `|c| ≥ 2` contributes `(−1)^{|c|−1} Π_{i∈c} q_{iσ(i)}`,
     because every factor is off-diagonal.
3. The map from `σ` to (the partition of `S` into cycle supports, a cyclic
   permutation on each support) is a bijection.
4. Hence `det(I + Q_SS) = Σ_{π∈Π(S)} Π_{B∈π} z_B`, where
   - `z_B = (−1)^{|B|−1} Σ_{γ cyclic on B} Π q` for `|B| ≥ 2`;
   - `z_{i} = 1 + q_ii`.
5. By the uniqueness part of C1, `chi_B = z_B`. ∎

**Corollary 3.2 (scalar coincidence with the infinitesimal cumulant of Q).** For
scalar ports and `|S| ≥ 2`,

    chi_S = [θ^S] log det(I + diag(θ) Q),

the multilinear coefficient.

*Proof.* `log det(I + ΘQ) = Σ_k (−1)^{k−1}/k tr((ΘQ)^k)`. The coefficient of
`θ^S` comes only from `k = |S|` and from closed walks that visit each vertex of
`S` exactly once. Each cyclic permutation contributes `k` rotations, which gives
`(−1)^{k−1} Σ_γ Π q`. Theorem 3.1 finishes the proof. ∎

(This is `Q`-based. The existing `cumulants.py` differentiates `log det(I + ΘK)`
with the *unnormalized* `K`, which still carries the local diagonal blocks. See §6.)

**Numerical check.** Relative residual ≤ 1.6e-14 (`Q_ii = 0`) and ≤ 6.2e-14
(`Q_ii ≠ 0`, `|S| ≥ 2`) over 100 random matrices with `n ≤ 6`. This is T4.

### 3.2 Block ports (p = 2): the scalar formula does NOT transfer

Channels are indexed by `ch(S)`, the union over `b ∈ S` of the channels of block
`b`. For `σ ∈ Sym(ch(S))`, join blocks `a` and `b` when some cycle of `σ`
contains channels of both. Call `σ` **block-connected** if this relation connects
all of `S`.

**Theorem 3.3 (block-connected expansion).** For every non-empty `S`,

    chi_S = Σ_{σ ∈ Sym(ch(S)) block-connected} sgn σ Π_{c∈ch(S)} (I + Q)_{c, σ(c)}.

*Proof.*

1. Leibniz gives `F(S) = Σ_σ sgn σ Π_c (I+Q)_{cσ(c)}`.
2. Every `σ` determines the partition `π(σ)` of `S` into the connected
   components of the block–cycle relation.
3. `σ` restricts to a block-connected permutation of `ch(B)` for each
   `B ∈ π(σ)`. Conversely, block-connected permutations chosen on the blocks of
   any `π` assemble into a unique `σ` with `π(σ) = π`.
4. Sign and weight are multiplicative over components with disjoint supports.
5. Hence `F(S) = Σ_π Π_B z_B`, with `z_B` the stated sum. By the uniqueness part
   of C1, `chi_B = z_B`. ∎

Within one block (`Q_bb = 0`), the identity gives weight 1 and the swap gives
`−(I+Q)_{12}(I+Q)_{21} = 0`, so `chi_{b} = 1`.

**Corollary 3.4 (two blocks, closed form).** With `X = Q_ab`, `Y = Q_ba` and
`Q_aa = Q_bb = 0`,

    chi_ab = F(ab) − 1 = det(I − XY) − 1 = −tr(XY) + det X · det Y.

The first equality uses the Schur complement. The last uses the 2×2 identity
`det(I − Z) = 1 − tr Z + det Z`. ∎

**Proposition 3.5 (the naive holonomy formula is refuted; what it is).**

- The "scalar formula with `q_ij` replaced by holonomy traces",
  `(−1)^{|S|−1} Σ_{cyclic orders} tr(Q_{b_1b_2} ⋯ Q_{b_kb_1})`, collects only the
  block-connected permutations made of **one** cycle that visits each block
  exactly once (every other channel is a fixed point of weight 1).
- That sum equals the multilinear coefficient of `log det(I + ΘQ)` with
  `Θ = blkdiag(θ_b I_2)`, by the argument of Corollary 3.2. It is the
  **infinitesimal** cumulant of `Q`.
- The finite-amplitude `chi_S` adds the remainder
  `R_S = chi_S − naive_S`. `R_S` collects block-connected permutations with
  several cycles, or with a cycle that visits a block twice. Already at `|S| = 2`,
  `R_ab = det Q_ab det Q_ba ≠ 0` generically.

Numerically, the relative gap `|chi_S − naive_S| / |chi_S|` is ≥ 0.045 over 30
random block draws. The theorem itself is confirmed by brute-force permutation
sums to ≤ 5.3e-15, and the closed form of Corollary 3.4 to ≤ 7.8e-16.

**Proposition 3.6 (gauge invariance).** Under an admissible port-basis change
`Q → S_b^{-1} Q S_b`, with `S_b` block diagonal:

- every `F(B)` is invariant, hence every `chi_B`;
- the single-cycle traces `tr(Q_{b_1b_2} ⋯ Q_{b_kb_1})` are invariant;
- hence the remainder `R_S` is also invariant.

Individual channel-level terms are not invariant. The numerical gauge drift is
≤ 7.7e-15.

## 4. Three orders that need not coincide (C4)

**Definitions.**

- `d_alg = max{|T| : μ_T ≢ 0}`, with `μ_T = Σ_{R⊆T} (−1)^{|T∖R|} F(R)` the Boolean
  Möbius term. This is the exact multilinear degree of `F` in the replacement
  indicators (`PRINCIPAL_MINOR_PORTFOLIO_STRUCTURE.md` §4).
- `d_conn = max{|T| ≥ 2 : chi_T ≢ 0}`, or 1 if there is no such `T`.
- `κ = min{|S| : S transversely unstable}`.

**Proposition 4.1 (Möbius terms from cumulants).**

    μ_S = Σ_{π∈Π(S)} Π_{B∈π} ψ_B,   with ψ_B = chi_B for |B| ≥ 2 and ψ_{i} = chi_{i} − 1.

In the normalized case (`chi_{i} = 1`), `μ_S` is the sum over partitions of `S`
**without singletons** of `Π chi_B`.

*Proof.*

1. `Σ_{T⊆S} Σ_{π∈Π(T)} Π ψ_B` is a sum over partial partitions of `S`.
2. Completing each partial partition with singletons (weight `ψ_{i} + 1 = chi_{i}`)
   gives `Σ_{π∈Π(S)} Π chi_B = F(S)`.
3. Möbius inversion is unique. ∎

Consequences:

- `μ_S` can be non-zero with `chi_S = 0`. For two independent pairs,
  `μ_1234 = chi_12 chi_34` (T1).
- `chi_S` can be non-zero with `μ_S = 0`, when the non-singleton partition sum
  cancels (example d below).

**Proposition 4.2 (logical independence).** None of the six inequalities between
the three orders holds in general. The counterexamples below all use the
lag-network toy `A_S = −aI − G_SS`, `Q(s) = G/(s+a)` with `a = 1`, `G` zero-diagonal
and scalar ports. For examples a–c, `G` is scaled to 1.02× its boundary value.

| example | κ | d_alg | d_conn | why (exact) |
|---|---|---|---|---|
| (a) directed ring 1→2→3→4→1 | 4 | 4 | 4 | Proper subsets carry no cycle, so `F(S) = 1` and they are stable. The single 4-cycle gives `chi_V = −Π q`. |
| (b) symmetric path 1–2–3–4 | 4 | 4 | **2** | The coupling graph is a tree: no cycle of length ≥ 3, so `chi_B = 0` for `|B| ≥ 3`. `μ_V = chi_12 chi_34 ≠ 0`. The Perron value of the path grows with length (`√2` for P3, golden ratio for P4), so only `V` crosses. |
| (c) symmetric star, centre 1 | 4 | **2** | 2 | A tree with singular adjacency, so `μ_V = det G/(s+a)^4 ≡ 0` and `μ_T = 0` for `|T| = 3`. Perron values `√2` (P3) against `√3` (star). |
| (d) pairs 1–2, 3–4 plus the ring 1→3→2→4→1, ring weight solved from `det G = 0` | 4 | 3 | **4** | `μ_V ≡ 0`, while the only Hamiltonian cycle gives `chi_V ≠ 0`. `d_alg = 3` comes from the directed triangle 1→3→2→1. This is the only tuned example, and the tuning is to a determinant zero. |
| (e) ring (w = 0.6) plus the reciprocal pair 1–3 (w = −1.3) | **2** | 4 | 4 | The pair has eigenvalue `−1 + 1.3 > 0`, while `det G ≠ 0` and the Hamiltonian cycle survives. |

The examples realize every non-implication:

- `d_alg > d_conn` (b) and `d_alg < d_conn` (d);
- `κ > d_conn` (b, c) and `κ > d_alg` (c);
- `κ < d_alg` and `κ < d_conn` (e).

**The one valid structural link is Theorem 2.3:** a minimal coalition is
cumulant-connected. It does **not** bound κ by `d_conn`. In (b), κ = 4 while
`d_conn = 2`.

**Same κ = 4, different connected structure (T3).** At the boundary:

| | `nu_V` | reading |
|---|---|---|
| ring (a) | 1/2 exactly | the zero is carried entirely by the 4-way cumulant; genuinely connected |
| path (b) | 0 exactly | the zero is carried entirely by the pair clusters `chi_12, chi_23, chi_34` and their product `chi_12 chi_34`; composite |
| star (c) | 0 exactly | composite, like the path |

Both kinds of κ = 4 coalition exist.

**Random toys (no tuning).** 500 random sparse 5-unit toys contain 184 minimal
coalitions with `|H| ≥ 2`:

- 0 violate Theorem 2.3;
- 5 are composite (`chi_H = 0` at the zero).

## 5. Boundary sensitivities (Phase 6 formulas)

**Jacobi's formula in adjugate form** is valid at a singular matrix:

    d F(B)/da = tr( adj(I + Q_BB) · dQ_BB/da ),

    d chi_S/da = Σ_π (−1)^{|π|−1}(|π|−1)! Σ_{B∈π} (dF(B)/da) Π_{B'∈π, B'≠B} F(B').

The adjugate is computed by cofactors (no inverse). `connected.cumulant_derivatives`
implements the second formula, which is the product rule applied to the
definition.

**Boundary motion.** At a simple zero `s*` of `F(H)(s, a)`, the implicit function
theorem gives

    ds*/da = − ∂_a F(H) / ∂_s F(H)
           = − ( ∂_a chi_H + ∂_a C_H ) / ∂_s F(H).

`sigma_a = ∂_a chi_H / ∂_a F(H)` is the share of the boundary motion that is
carried by the connected term.

The formula agrees with FC18's eigenvector form `−(p* F_a q)/(p* F_s q)`, for
two reasons:

- at a simple eigenvalue `−1` of `Q_HH`,
  `adj(I + Q_HH) = Π_{j≠crit}(1+λ_j) · q p*/(p* q)`;
- the local factor `Π det(I+M_ii)` is non-zero at `s*`, so it cancels in the
  ratio.

## 6. Relation to the objects already in the repository

| object | definition | relation to `chi` |
|---|---|---|
| Boolean Möbius `μ_T` (`PRINCIPAL_MINOR…` §2) | Möbius transform of `F` on the Boolean lattice | Prop. 4.1: `μ` sums products of cumulants over singleton-free partitions. `μ_V ≠ 0` does not imply `chi_V ≠ 0` (T1). |
| log-Möbius `c_T` (`PRINCIPAL_MINOR…` §3.5) | Möbius transform of `log F(S)` | Also vanishes under factorization, but is **undefined at a boundary** (`log 0`). `chi` is polynomial in `F` and stays finite there. |
| infinitesimal cumulants (`cumulants.py`) | mixed derivatives at `θ = 0` of `log det(I + ΘK)`, with the unnormalized `K` | A different normalization (it keeps the local blocks of `K`) and infinitesimal amplitude. Its `Q`-analogue equals `chi` for scalar ports (Cor. 3.2) and differs by the multi-cycle remainder for 2×2 ports (Prop. 3.5). |
| cycle / holonomy scores (E18) | spectral radius of `M_{a_1a_2} ⋯ M_{a_ka_1}` for one cycle | A single-cycle magnitude. `chi_S` is a *signed sum* over all block-connected permutations, so cancellations between cycles and the multi-cycle remainder are included. A large single cycle can coexist with a small `chi_S`, and conversely. |

## 7. What is and is not claimed

**Claimed (proved):**

- C1: inversion and uniqueness;
- C2 and its converse;
- the connectivity of minimal coalitions (Thm 2.3 / Cor 2.3′, under (A));
- the boundary identity (Cor 2.4);
- the scalar cycle formula for `|S| ≥ 2` and its failure at `|S| = 1` (Thm 3.1);
- the block-connected expansion (Thm 3.3);
- the two-block closed form (Cor 3.4);
- the refutation of the naive holonomy formula for 2×2 ports (Prop 3.5);
- gauge invariance (Prop 3.6);
- the logical independence of `d_alg`, `d_conn` and `κ` (Prop 4.2);
- the boundary sensitivity formulas (§5).

**Not claimed:**

- that `chi_H ≠ 0` is necessary for incompatibility (false: path, star);
- that `chi` "causes" anything. `chi_H` is an irreducible connected
  characteristic interaction, a term in an exact expansion;
- any block-port *cycle* theorem beyond Theorem 3.3 and Corollary 3.4;
- that `chi` is a better diagnostic than existing ones. That is an empirical
  question, tested in CC04.
