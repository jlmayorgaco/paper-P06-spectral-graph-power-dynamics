# Contextual marginals and discrete curvature (CDW hardening, H2)

Status: **PROVED** (elementary). Preregistered in `docs/CDW_HARDENING_PREREG_V1.md`
§1 (commit `05b507e3`).

**Novelty statement.** Every algebraic ingredient below is classical set-function
theory:
- increasing differences and lattice monotone comparative statics: Topkis
  (1978, 1998);
- submodular functions: Lovász (1983), Fujishige (2005), Bach (2013);
- the second-difference characterization of submodularity: textbook material in
  all of the above.

**Nothing in this note is claimed as new mathematics.** What the CDW work adds is
the *use* of these facts:
- as a certificate linking observed contextual sign reversals of inverter
  replacements to the sign of discrete curvature;
- as an explanation of why a context-independent (fixed) ranking of units must
  fail on this class.

## 1. Setting

Let V be a finite ground set: the candidate units whose synchronous generator can
be replaced by a converter. Let f : 2^V → ℝ be a set function. In CDW,
f(S) = α_⊥(S), the transverse spectral abscissa of portfolio S at a fixed policy.
On the V9 census every subset is solvable, so f is defined on the whole lattice.

For S ⊆ V and i, j ∈ V∖S, i ≠ j, define:

- the **marginal** (contextual effect of i in context S):

  $$\Delta_i f(S) = f(S\cup\{i\}) - f(S);$$

- the **one-step second difference** (discrete curvature):

  $$d_{ij}(S) = f(S\cup\{i,j\}) - f(S\cup\{i\}) - f(S\cup\{j\}) + f(S).$$

Two elementary facts follow directly from the definitions:

$$d_{ij}(S) = d_{ji}(S), \qquad d_{ij}(S) = \Delta_i f(S\cup\{j\}) - \Delta_i f(S). \tag{1}$$

**Sign convention in CDW.** Δ_i f(S) < 0 means replacing i in context S moves the
rightmost mode left, which is *stabilizing*. Δ_i f(S) > 0 is *destabilizing*.

**Definitions.**
- f is **submodular** if Δ_i f(S) ≥ Δ_i f(T) for all S ⊆ T ⊆ V and i ∉ T.
- f is **supermodular** if the reverse inequality holds.
- f is **modular** if it is both.

## 2. Chain identity

**Lemma 1 (chain identity).** Let S₀ ⊆ V, and let j₁, …, j_m be distinct elements of
V∖S₀. Put S_r = S_{r−1} ∪ {j_r}, r = 1, …, m. Then for every i ∉ S_m,

$$\Delta_i f(S_m) = \Delta_i f(S_0) + \sum_{r=1}^{m} d_{i j_r}(S_{r-1}). \tag{2}$$

*Proof.* By (1), d_{i j_r}(S_{r−1}) = Δ_i f(S_{r−1} ∪ {j_r}) − Δ_i f(S_{r−1})
= Δ_i f(S_r) − Δ_i f(S_{r−1}). Summing over r = 1, …, m telescopes to
Δ_i f(S_m) − Δ_i f(S₀). ∎

Every ordering of S_m∖S₀ gives a valid chain. Identity (2) therefore holds along
each of the m! maximal chains from S₀ to S_m. Individual terms depend on the
chain; their sum does not.

**Lemma 2 (classical characterization).** f is submodular if and only if
d_ij(S) ≤ 0 for all S and all distinct i, j ∉ S. The same holds for supermodular,
with d_ij(S) ≥ 0.

*Proof.*
- (⇒) Take T = S ∪ {j} in the definition; (1) gives d_ij(S) = Δ_i f(S∪{j}) −
  Δ_i f(S) ≤ 0.
- (⇐) For S ⊆ T and i ∉ T, list T∖S = {j₁, …, j_m} and apply (2): Δ_i f(T) −
  Δ_i f(S) is a sum of nonpositive terms.

The supermodular case is symmetric. ∎

This is the standard second-order characterization (e.g. Fujishige 2005; Bach
2013, Prop. 2.3). It is reproduced here only to make the certificates below
self-contained.

## 3. Nested sign reversal and curvature

A **nested reversal** of i is a pair of contexts S₀ ⊊ S_m with i ∉ S_m and
Δ_i f(S₀), Δ_i f(S_m) of opposite signs. With a materiality threshold τ > 0, it is
**material** if

- one marginal satisfies Δ_i f(·) ≤ −τ and the other satisfies Δ_i f(·) ≥ +τ;
- the direction is **stabilizing → destabilizing** (s→d) when the negative one is
  at the smaller context S₀, and **destabilizing → stabilizing** (d→s) otherwise.

**Proposition A (submodularity forbids s→d).** If f is submodular, then along every
growing chain the marginal of i is nonincreasing:
Δ_i f(S_m) ≤ Δ_i f(S₀). In particular no nested reversal of direction s→d exists,
material or not.

*Proof.* Lemma 2 gives d ≤ 0 everywhere. Lemma 1 then gives
Δ_i f(S_m) − Δ_i f(S₀) ≤ 0. An s→d reversal would need Δ_i f(S₀) < 0 < Δ_i f(S_m),
a strict increase. ∎

**Proposition B (supermodularity forbids d→s).** If f is supermodular, then
Δ_i f(S_m) ≥ Δ_i f(S₀) along every growing chain, and no nested reversal of
direction d→s exists. *Proof.* As for A, with the inequalities reversed. ∎

**Proposition C (both directions ⇒ neither property).**
- If f has a nested reversal of direction s→d, it is not submodular.
- If it has one of direction d→s, it is not supermodular.
- If both directions occur (possibly for different units, contexts or chains), f
  is neither submodular nor supermodular.

*Proof.* These are the contrapositives of A and B. ∎

**Proposition D (curvature witness with a magnitude bound).** Let S₀ ⊊ S_m, i ∉ S_m,
and D := Δ_i f(S_m) − Δ_i f(S₀), with m = |S_m∖S₀|. Then along **every** chain from
S₀ to S_m:
- if D > 0, some r satisfies d_{i j_r}(S_{r−1}) ≥ D/m > 0 (a submodularity
  violation);
- if D < 0, some r satisfies d_{i j_r}(S_{r−1}) ≤ D/m < 0 (a supermodularity
  violation).

For a material nested reversal, abs(D) ≥ 2τ, so the witness satisfies
abs(d) ≥ 2τ/m.

*Proof.* By (2), D is the sum of m curvature terms, so the largest term is at
least their mean D/m. The case D < 0 is symmetric. ∎

Each witness (S_{r−1}, i, j_r) is a one-step square: four portfolios differing by
at most two units. It is therefore a *minimal* counterexample of the form used in
the old E1b table. A nested reversal is certified by an explicit local
interaction: some unit j_r, added on the way from S₀ to S_m, changes the effect of
i in the direction of the reversal, by at least 2τ/m.

## 4. Why nestedness matters: arbitrary-context reversals are weaker

**Remark 1 (an arbitrary-context reversal is compatible with submodularity).** Let
V = {i, a, b} and f(∅) = f({a}) = f({b}) = f({a,b}) = 0, f({i}) = 2,
f({a,i}) = −1, f({b,i}) = 1, f({a,b,i}) = −2.

The six one-step second differences are:

$$d_{ab}(\emptyset)=0,\ d_{ai}(\emptyset)=-3,\ d_{bi}(\emptyset)=-1,\ d_{ab}(\{i\})=0,\ d_{ai}(\{b\})=-3,\ d_{bi}(\{a\})=-1.$$

All are ≤ 0, so f is submodular by Lemma 2. Yet Δ_i f({a}) = −1 < 0 < 1 =
Δ_i f({b}). An arbitrary-context sign reversal (the old CDW A1 statistic) is
therefore **not** by itself evidence against submodularity. A *nested* reversal
is, by Proposition A.

**Lemma 3 (meet lemma).** Let i ∉ S₁ ∪ S₂, M = S₁ ∩ S₂, and suppose
Δ_i f(S₁) < Δ_i f(S₂). Then at least one of the following holds:
- Δ_i f(M) < Δ_i f(S₂), a strict increase along M ⊆ S₂, so f is not submodular;
- Δ_i f(M) > Δ_i f(S₁), a strict decrease along M ⊆ S₁, so f is not
  supermodular.

The same statement holds with the join J = S₁ ∪ S₂ in place of M, whenever
i ∉ J, with the inclusions reversed.

*Proof.* If neither held, then Δ_i f(M) ≥ Δ_i f(S₂) > Δ_i f(S₁) ≥ Δ_i f(M), a
contradiction. The join case is identical. ∎

So an arbitrary reversal proves that f is not **modular**, but not which of the
two properties fails. A nested reversal identifies the failing property from its
direction (Proposition C). Neither lemma says the nested pair through M or J is
*material* or has *stable* contexts. That is an empirical question, and the
hardening census answers it (H3).

## 5. Curvature heterogeneity is what defeats fixed rankings

A fixed ranking orders units once, independent of context. Its pairwise decision
between i and j in context S is correct when the sign of
Δ_i f(S) − Δ_j f(S) agrees with the ranking.

**Proposition E (preference drift equals curvature heterogeneity).** For a chain
S₀ ⊂ … ⊂ S_m with i, j ∉ S_m, i ≠ j,

$$[\Delta_i f-\Delta_j f](S_m) - [\Delta_i f-\Delta_j f](S_0) = \sum_{r=1}^{m}\big(d_{i j_r}(S_{r-1}) - d_{j j_r}(S_{r-1})\big). \tag{3}$$

*Proof.* Apply Lemma 1 to i and to j along the same chain and subtract. ∎

**Corollary 1.** Suppose some pair (i, j) has a **material preference reversal**
along a chain: Δ_i f(S₀) ≤ Δ_j f(S₀) − τ and Δ_j f(S_m) ≤ Δ_i f(S_m) − τ. Then:
1. every fixed ranking of V errs, by at least τ, on at least one of the two
   contexts;
2. some r satisfies d_{i j_r}(S_{r−1}) − d_{j j_r}(S_{r−1}) ≥ 2τ/m.

*Proof.*
1. A fixed ranking places i before j or j before i. Either way it contradicts one
   of the two material inequalities.
2. The left side of (3) is ≥ 2τ; take the largest term. ∎

**Remark 2 (sub/supermodularity is neither necessary nor sufficient for
fixed-ranking adequacy).** Let f(S) = g(|S|) + Σ_{k∈S} w_k, with g arbitrary. Then
d_ik(S) = d_jk(S) for all i, j, k, so by (3) the preference between any two
units never drifts. The ranking by w is exact in every context, whether g is
concave (submodular), convex (supermodular) or neither. Conversely, a submodular
f can have material preference reversals. In the f of Remark 1, compare i and b
along ∅ ⊂ {a}:
- in context ∅, Δ_b f(∅) = 0 ≤ Δ_i f(∅) − 1 = 1, so b is preferred;
- in context {a}, Δ_i f({a}) = −1 ≤ Δ_b f({a}) − 1 = −1, so i is preferred.

With τ = 1 this is a material reversal, and (3) attributes it to
d_{ia}(∅) − d_{ba}(∅) = −3 − 0 = −3.

What obstructs a fixed ranking is the **heterogeneity** of curvature across units,
d_{ik} ≠ d_{jk}, not the sign of curvature. In CDW the two are measured
separately:
- the sign structure by nested reversals and E1b (H3);
- the ranking obstruction by p* and regret (H4).

## 6. Scope of the certificates in the spectral-abscissa setting

- **Algebraic scope.** The identities hold for any real set function, so they hold
  exactly for α_⊥ on the census lattice, whatever the dynamics.
  - Intermediate portfolios on a chain may be unstable, or may have their rightmost
    mode in a fast converter-control band.
  - The **terms** d_{i j_r}(S_{r−1}) can therefore be driven by a mode switch
    (α_⊥ is a maximum over modes, and a maximum of smooth branches can carry
    curvature of either sign where the argmax changes).
  - The hardening therefore records whether a witness exists whose four
    portfolios are all stable with an electromechanical rightmost mode (the
    "EM-clean witness", prereg H3).
- **Level-C scope.** Level-C eligibility restricts only the two **endpoint**
  marginals: both contexts are stable, and the same tracked EM mode is rightmost
  before and after the replacement. It does not make the chain's intermediate
  squares electromechanical.
- **Direction and curvature sign.** Along a growing chain, s→d means Δ_i
  increased. By Proposition D this requires positive curvature somewhere: a
  replacement that becomes *more* destabilizing as other units are replaced.
  d→s requires negative curvature somewhere.

## 7. What is and is not claimed

**Claimed.**
- Identities (2)–(3) and Propositions A–E: elementary and classical in substance.
- **Empirical claim (H3/H4).** On the tested inverter-rich IEEE-39 class, nested,
  material, electromechanical same-mode reversals of both directions occur. By
  Proposition C this certifies mixed-sign discrete curvature of α_⊥ on nested
  chains. By Corollary 1, preference reversals certify that no fixed ranking can
  be exact.

**Not claimed.**
- Any new identity.
- Any statement about greedy algorithms beyond these exact certificates.
- Any causal mechanism for the sign of a particular d_ij. The spectral abscissa
  is a maximum over modes, and the note does not attribute curvature to a
  physical interaction path.

## References (classical ingredients)

- D. M. Topkis, "Minimizing a submodular function on a lattice," *Operations
  Research*, 26(2):305–321, 1978.
- D. M. Topkis, *Supermodularity and Complementarity*. Princeton Univ. Press,
  1998.
- L. Lovász, "Submodular functions and convexity," in *Mathematical Programming:
  The State of the Art*, Springer, 1983, pp. 235–257.
- S. Fujishige, *Submodular Functions and Optimization*, 2nd ed., Annals of
  Discrete Mathematics 58, Elsevier, 2005.
- F. Bach, "Learning with submodular functions: a convex optimization
  perspective," *Foundations and Trends in Machine Learning*, 6(2–3):145–373,
  2013.

The DOIs are verified in `results/CDW_LITERATURE_GAP_MATRIX.csv`.
