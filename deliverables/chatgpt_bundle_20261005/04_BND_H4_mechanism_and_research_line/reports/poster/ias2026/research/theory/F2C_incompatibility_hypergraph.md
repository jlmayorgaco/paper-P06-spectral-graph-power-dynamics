# F2C — The minimal spectral incompatibility hypergraph

Companion to `theory/spectral_composability_order.md` (F2) and
`theory/F2B_composability_region_theorem.md` (F2B). Implemented in
`src/ibr_cycles/diagnosis/composability.py` (`incompatibility_hypergraph`,
`upward_closure_defect`, `hypergraph_move`, `maximal_free_sets`); the
combinatorial statements and the construction of §4 are
`tests/test_incompatibility_hypergraph.py`; Safeguard A is
`tests/test_reactive_policy_coordinate.py`. Numerical use: `experiments/F7_policy_hypergraph.py`,
reported in `docs/F7_POLICY_HYPERGRAPH.md`.

`kappa_Gamma` says *how many* replacements a planner may make before the
protected region is disturbed. It does not say *which*. Two parameter points
with the same `kappa` can be dangerous for different coalitions, and the
planner acts on the coalition. This note keeps the whole object.

## 1. Setting

As in F2: a finite candidate set `A` (vertices), a protected region `Gamma` with
piecewise-C1 boundary, and for each replacement set `S` the count
`N_Gamma(S; theta)` of eigenvalues of `A_red(S; theta)` in `Gamma`.

Two standing requirements.

- **Baseline-stable region.** `N_Gamma(empty; theta) = 0`. If the baseline
  itself disturbs `Gamma`, the literal definition below returns `H = {empty}`
  and says nothing about portfolios. The code accepts a nonzero baseline by
  using `Delta N != 0`, which reduces to `N > 0` here and matches F2 eq. (4).
- **Downward-closed tested family.** `Script_S` contains every subset of every
  tested set; in this project `Script_S` is the power set of `A`. Without it,
  "every proper subset" would range over untested sets, and a set could be
  called minimal only because its destabilizing subset was never evaluated.
  `incompatibility_hypergraph` refuses a family that is not downward closed.

## 2. Definitions

    U_Gamma(theta) = { S in Script_S : N_Gamma(S; theta) > 0 }       unsafe family
    C_Gamma(theta) = Script_S \ U_Gamma(theta)                         safe family

    H_Gamma(theta) = { S : N_Gamma(S; theta) > 0  and
                           N_Gamma(R; theta) = 0 for every proper subset R of S }
                   = min U_Gamma(theta)                               (minimal elements)

Vertices are candidate replacements. Hyperedges are **minimal destabilizing
portfolios**: sets that disturb `Gamma` while none of their proper subsets does.

## 3. Elementary properties

**Proposition 1 (antichain).** `H_Gamma(theta)` is an antichain under inclusion,
i.e. a clutter (Sperner family) on `A`.

*Proof.* Let `E, F` be in `H` with `E` a proper subset of `F`. Then `F` has a
proper subset, `E`, with `N_Gamma(E) > 0`, contradicting the minimality clause
in the definition of `F`. QED

**Proposition 2 (kappa is carried by H).**

    kappa_Gamma(theta) = min over S in H_Gamma(theta) of |S|,

with `min` of the empty family equal to `infinity`; `H_Gamma = empty` iff
`kappa_Gamma = infinity`.

*Proof.* `H` is contained in `U`, so `min_H |S| >= min_U |S|`. Conversely, `U` is
finite, so every `S` in `U` contains a minimal element `E` of `U` (descend through
unsafe proper subsets until none is left), and `|E| <= |S|`; hence
`min_H |S| <= min_U |S|`. With a zero baseline, `U = { S : Delta N(S) != 0 }`,
whose minimum size is F2 definition (4). QED

**Proposition 3 (what H does and does not determine).** Let
`up(H) = { S : S contains some E in H }`. Then

    U_Gamma(theta)  is a subset of  up(H_Gamma(theta)),

and equality holds **iff** `U_Gamma` is upward closed, **iff** `C_Gamma` is
downward closed, i.e. an abstract simplicial complex. In that case `H` is exactly
the family of minimal non-faces of `C_Gamma` and determines `C_Gamma` completely.

*Proof.* The inclusion is the descent argument of Proposition 2. If `U` is upward
closed, every superset of a hyperedge is unsafe, so `up(H)` is contained in `U`.
If `U = up(H)` then `U` is upward closed because `up` of anything is. `U` upward
closed is the same statement as its complement being downward closed. QED

**The stable-set family is not called a simplicial complex here.** Dynamic
stability is not monotone under adding replacements: a portfolio can contain a
destabilizing coalition and be safe because a further replacement re-stabilizes
it. The **closure defect**

    D_Gamma(theta) = up(H_Gamma(theta)) \ U_Gamma(theta)

(safe sets containing a hyperedge) measures this. `C_Gamma` is a simplicial
complex at `theta` exactly when `D_Gamma(theta)` is empty, and that is checked at
every F7 point (`upward_closure_defect`), never assumed. Without heredity `H`
does not determine `U`; `D` is additional information.

**Proposition 4 (H-free portfolios are safe, without heredity).** If `S` contains
no hyperedge of `H_Gamma(theta)`, then `N_Gamma(S; theta) = 0`.

*Proof.* If `S` were unsafe it would contain a minimal unsafe set, which is a
hyperedge (Proposition 2). QED

This is the planning content of `H`. A set `T` meeting every hyperedge (a
transversal) is a **retention set**: keeping the machines in `T` makes every
portfolio drawn from `A \ T` safe. The maximal `H`-free portfolios are the
complements of the minimal transversals (`maximal_free_sets`). The converse —
that a portfolio containing a hyperedge is unsafe — holds only where `D` is
empty.

**Proposition 5 (restriction to a sub-fleet is exact).** For `B` a subset of `A`,
the hypergraph computed on the power set of `B` is

    H_Gamma^B(theta) = { E in H_Gamma^A(theta) : E subset of B }.

*Proof.* Whether `E` is minimal unsafe depends only on the counts of `E` and of
its subsets, all of which lie in the power set of `B` when `E` does. QED

So the core-restricted hypergraph that F7 computes is exactly the trace of the
full-fleet hypergraph on the core. `kappa` restricts only as an inequality,
`kappa^B >= kappa^A`.

## 4. Local constancy of the entire hypergraph

**Theorem 6.** Assume (H0)–(H3) of F2B for **every** `S` in `Script_S` on
`Theta_reg`. Then `theta -> H_Gamma(theta)` is locally constant on `Theta_reg`,
hence constant on each connected component.

*Proof.* By F2B Theorem 1 each `N_Gamma(S; ·)` is locally constant on `Theta_reg`.
`Script_S` is finite, so at each `theta_0` there is one neighbourhood on which all
of them are constant. `U_Gamma(theta)` is determined by the finite vector of
indicators `1[N_Gamma(S; theta) > 0]`, and `H_Gamma = min U_Gamma` is a function of
`U_Gamma`, so `H_Gamma` is constant on that neighbourhood. Constancy on a connected
component follows because each level set of `H_Gamma` is open (local constancy)
and its complement is a union of other open level sets. QED

The hypothesis is (H3) for every `S`, not only for the witnesses. That is what
F2B's `Theta_reg` already requires. It is slightly more than needed: `H` depends
only on the indicators, so a crossing that takes some `N(S)` from 1 to 2 is a
spectral boundary that leaves `H` unchanged.

**Corollary 7 (the H-partition refines the kappa-partition).**

    P_k = union over { H : min_{E in H} |E| = k } of P_{k,H},
    P_{k,H} = { theta in Theta_reg : H_Gamma(theta) = H }.

Each `P_{k,H}` is a union of connected components of `Theta_reg`. Since `kappa` is
a function of `H` (Proposition 2) but not conversely, **two regions can have the
same `kappa` and different hypergraphs**: `H` separates regions that `kappa`
merges.

**Proposition 8 (the refinement is strict; an explicit construction).** Three
actions, one oscillatory mode at 0.6 Hz with real part

    sigma(S; theta) = -1 + 0.6 [1 in S] + (0.6 - 0.4 theta) [2 in S] + (0.3 + 0.4 theta) [3 in S].

For `theta` in `[0, 1]` every single and the pair `{2,3}` stay out of `Gamma`, and

| `theta` | `H_Gamma` | `kappa_Gamma` |
|---|---|---|
| `[0, 1/4)` | `{ {1,2} }` | 2 |
| `(1/4, 1/2)` | `{ {1,2}, {1,3} }` | 2 |
| `(1/2, 1]` | `{ {1,3} }` | 2 |

`kappa` is constant and `H` changes twice, at `theta = 1/4` (`{1,3}` reaches the
imaginary axis) and `theta = 1/2` (`{1,2}` leaves it). The system satisfies
(H0)–(H3) away from those two points. `tests/test_incompatibility_hypergraph.py::test_same_kappa_different_hypergraph_on_a_continuous_path`.
It is a construction, not evidence about power systems; F7 asks whether the
same thing happens on IEEE-39.

## 5. Where the hypergraph can change

**Corollary 9 (boundary characterisation).** If `theta_1`, `theta_2` lie in
`Theta_reg` and `H_Gamma(theta_1) != H_Gamma(theta_2)`, then every continuous path
between them leaves `Theta_reg`. If (H0)–(H2) hold along the path, the failure is
(H3): at some point of the path some `S` in `Script_S` has an eigenvalue on
`dGamma`.

*Proof.* Otherwise the path lies in one connected component of `Theta_reg`, on
which `H` is constant by Theorem 6. QED

**Proposition 10 (elementary moves).** Suppose that at `theta*` on a path exactly
one set `W` changes its indicator — the generic case of F2B Proposition 3, a
simple crossing by a single subset. Write `H-`, `U-` for the side before.

(a) `W` becomes unsafe. `H` changes iff every proper subset of `W` is safe, and then

    H+ = { E in H- : E does not contain W }  union  { W }.

(b) `W` becomes safe. `H` changes iff `W` is in `H-`, and then

    H+ = (H- \ {W})  union  { T in U- : T strictly contains W and
                              the only proper subset of T in U- is W }.

In both cases the witness `W` lies in `H- symmetric-difference H+`: every
hypergraph transition is carried by an identifiable hyperedge that enters or
leaves.

*Proof.* `H+ = min(U+)` with `U+ = U- union {W}` in (a) and `U- \ {W}` in (b).
(a) If some proper subset of `W` is unsafe, `W` is not minimal and no other set's
minimality changes, because every set containing `W` already contains that
unsafe subset. Otherwise `W` is minimal, the members of `H-` containing `W` stop
being minimal, and the others are untouched because `W` is not a subset of them.
(b) If `W` is not in `H-`, it has an unsafe proper subset `R`, which still blocks
every set that `W` blocked. If `W` is in `H-`, a set `T` becomes minimal exactly
when the unsafe proper subsets of `T` were `{W}` alone. QED

Naming used in F7 (`hypergraph_move`):

| move | meaning |
|---|---|
| `CONTRACTION` | a hyperedge is replaced by a proper subset of itself: (a) with `W` inside a lost hyperedge, e.g. `{30,33,35,37} -> {30,33,35}` |
| `EXPANSION` | the reverse |
| `INSERTION` | a new hyperedge incomparable with all others |
| `DELETION` | a hyperedge leaves and nothing replaces it |
| `REORGANISATION` | anything else, e.g. (b) with new supersets appearing |

A spectral boundary of a set that is not in `H- symmetric-difference H+` — for
instance the flagship crossing while one of its triples is already unsafe — is a
genuine eigenvalue crossing that moves neither `H` nor `kappa`. F7 records those
events separately.

## 6. Which part of dGamma is crossed — Safeguard B

For the band region of this project

    Gamma = { s : Re s > 0,  omega_lo <= Im s <= omega_hi },  omega_lo = 2 pi 0.3,  omega_hi = 2 pi 1.5

(one representative per conjugate pair), the boundary has three segments:

| segment | set | crossing means |
|---|---|---|
| `IMAGINARY_AXIS` | `Re s = 0`, `omega_lo <= Im s <= omega_hi` | an eigenvalue changes stability inside the band; the direct right-half-plane count changes |
| `LOWER_FREQUENCY_EDGE` | `Im s = omega_lo`, `Re s > 0` | an **already unstable** eigenvalue enters or leaves the band; the direct right-half-plane count does not change |
| `UPPER_FREQUENCY_EDGE` | `Im s = omega_hi`, `Re s > 0` | as above, at the top of the band |

plus `DAE/EQUILIBRIUM FAILURE` where (H1) or (H2) fails, and corners where two
segments are reached at once (codimension two, not generic). The magnitude floor
`|s| > 1e-3` is never reached inside the band because `|s| >= omega_lo` there.

**Only an imaginary-axis crossing may be called a dynamic-stability
composability transition.** A band-edge crossing is a **band-classification
transition**: nothing became stable or unstable, a mode was re-labelled as inside
or outside the protected window.

Classification rule (`_f7_common.locate_event`). Bisect the indicator of the
witness to a bracket, match each eigenvalue inside `Gamma` on the unsafe side to
the nearest eigenvalue on the safe side, and report which constraint of `Gamma`
the safe-side partner violates. Every event also records the direct
right-half-plane count on both sides, which gives an independent check: it
changes by one at an imaginary-axis crossing and by zero at a band edge.

## 7. Safeguard A — a regular reactive-policy path

F2 §5 asked for the reactive policy to be parameterised "inside a fixed
converter structure, with a gain going to zero". That is not enough.

**The naive coordinate is singular at zero.** Scale the gains of the plant
voltage regulator `q_cmd = kp_v e + x_v`, `x_v' = ki_v e` by `g`. At `g = 0`,
`x_v' = 0`: every converter carries a **decoupled marginal integrator**. The
reduced matrix gains one exact zero eigenvalue per converter and the equilibrium
is a continuum (`x_v` arbitrary), so (H1) holds only by choosing a branch and
the zero-mode count changes at `g = 0`. Measured on the flagship: **four exact
zero eigenvalues** at `g = 0`; at `g = 1e-4` the same four states sit at
`|lambda| = 4.7e-4 ... ` (real, `-4.7e-4`, `-5.8e-4`, ...) and at `g = 1e-3` at
`|lambda| = 0.0047 ... 0.0144`, collapsing onto the origin as `g -> 0`
(`results/F7/F7_safeguard_A_audit.json`). So `g = 0` is **excluded** from
`Theta_reg` for that coordinate.

**The coordinate F7 uses.** A leaky regulator with the same gains:

    q_cmd = q_ref + g kp_v e + x_v,        x_v' = g ki_v e - w (x_v - q_ref),     w = 0.05 rad/s.

Properties, each proved from the equations and each checked in
`tests/test_reactive_policy_coordinate.py`:

| check | status | reason |
|---|---|---|
| state dimension constant in `g`? | **yes** | `x_v` exists for every `g`; 86 states for the flagship at every `g` |
| equilibrium independent of `g`? | **yes** | at the matched operating point `e = 0` and `x_v = q_ref` solve the equations for every `g` |
| `g_z` regular? | **yes, and independent of `g`** | the network equations contain no regulator parameter; `cond(g_z) = 2.79e3` for the flagship at every `g` |
| artificial zero or marginal state absent? | **yes** | at `g = 0` the matrix is block triangular: spectrum = matched fixed-Q spectrum plus `-w` per converter; no eigenvalue below `1e-6` |
| `A_red` continuous? | **yes, affine** | `f` is affine in `g` at fixed `(x, z)` and the equilibrium does not move, so `A_red(g) = A_red(0) + g Delta` exactly (checked to `1e-7`, the central-difference level) |
| baseline count unchanged for structural reasons? | **yes** | the base case contains no converter, so `A_red(empty)` does not depend on `g` at all |

`g = 0` is therefore inside `Theta_reg` and **is** exact fixed-Q: its spectrum is
the matched-dispatch spectrum plus four copies of `-0.05`, which lie at frequency
zero, far from `dGamma`. So the `g = 0` edge of every F7 map reproduces the
matched policy exactly (to `2.4e-6`, central-difference level), and the F7A
`g = 0` edge must reproduce F2. At `g = 1` the regulator is the E30 voltage
regulator up to the leak (1.4 percent of the integral path at 0.57 Hz): the
flagship's worst band real part is `-0.0956` against `-0.0953` for the plain PI,
and the largest oscillatory-mode shift the leak causes is `0.064`, on a mode at
`-7.7 +- 2.5j` (damping ratio 0.95), far from `dGamma`.

Audit numbers (`experiments/F7_safeguard_A_audit.py`): `g = 0` spectrum against
matched-plus-leak-poles `2.4e-6`; equilibrium shift over `g` in `[0, 1]` at most
`8.1e-8` (solver tolerance); relative spread of `cond(g_z)` `4.5e-9`; affine
residual of `A_red` `5.4e-11`; base matrix bit-identical across `g`.

The excitation coordinates are regular for the same structural reason: the
reference `vref = V + efd/K` absorbs the gain, the time constant does not enter
the equilibrium, and neither appears in the network equations. Along **every**
F7 slice (H0)–(H2) therefore hold identically, and a `DAE/EQUILIBRIUM FAILURE`
boundary is structurally impossible. Every recorded transition must be an
imaginary-axis, band-edge or corner event; F7 checks this rather than trusting it.

## 8. Exceptions

| exception | consequence |
|---|---|
| baseline disturbs `Gamma` | `H = {empty}`; points excluded from every map |
| tested family not downward closed | minimality undefined; refused in code |
| closure defect `D` nonempty | `C_Gamma` is not a simplicial complex at that point; `H` does not determine `U` |
| multiple simultaneous flips | several hyperedges change at once; Proposition 10 does not apply and the event is flagged `simultaneous` |
| thin features | a region narrower than the finest refined lattice spacing, or a tangential contact of `dGamma`, can be missed by any finite lattice; points with a subset within `2e-3` of the axis are counted |
| finite tested family | `H = empty` means "no tested coalition", never "no coalition" |

## 9. What is and is not new

Minimal failing sets are old. In reliability theory they are minimal cut sets;
in power-system security analysis, minimal `N-k` contingency sets; in
combinatorics, the circuits of an independence system, a clutter and its blocker.
Proposition 4 is the standard duality between a clutter and its transversals.
Nothing in §3–§5 is claimed as mathematics.

What is specific here is the pairing of (i) a failure criterion defined by a
**spectral count in a protected region**, which is branch- and label-independent,
with (ii) the **local constancy of the whole clutter in continuous physical
parameters** and (iii) the classification of every change by the boundary segment
crossed. Together they make "the policy-dependent incompatibility hypergraph" a
well-posed map `theta -> H_Gamma(theta)`, so that a phase diagram labelled by
hyperedges is a mathematical object. Whether that map does anything interesting
on a real benchmark is an empirical question, and F7 is where it is asked.
