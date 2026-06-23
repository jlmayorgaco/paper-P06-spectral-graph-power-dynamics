# Mode-Resolved Substitutability Framework

This note records the added theory from the post-master-package update. It is a
planning interpretation layer, not a replacement for the controller-aware NEP
validation.

## 1. Structural lever space

The real symmetric action space decomposes as

```tex
Sym(n) = L_all \oplus D_sh,
```

where `L_all` is the span of all line Laplacian directions and `D_sh` is the
diagonal shunt/equipment space. The algebra is classical. The contribution is to
use it as a power-system robustness-lever taxonomy:

- existing-line reinforcement,
- candidate new-line topology,
- shunt/STATCOM/synchronous-condenser or virtual-inertia support,
- damping and controller tuning through the NEP self-energy.

## 2. Virtual-inertia support location

For the inertia-normalized Laplacian

```tex
L_tilde = M^{-1/2} L M^{-1/2},
```

a local inertia perturbation satisfies

```tex
d L_tilde / d M_i = -(E_i L_tilde + L_tilde E_i)/(2 M_i).
```

Its off-diagonal support is inherited from the existing network. Therefore the
virtual-inertia-induced matrix direction lives in existing-line directions plus
a diagonal residual:

```tex
d L_tilde / d M_i in L_ex \oplus D_sh.
```

It has no component in the new-line subspace. This does not mean that line
reinforcement and inertia are dynamically equivalent. It says where the action
lives in the matrix-action space.

## 3. Node-by-mode substitutability spectrum

Let

```tex
X_i = - d L_tilde / d M_i = X_i^ex + diag(r_i),
```

where `diag(r_i)` is the shunt residual after projecting onto existing-line
directions. For every retained network-family mode `q_k`,

```tex
s_{i,k} = 1 - |r_i^T q_k| / (||r_i|| ||q_k|| + eps).
```

Large `s_{i,k}` means that existing-line actions are a plausible substitute for
the local equipment/inertia action for that specific mode. Low `s_{i,k}` means
that the residual equipment direction is aligned with the mode, so line-only
repair is structurally incomplete.

The important point is mode dependence. A node can be copper-substitutable for
one mode and equipment-requiring for another. A scalar weak-node score collapses
this matrix and can lose the repair mechanism.

The user-provided Claude audit reports:

- 60 random graphs tested;
- mode-dependent node rankings disagreed in all graphs;
- the distinction persisted among candidate critical modes;
- in about 98% of graphs, the top node for at least one candidate critical mode
  differed from the top node under single-critical-mode scalarization.

These numbers should be treated as a structural audit until reproduced in the
repository and validated over ANDES.

## 4. Area and Y-bus projections

With a Fiedler area partition,

```tex
L_ex = L_intra \oplus L_inter.
```

The mode-dependent inter-area share is

```tex
alpha_{i,k}^inter =
|q_k^T X_i^inter q_k| / (|q_k^T X_i^ex q_k| + eps).
```

This separates within-area copper actions from corridor actions.

The same line-shunt algebra also has a complex admittance version,

```tex
Sym_C(n) = L_all^C \oplus D_sh^C.
```

This is a `Y`-bus taxonomy for series admittance versus shunt compensation. It
is related to, but not identical to, the inertia result because inertia enters
the dynamic normalization rather than the admittance matrix itself.

## 5. Damping-aware hosting capacity

For SG-to-IBR conversion at bus `i`, define a local repair set `A_i` and a
composite hidden margin `m_H`. A damping-aware hosting capacity can be written
as

```tex
H_i^damp = sup { rho in [0,1] :
                 exists u in A_i,
                 m_H(rho, i, u) >= m_star,
                 c^T |u| <= B }.
```

This is not a thermal or voltage hosting metric. It asks how much conversion can
be hosted while preserving damping and control-resonance margins. The possible
novelty is the copper-vs-equipment decision attached to each bus and mode.

## 6. Polytope caution

The feasible action set can be a polytope in lever coordinates, but the true
margin map is generally nonsmooth and nonconvex:

```tex
m_H(u) = min{m_phi(u), m_zeta(u), m_c(u)}.
```

The safe polytope statement is local only. Once the active modes and active
control-pole pairs are fixed, first-order sensitivity gives

```tex
max_{u in P} g_H^T u.
```

This linearized problem has an extreme-point optimum under ordinary LP
conditions. That is a local screening result, not a global optimization theorem.

