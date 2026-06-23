# Line-Shunt Lever Taxonomy And Virtual-Inertia Location

This note records the full derivation behind the compact theorem block added to
the paper.  The algebraic decomposition is standard linear algebra; the intended
contribution is its use as a robustness-lever taxonomy for power-system planning
and the structural location of virtual-inertia actions.

## 1. Network Model

Start from the linearized swing model

```math
M\ddot\theta + D\dot\theta + L\theta = p,
```

where `M=diag(M_i)>0`, `D=diag(D_i)>=0`, and `L` is the power-network stiffness
Laplacian at the operating point.  With `\vartheta=M^{1/2}\theta`,

```math
\ddot\vartheta + \widetilde D\dot\vartheta+\widetilde L\vartheta
=M^{-1/2}p,
\qquad
\widetilde L=M^{-1/2}LM^{-1/2},
\quad
\widetilde D=M^{-1/2}DM^{-1/2}.
```

The planning levers considered here are line changes, new-line additions, and
shunt-like actions such as virtual inertia or synchronous-condenser support.

## 2. Complete Line-Shunt Decomposition

Let `b_ij=e_i-e_j`.  Define

```math
\mathcal L_{\rm all}
=\operatorname{span}\{b_{ij}b_{ij}^T:1\le i<j\le n\},
\qquad
\mathcal D_{\rm sh}
=\{\operatorname{diag}(g):g\in\mathbb R^n\}.
```

The matrices `b_ij b_ij^T` span the symmetric zero-row-sum Laplacian subspace,
which has dimension `n(n-1)/2`.  The diagonal shunt space has dimension `n`.
Their intersection is trivial, because a diagonal matrix whose row sums are zero
must be the zero matrix.  Therefore

```math
\operatorname{Sym}(n)=\mathcal L_{\rm all}\oplus\mathcal D_{\rm sh}.
```

Equivalently, any symmetric perturbation `X` decomposes uniquely into a line
part and a shunt part:

```math
X=X_{\mathcal L}+X_{\mathcal D}.
```

The line part is the unique zero-row-sum matrix with the same off-diagonal
entries as `X`; the residual diagonal is the shunt part.  This gives an explicit
`O(n^2)` decomposition.

If the current network edge set is `\mathcal E`, then

```math
\mathcal L_{\rm all}=\mathcal L_{\rm ex}\oplus\mathcal L_{\rm new},
```

where `\mathcal L_{\rm ex}` uses existing-edge directions and
`\mathcal L_{\rm new}` uses candidate new-line directions.

## 3. Location Of Virtual Inertia

The real derivative of the normalized stiffness with respect to local inertia is

```math
\frac{\partial \widetilde L}{\partial M_i}
=-\frac{1}{2M_i}
\left(E_i\widetilde L+\widetilde L E_i\right),
\qquad
E_i=e_ie_i^T.
```

The negative sign is important dynamically: increasing inertia lowers normalized
modal stiffness.  The structural support statement, however, is independent of
the sign.

Multiplication by a diagonal matrix cannot create a new off-diagonal nonzero
entry.  Hence `E_i\widetilde L+\widetilde L E_i` has off-diagonal support only
on the current network edges.  By the line-shunt decomposition,

```math
\frac{\partial \widetilde L}{\partial M_i}
\in
\mathcal L_{\rm ex}\oplus\mathcal D_{\rm sh},
```

and it has no component in `\mathcal L_{\rm new}`.  In words: local virtual
inertia lives structurally in existing-line directions plus a shunt residual; it
does not require new topology.

This does not say that inertia and line reinforcement are dynamically
equivalent.  They have opposite signs in the conservative modal-stiffness
screen, and the shunt residual can be the part that actually couples to the
critical mode.

## 4. Modal Substitutability

Let

```math
X_i=-\frac{\partial\widetilde L}{\partial M_i}
=X_i^{\rm ex}+\operatorname{diag}(r_i),
\qquad
X_i^{\rm ex}\in\mathcal L_{\rm ex}.
```

For a normalized critical network-family mode `q_c`, define

```math
s_i=
1-\frac{|r_i^Tq_c|}{\|r_i\|+\varepsilon_{\rm reg}}.
```

High `s_i` means the shunt residual is nearly invisible to the critical mode;
existing-line actions are a plausible substitute.  Low `s_i` means the shunt
residual is aligned with the critical mode; inertia or synchronous-condenser
support is harder to replace by line-only actions.

## 5. Claim Discipline

Standard and should be cited, not claimed as new:

- edge Laplacians `b_ij b_ij^T`;
- the dimension count for zero-row-sum symmetric Laplacians;
- the direct-sum decomposition as linear algebra;
- the congruence derivative of `M^{-1/2}LM^{-1/2}`.

Potential contribution:

- using line/shunt/new-line decomposition as a complete robustness-lever
  taxonomy;
- identifying the structural location of virtual inertia as
  `existing lines + shunt`, not `new topology`;
- using modal substitutability `s_i` as a local planning diagnostic.

Scope:

- this is a structural statement for the network-inertia subsystem;
- controller self-energy `\Sigma(s)` can change the dominant pole mechanism, so
  full IBR claims still require validation through the controller-aware NEP.
