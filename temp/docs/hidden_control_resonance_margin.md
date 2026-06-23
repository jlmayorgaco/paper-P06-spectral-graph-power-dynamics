# Hidden Control-Resonance Margin

This note records the theoretical extension added to the manuscript after the
Phase-0D Schur-NEP bridge validation.  The contribution is intentionally framed
as a local diagnostic and repair formulation, not as a system-level collapse
certificate.

## 1. Rational Self-Energy Object

The retained-network nonlinear eigenvalue problem is

```math
T(s)=s^2M+s\Sigma(s)+\widetilde L,
\qquad
\Sigma(s)=D_0+C(sI-A_{cc})^{-1}B.
```

The controller term is frequency dependent.  A static scalar bridge loses this
dependence and therefore cannot see resonance between a network-family pole and
a condensed-control pole.

## 2. Distance To Condensed-Control Poles

For a NEP zero `s_k` and a condensed-control pole
`\mu_j in spec(A_cc)`, define

```math
d_{kj}^{ctrl}
=\frac{|s_k-\mu_j|}{|s_k|+|\mu_j|}.
```

This is dimensionless.  Small distance does not imply instability by itself; it
means a topology or controller perturbation may enter a high-resolvent region.

## 3. Resolvent Bound

Let

```math
\Delta_c(s)=dist(s,spec(A_cc)).
```

If `A_cc` is diagonalizable with eigenvector condition number `kappa_c`, then

```math
\|(sI-A_cc)^{-1}\|\le \frac{\kappa_c}{\Delta_c(s)},
```

and

```math
\|\Sigma'(s)\|
\le
\frac{\|C\|\|B\|\kappa_c^2}{\Delta_c(s)^2}.
```

This explains why a pole can look safe in damping ratio but still be fragile:
when it is close to a condensed-control pole, the sensitivity denominator
contains a large `s Sigma'(s)` contribution.

## 4. Generic NEP Sensitivity

For an admissible action parameter `p`,

```math
\frac{\partial s_k}{\partial p}
=
-\frac{y_k^*T_p(s_k)x_k}{y_k^*T_s(s_k)x_k},
```

with

```math
T_s(s)=2sM+\Sigma(s)+s\Sigma'(s),
\qquad
T_p(s)=s^2M_p+s\Sigma_p(s)+L_p.
```

The damping-ratio derivative is

```math
\frac{\partial\zeta_k}{\partial p}
=-\frac{Re(s_{k,p})}{|s_k|}
\frac{Re(s_k)Re(\overline{s_k}s_{k,p})}{|s_k|^3}.
```

This gives one language for line, shunt, inertia, damping, and controller
actions.

## 5. Hidden Action Margin

For a damping target `zeta_star`,

```math
m_p^\zeta(k)
=
\frac{\zeta_k-\zeta_\star}
{[-\partial\zeta_k/\partial p]_+}.
```

For a control-distance target `d_star`,

```math
m_p^d(k,j)
=
\frac{d_{kj}^{ctrl}-d_\star}
{[-\partial d_{kj}^{ctrl}/\partial p]_+}.
```

The combined hidden margin is

```math
m_p^{hid}(k,j)=min\{m_p^\zeta(k),m_p^d(k,j)\}.
```

This asks how much admissible action is needed before either the damping target
or the control-resonance distance target is violated.

## 6. Resolvent-Aware Weak Nodes And Links

For an action set at node `i`,

```math
W_i^{hid}
=max_{k,j,p in P_i}
\phi_k^\Sigma
[-\partial\zeta_k/\partial p]_+
[-\partial d_{kj}^{ctrl}/\partial p]_+.
```

For a line `e`,

```math
W_e^{hid}
=max_{k,j}
\phi_k^\Sigma
[-\partial\zeta_k/\partial w_e]_+
[-\partial d_{kj}^{ctrl}/\partial w_e]_+.
```

The product structure deliberately requires two things: the action must reduce
damping and it must move the pole toward a condensed-control resonance.

## 7. Planning Formulation

Let

```math
u=[u_ex,u_new,u_sh,u_M,u_D,u_kappa]^T.
```

The minimum-cost repair is

```math
min_u cost(u)
```

subject to

```math
\zeta_k+g_k^Tu \ge \zeta_\star,
\qquad
d_{kj}^{ctrl}+h_{kj}^Tu \ge d_\star.
```

This is a local linearized planner.  It becomes a contribution only after
comparison against full ANDES/NEP recomputation and against existing damping
ratio, QEP-sensitivity, impedance/SCR, and graph-topology baselines.
