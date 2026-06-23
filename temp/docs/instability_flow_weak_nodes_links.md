# Instability-Flow Weak Nodes and Weak Links

This note records the mathematical object added to the paper: weak nodes and
weak links as local instability-flow diagnostics derived from the controller-aware
rational NEP.

## 1. NEP pole sensitivity

The retained-network nonlinear eigenvalue problem is

```tex
T(s)x = 0,
T(s) = s^2 M + s Sigma(s) + L_tilde,
Sigma(s) = D_0 + C(sI - A_cc)^{-1}B.
```

For a simple zero `s_k` with right/left null vectors `x_k`, `y_k`,

```tex
ds_k/dp = - y_k^* T_p(s_k)x_k / (y_k^* T_s(s_k)x_k),
```

where

```tex
T_s(s) = 2sM + Sigma(s) + s Sigma'(s),
Sigma'(s) = -C(sI-A_cc)^{-2}B.
```

This is the exact local sensitivity of the pole for any scalar action `p`.

## 2. Damping-ratio sensitivity

With

```tex
zeta_k = -Re(s_k)/|s_k|,
```

and `s_{k,p}=ds_k/dp`,

```tex
d zeta_k/dp =
 - Re(s_{k,p})/|s_k|
 + Re(s_k) Re(conj(s_k)s_{k,p})/|s_k|^3.
```

An action is locally damaging for mode `k` if this derivative is negative.

## 3. Control-pole resonance expansion

If `A_cc` is diagonalizable with right/left eigenvectors `v_p`, `w_p` and
control pole `mu_p`, then

```tex
(sI-A_cc)^{-1} = sum_p v_p w_p^* / (s-mu_p),
Sigma(s) = D_0 + sum_p C v_p w_p^* B / (s-mu_p).
```

The sensitivity-relevant network/control resonance score is

```tex
R_{kp}^ctrl =
 |y_k^* C v_p| |w_p^* B x_k| / (|s_k-mu_p|^2 + eps).
```

The squared distance comes from `Sigma'(s)`, which contains `(s-mu_p)^(-2)`.

## 4. Line participation and modal propagation

For line `e=(i,j)`,

```tex
b_e = e_i - e_j,
G_e = b_e b_e^T.
```

In a symmetric graph-modal approximation,

```tex
q_k^T G_e q_k = (b_e^T q_k)^2,
chi_{e,k}=|b_e^T q_k|^2.
```

For two modes,

```tex
q_k^T G_e q_l = (b_e^T q_k)(b_e^T q_l).
```

Thus line perturbations can couple modes; they do not only affect the critical
mode directly.

For the full NEP, a normalized line participation is

```tex
chi_{e,k}^NEP =
 |y_k^* G_e x_k| / (||y_k|| ||G_e|| ||x_k|| + eps).
```

## 5. Instability current

The line-level instability current is

```tex
I_e^inst =
 sum_{k,p} chi_{e,k}^NEP R_{kp}^ctrl [- d zeta_k/dw_e]_+.
```

It is large when:

1. the line moves the retained mode,
2. the retained mode is exposed to a condensed-control pole,
3. reinforcing or perturbing that line lowers damping.

This score is a local first-order margin-loss diagnostic, not a global topology
theorem.

## 6. Node vulnerability and repairability

Node resonance exposure can be screened by

```tex
I_i^res = sum_{k,p} |q_{k,i}|^2 R_{kp}^ctrl.
```

For planning, define a local action set `A_i`. Vulnerability is

```tex
V_i = max_{a in A_i} [- d m_H/da]_+.
```

Repairability is

```tex
R_i^repair =
 max_{a in A_i} [d m_H/da]_+ / (cost(a)+eps).
```

A vulnerable node is not necessarily the best repair node.

## 7. Local certificate

If `m_H(a)` is differentiable and `dm_H/da < 0`, Taylor expansion gives

```tex
m_H(a+eta) = m_H(a) + eta dm_H/da + o(eta).
```

For sufficiently small positive `eta`, the margin decreases. Therefore the
negative derivative is a valid local weak-element certificate.

What remains empirical is whether the instability-current ranking outperforms
flow loading, betweenness, Fiedler participation, SCR/impedance, and classical
damping participation on full ANDES benchmarks.

