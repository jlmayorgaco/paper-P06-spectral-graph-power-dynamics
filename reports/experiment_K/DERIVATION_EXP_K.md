# Experiment K: fixed-support analytic design

## Physical architecture and cost

Let `epsilon_i=1-rho_i`. The support `S={i:epsilon_i>0}` is a discrete
architecture choice. There are 1024 supports for buses 30–39. At `epsilon_i=0`
the SG differential block is removed. At `epsilon_i=1` the GFL block is removed.
Neither endpoint is differentiated through. The initialized SG dispatch is
`P_i`, and the retained power is `J=P'epsilon`.

At a fixed interior architecture, each present SG contributes its full
differential equations and `epsilon_i` times its terminal-current Jacobian.
Each present GFL contributes its full equations and `(1-epsilon_i)` times its
terminal-current Jacobian. The exact algebraic reduction is

```
A_red = A - B G_y^{-1} C,
G_y = Y_static + D.
```

`Y_static` includes the initialized ZIP admittances reconstructed from the
frozen operating point, including colocated loads at buses 31 and 39. The
common-angle gauge vector `g` satisfies `A_red g≈0`; an orthonormal basis `Q`
for `g^perp` yields `A_q=Q' A_red Q`. Every spectral margin and reported
abscissa uses **all** eigenvalues of `A_q`.

## Derivatives inside one architecture

For a simple eigenvalue of `A_q`, with `A_q r=lambda r` and
`l^H A_q=lambda l^H`,

```
d lambda/dz = (l^H A_{q,z} r)/(l^H r).
```

For a share coordinate in an interior SG+GFL bus, `A` and `B` do not change,
while `C` and `D` do. Hence

```
A_{red,epsilon_i}
  = B G_y^{-1}(D_{epsilon_i} G_y^{-1} C - C_{epsilon_i}).
```

For a PLL gain, `C` and `D` are fixed while the GFL state and voltage-input
Jacobians change:

```
A_{red,K} = A_K - B_K G_y^{-1} C.
```

The implemented `gain_derivatives` are exact derivatives of the PLL frequency
and integrator rows. The gauge direction is constant within a fixed
architecture. At a mode collision, a single-mode derivative is not used as a
gradient of the spectral abscissa. For active poles
`A={j:Re(lambda_j)>=alpha-1e-6}`, its directional derivative is

```
D alpha[dZ] = max_{j in A} Re(d lambda_j[dZ]).
```

The derivative ledger compares these formulas against full quotient-spectrum
central differences. Its maximum absolute complex error at the selected
steps is `1.46e-7 s^-1` per unit parameter; this does **not** attain the
requested `1e-8` derivative target.

## One-sided SG re-entry

At `epsilon_i->0+`, the SG state block exists, but its network-current
feedback vanishes. Its limiting poles are those of the isolated SG state
matrix `A_SG,i`. For each branch the table reports a two-point estimate of

```
lambda_i,l(epsilon) = mu_i,l + b_i,l epsilon + O(epsilon^2).
```

The coefficient is measured at `1e-5` and `5e-6`, with nearest-pole matching.
It is a diagnostic, not a rigorous remainder bound. The near-zero collective
branch is separate and is not used as a complete stability constraint.
The single-bus interval table samples a fixed log+linear mesh and bisects
observed crossing brackets. Narrow unsampled feasible components are not
excluded.

## Gain box, multimode dual, and signs

For a local mode multiplier vector `y>=0`, let `v_K=G_K^T y`. The exact
support function of a gain movement box `[d_min,d_max]` is

```
sigma_K(v_K) = sum_k v_k * (d_max,k if v_k>=0 else d_min,k).
```

The maximizing coordinate is the corresponding bound (free when `v_k=0`).
For one modal feasibility inequality, the minimum weighted gain movement is
the active-set pseudoinverse solution, clipping violated box coordinates and
resolving on the remaining free coordinates. The code does not implement the
full multimode active-set KKT corrector required for a global claim.

Use constraints `g_j=Re(lambda_j)+sigma_req<=0` and, if an initial-RoCoF
retention requirement is imposed, `g_R=h_req-h'epsilon<=0`. For

```
L = P'epsilon + y'g + mu(h_req-h'epsilon)
    - lambda_lower'epsilon + lambda_upper'(epsilon-1),
```

the unboxed stationarity term is `P+G_epsilon^T y-mu h`. Define

```
Psi_i = -[G_epsilon^T y]_i + mu h_i - P_i.
```

Then `Psi_i<=0` at a continuous lower bound, `Psi_i>=0` at an upper bound,
and `Psi_i=0` in the interior, provided valid nonnegative KKT multipliers
exist. **A removed SG is a different architecture**; no derivative or Psi is
assigned across its removal boundary. The ledger's multipliers are diagnostic
pseudoinverse values, not certified KKT shadow prices. RoCoF was reported as
a post-design metric, so `mu=0` in that ledger.

## Direct normalized robustness

ExpG's frozen uncertainty is the additive full complex block
`A_q + (1 s^-1) Delta`, `||Delta||_2<=beta`. The margin-shifted resolvent is

```
M_sigma(jw) = [jw I - (A_q+sigma_req I)]^{-1}(1 s^-1 I),
beta_star = 1 / sup_w sigma_max(M_sigma(jw)).
```

The frequency search uses a deterministic log+linear mesh and golden-section
refinement around sampled local maxima. `beta` is normalized, not a physical
percentage. The frontier is the best of a declared finite analytical
candidate pool. It is **not** a solved robust Pareto optimum and no
frequency-domain global peak proof is asserted.

## BND explanation

At each active pole, the 20-coordinate exact port closure `T` is partitioned
around a two-coordinate pivot bus `k`:

```
T_eff,k = T_kk + Gamma_k,
Gamma_k = -T_kr T_rr^{-1} T_rk.
```

For any parameter `z`,

```
Gamma_z = -T_kr,z R T_rk - T_kr R T_rk,z
          +T_kr R T_rr,z R T_rk,  R=T_rr^{-1}.
```

The pole derivative splits into direct `T_kk,z` and collective `Gamma_z`
pieces using the same effective-matrix left/right null vectors and the full
`T_eff,s` denominator. The decomposition explains the candidate; it was not
used to replace the exact closure with a graph Laplacian law. The current
pairwise cancellation statistic is descriptive and has not been established
as an invariant graph authority measure.

## Certification boundary

Support enumeration is exhaustive; continuous active-set branch completion
is not. The deterministic correction searched 64 supports locally and left
the others at feasible seeds or unproven seed failures. No LICQ, SOSC,
critical-cone, or global continuous optimum certificate exists. The ExpE
provisional point, evaluated **after** ExpK freeze, has lower retained SG MW
and satisfies the same analytical spectral margin. This directly falsifies
global optimality of the frozen ExpK nominal candidate.
