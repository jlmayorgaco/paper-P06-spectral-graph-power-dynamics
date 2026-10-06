# Analytical SG→GFL co-design: derivations and claim limits

This document separates exact identities of the declared IEEE-39 fixed-interior model from a proposed modal approximation. All delay values are exogenous. The design variables are retained SG fractions `epsilon_i=1-rho_i` and positive PLL gains `Kp_i,Ki_i` within the frozen bounds. Network equilibrium, port power sharing, and ratings follow the repository model, not an independently invented replacement rule.

## 1. Nonlinear DAE with delayed PLL measurements

Let `x` be differential device states and `z` algebraic network voltages. The fixed-architecture model has the form

`x_dot(t)=f_0(x(t),z(t);rho,K)+sum_i b_i(K) e_i(x(t-tau_i),z(t-tau_i))`,

`0=g(x(t),z(t);rho)`.

The source PLL equations use a single scalar phase-detector error in both PI branches. In the repository notation, `b_i=Kp_i b_{p,i}+Ki_i b_{I,i}`. The delay applies to that measured scalar, not to current injection, current-control states, power setpoints, or the complete converter.

At an equilibrium `(x*,z*)`, delayed and present states coincide. Consequently fixed measurement delay does not change the equilibrium equations, though changing `rho` or gains may change equilibrium initialization and must be re-solved.

## 2. Exact algebraic elimination

Assume `g_z` is nonsingular at the equilibrium. Linearized KCL gives `delta z(t)=-g_z^{-1}g_x delta x(t)` and the same relation at `t-tau_i`. The exact network-reduced PLL measurement row is

`cbar_i^T=e_{i,x}-e_{i,z} g_z^{-1}g_x`.

The non-delayed differential Jacobian after network elimination is `f_x-f_z g_z^{-1}g_x`, with instantaneous terms collected consistently in `A0`. Thus the fixed-equilibrium linear DDE can be written

`delta x_dot(t)=A0 delta x(t)+sum_i b_i cbar_i^T delta x(t-tau_i)`.

This is exact for the linearization and the chosen state ordering where `g_z` is nonsingular. It is not a nonlinear model reduction or a statement about model variants with hard current limiting not represented in the frozen equations.

## 3. Exact low-rank characteristic closure

With `B=[b_1,...,b_m]`, `C=[cbar_1^T;...;cbar_m^T]`, and `E_tau(s)=diag(exp(-s tau_i))`,

`Delta(s)=sI-A0-B E_tau(s) C`.

For `Ac=A0+BC`, this becomes `Delta=sI-Ac-B(E_tau-I)C`. Where `sI-Ac` is nonsingular, the determinant lemma gives

`det Delta(s)=det(sI-Ac) det[I_m-(E_tau-I)C(sI-Ac)^{-1}B]`.

The individual delayed PLL contribution is `b_i cbar_i^T`, rank at most one. The small determinant is meromorphic; poles of the reference resolvent must be accounted for in any root count. The identity alone does not decide DDE stability.

## 4. Joint descriptor action space

Keep the augmented differential/algebraic descriptor ordering fixed at an interior SG/GFL architecture. Changing `epsilon_i` modifies the two real KCL terminal rows, so each replacement action has rank at most two in that descriptor. Changing both gains of one PLL uses the same phase-detector row, so its gain action has rank at most one. For ten terminals,

`T(s;y,tau)=T_ref(s,tau)+U Theta(s;y,tau) V^H`, with `rank(U Theta V^H)<=20+10=30`.

Where `T_ref` is nonsingular,

`det T=det T_ref det[I_r+Theta V^H T_ref^{-1}U]`, `r<=30`.

This is a fixed-coordinate identity under the source scaling contract. It does not cover architecture endpoints that change state dimension, singular reference pencils, or a claim that the Kron-reduced state matrix is affine in retention.

## 5. Conditional exact PI boundary law

Freeze every other device. Remove the `i`th variable PI update to form `Delta_{-i}(s)`. At a target `s_b=-sigma+j omega`, assume `Delta_{-i}(s_b)` nonsingular. Define `g_p=c_i^T Delta_{-i}^{-1}b_{p,i}` and `g_I=c_i^T Delta_{-i}^{-1}b_{I,i}`. The exact conditional target-root condition is

`Kp_i g_p(s_b)+Ki_i g_I(s_b)=exp(s_b tau_i)`.

Writing `g_p=a+jb`, `g_I=c+jd`, and `D_g=ad-bc != 0` gives

`Kp_i=exp(-sigma tau_i)[d cos(omega tau_i)-c sin(omega tau_i)]/D_g`,

`Ki_i=exp(-sigma tau_i)[a sin(omega tau_i)-b cos(omega tau_i)]/D_g`.

The formula places the chosen `s_b` on the characteristic equation. It does not establish that it is the rightmost pole, that the other roots satisfy the margin, or that the resulting gains lie in the frozen box.

## 6. Conditional scalar modal PI formula

Suppose a *validated* critical-mode reduction has scalar characteristic

`m(rho)s^2+d(rho)s+ell(rho)+g(s,rho)(Ki+sKp)exp(-s tau)=0`.

For `omega != 0`, `g(s_b,rho) != 0`, define

`Q=-exp(s_b tau)[m s_b^2+d s_b+ell]/g(s_b,rho)`.

Since `Ki+s_b Kp=(Ki-sigma Kp)+j omega Kp`, the unique real gains that make this scalar equation vanish are

`Kp=Im(Q)/omega`, `Ki=Re(Q)+sigma Im(Q)/omega`.

This is an exact algebraic solution **of that scalar model**, conditional on the model retaining the correct critical root and interactions. For multiple GFLs the projected determinant generally contains cross-coupling `I_cross(s)`; dropping it requires a measured approximation bound or a matrix-Rouché certificate. No IEEE-39 scalar-model accuracy is presumed here.

## 7. Exact simple-root sensitivities

For a simple root `T(lambda,p)v=0`, `w^H T(lambda,p)=0`, and nonzero `w^H T_s v`, implicit differentiation gives

`lambda_p=-(w^H T_p v)/(w^H T_s v)`.

For `Delta=sI-A0-sum_i A_i exp(-s tau_i)`, fixed-matrix partial derivatives are `Delta_s=I+sum_i tau_i A_i exp(-s tau_i)` and `Delta_{tau_i}=s A_i exp(-s tau_i)`. Design derivatives in `rho` must also differentiate the network closure, ratings, and equilibrium; they cannot be replaced by the PLL delay term. Near multiple roots the simple-root expression is ill-conditioned or invalid. Numerical use requires finite-difference checks.

## 8. Local controller authority

Let `alpha=Re(lambda_c)` for a verified simple active root, and let `u` denote a local admissible gain update. A first-order expansion gives `delta alpha=g_epsilon^T delta epsilon+g_K^T u`. With `a=-g_epsilon`, `c=-g_K`, write a local required inequality `a^T epsilon+c^T u>=b`. Under positive-definite `R_K` and `u^T R_K u<=kappa^2`, Cauchy-Schwarz yields

`A_K=kappa sqrt(c^T R_K^{-1} c)`,

`u*=kappa R_K^{-1}c/sqrt(c^T R_K^{-1}c)` when `c!=0`.

For an admissible box `lower_j<=u_j<=upper_j`, the exact support of this *linearized* stability functional is `A_K_box=sum_j max(c_j lower_j,c_j upper_j)`. These are exact convex support calculations for the local model, not exact full nonlinear stability bounds.

For a two-component real/imaginary modal correction `A_K k=b`, the gain box maps to the zonotope `{A_K k: lower<=k<=upper}`. Membership shows reachability of that linear correction only; full spectral and event constraints remain to be checked.

## 9. Reduced retention LP

If the validated local stability floor after optimal controller authority is `b_S=[b-A_K]_+`, and an independently derived instantaneous RoCoF necessary condition is `h^T epsilon>=b_R`, the proposed reduced problem is

`min_{0<=epsilon<=1} P^T epsilon` subject to `a^T epsilon>=b_S` and `h^T epsilon>=b_R`.

This finite LP has a global reduced-model optimum. The instantaneous RoCoF row requires the frozen GFL model to have no instantaneous synthetic-inertia feedthrough and must be derived from initialized physical inertia and the declared disturbances. If that derivation is unavailable, this LP is a conditional theorem rather than an evaluated physical bound.

## 10. Partial-retention theorem and anchors

Consider a bounded feasible LP in `n` retention variables and at most `m` independent system-wide active constraints. A linear objective attains an optimum at an extreme point. If `k` variables are strictly interior to their `[0,1]` bounds, at most `n-k` independent bound constraints are active. An extreme point needs `n` independent active constraints. Hence `n<=(n-k)+m`, so `k<=m`. For the one-stability/one-RoCoF LP, **there exists** an optimum with at most two partially retained SGs. This says nothing about a non-extreme alternative optimum or the full nonlinear multi-event optimum.

After fixing the `n-2` other variables at bounds, two partial anchors `i,j` satisfying both active rows solve

`[epsilon_i;epsilon_j]=[a_i a_j;h_i h_j]^{-1}[bar_b_S;bar_b_R]`,

provided this matrix is nonsingular and the resulting retentions lie in `(0,1)`. Cases with one/zero partial anchors, singular anchor pairs, and slack rows must also be enumerated.

## 11. Finite-frequency angular-port damping

The previously attempted Schur slope at `s=0` is undefined for its gauge-deflated partition: the hidden block condition is around `3.10e18`. For any alternative *homogeneous physically conjugate angular port* with a nonsingular hidden elimination at `s=j omega`, define

`D_G(omega,tau)=[T_theta(j omega,tau)-T_theta(j omega,tau)^H]/(2j omega)` for `omega!=0`.

On a synthetic quadratic pencil `T_theta(s)=s^2 M+sD+L` with real symmetric `M,D,L`, this returns exactly `D`. Outside those assumptions it is a Hermitian frequency-resolved operator, not automatically physical dissipation: asymmetric stiffness, nonnormality, or coordinate scaling can contribute. A delay penalty is conditionally `L_tau(omega)=D_G(omega,0)-D_G(omega,tau)`. No IEEE-39 graph-damping claim follows without port construction, gauge consistency, and pole-tracking validation.

## 12. Matrix-Rouché modal certificate

Let the full retained modal pencil be `T_modal(s)=D(s)+E(s)` with `D` analytic and nonsingular on a closed contour `Gamma`. If both terms are analytic in the enclosed domain and `sup_{s in Gamma} ||D(s)^{-1}E(s)||_2<1`, then the homotopy `D+tE` is nonsingular on `Gamma` for all `t in [0,1]`. Consequently `det(D+E)` and `det D` have the same number of zeros inside, counted with multiplicity. If the reduction is meromorphic, poles and removable singularities must be counted explicitly. A sampled maximum is not a certificate for the supremum.

## 13. Limitations and validation chain

The determinant and conditional PI formulas are exact at a fixed, regular linearization. Modal reduction, local pole sensitivities, ellipsoidal/box authority, and the retention LP become physical design tools only after error, mode-switching, and constraint validation. Full delayed feasibility requires a trustworthy DDE root oracle plus nonlinear delayed simulations of all frozen events. The zero-delay event witness alone cannot establish a positive-delay frontier or global optimum. A rigorous full-model upper bound additionally requires validated remainders or branch-and-bound; without one, report only the best fully validated design found.
