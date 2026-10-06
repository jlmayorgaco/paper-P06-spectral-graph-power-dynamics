# Graph expansions and feasible joint design

All exact statements below concern a regular fixed-interior model. A local
linearization is not an exact nonlinear transfer or a global stability theorem.

## Physical graph and controller signals

L=B diag(b_e |V_i||V_j| cos(theta_i-theta_j)) B^T is the declared lossless
synchronizing approximation. Eliminate passive nodes by a Schur complement to
obtain L_P on the ten device buses. Check symmetry, row sums, off-diagonal signs,
and positive semidefiniteness. The actual lossy AC Jacobian remains in the plant;
it is NOT replaced by L_P. Let L_P U=U Lambda, U orthogonal, u_0=1/sqrt(10).

Let Phi_q=[u_0,...,u_q]. Gains are nodal graph signals:
log Kp=log Kp_base + Phi_q a, log Ki=log Ki_base + Phi_q b.
The full constraint gradient pulls back by the chain rule, J_a=J_logKp Phi_q.
This is EXACT_IDENTITY for that parameterization; its usefulness is empirical.
Rho remains a free ten-dimensional nodal vector in all experiments.

For a linear score c and a Euclidean trust ball in log gains, projected authority
is kappa ||Pi_q c||, versus kappa ||c|| for free gains. Thus loss is characterized
by the discarded gradient component, not by low graph frequency alone. Different
gain metrics require the corresponding weighted projection/support function.
For the actual box used here, compute its support by LP; do not reuse the ball law.

## Frequency-resolved collective interaction and graph walks

At the physical SG torque-speed ports, Z(s)=M Schur(Delta(s)); D_H=(Z+Z^H)/2
on the imaginary axis. Define Zhat=U^T Z U=D+E, with D its diagonal. If D is
invertible and r=||D^-1 E||_2<1, then

    Zhat^-1 = sum_{k=0}^infinity (-D^-1 E)^k D^-1,
    ||remainder after k=p||_2 <= r^(p+1)/(1-r) ||D^-1||_2.

This is a standard Neumann-series theorem. Matrix powers represent weighted walks
among GRAPH MODES. It becomes a useful graph approximation only if the condition
holds and low-order truncation is accurate in the frequency band that matters.
It is not automatically valid near a limiting PLL root.

For any D_H, U^T[L_P,D_H]U has entry (lambda_a-lambda_b) Dhat_H,ab. Hence the
commutator measures graph-modal mixing with spectral-gap weights. A commuting
operator has no off-diagonal entries between distinct Laplacian eigenvalues;
degenerate eigenspaces require block statements. This does not prove stability.

## Exact parametric network closure

At the frozen initialized device operating points, the full-dispatch SG/GFL
terminal powers agree. Under rating sharing the same trim remains an equilibrium.
Build local zero-PLL-error-action differential Jacobians A_dev, B_v and terminal
current derivatives C_SG, C_GFL and D_SG. Then

    G(rho)=Y_Kron+sum_i (1-rho_i) D_SG,i,
    C_out(rho)=sum_i [(1-rho_i) C_SG,i+rho_i C_GFL,i],
    V_x=-G^-1 C_out,
    A0=A_dev+B_v V_x, C_PLL=E_theta+H_v V_x.

Thus Delta=sI-A0-B_PLL(K) exp(-s tau) C_PLL is an exact characteristic of the
linearized physical model. For a replacement direction,
V_x,rho=-G^-1(G_rho V_x+C_out,rho). This gives exact local root derivatives,
including the network effect of replacement, without assuming an affine reduced
state matrix. The exported model is checked against direct automatic derivatives.

## Nonlinear delayed discrete tangents

For delayed scalar error e_h with parameter derivative E_h, current error e(x,p),
and PLL injection B(p), use f_delay=f_ODE+B(e_h-e). Then

    A_current=f_ODE,x-B e_x,
    f_parameter=f_ODE,p+B(E_h-e_p)+B_p(e_h-e).

The past derivative is E_h=e_x(x_h) S_h+e_p(x_h), evaluated with the network
applicable at that history time. Implicit-stage derivatives solve the same
I-h gamma A_current factorization as the nonlinear stage. The derivative is exact
for the declared discrete method; convergence to the continuous DDE derivative
requires mesh verification. Initial trim sensitivity is zero for this sharing
contract, and pre-event error is zero for every design.

## Local constrained design and its limits

For p=(rho,log Kp,log Ki), linearized security constraints are g(p)+J dp<=0.
Pull them into each family's gain basis, keep identical nodal trust/gain bounds,
and maximize P0^T d rho. The LP optimum is global ONLY for this local linear model.
The adaptive nonlinear DDE and numerical spectral oracle accept/reject the
candidate. Trust fractions are fixed in PROTOCOL.md. No global convergence or
maximum nonlinear replacement follows. Failure of a graph subspace is retained.

## What graph expansion failure does and does not prove

The norm test r<1 is sufficient, not necessary. For a finite square matrix K,
the Neumann matrix series converges to (I+K)^-1 precisely when its spectral
radius is below one. TABLE_05B checks this stronger distinction independently.
Having r>1 alone is not a divergence proof. Having an eigenvalue of K outside
the unit disk makes the infinite matrix series diverge, even if one low-order
truncation happens to reduce error. None of these statements invalidate GSP as
a change of coordinates or as a design parameterization.

If graph-modal interactions cannot be treated perturbatively, split Zhat into
retained graph coordinates l and eliminated coordinates h. Where Z_hh is
invertible the exact dynamic closure is

    Z_l,eff = Z_ll - Z_lh Z_hh^-1 Z_hl.

This closure retains the feedback through the excluded graph coordinates. It
does not amount to discarding high graph frequencies, and need not be a static
Laplacian, symmetric operator or positive damping. It is an existing Schur
identity, offered as the correct algebraic alternative, not a new tested
compression method in this preregistered pilot.

## A possible next basis, distinguished from executed evidence

Low physical graph frequency need not coincide with controller authority for
the active constraints. Let Jp and Ji collect normalized gain derivatives of
declared constraints. For an orthonormal nodal basis Phi and fixed positive
semidefinite constraint weights W, define

    S = Jp^T W Jp + Ji^T W Ji,
    captured authority energy = tr(Phi^T S Phi).

For a fixed dimension and an imposed constant mode, the other eigenvectors of
the projected S maximize this Frobenius derivative-energy criterion. Adding a
declared graph roughness penalty beta*tr(Phi^T L_P Phi) replaces S by S-beta L_P
on the constant-orthogonal space (Rayleigh--Ritz/Ky Fan principle). This is a
standard variational result. It supplies a precise way to trade graph smoothness
against physical control authority, rather than assuming that the first q
Laplacian modes must be best. It does NOT directly maximize replacement or
provide a safety guarantee. Such a basis needs a separate, preregistered test
with held-out designs/events before making a comparative claim. It is not added
post hoc to the primary family comparison here.
