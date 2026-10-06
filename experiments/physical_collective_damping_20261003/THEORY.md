# Physical ports and finite-frequency collective damping

## 1. Physical normalization (EXACT_IDENTITY for this model)

Let nu_i=omega_i-1 (per-unit rotor-speed perturbation), and let
S_i=(1-rho_i) S_i^rated. The implemented swing equation has
2 H_i dot(omega_i)=tau_m,i-tau_e,i-D_i(omega_i-1).
Add an external torque perturbation d_i/S_i, where d_i is torque expressed as
MW at nominal speed. Then M_i=2 H_i S_i and the additive input is d_i/M_i.
The incremental mechanical kinetic-energy function is
E_inc=1/2 sum_i M_i nu_i^2, in MW s, and the input contribution to its derivative
is sum_i d_i nu_i. This is an incremental torque-speed supply, not total AC power
or a Lyapunov function for the entire controlled nonlinear network.

## 2. Exact finite-frequency port elimination (EXACT_IDENTITY)

After exact algebraic voltage elimination and linearization of the full nonlinear
device equations at a fixed equilibrium:

    Delta(s)=s I-A0-B diag(exp(-s tau_i)) C^T.

B injects the delayed PLL detector into PLL frequency/integrator derivatives;
C^T differentiates the actual terminal-voltage PLL detector. Partition states
as (nu,z), retaining all ten physical SG speeds. For invertible Delta_zz(s):

    S(s)=Delta_nn-Delta_nz Delta_zz^{-1} Delta_zn,
    Z(s)=M S(s),             d=Z(s) nu,
    Y(s)=C_n Delta(s)^{-1} B_d=Z(s)^{-1}.

Here B_d has entries 1/M_i at speed rows. This is an exact frequency-dependent
impedance of the full-order LINEARIZATION, including hidden SG, network and
converter dynamics. It is not an exact nonlinear transfer, a constant damping
matrix or a Laplacian. At hidden-block poles this representation
is undefined; no determinant equivalence across such poles is claimed.

For real omega>0 define D_H(omega)=(Z(j omega)+Z(j omega)^H)/2. For harmonic
physical phasors, average incremental mechanical input is

    <d(t)^T nu(t)>=1/2 Re(nu_hat^H d_hat)
                   =1/2 nu_hat^H D_H(omega) nu_hat.

The remaining skew-Hermitian part has zero real harmonic supply. This physical
normalization, not matrix density or symmetry alone, justifies the damping term.

## 3. Nodal and collective terms (EXACT_IDENTITY)

    D_node=diag(diag(D_H)), D_cross=D_H-D_node,
    W_node=1/2 nu_hat^H D_node nu_hat,
    W_cross=1/2 nu_hat^H D_cross nu_hat.

The real total is W_node+W_cross. A sign reversal between W_node and the total
means independent nodal coefficients fail for that specified motion and frequency.
It is not proof that removing off-diagonal terms is a physically possible design.
The physical nodal coordinate system and retained ports are fixed throughout.

For frozen rho and gains define the exact finite-frequency delay difference
D_tau(omega)=D_H(omega,0)-D_H(omega,tau). It need not be positive semidefinite:
latency can increase some projected damping and decrease other directions.
This is not the previously blocked derivative of a Schur operator at s=0.

## 4. Pole and design connection (SUPPORTED_LOCAL only when validated)

For a simple characteristic root s_c with Delta(s_c)v=0, and invertible
Delta_zz(s_c), Z(s_c) v_nu=0. At a Hopf point s_c=j omega_c the complete complex
equation implies v_nu^H D_H(omega_c) v_nu=0. The scalar equality alone is not
sufficient to locate a pole. Away from the imaginary axis evaluate both
Z(s_c) and D_H(Im(s_c)); do not identify harmonic damping with Re(s_c).

Local parameter sensitivity remains
ds_c/dp=-(w^H Delta_p v)/(w^H Delta_s v).
Equivalently ds_c/dp=-(ell^H Z_p r)/(ell^H Z_s r) where the port root is simple.
The latter can attribute design changes through physical interactions, but this
standard identity is not a new optimization algorithm. Any claim that the
mechanism raises feasible replacement requires a later constrained reoptimization
and the entire frozen nonlinear event set.

## 5. Scope of nonlinear test

For constant-delay retarded dynamics use the original nonlinear RHS and replace
only the PLL detector e_i(t) by e_i(t-tau_i). A method of steps evaluates past
states using completed dense ODE solutions, with step intervals no longer than
the smallest positive delay. History x*=x_eq+Re(a v exp(s_c t)) excites the
tracked branch. A matching linear DDE uses the same history. Agreement as a->0
tests the derivative and delayed implementation; it does not prove large-event
robustness, attraction-region size, or a complete spectral certificate.

## 6. Exploratory follow-up: collective action on an actual pole

Added after observing that positive harmonic damping and an unstable PLL root
coexist. This is a diagnostic follow-up, not a changed preregistered hypothesis.
Let X=Delta_zz^{-1} Delta_zn. For any parameter p with constant M,

    Z_p=M [Delta_nn,p-Delta_nz,p X
         -Delta_nz Delta_zz^{-1}(Delta_zn,p-Delta_zz,p X)].

Use the same expression with p=s. For simple port root with left/right vectors
ell,r, split Z_p into its physical nodal diagonal and off-diagonal remainder:

    ds/dp = -ell^H diag(diag(Z_p)) r / (ell^H Z_s r)
            -ell^H offdiag(Z_p) r / (ell^H Z_s r).

This is an EXACT_IDENTITY under the earlier regularity conditions. Validate
against the full characteristic sensitivity and centered finite differences.
It quantifies collective action on the actual pole, while retaining reactive
dynamics and left/right mode geometry. Neither summand alone is an independently
realizable intervention. The decomposition is coordinate/port specific, with
physical SG torque-speed coordinates fixed here. Its identity is standard
eigenvalue sensitivity calculus; no novelty claim is attached to the formula.
