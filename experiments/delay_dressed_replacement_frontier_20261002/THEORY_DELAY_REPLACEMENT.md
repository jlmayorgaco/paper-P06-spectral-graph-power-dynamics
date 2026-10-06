# Theory: exact delayed PLL channel, retained graph operator, and local replacement bound

Status labels in this file refer to mathematics derived from the frozen model form. They do not imply that the implementation or IEEE-39 result has passed D0/D1.

## 1. Exact linear retarded model for the declared channel

For an interior fixed support, let x collect all retained SG and GFL differential states after exact algebraic network elimination. The no-delay reduced model is dot(x)=F(x;rho,Kp,Ki). For GFL port i, define the scalar nonlinear detector output e_i(x) from the instantaneous network voltage and PLL angle. The only changed channels are
dot(DeltaOmega_i)=(xi_i+Kp_i e_i(x(t-tau_i))-DeltaOmega_i)/tau_lpf,
dot(xi_i)=Ki_i e_i(x(t-tau_i)).
The remaining vector field is unchanged.

At a locked equilibrium x*, e_i(x*)=0. Let A=F_x(x*) be the no-delay Jacobian, c_i^T=de_i/dx at x*, and b_i have two nonzero entries: Kp_i/tau_lpf in the DeltaOmega_i equation and Ki_i in the xi_i equation. Then A_i=b_i c_i^T is rank one and the exact linear DDE is
dot(y)=A_0 y(t)+sum_i A_i y(t-tau_i), with A_0=A-sum_i A_i.
For this explicit-state ODE reduction, E=I and the exact characteristic matrix is
Delta(s)=sI-A_0-sum_i A_i exp(-s tau_i).
This is an exact identity for the declared delayed-error model on the fixed support, provided A and c_i are derivatives of the same exact algebraically reduced vector field. The DDE has infinitely many roots for positive delays; D1 requires numerical root coverage, not merely low residual at a few roots.

At tau=0, Delta(s)=sI-A exactly. This gives an algebraic parity test. The no-delay nonlinear RHS is also recovered by replacing e_i(t-tau_i) with e_i(t).

For a simple root Delta(lambda)v=0 and left vector w^H Delta(lambda)=0,
d lambda/dp = -(w^H Delta_p v)/(w^H Delta_s v),
Delta_s=I+sum_i tau_i A_i exp(-s tau_i),
Delta_tau_i=s A_i exp(-s tau_i).
For the declared PI channels, Delta_Kp_i=-(e_DeltaOmega_i/tau_lpf)c_i^T exp(-s tau_i) and Delta_Ki_i=-e_xi_i c_i^T exp(-s tau_i). Derivatives with respect to rho must include all changes in the network Schur closure and fixed-support Jacobian. The formula is EXACT_IDENTITY for simple roots; computed derivatives remain NUMERICALLY_VALIDATED only after withheld centered finite differences pass. Near multiple roots, report conditioning and do not apply the simple-root formula.

## 2. Exact retained operator

For the executed prototype at the frozen interior support rho_i=0.875, choose q as the ten surviving SG rotor-angle state directions projected into the physical rotational quotient; all other dynamic device states are hidden coordinates c. Partition the exact DDE characteristic matrix:
Delta(s)=[Delta_qq Delta_qc; Delta_cq Delta_cc].
Where Delta_cc(s) is nonsingular, the exact Schur operator is
T_r(s)=Delta_qq-Delta_qc*inverse(Delta_cc)*Delta_cq.
Zeros of det(Delta) and det(T_r) agree locally only where the eliminated block is nonsingular. Record its condition number and verify determinant/pole consistency on test points. Do not claim global equivalence through singular eliminated blocks.

Define D_G(tau)=Sym(dT_r/ds at s=0). Define D_base as the declared SG/network base contribution under the same retained-coordinate convention, L_D=D_G(0)-D_base, and L_tau=D_G(0)-D_G(tau), so D_G(tau)=D_base+L_D-L_tau by definition. L_D and L_tau are graph operators. They are Laplacians only if symmetry, zero row sums, and the required off-diagonal sign conditions are separately verified. The skew part of dT_r/ds is modal rotation/non-normal coupling, not direct dissipation; only x^T Sym(dT_r/ds) x is a quadratic dissipative contribution.

The simple commutator approximation is conditional. In a symmetric special case with a common controller gain and symmetric C, expanding exp(-sT_tau) gives a first-order penalty proportional to Sym(C*T_tau); its skew part is proportional to the commutator [T_tau,C]. The requested formulas L_tau=(Ki/2){C,T_tau} and D_A=-(Ki/2)[T_tau,C] apply only after deriving that exact reduced pencil and checking the chosen signs/coordinates. They must not be imposed on the full heterogeneous IEEE-39 model.

## 3. Physical graph and GSP quantities

Build a physical angle Jacobian from the solved AC operating point and separately build the declared lossless branch approximation B*diag(b_e V_i V_j cos(eta_e))*B^T. For the exact lossy Jacobian report symmetry, row-sum and off-diagonal-sign residuals; if it fails the Laplacian tests, call it the exact angle Jacobian, not a Laplacian. Deflate the rotational mode.

On interior supports with positive retained inertia M, use the mass-normalized lossless physical operator M^(-1/2)L_P M^(-1/2)=U Lambda U^T. If any retained inertia is zero, this normalization is singular and the endpoint must be treated separately; do not add silent regularization. Project the symmetric delay-dressed operator to obtain Dhat_G=U^T M^(-1/2)D_G M^(-1/2)U and eta_cross=norm(offdiag(Dhat_G),F)/norm(Dhat_G,F). Treat the frozen tau vector as a graph signal, tauhat=U^T tau, and graph roughness r_tau=tau^T L_P tau. The scale-normalized commutator metric is chi_tau=norm([L_P,diag(tau)],F)/norm(L_P,F). For symmetric L_P, its squared numerator equals 2*sum over undirected edges w_ij^2*(tau_i-tau_j)^2. This is an EXACT_IDENTITY, not evidence that chi_tau predicts capacity.

## 4. Lossless nonlinear power Taylor identity

For branch incidence matrix B and declared branch orientation, P(theta)=B*diag(gamma_e)*sin(B^T theta), eta*=B^T theta*, and q a phase perturbation:
P(theta*+q)=P(theta*)+L_P q + 1/2 Q(q,q)+1/6 C(q,q,q)+O(norm(q)^4),
L_P=B*diag(gamma_e cos(eta_e*))*B^T,
Q(q,q)=-B*diag(gamma_e sin(eta_e*))*(B^T q)^(elementwise 2),
C(q,q,q)=-B*diag(gamma_e cos(eta_e*))*(B^T q)^(elementwise 3).
This is an EXACT Taylor identity for the lossless sinusoidal branch model. It is not an exact Taylor model of lossy AC power, nor a nonlinear stability proof. TABLE_D03 must compare direct branch-power differences against this expansion over preregistered directions and amplitudes.

## 5. Local controller authority and reduced retention LP

Let alpha be a simple active spectral abscissa root at a frozen operating point. Linearize the active stability constraint in epsilon and controller coordinates u=log(K/K_ref). With a=-grad_epsilon(alpha), c=-grad_u(alpha), and a consistently computed intercept b, the local inequality has form a^T epsilon+c^T u>=b. Under the declared controller trust budget u^T R_K u<=kappa^2, Cauchy-Schwarz gives the local maximum controller authority
A_K=kappa*sqrt(c^T inverse(R_K)c),
u*=kappa*inverse(R_K)c/sqrt(c^T inverse(R_K)c).
This is a reduced-model identity for the linearized subproblem, not a nonlinear design guarantee. The remaining local synchronous-support constraint is a^T epsilon>=b_S=[b-A_K]_+.

A separate instantaneous RoCoF condition, if the exact GFL model has no synthetic-inertia feedthrough, is h^T epsilon>=b_R. Derive h,b_R from initialized inertia and the event injection, then validate against a direct short-time simulation. Do not assert this necessary condition before that check passes.

The one-mode plus one-RoCoF reduced LP is min P^T epsilon subject to a^T epsilon>=b_S, h^T epsilon>=b_R, and 0<=epsilon<=1. PROVED_REDUCED_MODEL: every extreme-point solution has at most two partially retained generators. Proof: if k coordinates are strictly interior, the remaining n-k active box constraints have rank at most n-k. An extreme point in n dimensions needs n linearly independent active constraints, while at most two additional non-box inequalities exist; therefore k<=2. This statement applies only to this LP with its declared active-constraint structure. It does not apply to the full multi-event nonlinear/DDE design.

For a proposed partial pair i,j with nonsingular matrix, solve [a_i a_j; h_i h_j]*[epsilon_i;epsilon_j]=[b_S_bar;b_R_bar], then check 0<=epsilon_i,epsilon_j<=1 and every other constraint. Enumerate all one- and two-anchor candidates and compare with the LP solver. Singular pairs and boundary anchors must be handled explicitly.

## 6. Status boundary

EXACT_IDENTITY: delayed-channel characteristic matrix, simple-root sensitivity formula, lossless branch Taylor coefficients, symmetric-graph commutator identity, and the reduced-LP two-anchor theorem under their stated assumptions.

PROVED_REDUCED_MODEL: only the local authority formulas and two-constraint LP theorem after executable checks of assumptions and solver agreement.

NUMERICALLY_VALIDATED: reserved for passing source-matched D0/D1 tests, finite-difference sensitivities, and reduced/full-model comparisons.

SUPPORTED_LOCAL: local optimizer or local sensitivity evidence only.

SUPPORTED_EMPIRICAL: finite fixed delay-pattern/event observations only.

NEGATIVE_RESULT: registered hypothesis fails its preregistered metric/gate.

BLOCKED: model parity, DDE root coverage, nonlinear DDE integration, or required inputs are unavailable. No optimality/near-optimality claim without a validated upper bound.

## 7. Executed outcome and limits

- D0 passed no-delay parity for three stored designs: the maximum matched finite-pole error is 5.90e-10 and the maximum rightmost-pole error is 3.22e-11. This confirms implementation parity only.
- The historical 90.047% GFL co-design is not a fully feasible lower point under the frozen five holdouts: its reproduced frequency peaks are 0.50569 Hz at bus 8, -100 MW and 0.51772 Hz at bus 16, -100 MW, above the 0.5-Hz cap. Holdout actuator slack was not recomputed here.
- D1 is BLOCKED_EXACT_DDE_SPECTRUM. Tau=0 identity and local tracked roots pass residual checks, but the positive-delay rightmost-root count is absent.
- D2 finite-difference checks validate only tau/Kp/Ki sensitivities of the locally tracked pair; rho sensitivity and root completeness are absent.
- The lossless branch Taylor expansion passed 60 numerical checks, with median relative third-order error 5.30e-8 at the maximum 0.1-rad amplitude. No lossy AC tensors or nonlinear stability conclusion were validated.
- D4 is BLOCKED: after gauge deflation, the hidden Schur block at s=0 is singular (condition estimate about 3.10e18). Therefore D_G, L_D and L_tau are not defined by the attempted partition. The pre-gauge values are retained only in files explicitly named PRE_GAUGE_DIAGNOSTIC and must not support a scientific claim. See D04_GATE_REPORT.md.
- Consequently no delayed co-design, fully validated lower frontier R_F, delay shadow price, rigorous upper frontier R_U, or optimality gap is available. The reduced two-anchor LP theorem remains a conditional mathematical result, not a description of the full optimum.
