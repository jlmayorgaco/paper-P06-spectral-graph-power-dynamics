# Experiment G derivation

## Isolated PLL seed

Linearizing the implemented phase detector at aligned terminal voltage gives e = -V delta-theta. Together with dot(theta)=delta-omega, tau dot(delta-omega)=delta-omega-i + Kp e - delta-omega, and dot(delta-omega-i)=Ki e, elimination yields tau s^3+s^2+V Kp s+V Ki=0. Matching tau(s+r)^2(s+q) and tau(2r+q)=1 gives q=1/tau-2r, Kp=(2r-3 tau r^2)/V, and Ki=(r^2-2 tau r^3)/V. Since the sum of the three roots is -1/tau, a common decay r obeys 3r<=1/tau; the triple-root seed is r=1/(3tau). It is only a local reference.

## Collective authority and endpoint singularity

The all-GFL state matrix has an exact global-angle gauge vector g, checked by ||Ag||/(||A||||g||). The quotient basis Q=null(g') defines Aq=Q' A Q, whose spectrum removes precisely one gauge root. ExpG estimates gamma_i=-d Re(lambda_c)/d epsilon_i by centered full-state differences and compares steps 1e-4 and 5e-5. Because the endpoint also has a physical zero and the unreduced zero is defective, this derivative is a local quotient-branch seed, not a full design law.

## Finite surrogate enumeration

For frozen gains and a frozen robustness radius, the surrogate is min P'epsilon subject to gamma'epsilon>=bS, h'epsilon>=bR, 0<=epsilon<=1. Every upper-bound subset is enumerated; with two active scalar inequalities, a vertex has at most two remaining fractional coordinates. Single-anchor and pair intersections are solved explicitly. This certifies the linearized authority surrogate only.

## Robust margin

For the reduced, gauge-free matrix Aq, the full-block uncertainty is Aq+scale*Delta, ||Delta||2<=beta. On the shifted boundary s=-sigma+j omega, ExpG adaptively refines the largest singular value of (sI-(Aq+sigma I))^-1; the numerical small-gain test is beta||M_sigma||infinity<1. It also computes R_delta=sum ||r_j||||l_j||/|l_j'r_j| and imposes the conservative nominal separation alpha<=-sigma-beta R_delta. The complex full block is an outer uncertainty model, not an exact statement about structured real parameter uncertainty.

## RoCoF and transients

With h_i=H_i S_i, the initial COI bound is |df/dt(0+)|=f0|DeltaP|/(2 h'epsilon). The chosen illustrative scenario is recorded in DESIGN_SCENARIO_FROZEN.toml; it is not a universal grid standard. The final transient uses the same descriptor-reduced A, a one-MW current disturbance at the declared load bus, and an inertia-weighted retained-SG COI frequency output. Modal residues and direct matrix exponentials are compared on the same frozen linear model.

## Exactness and scope

The final candidate is checked using the complete finite spectrum of Ared, not only det(C). The closure determinant identity is tested at a regular complex point, and the low-rank determinant lemma is kept as a boundary predictor. The deterministic coordinate candidate is not a globally certified full-closure optimum; no KKT/SOSC certificate is claimed.
