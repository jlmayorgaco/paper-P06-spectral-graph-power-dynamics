# ExpH derivation: Woodbury gain authority and collective self-energy

## One-scalar gain dependence per inverter port

For the frozen analytical GFL port, the gain-dependent state-matrix terms share one PLL phase-error covector. Write the PLL contribution as

\[
A(K)=A_0+h(K)e_x^T,\qquad B(K)=B_0+h(K)e_u^T,
\quad h(K)=K_p h_p+K_i h_i.
\]

For resolvent (R_0(s)=(sI-A_0)^{-1}), the matrix inversion lemma gives a single scalar denominator

\[
d(s,K)=1-e_x^TR_0(s)h(K),
\]

and the port transfer

\[
Y_{GFL}(s,K)=D+C\left[R_0B_0+R_0h(K)
\frac{e_x^TR_0B_0+e_u^T}{d(s,K)}\right].
\]

This is an exact identity for the model, not a fitted polynomial. `src/bnd_design_h/BNDDesignH.jl` implements it and audits the one-scalar factorization per bus in `tables/TABLE_H06_woodbury_gain_structure.csv`. The representation is low-rank per device; it does not by itself prove a globally low-dimensional controller optimization.

## Active gain direction

ExpH forms the normalized mixed-authority Jacobian (J_K=\partial A/\partial (\log K_p,\log K_i)), then takes its right singular vectors. The first singular direction captures more than 99.9% of the Frobenius energy, and the numerical rank is one at all three declared tolerances. This is the exact first-order active basis (V_{K,active}) used by the diagnostic; it is not selected from a graph heuristic.

The gain dependence is nonzero, while the final robust candidate's (K_p,K_i) differ from nominal only at roughly floating-point / (10^{-8})-relative scale. Thus the data support mixed authority, but not a meaningful realized PLL retuning in the final candidate.

## Direct term and Feshbach self-energy

Partition the exact port operator into the critical and remaining coordinates. Eliminating the latter gives

\[
F=T_{cc}-T_{cr}T_{rr}^{-1}T_{rc},\qquad
\Gamma=-T_{cr}T_{rr}^{-1}T_{rc}.
\]

The reported retained-SG authority splits into a direct term and the derivative of the collective self-energy. Across buses the magnitudes of these two terms are often hundreds or thousands of times the net authority, with opposite signs. The tabulated `self_energy_share` is therefore a cancellation ratio, not a literal percentage of damping. For example, bus 30 has direct (-1256.61), self-energy (+1256.83), and net (0.214). This is strong evidence that nodal-only damping is inadequate for this coefficient. Path-by-path (l,m) contributions and cancellation ratios are in `tables/TABLE_H08_self_energy_pathways.csv`.

## Noncommutative graph reduction

The physical conductance and susceptance port operators do not commute. ExpH uses the word span \(I,L_G,L_B,L_G^2,L_GL_B,L_BL_G,L_B^2\) without simultaneous diagonalization. Its degree-two graph basis captures 99.24% of the active mixed-gain subspace energy, with principal angle about 5.00 degrees. This is a useful empirical reduction for this case, not evidence for an exact single-Laplacian modal law.

## Design limits

The exact local coefficient predicts bus 38 as the best retained-SG anchor at 7.38 MW. Exact all-pole checking rejects that branch; the best robust candidate found by the declared deterministic \(0.02\)-epsilon grid retains bus 36 at 11.2 MW. This is an exact feasible grid point for the analytical spectral and direct small-gain constraints, not a continuum optimum. The active transient frequency constraint was not incorporated into a KKT solve, and the independent PowerDynamics spectral abscissa misses the required margin. Accordingly, this derivation supports structural compression and self-energy mechanisms, not a completed globally certified co-design.
