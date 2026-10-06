# ExpH derivation: gauge quotient and retained-SG root law

## Question and conclusion

The all-GFL endpoint has a rotational gauge mode and a second zero algebraic direction. Removing the gauge analytically leaves one simple physical zero. The retained-SG perturbation moves that quotient root linearly in retained fraction; it does not produce a square-root Puiseux branch. The first-order coefficient depends on PLL gains.

## Exact gauge quotient

Let (A_{red}) be the finite reduced state matrix obtained from the frozen collective model, and let (g) be the explicitly constructed uniform-angle vector. ExpH verifies (A_{red}g=0) to numerical precision and constructs an orthogonal complement (Q\) before computing quotient nullity or spectrum. It does not remove a numerically selected eigenvalue. At the all-GFL endpoint, the full matrix has nullity one and algebraic multiplicity two; the quotient has nullity one and algebraic multiplicity one for nominal gains, three independent gain patterns, and both gain-box corners. The quotient classification is `SIMPLE_ZERO_AFTER_QUOTIENT`; the full gauge chain has length two.

The residual physical zero across six gain patterns is tabulated in `tables/TABLE_H01_gauge_quotient.csv`. The quotient Jordan and chain residual audit is in `tables/TABLE_H02_jordan_chain.csv`.

## Linear scaling from the quotient determinant

For one retained-SG support (i), factor the common gauge root from the exact port-closure determinant. Around ((s,\epsilon_i)=(0,0)), the remaining scalar quotient determinant has the local form

\[
F_i(s,\epsilon_i,K)=a_2(K)s+a_{1,i}(K)\epsilon_i+
O(s^2,s\epsilon_i,\epsilon_i^2).
\]

The simple quotient root therefore satisfies

\[
\lambda_{phys,i}(\epsilon_i,K)=-A_i(K)\epsilon_i+O(\epsilon_i^2),
\qquad A_i(K)=\operatorname{Re}\!\left(\frac{a_{1,i}(K)}{a_2(K)}\right).
\]

The coefficients are computed from the exact collective closure at symmetric (s=\pm h), using the determinant derivative identity \(d\det C=\det(C)\operatorname{tr}(C^{-1}dC)\). No ExpG scalar-authority formula is assumed. The fitted exponent over \(10^{-8}\) to \(10^{-2}\) averages 1.01088, within 5% of one for every support. The tighter small-epsilon exponents are closer to one. See `tables/TABLE_H03_scaling_exponent.csv` and `tables/TABLE_H03_root_law_coefficient_validation.csv`.

The large coefficient-vs-state-root error at some supports is expected when the leading coefficient is small and finite-epsilon higher-order terms matter; those values remain explicit in the validation table rather than being used to retune the fit.

## Gain dependence

At \(\epsilon=0\), changing PLL gains does not remove the quotient physical zero: the six deterministic gain patterns retain the same quotient class. Gains instead change the retained-SG coefficient (A_i(K)). The mixed derivative matrix is computed analytically from determinant-trace derivatives, including the mixed port term, and checked independently by complex-step perturbations of the exact scalar Woodbury port closure. For entries above the predeclared near-zero cutoff, the relative-error median is \(1.28\times10^{-6}\) and the 95th percentile is \(1.78\times10^{-6}\), passing the fixed \(10^{-5}\) and \(10^{-4}\) targets. Values excluded by the cutoff are listed rather than assigned an unstable relative error.

The normalized mixed-authority Jacobian has numerical rank one at relative singular-value thresholds \(10^{-2},10^{-3},10^{-4}\). The root law is therefore linearly mixed in \(\epsilon\) and \(K\) through (A_i(K)), even though the all-GFL zero itself is gain-invariant.

## Scope and limit

This derivation establishes the local analytical root law and the gain dependence of its coefficient. It does not prove a globally optimal \((\rho,K)\) co-design. The local predictor materially misranks some supports, exact all-pole correction is required, the final controller gains remain effectively nominal, the frequency excursion diagnostic is active at the stated 100 MW pulse, and the independent PowerDynamics spectral margin does not agree with the analytical margin. These findings make the overall theory result partial.
