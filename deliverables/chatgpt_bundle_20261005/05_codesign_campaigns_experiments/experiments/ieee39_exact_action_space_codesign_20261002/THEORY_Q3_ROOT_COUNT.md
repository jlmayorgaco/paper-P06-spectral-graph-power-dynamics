# Q3 — finite contour for the delay margin

For the Q2 factorization, set \(A_c=A_0+BC\), using \(C\in\mathbb R^{m\times n}\) in this derivation. The low-rank factor is

\[
R_\tau(s)=(E_\tau(s)-I)C(sI-A_c)^{-1}B.
\]

On the half-plane \(\Re s\ge\gamma\),

\[
\|E_\tau(s)-I\|_2\le \max_i(1+e^{-\gamma\tau_i})=:d_\gamma.
\]

For \(r=|s|>\|A_c\|_2\), the Neumann resolvent bound and resolvent identity give

\[
\|C(sI-A_c)^{-1}B\|_2
\le \frac{\|CB\|_2}{r}
+\frac{\|CA_c\|_2\|B\|_2}{r(r-\|A_c\|_2)}.
\]

Choose \(r_\star>\|A_c\|_2\) so that the right side times \(d_\gamma\) is strictly below one. Then \(sI-A_c\) is nonsingular and \(I-R_\tau(s)\) is nonsingular whenever \(\Re s\ge\gamma\) and \(|s|\ge r_\star\). Thus all roots in the tested half-plane lie inside a finite contour with imaginary and positive-real extent exceeding \(r_\star\). This derives a frequency extent from the model; it is not a selected cutoff.

The contour winding is evaluated in Float64 using a Schur factorization and the low-rank determinant. The norm inequality is analytic, but the current arithmetic is not interval/rounding-directed, so this implementation can at best return a numerical contour count, never a rigorous certificate. A zero count still requires phase-resolution refinement and independent root continuation before it can be labeled numerically safe.
