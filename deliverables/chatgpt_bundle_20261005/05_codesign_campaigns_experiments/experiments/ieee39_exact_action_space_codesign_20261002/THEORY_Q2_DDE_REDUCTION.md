# Q2 — exact low-rank characteristic reduction

Let the gauge-quotiented zero-delay matrix be

\[
A_c=A_0+BC,\qquad B=[b_1,\ldots,b_m],\quad C=[c_1^T;\ldots;c_m^T],
\]

with fixed delays and \(E_\tau(s)=\operatorname{diag}(e^{-s\tau_i})\). Because each PLL channel contributes \(b_i c_i^T e^{-s\tau_i}\),

\[
\Delta(s)=sI-A_0-BE_\tau(s)C
=sI-A_c-B(E_\tau(s)-I)C.
\]

At any \(s\) for which \(sI-A_c\) is nonsingular, the matrix determinant lemma gives

\[
\det\Delta(s)=\det(sI-A_c)\det\!\left[I_m-(E_\tau(s)-I)C(sI-A_c)^{-1}B\right].
\]

Thus the dynamic correction lives in an \(m\times m\) matrix (at most ten here), while the full finite state dimension is much larger. The factorization is an exact meromorphic identity. At poles of \((sI-A_c)^{-1}\), the product is interpreted by analytic continuation; one must not evaluate the two factors independently there or count uncancelled poles as DDE roots.

Numerical validation compares the characteristic operators before determinants and, away from reference poles, compares complex log-determinants. A delay-dependent root count additionally has to account for poles of the reference resolvent and certify the contour truncation. That is a separate gate; a small matrix alone is not a root-count certificate.

**Claim class:** exact identity for the fixed linearized DDE; numerical operator checks are in `TABLE_Q02_DDE_REDUCTION_VALIDATION.csv`.
