# D1 exact DDE spectrum gate

**Status: BLOCKED_EXACT_DDE_SPECTRUM.**

The delayed-error characteristic matrix is assembled without Padé. Its tau=0 limit reproduces the full 203-state finite ODE spectrum for the baseline and joint design (matched max errors 1.97e-11 and 3.76e-11; nonlinear-eigenvalue residuals below 3.4e-17).

For three frozen patterns, a bordered-Newton continuation tracked the dominant zero-delay conjugate pair through four delay steps. Each returned root has nonlinear-eigenvalue residual below 3.4e-18. These are local tracked-root results only; the CSV labels the result SUPPORTED_LOCAL.

The required rightmost root count was not completed. A norm-derived finite contour combined with a low-rank determinant-lemma evaluation was implemented, but the bound-driven contour exceeded the practical compute budget before a stable argument-principle count could be returned. No count, completeness, or exact rightmost-root claim is accepted. The DDE stability frontier and optimization therefore remain blocked; Padé was not substituted.

The dominant tracked pair shifts by only the values in TABLE_D01_LOCAL_DDE_ROOTS.csv across the tested patterns. Because other characteristic roots were not counted, this does not establish the full spectral abscissa or falsify the spatial-delay hypothesis.
