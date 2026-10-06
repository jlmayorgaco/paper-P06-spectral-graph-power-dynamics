# TX3 E03/E04 preregistration

Freeze date: 2026-08-27. Git SHA: `14dd51ce1528533843508b9f5af08284528bc728`. E02C aggregate model hash: `0f40a8b33e7b93e16a7d8a9f5d1b7ed9f065a64f52c074591890fae9cea6e072`.

The principal system is the accepted `ieee39_tx3_gfl3` case with GFLs at buses 36, 37, and 38 and the original 31–6 tap 0.97143. E02C is not reopened. No new ParaEMT result is permitted. E05 and later stages are outside scope.

Eight actions are frozen in `ACTION_SET_FREEZE.json`; all 28 pairs and 56 triples are evaluated at 24 preregistered operating points (16 discovery, 8 holdout). The primary functional is the baseline-relative real log-absolute descriptor determinant over 0.1–30 Hz with normalized logarithmic-frequency weight. The secondary 0.1–50 Hz band is sensitivity only.

The numerical pilot is restricted to OP00/OP01 and coalitions (A1,A2), (A7,A8), (A1,A2,A3), and (A1,A7,A8). It may select only the first frequency, contour, derivative-threshold, and hypercube rule satisfying the predeclared convergence rules. It may not alter actions, amplitudes, operating points, splits, bands, K=5, or the C1/C2 decision rules. Holdout remains unread until `NUMERICAL_METHOD_FREEZE.json` exists.

Direct externalities use independently re-equilibrated coalition vertices. Connected reconstruction uses sparse LU, exact active-column resolvent traces, Gauss action-cube integration, and total mixed derivatives of fully re-equilibrated operators. Numerical errors are combined conservatively by summation; no independence assumption is made.
