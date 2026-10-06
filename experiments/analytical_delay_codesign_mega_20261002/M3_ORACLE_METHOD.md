# M3: finite-contour DDE root counting

The frozen PLL delay enters only the measured error rows of the exact 203-state linearized DDE. Write `Ac=A0+BC` for the zero-delay Jacobian and `E_tau(s)=diag(exp(-s*tau_i))`. The characteristic determinant is

`det Delta(s)=det(sI-Ac) det F(s)`, where

`F(s)=I-(E_tau(s)-I)C(sI-Ac)^(-1)B`.

The first determinant has zero roots to the right of the declared margin `gamma=-0.05 s^-1` for both tested zero-delay designs. No Padé approximation is used. The second determinant is 10 by 10. Its logarithmic derivative uses the analytic derivative of `F` and is integrated around a rectangle. An independent discrete determinant-phase winding is computed from the adaptively evaluated contour points. A diagonal similarity transformation balances `Ac`; it leaves the characteristic roots unchanged.

The outer contour radius follows a model-specific norm inequality. For `|s|=r>||Ac||_2` and `Re(s)>=gamma`, the identity `C(sI-Ac)^(-1)B=CB/s+CAc(sI-Ac)^(-1)B/s` gives

`||(E_tau-I)C(sI-Ac)^(-1)B||_2 <= max_i(1+exp(-gamma*tau_i)) [||CB||_2/r + ||CAc||_2 ||B||_2/{r(r-||Ac||_2)}]`.

The chosen radius makes this bound less than one, excluding further characteristic zeros in the exterior right-hand region in exact arithmetic. Matrix norms and contour integrals here are evaluated in Float64; the inequality is analytically valid but no outward-rounded proof of the numeric radius or winding has been produced. The adaptive Simpson error is an estimate, **not** a rigorous quadrature bound. It is labeled `quadrature_error_estimate` in the tables.

The contour phase and trace integral agree near integers for the executed 20, 30 and 40 ms cases. At 40 ms, exact-characteristic nonlinear-root refinement found three distinct positive-imaginary roots for the 87.5% seed and six for the 88.455% zero-delay design. Their conjugates match the contour counts of 6 and 12 respectively. The small-factor singular values at the refined roots are at most `4.72e-10`; the corresponding full characteristic normalized minimum singular values are near `4.5e-16`. These are numerical root validations, not a proof covering arbitrary parameter values or delayed nonlinear events.

An earlier phase-only implementation in Experiment Q reported unverified winding indices 64 and 62 at 20 and 40 ms. The present independently evaluated trace integral, adaptive phase calculation, and 40 ms root refinement reject those old indices. The old outputs remain untouched and are not used as root counts.

This M3 gate is **NUMERICALLY_VALIDATED for the six tested uniform-delay/design combinations only**. It is not a directed-rounding certificate and does not validate a maximum replacement frontier. There is no nonlinear method-of-steps event simulation in this experiment.
