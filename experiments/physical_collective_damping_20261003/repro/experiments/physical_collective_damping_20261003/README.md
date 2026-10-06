# Physical collective damping on IEEE-39

Start with **REPORT_ES.md** for the executed result, **STATUS.md** for claim
limits, **THEORY.md** for the physical derivation, and **PROTOCOL.md** for the
preregistered primary test. **NUMERICAL_NOTES.md** records rejected numerical
implementations and corrections. Existing frozen experiments were not modified.

This experiment tests a real Beyond Nodal Damping mechanism at fixed replacement.
It does not calculate a new maximum SG-to-GFL replacement or a global optimum.
The graph operator is frequency dependent; the nonlinear simulator uses pure
delay on the PLL measurement error, not Pade or an added first-order filter.

The bundle contains the code, minimal original repository inputs, matrices,
complete CSV tables, PNG/PDF figures, report, exact source hashes, Project and
Manifest files. Extract it and follow REPRODUCE.md. No commit or push was made.
