# F4 — Steady-frequency structure

Tested 200 deterministic mixed/support designs. For each, `H0=-C A⁻¹B+D` was computed from the exact reduced model at all generator-bus frequency outputs. Four cases had a nonfinite steady response; the remaining 186 use 100 training cases and 86 holdout cases for the scalar stiffness fit `Keff=-1/mean(H0)`. Input gains were sampled over the existing frozen bounds. This is a diagnostic affine fit, not a global certificate.

- Classification: **NONAFFINE**.
- Holdout max relative error: 1.7802360964952553; max absolute error 3298.2702948913593 MW/Hz.
- Fitted Kload: 302.1132202618006 MW/Hz.
- Fitted c_i: 283.78342455317465, 173.7689317164427, 339.4620809517033, 167.13284674368606, 164.71522992742157, 122.9744698936073, 241.8750949384369, 157.65035947471702, 550.9132906323591, -23.930601024850553.
