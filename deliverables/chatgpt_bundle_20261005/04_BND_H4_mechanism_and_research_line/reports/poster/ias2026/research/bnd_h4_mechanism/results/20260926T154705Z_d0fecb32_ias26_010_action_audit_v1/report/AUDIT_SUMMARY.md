# IAS26-010 action-space representation audit

- RUN_ID: `20260926T154705Z_d0fecb32_ias26_010_action_audit_v1`
- Exact critical value: `0.12700646777028055 +3.9098984764449018j`; matrices/eigenvectors are frozen in binary `.npy` files.
- Low-rank identity `TH = T0 + E D E^T`: **PASS**; absolute Frobenius residual `7.1401117121076822e-15`, relative `2.7201720726061947e-18`, maximum entry `7.1054273576010019e-15`.
- Action-space relative Frobenius residual: `2.3407240752625231e-11`; frozen gate `1e-12`: **FAIL**.
- Same-vector 80-digit relative residual: `2.34072191872e-11`.
- Action-to-terminal mapping backward residual: `1.5039737e-17` double; `5.9224257e-19` at 80 digits (scaled by `||TH||_F ||y||_2`).
- `cond2(T0)`: `3893.8353`; `cond2(I + D K)`: `3.12876504e+10`.
- Independent candidate-port/DAE-Schur mismatch: `1.3083706e-11` relative, same order as the action residual.
- Classification: **NUMERICAL_FLOOR_OF_FROZEN_DAE_PORT_REPRESENTATION**. The frozen low-rank/action-to-terminal algebra is consistent; raising arithmetic precision does not remove the residual. No gate, physics, or source code changed in this audit.
