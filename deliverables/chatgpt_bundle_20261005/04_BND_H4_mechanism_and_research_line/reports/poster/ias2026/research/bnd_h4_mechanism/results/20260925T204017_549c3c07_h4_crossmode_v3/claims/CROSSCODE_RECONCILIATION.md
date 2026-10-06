# Cross-code reconciliation

Status: **BLOCKED_CROSSCODE_MISMATCH**

The active Python TX4 run passes G0--G11. The independent Julia reference
also solves its own frozen snapshot with small equilibrium residuals, but it
does not reproduce the active H4 contract:

| quantity | active Python TX4 | Julia reference |
|---|---:|---:|
| dynamic states | 86 | 82 |
| transverse alpha (s^-1) | 0.1270064678 | 0.1446702207 |
| critical frequency (Hz) | 0.6222796695 | 0.5746810861 |
| max algebraic residual | 9.38e-13 | 1.63e-13 |

The Julia output is retained under `raw/reconciliation/` as historical
cross-code evidence. It is not combined with the Python science claim, and
no tuning or model substitution was performed.
