# Exact GFL11 cross-code reconciliation

Status: **PASS_GFL11_EXACT_REPRODUCTION**

The correct Julia runner was `code/tx4/run_tx4_exact_p4_ieee39.jl`, not the
historical 82-state GFL10 runner. For H4 it reports 86 states, matching the
active Python contract:

| quantity | Python | Julia GFL11 |
|---|---:|---:|
| states | 86 | 86 |
| transverse alpha (s^-1) | 0.1270064678283 | 0.1270064683223 |
| frequency (Hz) | 0.6222796695037 | 0.6222796695187 |

The exact runner executed the complete 16-portfolio census. The old 82-state
output was not used as a gate.
