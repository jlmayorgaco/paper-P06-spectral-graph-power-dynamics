# M1A DAE-to-retained-operator reconciliation

This run diagnoses the strict M1 failure without changing the frozen model,
P4, mode, threshold, or previous run. It executes only R0--R7 and stops at
the first material degradation in the exact algebraic ladder.

The same in-memory H4 case and Jacobian blocks are used throughout:

```text
A_red = fx - fz solve(gz, gx)
T_raw(s) = gz + gx solve(sI-fx, fz)
```

The retained candidate-port and action-space routes are compared matrix by
matrix. No determinant-only conclusion is accepted. All solves use linear
systems; explicit inverse calls are forbidden in this diagnostic.

The original M1 absolute threshold remains `|h_k(lambda_c)| < 1e-8`. Relative
backward errors are recorded additionally but cannot promote M1 to PASS.

Outputs are limited to one immutable run root:

- `tables/M1A_ERROR_LADDER.csv`
- raw and detailed reconciliation JSON/CSV diagnostics
- a run report and claim-status JSON

No eta continuation, graph modes, damping budget, poster figures, tuning, or
modification of `H4_MODAL_SCHUR_V1` is allowed.
