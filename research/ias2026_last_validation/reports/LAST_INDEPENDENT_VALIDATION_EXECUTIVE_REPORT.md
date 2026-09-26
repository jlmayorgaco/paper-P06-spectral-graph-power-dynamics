# IAS2026 last independent validation

## Verdict

The exact frozen custom IEEE-39 SG/AVR/PSS/load/network model is now implemented
in Julia and reconciled against the canonical Python implementation. The base,
the all-target V4 portfolio, and all 16 subsets pass the cross-code gate. The
all-target 30+33+35+37 portfolio is transversely unstable with
`alpha=0.1446702204 s^-1` at `0.574681 Hz`; its 15 proper subsets are stable.

The Julia terminal-operator audit supports a collective mechanism: the
four-target update is represented in an 8-dimensional action space, determinant
factorization residual is `7.22e-16`, maximum solve residual is `1.91e-17`, and
the minimum local factor singular value is `0.294`. The closest collective Q
eigenvalue on the 0.3–1.5 Hz grid has distance `0.0881` from `-1`.

Fresh nonlinear Julia traces solve the base, proper, repaired, physical-pulse
blocker, and mode-seeded blocker cases with algebraic residual at or below
`1e-9`. The common short pulse does not visibly excite the RHP mode, so the
nonlinear trace is recorded with an observability limitation rather than used
as a second instability proof.

The corrected official SimpleGFLDC scope is matched scheduled P/Q on the
common 100-MVA system base. Rating equivalence is not established. Its full
16-case V4 census is stable, and its representative PLL/current bandwidth
search (0.125, 0.5, 1, 2, 8 multipliers on proper and flagship portfolios)
found 20/20 stable cases and no mixed policy region. Consequently a >=12-point
mixed-class holdout was not honestly constructible and is recorded as a
negative result rather than fabricated.

## Scope

This bundle executes L1–L10 from the external closure audit as far as the
observed model evidence permits. It does not start robustness, scaling, paper
rewriting, or poster rewriting. The claim ledger is the controlling statement
of what may be promoted.

## Reproduction entry points

- `code/julia/run_true_same_model_ieee39.jl`
- `code/julia/run_true_same_model_mechanism.jl`
- `code/julia/run_true_same_model_tds.jl`
- `code/julia/run_p5_simplegfldc_ieee39.jl`
- `code/julia/run_second_model_discovery.jl`
- `code/python/run_true_same_model_oracle.py`
- `code/python/reconcile_true_same_model.py`
- `code/python/run_true_same_model_gate.py`
- `code/python/freeze_mixed_holdout.py`

The package excludes the prohibited generated `FINAL_PAPER.pdf` and
`FINAL_POSTER_DRAFT.pdf` artifacts.
