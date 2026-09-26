# PD39 statistical and reporting plan

This plan is paired with `PD39_ROBUST_TRANSITION_PREREG.md`.

## Estimands

The primary estimands are deterministic because the benchmark and portfolio set are finite:

1. nominal baseline non-gauge dynamic margin;
2. nominal single-replacement margin loss for each candidate bus;
3. robust margin of a portfolio, defined as the minimum `m_dyn` over the nine fixed controller scenarios;
4. minimum feasible intervention cardinality;
5. rank association between static and dynamic node scores.

No p-value is required for the exhaustive primary design. Uncertainty is represented by the preregistered controller box, not by post-hoc resampling.

## Numerical rules

- Preserve raw complex eigenvalues and all solver/error statuses.
- Use `gauge_tol=1e-8` and `margin_tol=1e-8` exactly.
- Do not round values before ranking or feasibility decisions.
- If a fixed-point or finiteness check fails, label the case `equilibrium_failed`; do not label it unstable.
- If eigenvalue computation fails after a qualified equilibrium, label the case `spectrum_failed`; do not infer a margin.
- Report counts of attempted, equilibrium-failed, spectrum-failed, unstable, nominal-safe, and robust-feasible cases.

## Ranking and ties

Static and dynamic rankings use the same bus universe for node comparisons. Ties are retained in the score table and handled by the standard Spearman tie correction. The code must report the score values, not only the rank coefficient. A rank reversal is descriptive evidence of context dependence; it is not a proof that the static proxy is invalid in all systems.

## Multiple comparisons and exploratory analysis

The 256 portfolios are an exhaustive design space, not 256 independent hypothesis tests. Exploratory pairwise or subgroup claims must be labeled exploratory and must not be presented as preregistered primary evidence. Any new uncertainty scenario, threshold, candidate, or metric is a deviation and is analyzed separately.

## Reporting template

Every table must expose:

```text
portfolio, replaced_buses, scenario, equilibrium_status,
dynamic_margin, max_real, stable, robust_feasible_for_scenario,
error_type, julia_version, package provenance
```

The final design row must include its full scenario vector and worst-case margin. If no robustly feasible row exists, the report must state that result directly and retain the complete attempted table.
