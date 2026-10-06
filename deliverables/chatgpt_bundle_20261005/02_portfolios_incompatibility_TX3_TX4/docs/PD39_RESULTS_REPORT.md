# PD39 executed experiments report

Execution date: 2026-09-13. Environment: Julia 1.11.9, PowerDynamics 5.0.0, isolated worktree `C:\tmp\pd39`.

## Coverage

The primary campaign completed all 256 candidate portfolios across all 9 preregistered controller scenarios: 2,304/2,304 rows were generated, all with `equilibrium_status=ok`, and no portfolio/scenario row was dropped. The 46 one-at-a-time branch-removal diagnostics also completed. Of those, 35 had a qualified equilibrium and 11 were retained as `equilibrium_failed`.

## Primary design result

The preregistered lexicographic objective selected the empty portfolio `none` with robust margin `0.09806540933624461 s^-1`. This is a valid null-design result but not a useful nonzero transition plan: the preregistration did not impose a minimum IBR penetration. It must not be described as evidence that the transition is solved.

For a descriptive nonempty frontier, the best one-bus portfolio is bus 37 with robust margin `0.1379049049215047 s^-1`. All 8 one-bus portfolios are robustly feasible at the preregistered `0.05 s^-1` threshold. The full 8-bus replacement is not robustly feasible; its worst-case margin is negative in the high-PLL scenarios. This nonempty comparison is secondary/descriptive because the primary optimizer was allowed to choose the null portfolio.

The cardinality frontier is recorded in `results/pd39/summary/portfolio_by_cardinality.csv`. It contains robustly feasible portfolios at cardinalities 0 through 7, while the all-eight portfolio has zero robustly feasible designs.

## Weakness results

The single-replacement dynamic score is margin loss `m_baseline - m_i`. Its Spearman association with the preregistered static node proxy is `rho=0.11904761904761904` over the 8 candidates. The largest positive losses are buses 34, 36, and 32; bus 37 has the largest margin improvement. This is context-specific evidence that the static proxy and dynamic response are weakly aligned in this benchmark, not a universal weakness theorem.

For branch removal, the static-vs-dynamic Spearman association over the 35 qualified outages is `rho=-0.050577936552638356`. The remaining 11 branch removals are not estimable under the no-redispatch/no-topology-repair rule because their equilibrium qualification failed. The largest qualified dynamic margin losses are branch 42 (26→27, `0.0081807976586932 s^-1`) and branch 6 (3→4, `0.005476589299449219 s^-1`).

## Interpretation boundary

The campaign supports a narrow benchmark statement: under the frozen PowerDynamics IEEE-39 model, frozen `SimpleGFLDC` controller, fixed candidate set, and fixed ±20% controller-factor box, dynamic margins depend strongly on which buses are replaced and are not well predicted by the static proxies used here. It does not support a universal weak-bus/link index, a universal causal mechanism, or a nonzero minimum-transition design without adding a preregistered penetration requirement.

## Reproducible artifacts

- Full portfolio/scenario table: `results/pd39/portfolio_campaign/portfolio_scenario_results.csv`
- Full branch diagnostic: `results/pd39/link_outage/link_outage_results.csv`
- Summary metrics: `results/pd39/summary/summary.toml`
- Cardinality frontier: `results/pd39/summary/portfolio_by_cardinality.csv`
- Dynamic node scores: `results/pd39/summary/dynamic_node_scores.csv`
- Qualified dynamic link scores: `results/pd39/summary/dynamic_link_scores_qualified.csv`
- Claims matrix: `docs/PD39_ROBUST_TRANSITION_CLAIMS.csv`
