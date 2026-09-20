# TX4 Robustness ChatGPT Handoff

## Delivery state

Branch: `research/tx4-final-robustness-statistics`
Parent: `f64db0004026ceafdb08dd13b5e2ff59d6060742`
No push: YES

## Main result

The complete master table covers 20,590 conditions and all 16 target-family portfolios. H4 rows are exact reduced-DAE evaluations. Exact all-portfolio calibration covers 426 conditions. Proper-subset rows outside calibration are ExtraTrees predictions and carry `SURROGATE_CALIBRATED`.

## Probabilistic outputs

QMC H4 presence: 0.5779. MC H4 presence: 0.5814. QMC exact-H4 minimality: 0.0735. MC exact-H4 minimality: 0.0712. Interpret these as bounded-box engineering probabilities, not calibrated physical probabilities.

## Key limitations

The prior PowerDynamics same-model network gate stopped before full parity. Random Julia cross-code conditions were not executed. The TDS file uses a spectral trace proxy and is not a new nonlinear TDS run. No physical robust radius, EMT, current-limit, DC-link, protection, or hardware claim is made.

## Files

See `results/TX4_ROBUSTNESS_MASTER_CONDITIONS.parquet`, compressed CSV, summary tables, claim matrix, headline JSON, and the two ZIP bundles.

## Final case

Case C: the poster and six-page paper can be strengthened with scoped tiered robustness language; TPWRS-ready universal robustness is not supported.
