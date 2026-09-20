# Reviewer attack matrix

| Likely attack | Evidence | Safe response |
|---|---|---|
| Python and Julia implement different synchronous-machine models | Frozen source manifest, full census reconciliation, and mode MAC table | State the exact frozen custom model and report the numerical tolerances |
| Stable alternative-model census is presented as a refutation | P2/P5 reports and raw CSVs | Concede scope: both are negative results for their named models only |
| TDS trace is used as modal proof | P4 trace and solver table | Do not claim modal-frequency agreement |
| No mixed second-model holdout exists | Second-model policy search and negative holdout status | Report the negative search; do not fabricate a mixed holdout |
| Robustness is inferred from nominal data | robustness status file | Concede `NOT_TESTED`; no physical uncertainty claim |
| Scaling table is called a new asymptotic benchmark | scaling status file | Label the old result retrospective |
| Legacy PDFs are mistaken for final submission artifacts | ZIP exclusion manifest | The final bundle contains only the executive audit PDF |

The attack matrix is a scope-control document, not a substitute for missing
experiments.
