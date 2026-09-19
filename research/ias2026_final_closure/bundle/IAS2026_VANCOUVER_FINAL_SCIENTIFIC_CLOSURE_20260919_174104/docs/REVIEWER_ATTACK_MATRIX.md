# Reviewer attack matrix

| Likely attack | Evidence | Safe response |
|---|---|---|
| Python and Julia implement different synchronous-machine models | Frozen source manifest and stopped gate report | Concede: same-model parity is not established |
| Stable alternative-model census is presented as a refutation | P2/P5 reports and raw CSVs | Concede scope: both are negative results for their named models only |
| TDS trace is used as modal proof | P4 trace and solver table | Do not claim modal-frequency agreement |
| Four stable holdout cases prove prediction | P6 metrics | Report single-class limitation; require second mixed holdout |
| Robustness is inferred from nominal data | robustness status file | Concede `NOT_TESTED`; no physical uncertainty claim |
| Scaling table is called a new asymptotic benchmark | scaling status file | Label the old result retrospective |
| Legacy PDFs are mistaken for final submission artifacts | ZIP exclusion manifest | The final bundle contains only the executive audit PDF |

The attack matrix is a scope-control document, not a substitute for missing
experiments.
