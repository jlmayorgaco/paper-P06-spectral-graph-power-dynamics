# TX4: Policy-dependent minimal incompatibility of SG-to-inverter replacement portfolios

`main.pdf` is the 12-page IEEEtran journal manuscript. It is not the IAS poster
and not the TX3 paper; neither of those was modified.

## Build

`latexmk` needs Perl, which is not installed here. Build with:

```
pdflatex main && bibtex main && pdflatex main && pdflatex main
```

`p2_numbers.tex` holds the IEEE-39 P2 Monte Carlo numbers used in Table III and
Section VI.

## Figures

Every figure in `figures/` is produced by
`reports/poster/ias2026/research/experiments/paper_tx4/MC03_paper_figures.py` and
comes with a `*_source.csv` holding the plotted data.

## Provenance of every number

| paper element | source |
|---|---|
| Theorems 1–5, Props 1–4 | `research/theory/` notes (TRANSVERSE_STABILITY_QUOTIENT, PRINCIPAL_MINOR_PORTFOLIO_STRUCTURE, FINAL_COMBINATORIAL_THEOREMS, SYMMETRY_DEFLECTED_PORT_CLOSURE) |
| Table III, Fig. 1 | Monte Carlo run `outputs/ias2026/paper_mc_20260911T083722` (MC01 synthetic, MC02 IEEE-39, MC04 post-hoc); summaries copied to `research/results/paper_mc/` |
| preregistration | `research/configs/ias2026/paper_mc_validation_v1.yaml`, committed as 7a772808; amendment v1.1 as f2946257 |
| Tables IV–V, Figs. 4–6 | FC18 targeted port checks (`research/results/FC18/`, `docs/FC18_TARGETED_PORT_CHECKS.md`) |
| Fig. 3 and the F7 statements | FC01 transverse re-audit (`research/results/TSQ_ieee39_reaudit.csv`) |
| Table VI, Fig. 7 | FC03 governed replication |
| Figs. 8–9, Table VII | FC05/FC06 nonlinear thresholds, FC07 v2 curvature, FC13 Hopf |
| Table VIII | FC10 census, FC12 planning |
| Fig. 10 | FC02 zero-frequency port holdout |
| claim status and wording | `research/docs/FINAL_TRANSACTION_THEORY_AND_EVIDENCE.md`, `FINAL_TPWRS_CLAIMS.md`, `FINAL_REJECTED_WORDING.md` |

## Deviations disclosed in the paper (Section VI-C)

1. **P1 amendment.** The spectral-identity rule matched the numerically split
   Jordan pair. It was amended before the run, with the threshold unchanged.
2. **P2 failed as implemented.** The bisection stopped at the classifier band.
   The post-hoc re-bisection (MC04) is reported beside the failure.
3. **S7 generator defect.** Three coalescence-class draws contained a genuine
   zero crossing. In addition, 183 S7 draws were left undecided by the
   device-flip rule.
4. **P1 chart redraws.** Two draws needed a redrawn off-equilibrium state.
5. **FC07 v1.** The first-order-hold coefficient was wrong. The v1 run is
   archived; only v2 is reported.

## Before submission (author action)

- Confirm the affiliation and e-mail line.
- Add funding and acknowledgments, if any.
- Choose the repository or DOI to cite in the Reproducibility statement.
