# TX3 final manuscript

This directory contains the terminal IEEE Transactions manuscript and its two
independent reviewer audits.

## Build

From this directory, run pdflatex, bibtex, and two final pdflatex passes on
main.tex. The figures directory makes the manuscript source self-contained.
The final review PDF is main.pdf (10 pages, IEEEtran journal format).

## Terminal scientific status

- C1, C2, D1, D2, and conditional D3 are supported.
- Nominal modal materiality (N1/C3a) and load-conditioned materiality
  (N2/C3b) are rejected under their frozen gates.
- Weak-grid materiality (N3/C3c) is unresolved because E05D was baseline
  ineligible; no coalition outcome was inspected.
- Cycle/SCC localization, mechanism surgery, industrial mitigation, nonlinear
  consequence, and blind EMT transfer remain unassessed.
- E02C is benchmark infrastructure validation only.

The controlling wording is frozen in TX3_FINAL_CLAIM_FREEZE.md. The
descriptive post-mortem is in TX3_CRITICAL_MODE_BASELINE_AUDIT.md; it cannot
change any registered claim status.
