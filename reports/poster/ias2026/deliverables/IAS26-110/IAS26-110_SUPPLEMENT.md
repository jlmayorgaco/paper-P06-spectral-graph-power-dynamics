# IAS26-110 poster package and reproducibility supplement

Status: **frozen poster package; IAS26-120 audited with explicit limitations**.  
Poster source: `main.tex` in this directory.  
Frozen evidence run: `20260927T144629Z_d0fecb32_ias26_060_operating_v1`.

Final poster: `IAS2026_IAS26-110_FINAL.pdf` (one page, 36 × 48 in). SHA-256: `9F0A4957162BA96F627A19ADEE714347006A8C68BE82853D8B1A39B6EC03930C`. The final source was rebuilt twice from this directory; the generated PDF hash matches this deliverable exactly. Its 30-dpi page render is pixel-identical to the prior reviewed export. The rendered page was visually inspected; no TeX overfull/underfull warnings were present. The final source SHA-256 is `7EABD3E71E2F9861C192508134A633982881C9B68E08CA26BAAF6AABAD613F91`.

## Versioned visual mapping

This assembly follows the versioned mapping `config/IAS2026_POSTER_ARCHITECTURE_FINAL_V1.json` in the frozen evidence run. It changes only the poster mapping, not scientific protocols, models, thresholds, or historical results:

| Poster panel | Evidence | Source run / table |
|---|---|---|
| P1 - Stability cliff | Exact 16-portfolio H4 lattice, full transverse `alpha_perp` | IAS26-020 `20260926T161116Z_d0fecb32_ias26_020_f1_audit_v1`; `tables/F1_FINAL_AUDIT.csv`; final P1 PDF |
| P2 - Local/collective closure | Physical pre-normalized `I+M_ii`, collective `I+Q_H`, common boundary audit | IAS26-030 `20260926T162252Z_d0fecb32_ias26_030_f2_closure_v1`; `tables/F2_BOUNDARY_PORT_AUDIT.csv`; final P2 PDF |
| P3 - Operating-point ensemble | 1,000 frozen IDs, paired V4 transverse spectra, valid/infeasible status, Wilson intervals | IAS26-060; `tables/MC_EVENT_RATES.csv`, `derived/SCENARIO_METRICS.csv`, `derived/MC_CASES.csv`; final P3 PDF |
| P4 - Fixed intervention and TDS | Paired fixed-`g` results and the three completed canonical phasor-DAE traces | IAS26-060/080; `tables/IAS26-080_TDS_RESULTS.csv`, `raw/tds/traces/`, final P4 PDF |
| Evidence strip | Exact GFL11 Python/Julia verdict parity, same model | `tables/IAS26-FINAL_EVIDENCE_STRIP.csv`; final evidence-strip PDF |

F3 is omitted because IAS26-010 remains `BLOCKED_M1_STRICT`; F4 is omitted because the GFL/GFM lattice was not run. The four principal panels are the four pre-existing audited outputs P1-P4; the strip is supporting evidence, not a fifth headline figure. No new scientific solve or figure-data calculation is part of this layout step.

## V4 atlas supplement

The descriptive, read-only exact-mask analysis is recorded in [`IAS26-060_V4_BLOCKER_ATLAS_DESCRIPTIVE.md`](../../research/bnd_h4_mechanism/docs/IAS26-060_V4_BLOCKER_ATLAS_DESCRIPTIVE.md). It recomputed all inclusion-minimal unstable masks from the 16 original-treatment `alpha_perp` rows for each of 967 valid IDs and matched the persisted classifications exactly (0/967 mismatches). The exhaustive partition is:

| Exact minimal-blocker set | IDs |
|---|---:|
| None (all 16 V4 portfolios stable) | 412 |
| `{30,33,35}` only | 77 |
| H4 only | 190 |
| `{30,33,35}` and `{30,33,37}` | 151 |
| `{30,33,35}`, `{30,33,37}`, and `{33,35,37}` | 137 |

Across the 365 alternative-witness scenarios, H4 itself remains unstable but loses minimality to one or more triples. This is evidence for witness switching within V4 only. The collective `I+Q_B` mechanism was not re-evaluated for those changing witnesses.

## Claim-to-evidence ledger

| Claim | Status and safe wording | Run/table row or figure | Limitation carried on poster |
|---|---|---|---|
| C1 | `NUMERICALLY_VALIDATED`, nominal H4 is a minimal transverse blocker; all 15 proper subsets stable | IAS26-020 `F1_FINAL_AUDIT.csv`, H4 and proper-subset rows; P1 | Canonical nominal witness only; not universal across operating points |
| C2 | `NUMERICALLY_VALIDATED`, physical local factors are nonsingular while collective closure is near singularity at the audited boundary | IAS26-030 `F2_BOUNDARY_PORT_AUDIT.csv`, common H4 boundary rows; P2 | Nonsingularity is not positive local damping; M1 strict reduction is blocked |
| C3 | `EMPIRICAL_HOLDOUT`/descriptive synthetic ensemble, H4 minimal in 190/967 valid IDs; any V4 collective blocker in 555/967 | IAS26-060 `MC_EVENT_RATES.csv`, E1/E2; P3 | 33 infeasible IDs retained, no replacements; not real-world probabilities |
| C3a | `EXPLORATORY_DESCRIPTIVE`, alternative minimal triples in 365 valid IDs; exact V4 partition above | IAS26-060 `MC_CASES.csv` + `SCENARIO_METRICS.csv`; V4 atlas memo | Post hoc description, V4 only; no invariant-mechanism claim |
| C4 | `EMPIRICAL_HOLDOUT`/descriptive synthetic ensemble, fixed `g=0.25` improves 594/967, rescues 186/967, deteriorates 373/967 | IAS26-060 `MC_EVENT_RATES.csv`, E3/E4/E5; P4 | Fixed candidate only; not a universal cure or optimum |
| C5 | `NUMERICALLY_VALIDATED`, Python/Julia verdict agreement 100/100 for exact GFL11 same-model reproduction | IAS26-060 `CROSSCODE_SUMMARY.json`; evidence strip | Mode identity unresolved; not independent-model or field validation |
| C6 | `CONDITIONAL`, TDS primary-extractor signs agree in 56/60 completed finite fits | IAS26-080 `IAS26-080_TDS_RESULTS.csv`; canonical traces in P4 | Four disagreements are at 0 Hz with unavailable `R2`; 30/90 solver failures retained; mode unresolved; phasor DAE, not EMT |

## Build command

From this directory, with the repository's TeX distribution and the frozen run outputs present. Build intermediates stay under the ignored repository build tree, not beside deliverables:

```powershell
New-Item -ItemType Directory -Force -Path ../../build/IAS26-110/tmp/pdfs | Out-Null
xelatex -interaction=nonstopmode -halt-on-error -output-directory ../../build/IAS26-110/tmp/pdfs main.tex
xelatex -interaction=nonstopmode -halt-on-error -output-directory ../../build/IAS26-110/tmp/pdfs main.tex
```

The expected output is `../../build/IAS26-110/tmp/pdfs/main.pdf`. The two passes resolve page labels and outlines. `latexmk` is unavailable on the current MiKTeX installation because its Perl script engine is missing; direct `xelatex` is the verified build path. Auxiliary TeX outputs are build intermediates and should remain under the ignored build tree or be removed after QA.
