# CDW review bundle — how to read this

Campaign: Contextual Dynamic Weakness in Inverter-Rich Power Networks.
Branch `research/contextual-dynamic-weakness`, from tag
`TX4_FINAL_MANUSCRIPT_FREEZE` (commit `69f200df`). Not pushed. TX4 was read
only and never modified.

## Read first

1. `docs/CDW_FULL_CAMPAIGN_REPORT.pdf` — the final report (start with its
   Executive Summary, section 1).
2. `docs/CDW_PREREG_V1.md` — everything was frozen here before any numerical
   result.
3. `results/CDW_MASTER_CLAIM_MATRIX.csv` — one row per claim, its status, and
   a direct link to the code/data/figure that supports it.

## Directory map

- `docs/` — theory, preregistration, experiment plan, claim ledgers, the
  Africano gate record, deviation log, and the final report (md/tex/pdf).
- `theory/` — the graph-DAE model and notation.
- `experiments/cdw/` — all code: model library (`_cdw.py`), sensitivity
  engine (`_sens.py`), analysis helpers (`_analysis.py`), one `E*.py` module
  per experiment, the orchestrator (`CDW_MASTER_RUN.py`), figures
  (`CDW_FIGURES.py`), tables (`CDW_MASTER_TABLES.py`), the determinism check
  (`CDW_DETERMINISM.py`), report assembly (`CDW_ASSEMBLE_REPORT.py`,
  `md2tex.py`), and unit tests (`tests/`).
- `results/` — every summary CSV/JSON/parquet, one file per experiment, plus
  the three master tables and the determinism result.
- `figures/` — F1–F13 (pdf/svg/png each); F14–F17 do not exist (Africano
  blocked).
- `logs/` — orchestrator and per-run logs.
- `raw/` — per-task JSON checkpoints (one file per case; this is what the
  raw-data bundle archives).

## Strongest positive result

**GOLD-B** (§10 of the report): total, re-equilibrated dynamic line
sensitivity predicts held-out finite branch-strengthening effects far better
than any static baseline (median Spearman 0.987 vs 0.50, on 32 held-out
conditions), cross-validated by an independent port-form derivative.

## Strongest negative result

**GOLD-C** (§18–19): no reduced model achieves a genuine end-to-end speedup
(all ≈ 1.00×) together with a materially smaller state count — the
bottleneck is the nonlinear equilibrium solve, not the eigenanalysis these
reductions target. The certificate itself is accurate (0 false
certifications) whenever it fires.

## Blocked phases

Africano/PV (Phases 18–21, GOLD-E/F): the material gate (E17) found no source
document, feeder data, metric definition or reported result to reproduce. See
`docs/CDW_AFRICANO_MISSING_INPUTS.md`. No generic feeder was substituted.

E22 (nonlinear recovery vs α_⊥) was deferred by design (optional, listed as a
next step in the report, §27).

## Exact reproduction commands

```
# environment
.venv/tx3-analysis/Scripts/python.exe   (numpy 2.5.2, scipy 1.18.1, pandas 3.0.5, cvxpy 1.9.2)

# unit tests
python -m pytest experiments/cdw/tests -q

# full campaign (resumes from experiments/cdw/../raw/<phase>/*.json if present)
python experiments/cdw/CDW_MASTER_RUN.py --max-workers 12

# one phase only
python experiments/cdw/CDW_MASTER_RUN.py --only E01 --max-workers 12

# determinism check (reruns headline deterministic cases into DET_* stores)
python experiments/cdw/CDW_DETERMINISM.py

# figures, master tables, report
python experiments/cdw/CDW_FIGURES.py
python experiments/cdw/CDW_MASTER_TABLES.py
python experiments/cdw/CDW_ASSEMBLE_REPORT.py
python experiments/cdw/md2tex.py docs/CDW_FULL_CAMPAIGN_REPORT.md docs/CDW_FULL_CAMPAIGN_REPORT.tex
lualatex docs/CDW_FULL_CAMPAIGN_REPORT.tex   # twice, from docs/
```

## Commit hashes (this campaign, on `research/contextual-dynamic-weakness`)

- `b8ae3082` — preregistration (theory, prereg, experiment plan, claim
  ledger, graph-DAE model), before any CDW numerical result.
- `3ddf2bc2` — orchestrator, model library, sensitivity engine, experiment
  modules, unit tests.
- `0b53b136` — E1/E1b/E10/E3/E4 results (GOLD-B).
- `66294b44` — E2/E6/E8/E12 results (GOLD-A completed), plus two
  infrastructure bug fixes (documented in
  `docs/CDW_PREREG_V1_DEVIATIONS.md`).
- `d73b09a2` — E9/E13/E14 results (GOLD-D, GOLD-C).
- `7c1c15aa` — E7/E16/E17 results.
- (final commit hash: see `git log` / the closing message of this session)
