# CDW hardening campaign: guide for reviewers

The campaign runs on branch `research/contextual-dynamic-weakness-hardening`, created from `e73dd355`. It was preregistered at commit `05b507e3` before any new numerics. TX4 is frozen at tag `TX4_FINAL_MANUSCRIPT_FREEZE` and is untouched.

## Read in this order

1. `docs/CDW_HARDENING_FINAL_REPORT.md` (also `.pdf`): verdict, phase status, the 22 final answers, claim matrix, negatives.
2. `paper/main.pdf` and `paper/supplement.pdf`: the IEEE TPWRS manuscript, revised after two internal reviews.
3. `docs/CDW_REVIEWER1_FINAL.md`, `docs/CDW_REVIEWER2_FINAL.md`, `docs/CDW_AUTHOR_RESPONSE_TO_INTERNAL_REVIEW.md`: the reviews and what changed.
4. `docs/CDW_HARDENING_PREREG_V1.md`, `docs/CDW_HARDENING_STATISTICAL_PLAN.md`, `docs/CDW_HARDENING_EXPERIMENT_MATRIX.md`: the frozen preregistration.
5. `docs/CDW_HARDENING_DEVIATIONS.md`: every deviation, including:
   - the post-hoc analyses;
   - the corrected timestamps;
   - the H18 defect;
   - the pinned ALT pole.
6. `results/CDW_HARDENED_CLAIM_MATRIX.csv`: allowed and prohibited wording per claim.

## Where the numbers come from

- **Gates and statistics:** `results/hardening/H*_gate.json`, `H05_stats.json`, `H19_evidence.json`.
- **Post-hoc reviewer analyses:** `H31_revision.json`, `H31_explore.json` and `H31_*.csv`.
- **Macros:**
  - The paper's numbers are macros in `paper/cdw_numbers.tex`, written by `experiments/cdw_hardening/CDWH_NUMBERS.py`.
  - `paper/tab_*.tex` and `paper/sec_corridors.tex` are written by `CDWH_TABLES.py`.
- **Cross-model data:** `results/hardening/alt/`.
  - ALT-WECC case files (`cases/`, 488).
  - The QFLAG=1 variant (`cases_qflag1/`).
  - Qualification (`H17_qualification_andes.json`).
- **Raw checkpoints:** per-task records under `raw/H_*` are not included because of their size. `results/hardening/CDW_HARDENING_RAW_MANIFEST.csv` lists each store with its file count, byte count and a digest of the sorted file hashes.

## Reproduce

Use `.venv/tx3-analysis`, with OPENBLAS/OMP/MKL threads set to 1.

1. Run the compute phases with `experiments/cdw_hardening/CDWH_MASTER_RUN.py`. It checkpoints each phase and resumes after an interruption.
2. Run the analyses `H03` to `H19` and `H31_*`.
3. Run the generators: `CDWH_NUMBERS.py`, `CDWH_TABLES.py`, `CDWH_FIGURES.py`, `CDWH_SUPPLEMENT.py`, `CDWH_MASTER_TABLES.py`, `CDWH_REPORT.py`.

The ANDES cross-model runs (`H17_andes_alt.py`) use `.venv/xtool-andes-gfl` with `USERPROFILE`/`HOME` set to `.venv/xtool-andes-gfl/home`.

## Scope limits

- Phasor-domain small-signal analysis on the IEEE 39-bus system with one custom GFL.
- No EMT claim.
- No GFM model.
- No Africano/PV result (source material missing).
