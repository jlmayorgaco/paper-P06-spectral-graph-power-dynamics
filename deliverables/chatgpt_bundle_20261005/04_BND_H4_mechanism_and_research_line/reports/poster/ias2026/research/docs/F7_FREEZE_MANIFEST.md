# Freeze — `IAS2026_TRACKA_F7_POLICY_HYPERGRAPH_FREEZE`

Frozen 2026-09-10, after F7 and before any F8 work. **Nothing under
`results/F7/`, `figures/F7*`, `theory/F2*`, or the F7 sections of `docs/CLAIMS.md`
may be modified after this tag.** Later phases write to new directories
(`results/F8/`, ...) and add new ledger rows; a correction to a frozen number is a
new row that cites the frozen one, never an edit.

## Scope of the tag

Everything under `reports/poster/ias2026/research/`:

| content | location |
|---|---|
| theory F2, F2B, F2C, F2D | `theory/` |
| code: library, experiments E00–E41, F1, F1B, F2, F7 | `src/ibr_cycles/`, `experiments/` |
| tests (508 passing at the tag) | `tests/` |
| configs, preregistrations, protocols, environment lock | `configs/`, `configs/ias2026/ENVIRONMENT_FREEZE_tx3-analysis.txt` |
| run manifests and progress log | `outputs/ias2026/final_validation_overnight_20260910T003225/` |
| source data, all phases | `results/` (F1, F1B, F2, F7, tables, manifests, discovery log) |
| F7 point tables | `results/F7/F7{A,B,C}_points.csv.gz`, gzip `-9 -n`, lossless |
| validation scripts | `experiments/F7_safeguard_A_audit.py`, `F7_safeguard_A_maps.py`, `F7_leak_sensitivity.py`, `F7_leak_qualitative.py`, the spot-check inside `F7_policy_hypergraph.py` |
| figures and their exact source data | `figures/` (F1B, F7) + the CSVs they are drawn from in `results/F1B_*`, `results/F7/` |
| ledgers and logs | `docs/CLAIMS.md`, `docs/FAILED_EXPERIMENTS.md`, `results/discovery_log.csv`, `results/F7_run.log` |

Unrelated working-tree changes elsewhere in the repository (TX3 manuscript,
thesis, poster sources) are deliberately **not** in this commit.

## Integrity

- `results/F7/SHA256SUMS_points.txt`: SHA-256 of the three **uncompressed** point
  tables as produced by the run, and of their `.csv.gz` files.
  `gzip -dc F7A_points.csv.gz | sha256sum` reproduces the first line.
- `FREEZE_SHA256SUMS.txt` at the research root: SHA-256 of every file in the tag.
- Reproducibility check performed before tagging: with the uncompressed tables
  removed, `experiments/F7_report.py` regenerated `results/F7/F7_summary.json`
  and `F7_report.txt` **identically** from the `.csv.gz` files.

## Environment

Python 3.13.14, numpy 2.5.2, scipy 1.18.1, pandas 3.0.5, matplotlib 3.11.1,
Windows 11; full list in `configs/ias2026/ENVIRONMENT_FREEZE_tx3-analysis.txt`.
Interpreter: `.venv/tx3-analysis/Scripts/python.exe` at the repository root.

## Reproducing F7

From `reports/poster/ias2026/research/experiments/`:

    python F7_safeguard_A_audit.py          # Safeguard A audit        (~1 min)
    python F7_policy_hypergraph.py          # three maps + spot-check  (~53 min, 16 workers)
    python F7_tongue.py                     # Safeguard D              (~10 min)
    python F7_report.py                     # tables, summary          (~1 min)
    python F7_figures.py                    # figures                  (~1 min)
    python F7_safeguard_A_maps.py           # direct structure checks  (~5 min)
    python F7_leak_sensitivity.py           # leak label agreement     (~5 min)
    cd .. && python experiments/F7_leak_qualitative.py

The adaptive refinement is deterministic; the spot-check and leak samples are
seeded (`20260911`, `20260910`).

## Frozen headline numbers (for later audit against drift)

| quantity | value |
|---|---|
| distinct hypergraphs over F7A/B/C | 36 (30 / 33 / 16) |
| located subset crossings | 29 979 (11 586 / 10 928 / 7 465) |
| imaginary-axis crossings port-detected | 29 851 of 29 851 |
| band-edge re-classifications | 128 (F7B only) |
| spot-check errors in homogeneous leaves | 0 of 1 153 |
| tongue tip `(g, k)` | `(0.10257, 1.24954)`, 0.681 Hz |
| gap-closure fold `(g, k)` | `(0.04092, 1.30710)`, 0.622 Hz |
| F2 regression on the F7A `g = 0` edge | `k* = 1.1488` (4 -> 3), `1.8168` (3 -> 2) |
