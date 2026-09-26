# Freeze — `IAS2026_TRACKA_F8_F12_POST_F7_FREEZE`

Frozen 2026-09-10, after the post-F7 phases were accepted and before the final
journal gates. Separate from, and later than, the F7 freeze
(`IAS2026_TRACKA_F7_POLICY_HYPERGRAPH_FREEZE`, commit `f12ae3a0`) and the F12
preregistration (commit `34a64bdc`). Nothing produced under this tag may be
modified afterwards; the journal gates write to new directories.

## Contents

| phase | code | data | report |
|---|---|---|---|
| F8 service attribution | `experiments/F8_service_attribution.py`, `F8_report.py`, `F8_tables.py`, `F8D_avr_blend.py` | `results/F8/` | `docs/F8_SERVICE_ATTRIBUTION.md` |
| F8B mean vs heterogeneity | `experiments/F8B_mean_vs_heterogeneity.py` | `results/F8B/` | same |
| F8C service-induced closure | `experiments/F8C_service_closure.py` | `results/F8C/` | same |
| F10 baselines | `experiments/F10_baselines.py` | `results/F10/` | `docs/F10_F11_BASELINES_AND_SCALING.md` |
| F11 scaling | `experiments/F11_scaling.py`, `src/ibr_cycles/models/port_core.py` | `results/F11/` | same |
| F12 Kundur | `experiments/F12_import_kundur.py`, `F12_kundur.py`, `F12_outside_band.py`, `F12_full_rhp_sensitivity.py`, `F12_figure.py` | `results/F12/`, `configs/kundur/` | `docs/F12_SECOND_BENCHMARK.md` |
| regression of F7 after library changes | `experiments/F7_regression_after_freeze.py` | `results/F7_regression_after_freeze.json` (120/120) | — |
| final decision | — | — | `docs/TRANSACTIONS_FINAL_RESULT_AUDIT.md` |

Library changes since the F7 tag, all default-off (IEEE-39 results unchanged,
F7 labels 120/120 reproduced): service switches (`avr_blend`, `flux_blend`,
condenser damping and reactive share, converter synthetic inertia), the SEXS
lead-lag exciter state, per-network controller data. Tests:
`tests/test_service_interventions.py`, `tests/test_port_core.py`; 519 passing.
Manifests: `outputs/ias2026_post_f7/`. Ledgers: `docs/CLAIMS.md` rows O72–O81,
N23–N24; `docs/FAILED_EXPERIMENTS.md` F15–F17.

Integrity: `POST_F7_SHA256SUMS.txt` at the research root lists every file added
or changed by this commit.
