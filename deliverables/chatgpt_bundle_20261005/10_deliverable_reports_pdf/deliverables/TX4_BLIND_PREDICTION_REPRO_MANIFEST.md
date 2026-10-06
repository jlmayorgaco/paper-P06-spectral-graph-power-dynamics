# TX4 blind prediction reproducibility bundle

This bundle is tied to parent freeze
`69f200dfe1dfb74cc4f678ad25a0c3b1751d62a6` and the campaign branch
`research/tx4-blind-portfolio-prediction-final`. It contains the blind
predictor, post-reveal comparison, contextual-return audit, TDS wrapper,
figure generator, frozen preregistration, and the raw/key tables needed to
inspect the reported finite results.

The bundle is not a standalone environment image. It assumes the repository
dependencies and the frozen TX4 model code at the parent commit. No push was
made.

## Included source

- `reports/poster/ias2026/research/experiments/tx4_blind_predict.py`
- `reports/poster/ias2026/research/experiments/tx4_reveal_runtime.py`
- `reports/poster/ias2026/research/experiments/tx4_contextual_return.py`
- `reports/poster/ias2026/research/experiments/tx4_g_boundary_reduced.py`
- `reports/poster/ias2026/research/experiments/tx4_modal_mechanism.py`
- `reports/poster/ias2026/research/experiments/tx4_tds_final.py`
- `reports/poster/ias2026/research/experiments/tx4_postprocess.py`
- `reports/poster/ias2026/research/experiments/tx4_figures.py`

## Included protocol and audit documents

- `docs/TX4_BLIND_PREDICTION_CURRENT_STATE.md`
- `docs/TX4_BLIND_PREDICTION_PREREG.md`
- `docs/TX4_BLIND_PREDICTION_CLAIMS_V0.csv`
- `docs/TX4_BLIND_PREDICTION_EXECUTION_PLAN.md`
- `docs/TX4_BLIND_PREDICTION_DEVIATIONS.md`
- `docs/TX4_BLIND_FIREWALL.md`
- `docs/TX4_BLIND_PREDICTION_CHATGPT_HANDOFF.md`
- `docs/TX4_CONTEXTUAL_RETURN_THEOREM.md`

## Included key outputs and raw traces

- all `results/TX4_*` blind, reveal, runtime, contextual-return, modal,
  remediation, and claim tables;
- `results/TX4_TDS/` including the three compressed raw traces;
- `reports/poster/ias2026/research/results/PCV/PCV02/` frozen answer tables;
- `reports/poster/ias2026/research/results/FINAL_CLOSURE/` frozen closure
  references;
- `reports/poster/ias2026/research/results/G2/` frozen TDS summaries and
  supporting traces;
- the existing E31 ANDES validation CSV/manifest.

The curated upload ZIP is intentionally smaller and contains only the final
report, handoff, headline/claim files, contextual-return theorem, key CSVs,
and poster PNGs.
