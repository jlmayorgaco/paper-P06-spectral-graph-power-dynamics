# TX3 E05D — weak-grid materiality challenge

This stage executes the final authorized TX3 stress experiment. It weakens the
34 frozen external AC transmission lines by scaling series R and X together,
calibrates eight endpoints with the empty coalition only, freezes each seed's
actual margin-setting oscillatory mode, and evaluates all 15 E05B-reproducible
pair/triple coalitions on the nine-point stress grid.

Run, in order:

```text
python experiments/tx3/E05D_weak_grid/preregister_e05d.py
python experiments/tx3/E05D_weak_grid/calibrate_e05d_stress.py --workers 8
python experiments/tx3/E05D_weak_grid/freeze_e05d_endpoints_modes.py
python experiments/tx3/E05D_weak_grid/run_e05d_campaign.py --workers 8
python experiments/tx3/E05D_weak_grid/finalize_e05d.py
python experiments/tx3/E05D_weak_grid/validate_e05d.py
python experiments/tx3/E05D_weak_grid/package_e05d.py
```

No E05E, E06, or new ParaEMT execution is part of this stage.
