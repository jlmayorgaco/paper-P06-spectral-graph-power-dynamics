# TX3 critical-mode baseline audit

This is a descriptive post-mortem at `kappa_grid=1`. It uses the frozen E05D
holdout, the complete frozen E05B coalition population, and the original A1–A8
endpoints. It cannot modify any historical claim status or authorize E06.

```text
python experiments/tx3/critical_mode_baseline_audit/freeze_protocol.py
python experiments/tx3/critical_mode_baseline_audit/run_audit.py --workers 8
python experiments/tx3/critical_mode_baseline_audit/finalize_audit.py
python experiments/tx3/critical_mode_baseline_audit/validate_audit.py
```
