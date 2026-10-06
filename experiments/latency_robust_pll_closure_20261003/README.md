# Fixed-replacement latency-robust PLL co-design closure

This folder continues the frozen experiment in
`../latency_robust_pll_codesign_20261003/`. It does not change the ten
SG→GFL replacement fractions, the five `±100 MW` events, gain bounds, or
the `−0.05 s⁻¹` spectral margin. All generated data live here; the parent
experiment is read-only.

The run begins by reconciling the parent Z, N, and old-best results in
`F00_PARENT_REPRODUCTION.csv`, then validates the stored 41 ms predictor in
`F01_PENDING_CANDIDATE_VALIDATION.csv`. Each subsequent design has a TOML
file under `designs/`, exact exponential-characteristic results under
`evaluations/`, and zero-delay nonlinear event results under
`event_validations/`. A candidate is accepted only after root coverage,
contour bracketing, and all five frozen events pass.

`TABLE_F8_POSITIVE_DELAY_VALIDATION.csv` and matching `TIME_DOMAIN_V2_*.csv`
are **linear exact-DDE method-of-steps evidence**. They do not constitute
nonlinear delayed event validation. The positive-delay nonlinear model remains
an open validation gate. `TABLE_F9_UNCERTAINTY_AUDIT.csv` separates actuator
limit failures from solver noncompletion.

The local optimizer uses the exact simple-root latency derivative and a
finite-difference actuator sensitivity measured on the nonlinear bus-16
`+100 MW` trajectory. Its LP KKT multipliers apply to that local predictor,
not to a globally optimal full-model design. See `THEORY_CLOSURE.md`.

`tau_crit` is the first uniform delay at which a tracked full-characteristic
root reaches the declared security line `Re(s) = -0.05 s^-1`; it is not a
certified global stability threshold. The first crossing is determined by
numerical contour/root coverage rather than interval arithmetic. The largest
fully checked design is listed in `BEST_VALIDATED_DESIGN_ID.txt`, and the
evidence, limitations, and exact final numbers are in `FINAL_REPORT.md` and
`POSTER_CLAIM_LEDGER.md`.

`TABLE_F5_MULTISTART.csv`, `TABLE_F6_GENERIC_OPTIMIZER_COMPARISON.csv`, and
`TABLE_F7_PARETO_FRONTIER.csv` are explicit NOT_EXECUTED records. The figures
with “Pareto” in their filenames plot observed designs only. See
`VALIDATION_GATES.csv` and `REPRODUCE.md` before using the plots in a poster.
`TABLE_F13_EVENT_SIZE_HEADROOM.csv` and `FIG_F7_EVENT_SIZE_HEADROOM.png` record
the narrow one-event actuator bracket measured at the accepted step-13 design.

To regenerate tables and figures after all candidate runs complete:

```powershell
python experiments/latency_robust_pll_closure_20261003/audit_parent.py
python experiments/latency_robust_pll_closure_20261003/finalize_closure.py
```

Julia scripts use the repository's existing `--project=.` environment.
