# Phase E0 Inertia-Control Cross-Term Gate

This report tests whether the joint inertia/control perturbation is measurably non-additive.
The full ANDES eigensolve is the pole reference.  The Schur partition is used only to
report control participation of the tracked mode.

- Status: **NOT_MEASURABLE**
- Base case: `C:\Users\walla\Documents\Github\paper-P06-spectral-graph-power-dynamics\validation\ieee39_rational_filter\cases\phase1_calibrated\no_pss\ieee39_ibr_mix60.xlsx`
- Mode selector: `critical`
- eps: `0.01`, half-step h: `0.005`
- Commit: `46f48b23d2ca168678e843f32d844eae5f4d93e5`

## Base Tracked Mode

- pole = -10.7514 + j86.3419, freq = 13.7417 Hz, zeta = 0.123566, pi_c = 0.8578, family = II_control

## Cross-Term Metrics

- |true joint shift| = `4.369660e-01`
- |estimated cross shift at eps| = `1.675444e-09`
- cross/true shift fraction = `3.83427e-09`
- additive complex error = `5.504618e-04`
- cross complex error = `5.504610e-04`
- complex-error improvement = `1.31993e-06`
- additive zeta error = `7.608837e-07`
- cross zeta error = `7.608997e-07`
- zeta-error improvement = `-2.10311e-05`

## Interpretation

The cross term is not measurable under this registered local test.  This does not invalidate the attribution framework, but it means the inertia-control cross term should not be used as a headline empirical claim for this benchmark without a stronger case.