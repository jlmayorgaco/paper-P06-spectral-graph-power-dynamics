# Phase E0 Inertia-Control Cross-Term Gate

This report tests whether the joint inertia/control perturbation is measurably non-additive.
The full ANDES eigensolve is the pole reference.  The Schur partition is used only to
report control participation of the tracked mode.

- Status: **NOT_MEASURABLE**
- Base case: `C:\Users\walla\Documents\Github\paper-P06-spectral-graph-power-dynamics\validation\ieee39_rational_filter\cases\phase1_calibrated\no_pss\ieee39_ibr_mix60.xlsx`
- Mode selector: `network`
- eps: `0.01`, half-step h: `0.005`
- Commit: `46f48b23d2ca168678e843f32d844eae5f4d93e5`

## Base Tracked Mode

- pole = -0.656641 + j4.88482, freq = 0.777443 Hz, zeta = 0.133226, pi_c = 0.5479, family = I_network

## Cross-Term Metrics

- |true joint shift| = `1.556843e-02`
- |estimated cross shift at eps| = `3.485490e-06`
- cross/true shift fraction = `0.000223882`
- additive complex error = `2.804225e-05`
- cross complex error = `3.152763e-05`
- complex-error improvement = `-0.124291`
- additive zeta error = `4.679823e-07`
- cross zeta error = `5.318994e-07`
- zeta-error improvement = `-0.13658`

## Interpretation

The cross term is not measurable under this registered local test.  This does not invalidate the attribution framework, but it means the inertia-control cross term should not be used as a headline empirical claim for this benchmark without a stronger case.