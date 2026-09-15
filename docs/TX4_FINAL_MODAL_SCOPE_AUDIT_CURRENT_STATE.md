# TX4 final modal-scope audit — current state

## Provenance

- Parent freeze: `69f200dfe1dfb74cc4f678ad25a0c3b1751d62a6`
- Starting branch: `research/tx4-blind-portfolio-prediction-final`
- Starting commit: `2cbb860eb4ba053edee99c1a040c8024f3f6828e`
- New branch: `research/tx4-final-modal-scope-audit`
- Push: prohibited

This is a clean descendant. Earlier commits and earlier result files are not
rewritten.

## Audited frozen material

The raw boundary table reports, at the H4 controller boundary, a nearest
eigenvalue of `Q_H` near `-1` and a nearest eigenvalue of the contextual
return near `+1`. The previous theorem/handoff wording incorrectly described
both objects with the same sign. The correction is mathematical wording only;
the numerical values are not changed.

The previous `TX4_LOCAL_VS_COLLECTIVE.csv` table reports local singular values
near one because it was computed from diagonal blocks of the normalized
matrix `I+Q`. Those blocks are identity by construction after normalization.
That table is not evidence that the original physical factors `I+M_ii` are
regular. This audit recomputes the physical factors from `M=D K`.

The frozen reduced predictor searches the positive-frequency band
`0.3--1.5 Hz`. The frozen V9 comparison used the global full-order status,
which includes modes outside that band. This audit will report both targets:

```text
alpha_EM     = max Re(lambda) over transverse eigenvalues in 0.3--1.5 Hz
alpha_global = max Re(lambda) over the complete admissible transverse spectrum
```

Neither target will be substituted for the other.

## No-new-campaign boundary

This branch will only correct the return-sign convention, replace the invalid
local-factor evidence, reclassify the already frozen V9 spectra by modal
scope, optionally inspect at most three clean slow false-safe cases with the
same frozen TDS, and rebuild the final figures/report/package. It will not
run a parameter search, repair optimization, new benchmark, structured
radius, weak-node/link campaign, or new theory.
