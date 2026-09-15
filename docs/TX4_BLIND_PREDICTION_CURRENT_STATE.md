# TX4 blind-prediction current state

This is the read-only archaeological record immediately before the blind
prediction campaign. It is not a new numerical result.

## Repository state

- Branch: `research/tx4-blind-portfolio-prediction-final`
- Parent/tag: `69f200dfe1dfb74cc4f678ad25a0c3b1751d62a6`
  (`TX4_FINAL_MANUSCRIPT_FREEZE`)
- Clean worktree: `C:\\w\\tx4blind`
- No push.

The main user checkout is a separate dirty branch. This campaign does not
modify it, and it does not merge CDW, PD39, or ParaEMT research branches.

## Frozen candidate sets and model

V4 is the four candidate set `{30,33,35,37}`. V9 is the existing nine-
candidate set `{30,31,32,33,34,35,36,37,38}` from the frozen TX4 census. The
model is the full-order IEEE-39 phasor-domain DAE with the frozen fourth-order
SGs, first-order AVR policy, no governor in the transverse reference model,
and the ten-state custom GFL model used by TX4. Equilibria are re-solved for
each direct full-order case.

The flagship reveal retained from the frozen record is

```
H4 = {30,33,35,37}
P4 = (g=0.03625, k=1.425, t=1.5, h=1)
alpha_perp approximately +0.127006467 s^-1
critical frequency approximately 0.62228 Hz
```

All 15 proper H4 subsets are historically recorded as stable. The historical
clean g-only path fixes `(k,t,h)=(1.425,1.5,1)` and places the full-order root
near `g=0.2076814045`. These values are an answer-key reveal, not inputs to
the blind predictor.

## Exact decomposition

From the all-SG baseline, the frozen port implementation constructs

```
T0(s), K(s)=E^T T0(s)^(-1) E, Delta Y_i(s),
D_S(s)=blkdiag(Delta Y_i), M_S=D_S K_SS,
Q_S=(I+D_S K_d,SS)^(-1) D_S K_o,SS.
```

The exact determinant factors are checked with the audited fixed rectangular
port convention. Local changes are device-port models only; no full-order
portfolio spectrum, alpha label, boundary table, F4/V4 answer table, or V9
answer table may be read by the prediction code or the prediction run.

## Known limitations

The reduced closure is an exact port determinant representation where its
baseline resolvent and local device realization are regular. It is not claimed
to be a uniform approximation to spectral abscissa away from a root. The
historical ANDES result is independent phasor evidence with different SG,
exciter, stabilizer, and network assembly; it cannot be upgraded to custom-GFL
EMT validation. ParaEMT portfolio validation is unresolved. No universal
cross-model/network claim is allowed.

## Campaign boundary

This branch adds only blind V4/V9 prediction and reveal, contextual-return
validation, aggregate-matched comparison, runtime measurement, lower-order
screen comparison, fixed g-only remediation, frozen phasor TDS, read-only
ANDES reconciliation, literature audit, report, figures, and handoff. No new
phenomenon, controller search, weak-element/radius campaign, planner,
co-design, PD39, CDW, IEEE-68, or EMT portfolio run is authorized.
