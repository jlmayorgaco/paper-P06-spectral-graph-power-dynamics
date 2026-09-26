# TX4 blind portfolio prediction — ChatGPT handoff

## Git

- Branch: `research/tx4-blind-portfolio-prediction-final`
- Parent freeze: `69f200dfe1dfb74cc4f678ad25a0c3b1751d62a6`
- Artifact HEAD before final packaging: `96052a3e69964ff76a0fb47aa4767b2d8944f7ca`
- Freeze tag: `TX4_FINAL_MANUSCRIPT_FREEZE`
- Push: none
- Working tree: the campaign artifacts are committed; LaTeX auxiliary logs
  are intentionally not part of the package.

The initial blind campaign commits are `ae1ed9ac`, `9b86b0ef`, `034e2d49`,
`7c2a8fd6`, `02d2d899`, `13edf644`, and `4f42ca0e`. The final report and
validation artifacts were added in `96052a3e`. A later packaging commit may
move HEAD without changing the numerical results.

## What was tested

This is the TX4 blind-prediction and contextual-return campaign. It uses the
frozen IEEE-39 phasor-domain DAE, P4
`(g,k,t,h)=(0.03625,1.425,1.5,1)`, and the positive-frequency band
`0.3--1.5 Hz`.

The blind stage built one all-SG network kernel and one local SG-to-GFL
replacement model for each candidate. It enumerated all 16 V4 subsets and all
512 V9 subsets, wrote predictions, hashed them, and committed them before the
full-order answer tables were loaded. The V4 answer is therefore the primary
blind result; V9 is a deliberately negative transfer benchmark.

## Headline results

### V4 flagship

- Candidates: `{30,33,35,37}`; all 16 subsets.
- Full-order reveal: 15 stable, H4=`30+33+35+37` unstable.
- Blind accuracy: `16/16`.
- False-safe: `0`; false-unstable: `0`.
- Exact blocker family: yes.
- Exact minimum cardinality: `kappa=4`.
- H4 full-order alpha: `+0.1270064671 s^-1` at `0.6222796696 Hz`.
- Reduced H4 root: `+0.1270063376 s^-1` at `0.6222796710 Hz`.

### Contextual return

The standard block determinant identity used is

```text
det(I + Q_H) = det(I + Q_RR) det(I - R_i|R)
R_i|R = Q_iR (I + Q_RR)^(-1) Q_Ri.
```

At the controller boundary, the corrected sign audit finds a nearest `Q_H`
eigenvalue at `-1` within `3.51e-8`; the contextual-return eigenvalues reach
`+1` within at most `2.28e-7`. The maximum Schur residual is `3.63e-16`.
The fresh physical local factors `I+M_ii` remain regular (minimum singular
value `0.2973` over the sweep), while the collective factor approaches
singularity (`1.26e-8`). This is standard Schur/block algebra applied to this
model, not a new theorem.

The full-order g boundary is `0.20768140519`; the reduced reproduction is
`0.20768138607` at `0.706424783 Hz`, absolute difference `1.91e-8`.
Important: the reduced boundary search was executed after reveal, so it is
reported as a reduced reproduction, not as a strictly blind boundary claim.

### V9 transfer benchmark and modal-scope correction

- Candidates: `{30,31,32,33,34,35,36,37,38}`; all 512 subsets.
- The archived FC10 table reported 327 stable and 185 unstable, with 395/512
  correct, 116 false-safe, and 1 false-unstable. That table is retained as a
  historical frozen artifact, but it is not on the same P4 policy as the
  blind predictor.
- Same-policy P4 recomputation for this audit gives 332 stable and 180
  unstable. The global comparison is 402/512 correct, 110 false-safe, and 0
  false-unstable.
- The targeted `0.3--1.5 Hz` comparison is 511/512 correct, with 1
  false-safe and 0 false-unstable. It is not a global safety classifier.
- All 110 same-policy global false-safe cases have an aperiodic real global
  critical mode; there are no slow false-safe cases and no target-band global
  false-safe cases.
- Exact blocker-family recovery: no.
- Minimum blocker cardinality recovery: yes (`kappa=4`).

Do not describe the method as a general exact portfolio predictor.

### Runtime

| set | full total (s) | reduced total (s) | full/reduced |
|---|---:|---:|---:|
| V4, 16 cases | 2.3965 | 4.3599 | 0.5497 |
| V9, 512 cases | 72.4899 | 309.6489 | 0.2341 |

The reduced method was slower on this hardware. Its value is finite blind
prediction and mechanism decomposition, not speed.

### Nonlinear TDS

The frozen nonlinear phasor-domain disturbance is the G2 bus-20 active-load
pulse: `+2%` for `0.2 s`. Three completed cases are in the package:

1. proper triple `30+33+35`: stable/decay;
2. H4 at `g=0.03625`: unstable/growth;
3. H4 at fixed `g=0.25`: stable/decay.

Growth-sign agreement is `3/3`; maximum primary exponent error is
`3.51e-4 s^-1`. This is phasor-domain TDS, not EMT validation.

### Aggregate composition counterexample

The fixed-rule nearest same-cardinality V9 pair is:

```text
A = 30+31+32+34+36+38, MW=3659.869758, MVA=6509.2, alpha=+324.2486, unstable
B = 31+32+34+35+37+38, MW=3652.304711, MVA=6499.9, alpha=-0.166389, stable
```

Normalized MW/MVA distance is `0.2513%`. Exact aggregate matching was not
available, so this is supporting near-match evidence, not an exact equal-MW
theorem.

## Interpretation and claim limits

The critical mode is classified as a physical inter-area electromechanical
mode with network-mediated SG-to-GFL participation. Across the audited cases,
maximum equilibrium residual is `1.30e-12` and the independent eigensolver
error is `0` to reported precision. No DAE/index or numerical pathology is
indicated. The frozen model is one IEEE-39 phasor-domain implementation with
one custom ten-state GFL model.

The campaign supports the narrow statement:

> In the frozen TX4 model at P4, one reusable network kernel plus four local
> replacement models predicted the minimal H4 dynamic instability before
> full-order reveal, and the contextual return identifies a collective rather
> than local closure failure.

It does not support universal minimality, probabilistic generalization,
global stability-radius claims, EMT validation, or TPWRS readiness. Existing
ANDES evidence is limited to supported phasor comparisons; the custom-GFL
g-boundary equivalence was not tested in ANDES. The reduced g-boundary timing
deviation and V9 failures must remain visible on the poster.

## Recommended ChatGPT review

Please audit the final report and claim matrix against the CSV evidence, check
that the poster wording does not upgrade the V4 result into universal
generalization, and assess whether the finite blind prediction plus
contextual-return mechanism is sufficiently distinct from classical return
ratio, impedance, Schur-complement, and eigenvalue-sensitivity methods.

The modal-scope audit supports `CASE A` for the corrected, targeted claim:
the V4 mechanism and the `0.3--1.5 Hz` transfer task are clean, while the
global V9 result remains explicitly limited by a modal-scope correction and a
policy-confounded legacy table. This is still not a fully generalized journal
result.

The optional IEEEtran note in
`reports/papers/ias_minimal_dynamic_incompatibility/main.pdf` is a verified
3-page technical note. It is not represented as a six-page IAS paper.

## Main files to inspect

- `docs/TX4_BLIND_PREDICTION_FINAL_REPORT.pdf`
- `docs/TX4_BLIND_PREDICTION_FINAL_REPORT.md`
- `docs/TX4_BLIND_PREDICTION_NOVELTY_AUDIT.md`
- `docs/TX4_CONTEXTUAL_RETURN_THEOREM.md`
- `docs/TX4_IAS_FINAL_POSTER_CONTENT.md`
- `results/TX4_BLIND_PREDICTION_HEADLINE.json`
- `results/TX4_FINAL_CLAIM_MATRIX.csv`
- `results/TX4_V4_BLIND_VS_FULL.csv`
- `results/TX4_V9_BLIND_VS_FULL.csv`
- `results/TX4_G_BOUNDARY_BLIND_PREDICTION.csv`
- `results/TX4_CONTEXTUAL_RETURN_BOUNDARY.csv`
- `results/TX4_MINIMALITY_SEPARATION.csv`
- `results/TX4_LOCAL_VS_COLLECTIVE.csv`
- `results/TX4_RETURN_DERIVATIVE_CHECK.csv`
- `results/TX4_TDS_FINAL.csv`
- `results/TX4_RUNTIME_SUMMARY.csv`
- `figures/tx4_blind_prediction/`

The upload ZIP is the curated analysis package. The repro ZIP contains the
source and raw artifacts needed to inspect or rerun the campaign.
