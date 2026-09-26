# TX4 final modal-scope audit - ChatGPT handoff

## Provenance

- Branch: `research/tx4-final-modal-scope-audit`
- Starting commit: `2cbb860eb4ba053edee99c1a040c8024f3f6828e`
- Parent freeze: `69f200dfe1dfb74cc4f678ad25a0c3b1751d62a6`
- Content artifact commit: `1380413d` (packaging-only changes may advance HEAD)
- Push: none
- Scope: sign correction, physical local-factor audit, V9 global versus EM
  modal truth, P7 decision, and equal-horizon replot only.

## Executive answer

`FINAL CASE: A`.

The V4 mechanism survives the correction. The exact closure signs are
`Q_H -> -1` and contextual return `R_i|R -> +1`. The physical local factors
`I+M_ii` remain regular while the collective factor `I+Q_H` becomes singular.

The same-policy P4 V9 census is 402/512 correct globally, with 110
false-safe and 0 false-unstable. In the targeted 0.3--1.5 Hz band it is
511/512 correct, with 1 false-safe and 0 false-unstable. All 110 global
false-safe cases are aperiodic real-mode misses; there are zero slow
false-safe cases. The archived 395/512 and 116 false-safe table used a
different FC10/default policy and is preserved as a policy-confounded legacy
artifact, not merged with the new truth.

## Key numerical checks

| quantity | result |
|---|---:|
| nearest Q distance to -1 | `3.509875e-08` |
| maximum return distance to +1 | `2.280538e-07` |
| maximum Schur residual | `3.625891e-16` |
| minimum physical local sigma over sweep | `0.2973268809593455` |
| collective sigma minimum | `1.261484530209848e-08` |
| minimum proper-subset collective sigma at boundary frequency | `0.2945460505423897` |
| V4 H4 alpha | `+0.1270064671 s^-1` |
| V4 H4 frequency | `0.6222796696 Hz` |

## Modal scope

Global V9 counts under the same P4 policy: 332 stable, 180 unstable, 110
false-safe, 0 false-unstable. Target-band counts: 442 stable, 70 unstable,
1 false-safe, 0 false-unstable. The global mode census has 271 aperiodic, 8
slow-outside-band, 233 target-band, and 0 fast cases. The global false-safe
set is 110 aperiodic cases. The targeted EM false-safe is a separate case and
must not be added to the global taxonomy.

## P7 and TDS

P7 is `SUPPLEMENT ONLY`. Its unstable member has a global aperiodic alpha of
`+422.8355 s^-1` but is stable in the target band at `-0.5313 s^-1`, while
the comparison member is stable at `-0.1511 s^-1` and `1.3203 Hz`. It is not a
validated approximately 0.62 Hz mechanism pair.

No optional slow false-safe TDS was run because the preregistered slow set is
empty. Existing phasor-domain TDS traces were re-plotted on equal 0--30 s
axes. This is not EMT validation.

## Claims to retain

- Finite V4 H4 collective closure blocker at P4.
- Physical local factors remain regular in the corrected H4 boundary sweep.
- Correct opposite signs for Q closure and contextual return.
- Targeted EM-band transfer is nearly exact in this frozen V9 census.
- Global V9 transfer is limited by 110 same-policy aperiodic false-safe cases.

## Claims to remove or downgrade

- Do not say the contextual return approaches `-1`.
- Do not use normalized `I+Q` diagonal blocks as physical local evidence.
- Do not cite 395/512 or 116 as a same-policy P4 result.
- Do not call the targeted EM classifier a global safety classifier.
- Remove P7 from the main poster; keep it supplement-only.
- Do not claim universal generalization, EMT validation, six-page-paper
  readiness, or TPWRS readiness.

## Review request

Please inspect the final report, claim matrix, sign audit, physical
local-versus-collective table, V9 confusion/taxonomy tables, P7 decision, and
figures. Check that the policy-confounded legacy table remains visible but is
not used as same-policy evidence. The correct recommendation is an IAS
poster redesign with a narrow, scope-aware claim; broader journal work should
be narrowed or deferred.
