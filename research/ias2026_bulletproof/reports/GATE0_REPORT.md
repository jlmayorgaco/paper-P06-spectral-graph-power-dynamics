# Gate 0 Report

status: PASS

Gate 0 uses the corrected same-policy modal-scope source staged under `raw/modal_scope`. The stale FC10/default-policy V9 reveal is not used for current P4 scope conclusions.

## V4 and modal closure

- `P4 H4 alpha [s^-1]`: expected `0.1270065`, observed `0.1270064671440848`, status `PASS`, label `IEEE39_VALIDATED`
- `P4 H4 frequency [Hz]`: expected `0.6222797`, observed `0.6222796695779256`, status `PASS`, label `IEEE39_VALIDATED`
- `P4 governed H4 alpha [s^-1]`: expected `-0.0745`, observed `-0.07454098428813066`, status `PASS`, label `IEEE39_VALIDATED`
- `P4 proper subsets stable`: expected `15`, observed `15`, status `PASS`, label `IEEE39_VALIDATED`
- `target-family minimal blocker antichain exact`: expected `True`, observed `True`, status `PASS`, label `NUMERICALLY_VERIFIED`
- `target-family minimal blockers`: expected `14`, observed `14`, status `PASS`, label `NUMERICALLY_VERIFIED`
- `target-family blocker orders`: expected `5x4, 6x5, 3x6`, observed `5x4, 6x5, 3x6`, status `PASS`, label `NUMERICALLY_VERIFIED`
- `target-family kappa`: expected `4`, observed `4`, status `PASS`, label `NUMERICALLY_VERIFIED`
- `H4 boundary g*`: expected `0.2076814`, observed `0.2076814051903784`, status `PASS`, label `RETROSPECTIVE`
- `H4 boundary frequency [Hz]`: expected `0.7064248`, observed `0.7064247880117632`, status `PASS`, label `RETROSPECTIVE`
- `Q_H eigenvalue sign near -1`: expected `-1`, observed `-1`, status `PASS`, label `NUMERICALLY_VERIFIED`
- `contextual return eigenvalue sign near +1`: expected `1`, observed `1`, status `PASS`, label `NUMERICALLY_VERIFIED`
- `distance of Q_H to -1`: expected `0.0`, observed `3.5098754538943794e-08`, status `PASS`, label `NUMERICALLY_VERIFIED`
- `max distance of contextual return to +1`: expected `0.0`, observed `2.2805379282384656e-07`, status `PASS`, label `NUMERICALLY_VERIFIED`
- `maximum Schur residual`: expected `0.0`, observed `3.6258911332571516e-16`, status `PASS`, label `NUMERICALLY_VERIFIED`
- `minimum physical local sigma_min(I+M_ii)`: expected `0.2973`, observed `0.2973268809593455`, status `PASS`, label `NUMERICALLY_VERIFIED`
- `collective sigma_min(I+Q_H)`: expected `1.26e-08`, observed `1.261484530209848e-08`, status `PASS`, label `NUMERICALLY_VERIFIED`
- `proper-subset collective minimum`: expected `0.2945`, observed `0.2945460505423897`, status `PASS`, label `NUMERICALLY_VERIFIED`

## Same-policy V9

The global spectrum and targeted 0.3-1.5 Hz task are reported separately; the targeted result is not a global safety certificate.
- `V9 same-policy global portfolios`: expected `512`, observed `512`, status `PASS`, label `NUMERICALLY_VERIFIED`
- `V9 same-policy global correct`: expected `402`, observed `402`, status `PASS`, label `NUMERICALLY_VERIFIED`
- `V9 same-policy global false-safe`: expected `110`, observed `110`, status `PASS`, label `NUMERICALLY_VERIFIED`
- `V9 same-policy global false-unstable`: expected `0`, observed `0`, status `PASS`, label `NUMERICALLY_VERIFIED`
- `V9 same-policy global false-safe class`: expected `APERIODIC_REAL`, observed `APERIODIC_REAL`, status `PASS`, label `NUMERICALLY_VERIFIED`
- `V9 target 0.3-1.5 Hz correct`: expected `511`, observed `511`, status `PASS`, label `NUMERICALLY_VERIFIED`
- `V9 target 0.3-1.5 Hz false-safe`: expected `1`, observed `1`, status `PASS`, label `NUMERICALLY_VERIFIED`
- `V9 target 0.3-1.5 Hz false-unstable`: expected `0`, observed `0`, status `PASS`, label `NUMERICALLY_VERIFIED`

## Retrospective and limitations

Frozen G2 TDS and topology records remain retrospective. PowerDynamics, second-model, robust, blind, EMT, and full repair-path gates are not promoted by this Gate 0 run.

Sources are bundled under `raw/`; derived tables are under `derived/tables`.
