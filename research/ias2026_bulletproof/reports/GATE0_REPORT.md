# Gate 0 Report

status: FAIL (negative evidence preserved)

The canonical frozen P4 result is reproduced from immutable source tables. The historical TX4 V9 reveal is audited as retrospective evidence; its plan-brief expectation of 511/512 is contradicted by the archived 395/512 result and is therefore not promoted.

## Results

### `gate0_v4.csv`

- `P4 H4 alpha [s^-1]`: expected `0.1270065`, observed `0.1270064666836382`, status `PASS`, label `IEEE39_VALIDATED`
- `P4 governed H4 alpha [s^-1]`: expected `-0.0745`, observed `-0.07454098428813066`, status `PASS`, label `IEEE39_VALIDATED`
- `P4 proper subsets stable`: expected `15`, observed `15`, status `PASS`, label `IEEE39_VALIDATED`

### `gate0_boundary.csv`

- `H4 boundary g*`: expected `0.20768`, observed `0.20768140441926236`, status `PASS`, label `RETROSPECTIVE`
- `H4 boundary frequency [Hz]`: expected `0.706`, observed `0.7064247878829905`, status `PASS`, label `RETROSPECTIVE`
- `boundary Q sign`: expected `-1`, observed `-1`, status `PASS`, label `RETROSPECTIVE`
- `contextual return sign`: expected `1`, observed `1`, status `PASS`, label `RETROSPECTIVE`

### `gate0_v9.csv`

- `V9 portfolios`: expected `512`, observed `512`, status `PASS`, label `RETROSPECTIVE`
- `V9 correct classifications`: expected `511`, observed `395`, status `FAIL`, label `REFUTED`
- `V9 false-safe`: expected `0`, observed `116`, status `FAIL`, label `REFUTED`
- `V9 false-unstable`: expected `0`, observed `1`, status `FAIL`, label `REFUTED`
- `V9 blocker antichain exact`: expected `True`, observed `False`, status `FAIL`, label `REFUTED`
- `V9 kappa exact`: expected `True`, observed `True`, status `PASS`, label `RETROSPECTIVE`

### `gate0_tds.csv`

- `declared TDS verdict agreement`: expected `32`, observed `32`, status `PASS`, label `NONLINEAR_TDS_VALIDATED`
- `P4 H4 TDS unstable rows`: expected `2`, observed `2`, status `PASS`, label `NONLINEAR_TDS_VALIDATED`

### `gate0_topology.csv`

- `F7A boundary records`: expected `None`, observed `12195`, status `PASS`, label `RETROSPECTIVE`
- `F7A distinct witness records`: expected `None`, observed `16`, status `INFO`, label `RETROSPECTIVE`

## Gate decision

G0 nominal P4 reproduction: PASS for the available frozen IEEE-39 records.
G0 V9 expected-number check: FAIL/REFUTED as a historical expectation; the actual archive is retained.
Dependent new-science phases must not treat V9 exact antichain transfer as established.

Sources are copied into `raw/gate0`; derived tables are under `derived/tables`.
