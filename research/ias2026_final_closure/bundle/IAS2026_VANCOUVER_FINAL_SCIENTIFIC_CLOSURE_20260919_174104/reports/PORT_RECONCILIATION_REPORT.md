# Port and Mechanism Reconciliation

## Executed

The algebraic port identities were tested on randomized cases in four pytest
modules. The tests cover local/collective factorization, determinant identity,
the bridge identity without inverting the disturbance matrix, the contextual
Schur identity under its explicit zero-diagonal assumption, and the transverse
rotation/drift construction.

`logs/theory_tests.log` records four passing modules. The contextual-return test
was corrected to state its theorem assumption (`Q_ii = 0`) explicitly; random
nonzero diagonal terms are outside that identity and were not hidden as a
passing case.

## Frozen IEEE-39 evidence

Gate 0 reproduces the nominal P4 H4 value `0.1270064666836382 s^-1` against
the frozen expected `0.1270065 s^-1`, and the P4 proper-subset count `15`.
The documented P4 governor value is `-0.07454098428813066 s^-1`.

## Boundary

The campaign establishes algebraic and frozen-artifact consistency. It does
not claim that the PowerDynamics tutorial has the same converter port, state
ordering, controller policy, or mechanism. Same-model parity is therefore
`STOPPED_BY_GATE`, not a silent success.
