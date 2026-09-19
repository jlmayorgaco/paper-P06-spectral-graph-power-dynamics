# PowerDynamics Reconciliation

This report is populated by the Julia gate. The required distinction is:

- `POWERDYNAMICS_VALIDATED`: an isolated PowerDynamics model reaches the declared equilibrium and supports the declared test.
- `STOPPED_BY_GATE`: package, model, or equilibrium prerequisites fail.
- `RETROSPECTIVE`: existing ANDES or frozen same-code evidence, not PowerDynamics evidence.

No PowerDynamics result may be described as same-model mechanism validation until the equilibrium and port/mode reconciliation gates pass.

## Executed result

The isolated Julia project completed Gate A for the official PowerDynamics
IEEE-39 tutorial equilibrium: three nonlinear tolerances, three deterministic
initial guesses, non-mutating and mutating initialization paths, residual
checks, spectrum consistency, and Jacobian conditioning.

| Field | Observed |
|---|---|
| Package | PowerDynamics 5.0.0 |
| Julia | 1.11.9 |
| Network | 39 buses, 46 branches |
| State | `NWState` |
| Power-flow interface entries | recorded by the official tutorial run |
| Result label | `POWERDYNAMICS_VALIDATED` |

The executable evidence is in `raw/powerdynamics/pd39_equilibrium_gate.md` and
`raw/powerdynamics/gate_a_paths.csv`. The source is the package's official
`docs/examples/ieee39_part1.jl` tutorial.

The subsequent P1 device-level check independently transcribes the frozen
GFL11 equations in Julia and compares them with a separate Python oracle over
64 cases. It passes with maximum relative error `5.14e-16`, and the exact
injector compiles with `current_source=true`. The minimal one-device
current-source/infinite-bus network harness is nevertheless
`STOPPED_BY_GATE`: its returned state retained zero dynamic states. This is a
real integration negative result and does not establish full network parity.

## Scope boundary

This is an independent package/tutorial equilibrium check, not a same-model
reimplementation of the frozen L0 GFL network. The tutorial model contains
synchronous-machine/governor components and does not expose the campaign's
converter state, port convention, contextual-return construction, or blocker
mechanism. Therefore the same-model parity and same-mechanism gates remain
`STOPPED_BY_GATE`/`NOT_TESTED`; no PowerDynamics result is used to promote a
claim about the frozen GFL mechanism. The independent modal assignment audit
for Gate A is in `raw/powerdynamics/gate_a_assignment_report.md`; unreduced
Jacobian conditioning is recorded as a diagnostic, not a passed condition.
