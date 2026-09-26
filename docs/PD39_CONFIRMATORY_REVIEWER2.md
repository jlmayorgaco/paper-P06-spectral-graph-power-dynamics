# PD39 confirmatory campaign — Reviewer 2

## Model and physical credibility

The stock PowerDynamics IEEE-39 network uses documented synchronous-machine,
AVR/governor, and `SimpleGFLDC` converter components. The disturbance API was
implemented faithfully as a deterministic +1% or +5% P/Q pulse at the largest
active-load bus, from 1.0 to 1.1 s, with simulation to 20 s. All ten rerun
traces have finite voltage and frequency observables.

The physical interpretation remains qualified. V8's critical mode is a
complex electromechanical/control mode dominated by machine states at bus 39
and control states at bus 31; it is not a clean converter-only pole. The modal
rerun is reproducible, but reduced Jacobian conditioning is poor: approximately
`10^18–10^22`, with nonzero smallest singular values down to about `4.7e-16`.
That conditioning prevents a strong claim that the V8 cliff is independent of
DAE/gauge numerical sensitivity.

The model does not expose a frozen converter fault-current/short-circuit model
for a mathematically valid SCR/gSCR comparison, so SCR/gSCR were not fabricated.
Controller limits, thermal limits, and detailed protection actions are also not
available as validated constraints in this campaign. Voltage acceptance was
fixed at 0.90–1.10 pu for confirmatory feasibility checks.

## Scores (0–10)

| Criterion | Score | Reason |
|---|---:|---|
| Novelty | 5 | The transition question is useful; the mode itself is not uniquely novel. |
| Correctness | 6 | TDS and outage diagnostics are now instrumented; several radius cases remain unresolved. |
| Model adequacy | 4 | One simplified benchmark model, no cross-model check, severe conditioning. |
| Evidence | 5 | TDS supports ordering, but repairs fail the fresh holdout. |
| Usefulness | 4 | P1 and static P2 sometimes outperform P3 on holdout. |
| IAS impact | 4 | Negative results are valuable but the method is not yet validated. |
| TPWRS readiness | 3 | Needs credible limits, topology/co-design, and independent model validation. |

## Recommendation

Do not claim that V8 is a universal physical incompatibility. State exactly:
V8 is nominally stable but becomes truly unstable in four high-PLL discovery
conditions and fails the 0.05 s^-1 screen in all nine frozen conditions. TDS
supports the relative ordering under the common pulse, while the numerical
conditioning and failed 24-condition repair validation require the result to
be presented as a bounded case study.
