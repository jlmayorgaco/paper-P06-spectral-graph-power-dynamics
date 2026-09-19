# P1–P6 execution status

Release status: `NOT_FINAL`.

The required sequence was started at P1. P1’s independent device-level
equation parity passed, but the complete P1 gate was stopped at the
PowerDynamics one-device/infinite-bus integration layer. The minimal
current-source loopback topology returned zero retained network states, so it
cannot certify dynamic network parity. This is recorded as a negative
integration result rather than promoted to a fresh IEEE-39 claim.

| gate | status | result |
|---|---|---|
| P1 exact frozen GFL11 | `STOPPED_BY_GATE` | Julia/Python device parity passes (`max=5.14e-16`); injector compilation passes; network integration is unresolved. |
| P2 PowerDynamics V4 census | `STOPPED_BY_GATE` | Not run; no fresh 16-portfolio census is substituted. |
| P3 collective mechanism | `STOPPED_BY_GATE` | Not run; no fresh Julia blocker exists to analyze. |
| P4 Julia TDS | `STOPPED_BY_GATE` | Not run; historical frozen TDS is not relabeled as fresh. |
| P5 second GFL model | `STOPPED_BY_GATE` | Not run; no materially different model was frozen. |
| P6 blind holdout | `STOPPED_BY_GATE` | Not run; no prediction hash/reveal exists for this branch. |

P7 robustness is not authorized by this status. See
`raw/p1_p6_gate_matrix.csv`, `reports/P1_GFL_PARITY.md`, and
`raw/gfl11/P1_FAILURE_MODES.md`.
