# Evidence-backed poster claim ledger

| ID | Statement | Status | Direct evidence and boundary |
|---|---|---|---|
| C1 | At fixed 88.4551408% GFL, gain retuning raised the numerical delay-to-`Re(s)=−0.05 s⁻¹` margin from 37.38719 to 43.79719 ms (+17.1449%). | FULL-MODEL NUMERICAL RESULT | `F00_PARENT_REPRODUCTION.csv`, `TABLE_F2_FINAL_OPTIMIZATION_TRACE.csv`, `evaluations/full20_step15_medium/ROOTS.csv`, five-event `Q0_RESULT.toml`. Not a proven maximum. |
| C2 | At fixed replacement, twenty PLL gains modify the 203-state linearization through an action of rank at most ten. | EXACT STRUCTURAL IDENTITY, NUMERICALLY CHECKED | `THEORY_CLOSURE.md`, `TABLE_F11_GAIN_ACTION_RANK.csv`; largest relative reconstruction error `2.98e−14`. Individual differences may have rank below ten. |
| C3 | The best tested design has five delayed modal crossings within 0.001194 ms; the local multimode LP balances all five and the actuator surrogate. | FULL-MODEL NUMERICAL ROOTS + LOCAL REDUCED-MODEL KKT | `ROOTS.csv`, `KKT_LP_full20_step15_medium.csv`. The LP multipliers do not prove a full-model optimum. |
| C4 | A nearly orthogonal limiting eigendirection appeared between the old best and the accepted 41 ms design (right MAC 0.006457). | SUPPORTED NUMERICALLY | `TABLE_F3_ACTIVE_MODAL_FAMILIES.csv`; reference-overlap labels are provisional, not a certified branch continuation. |
| C5 | Linear positive-delay DDE time histories decay at 0.90/0.98 and grow at 1.02 times the measured security threshold. | NUMERICALLY VALIDATED LINEAR V2 | `TABLE_F8_POSITIVE_DELAY_VALIDATION.csv` and `TIME_DOMAIN_V2_full20_step15_medium.csv`; nonlinear delayed event safety is untested. |
| C6 | For the step-13 bus-16 event, 100.01 MW passes the actuator guard but 100.02 MW fails. | NEGATIVE ONE-EVENT RESULT | `TABLE_F13_EVENT_SIZE_HEADROOM.csv`, one-event output records. This is specific to one tuning, not a universal architecture ceiling. |
| C7 | The final tuning passes all five frozen zero-delay events with 0.490019 Hz maximum excursion, 0.205474 Hz/s maximum RoCoF, and 0.002020306 minimum SG actuator slack. | FULL-MODEL NUMERICAL RESULT | `event_validations/full20_step15_medium/Q0_EVENT_METRICS.csv`; the GFL current ratio is observed but no current-limiter safety limit was enforced. |

Do **not** claim global/local optimality, a certified maximum delay, full
nonlinear positive-delay security, robustness to the +1 MW event, or an
optimized Pareto frontier. `VALIDATION_GATES.csv` and `FINAL_REPORT.md` list
the missing gates.
