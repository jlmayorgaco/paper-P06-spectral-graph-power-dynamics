# Claim audit

All experimental numbers in `main.tex` are referenced through macros emitted by `generate_poster_assets.py`; the poster source contains no hand-entered experimental result values. Figure rows are regenerated from the cited source tables. See `generated/results.tex` and `generated/claims.tex` for the editable values and release states.

| Poster content | Value source | Release state |
|---|---|---|
| IEEE-39 witness, stable proper subsets, nominal growth rate and frequency, P4 parameters, displaced MW/MVA | Frozen `POSTER_NUMBERS_FINAL.csv`, `IAS26-050_INTERVENTION_V1.json`, and audited canonical regression | VERIFIED |
| Four-port boundary gain, zero real part, frequency, local device singular values, collective closure | `F2_BOUNDARY_PORT_AUDIT.csv` | VERIFIED for the audited boundary; strict M1 reduction remains BLOCKED |
| M1 residual and threshold | `STATUS.json` ticket IAS26-010 | BLOCKED; values shown with that state |
| F7A/F7B/F7C distinct blocker counts and policy atlas | Frozen policy-map table and claim matrix | CONDITIONAL / historical atlas |
| Valid and infeasible scenario totals; H4-minimal, any-collective, alternative-witness counts | `SCENARIO_METRICS.csv` and `MC_EVENT_RATES.csv` | PENDING USER REVIEW; labeled as a descriptive synthetic ensemble |
| Control-retuning outcomes and editable scatter | `IAS26-FINAL_P4_MC_DATA.csv` | PENDING USER REVIEW |
| Phasor-DAE traces | `IAS26-FINAL_P4_TDS_TRACES.csv` | Same-model phasor DAE; finite-TDS agreement remains limited to the audited 56/60 set |
| Synchronous support threshold | Claim-matrix B05 | CONDITIONAL |
| 9-candidate safe-set census, topology blind test, finite-disturbance envelope | No released result used | PENDING; no favorable number is printed |
| Alternative converter and IEEE-68 | Readiness report | NOT RUN |

The policy image, scenario progress bar, intervention scatter, TDS traces, and result strip are generated from the listed files. The scenario ratios describe this frozen synthetic ensemble; they are not presented as probabilities. No EMT validation, physical switching transient, current-limit, or DC-link result is claimed.
