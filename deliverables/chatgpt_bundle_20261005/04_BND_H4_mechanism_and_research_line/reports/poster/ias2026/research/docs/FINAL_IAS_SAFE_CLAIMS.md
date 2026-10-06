# Final IAS poster: safe claims

This file does **not** edit the poster. It lists wording the poster may use
after the closure campaign. It supersedes `POSTER_FINAL_SAFE_CLAIMS.md` and v2
where they differ. The prohibited forms are in `FINAL_REJECTED_WORDING.md`.

## Required on the poster

- The scope sentence (one line): *"IEEE-39 as frozen: no primary frequency
  control (D = 0, no governor); stability is relative / transverse."*
- Units: *"2096.6 MW of active dispatch on 4270.7 MVA of converter rating."*
- Simulation type: *"nonlinear phasor-domain simulation (not EMT)."*

## Safe claims (with the number to show)

| # | claim | number | source |
|---|---|---|---|
| S1 | The converter reactive policy changes **which** minimal replacement coalitions fail, not only whether. | 30 / 36 / 16 distinct hypergraphs on three policy planes; 0 label changes in the transverse re-audit of 345 229 points | FC01, F7 |
| S2 | Every proper subset of the four-bus flagship is stable, while the full portfolio is unstable under one policy and composable under full voltage-support gain. | `kappa = 4` at P4; empty at `P_inf` | FC03 (frozen columns) |
| S3 | The failing mode is the 0.64–0.71 Hz inter-area family, not the amount of dispatch replaced. | 25/25 Pg-matched controls leave it stable | FC04 |
| S4 | A damped condenser restores composability; an undamped one creates its own unstable mode. | 106.0 MVA (2.48 %) at P4; TDS 2.47 % | G1, FC12 |
| S5 | Any-order safety = containing no minimal incompatible coalition. Implementation order matters for a few targets. | 3 of 327 stable targets | FC10 |
| S6 | Time-domain simulation confirms the small-signal verdicts, and within model validity the finite-disturbance structure equals the spectral one. | 32/32; 6/6 (point, family) | G2, FC05/06 |
| S7 | *Caveat box:* with the source's documented turbine governors, the four-bus coalition at P4 is stabilized, and the structure persists in a smaller region. | `alpha` +0.127 → −0.075 | FC03 |

## Safe only with the qualifier shown

| claim | required qualifier |
|---|---|
| subcritical boundaries | "at the three examined boundaries" |
| port closure | "on a preregistered Kundur holdout (28/29)" |
| planning rating 684–792 MVA | "to keep a 200 MW disturbance within the declared converter voltage envelope, in this model; not a stability requirement" |
| Kundur recurrence | "on the Kundur two-area system"; never "universal" |
| ANDES | "reproduces the network and base inter-area mode"; never "validates the converter results" |

## Not for the poster

Nonlinear composability as a headline, the resilience-complex topology, the
NP-completeness result, and any certificate are excluded (see
`FINAL_REJECTED_WORDING.md`).
