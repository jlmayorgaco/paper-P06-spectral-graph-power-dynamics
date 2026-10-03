# Claim and value audit

Experimental values are generated from frozen evidence inputs into `generated/results.tex`, `generated/claims.tex`, `generated/tables/` and `generated/figures/`. Result values are not manually duplicated in panel prose.

## IEEE-39 flagship

- `H_4 = {30, 33, 35, 37}`, 15/15 stable proper subsets, nominal `+0.127 s^-1`, and the 0.6223 Hz F1 target-band frequency are generated from the frozen evidence inputs. Mode identity remains explicitly unresolved.
- The 2096.6 MW displaced-generation and 4270.7 MVA converter values are prior-audited IEEE-39 values (claim B02); the poster calls this a canonical witness, not a universal four-unit rule.
- Strict M1 reduction is **BLOCKED**. `research/bnd_h4_mechanism/tickets/STATUS.json` records action-space residual 2.340724075262523e-11 above the frozen 1e-12 gate. No reduced-pole equivalence is claimed.
- Port-ratio residual and missing per-subset `Q_{SS}` traces remain **PENDING**.

## Policy atlas and scenario campaign

- Historical F7A/F7B/F7C distinct-hypergraph counts (30/36/16) and F7B path `4 → 3 → 2 → 3 → 4 → ∅` come from the audited claim matrix and policy-atlas CSV; the plot is kept separate from the new campaign.
- The frozen campaign evaluated 1,000 IDs: 967 feasible and 33 infeasible. Event numerators 190, 555, 594, 186 and 373 display **PENDING USER REVIEW** and are described as a synthetic ensemble, not a probability sample. Sources: `generated/tables/MC_EVENT_RATES.csv` and campaign state in `research/bnd_h4_mechanism/tickets/STATUS.json`.
- Pending scenario results do not replace or blank the historical F7A policy atlas.

## Validation and actionability

- 100/100 Python–Julia verdict parity is scoped to same-model agreement, not independent-model validation.
- TDS sign agreement is 56/60 among completed finite records. Four mismatches are zero-frequency cases with missing `R_2`; mode identity remains unresolved. The 30 solver failures of 90 planned records remain represented.
- Alternative-model blind holdout is 4/4 but single-class/non-discriminative; a second preregistered holdout is required. Source: `research/ias2026_final_closure/reports/P6_GENUINE_HOLDOUT_STATUS.md`.
- 106.0 MVA P4 condenser threshold is prior-audited/conditional. Topology modification and safe-plan census remain pending.

## Display policy

No PENDING or BLOCKED result is marked PASS or given a check mark. Reserved panels remain present, historical evidence is scoped, and no favorable result is fabricated from the unreviewed campaign.
