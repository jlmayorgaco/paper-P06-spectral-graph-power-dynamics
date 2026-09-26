# CDW claim ledger V1 (initial, before any CDW experiment)

- **Machine-readable matrix:** `results/CDW_CLAIM_MATRIX_V1.csv`, written by
  `experiments/cdw/CDW00_claim_matrix.py`.
- **Allowed statuses:**
  - INHERITED_PROVED
  - INHERITED_VALIDATED
  - HYPOTHESIS
  - NOT_YET_TESTED
  - REFUTED
- No new CDW claim is VALIDATED before its experiment runs.
- Final statuses are written to `results/CDW_MASTER_CLAIM_MATRIX.csv` at the end
  of the campaign. This V1 file is never edited.

## Inherited, proved (TX4)

| id | claim |
|---|---|
| I-P01 | transverse quotient (V01) |
| I-P02 | antichain / any-order safety (V02) |
| I-P03 | boundary localization |
| I-P04 | descriptor-affine closure factorization (V08) |
| I-P05 | boundary sensitivity (V10) |
| I-P06 | deflated zero-frequency port (T08) |
| I-P07 | NP-completeness (V26) |
| I-P08 | connected-cumulant identities (V22), secondary |

## Inherited, validated (IEEE-39, benchmark-specific)

| id | claim | limit |
|---|---|---|
| I-V01 | P4 witness H4 (V03) | model-specific witness |
| I-V02 | policy-dependent H at fixed statics (V06, V07) | — |
| I-V03 | port-derivative line ranking at one condition (V11) | — |
| I-V04 | non-composability robust, witness not (V15, V16) | — |
| I-V05 | frozen vs re-equilibrated line-derivative gap | 12 lines, one condition |
| I-V06 | subcritical Hopf / TDS confirmation (V30) | — |

## Inherited, refuted (kept)

| id | refuted claim |
|---|---|
| I-R01 | static quantities identify failing coalitions |
| I-R02 | GFL line ranking transfers across converter models |
| I-R03 | κ = 4 robust to governors |
| I-R04 | cycle, path or rank explanations |
| I-R05 | \|χ\| design target |

## New CDW claims: all HYPOTHESIS or NOT_YET_TESTED

| id | claim | test |
|---|---|---|
| C-H1 | reversal recurrence | A1 |
| C-H1t | mode-tracked reversal | A2 |
| C-H1b | node-only ranking insufficiency | A3 |
| C-H1c | robust contextuality | A4 |
| C-HS | submodularity | — |
| C-IV | IFT derivative validity | — |
| C-R1 | structural frozen = total for dynamic-only controls | — |
| C-H2n / C-H2l | frozen ≠ total | — |
| C-GB | dynamic links beat static | GOLD-B |
| C-H3 | robust weak corridors | — |
| C-H4 | topology-only removal / creation | — |
| C-EX | exchange rates | — |
| C-H5 | plan-level design | GOLD-D |
| C-SH | Shapley hides reversals | — |
| C-H6s / C-H6m / C-H6r | spectral layer | — |
| C-H7 | certified reduction | GOLD-C |
| C-H8 | modal energy vs limiting mode | — |
| C-H9 | cross-model transfer | — |
| C-AF | Africano / PV | gated |
| C-NL | nonlinear recovery | optional |

**Never claimed, whatever the results:**
- universal weak buses or lines;
- probability from envelope coverage;
- static topology sufficiency;
- cycle or cumulant causality;
- universal converter-model transfer;
- that an α_⊥ improvement implies a larger region of attraction;
- EMT validity.
