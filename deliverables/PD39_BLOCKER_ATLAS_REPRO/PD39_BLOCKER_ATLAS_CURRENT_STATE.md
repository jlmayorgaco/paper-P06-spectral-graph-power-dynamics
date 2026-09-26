# PD39 Physical Blocker Validation + First Compatibility Atlas

## Frozen parent and provenance

- Parent branch: `research/pd39-255plus1-mechanism-validation`
- Parent HEAD: `e81d18ab945bfcf90dd2f4d3df0d0aa76ec2aef3`
- Parent campaign base: `902403cfb0e8dc38323a5623caa51dd663957764`
- New branch: `research/pd39-physical-blockers-atlas-v1`
- No push is authorized.

The authoritative inputs are the committed 255+1 mechanism-validation tables,
especially `results/PD39_255PLUS1_HOLDOUT_CENSUS.csv`,
`results/PD39_255PLUS1_HOLDOUT_CONDITION_SUMMARY.csv`,
`results/PD39_255PLUS1_NUMERICAL_TRUTH_AUDIT.csv`,
`results/PD39_255PLUS1_MODE_TRACKING.csv`,
`results/PD39_255PLUS1_HIGH_PLL_TDS.csv`, and the corresponding final report
and handoff. The prior discovery and confirmatory outputs remain unchanged.

## Model and classifier

- PowerDynamics: 5.0.0, from the locked Julia manifest.
- Julia: 1.11.9, from `Manifest.toml`.
- Network: PowerDynamics IEEE-39 implementation in `src/pd39/model.jl` and
  `src/pd39/confirmatory.jl`.
- GFL replacement candidates: buses `30, 32, 33, 34, 35, 36, 37, 38`.
- Portfolio notation: semicolon-separated candidate buses; `none` is empty.
- True instability classifier: `H0`, with `alpha >= 0`.
- Engineering robustness classifier: `H0.05`, with `alpha > -0.05`.
- Minimality means every proper subset satisfies the relevant strict stable
  inequality. Numerical equality at a boundary is not silently rounded away.

## Frozen operating/control conditions

The 24 fresh conditions are `C01`--`C24` in
`results/PD39_HOLDOUT_CONDITIONS.csv`. They are deterministic maximin-LHS
conditions from the prior campaign and are not resampled or screened. The
condition coordinates are PLL, filter, current-control, load P/Q, IBR P, and
46 line-strength perturbations. The nine discovery scenarios remain the
original nominal and frozen controller corners used in the prior campaign.

## Complete fresh-holdout blocker census

The committed census has 6144 portfolio-condition rows and 17 retained solver
failures. The smallest H0 blockers are all in C12:

| cardinality | condition | portfolio | alpha (s^-1) | converted MW |
|---:|:---:|:---|---:|---:|
| 2 | C12 | `35;36` | 10.187526726660721 | 1210.0 |
| 2 | C12 | `37;38` | 22.52803405237259 | 1370.0 |
| 3 | C12 | `30;32;33` | 21.153510391046698 | 1532.0 |
| 3 | C12 | `33;36;37` | 10.670269367512065 | 1732.0 |
| 3 | C12 | `34;36;37` | 16.322741202239445 | 1608.0 |
| 4 | C12 | `30;33;34;36` | 30.331163918105393 | 1950.0 |
| 4 | C12 | `32;33;36;38` | 20.488951482685103 | 2672.0 |

The same low-order list is obtained under the H0.05 minimal-blocker rule for
cardinalities 2--4. The current census stores alpha but not critical
frequency, so no frequency is inferred here; Gate A must independently record
the frequency, eigenvalue, equilibrium, algebraic conditioning, and physical
classification for every selected blocker and required predecessor.

## Deterministic blocker/control selection

The two cardinality-2 blockers are selected by ascending `(cardinality,
condition, portfolio)` order. The next larger cardinality is represented by
the first cardinality-3 blocker in that same order. Stable controls are chosen
at the same condition and cardinality, first minimizing converted-MW
difference, then IBR-MVA difference, then remaining-inertia difference, then
portfolio lexicographically:

| blocker | control | blocker/control MW | blocker/control IBR MVA | reason |
|:---|:---|---:|---:|:---|
| `35;36` at C12 | `32;36` | 1210 / 1210 | 1500 / 1500 | exact MW/MVA; closest inertia |
| `37;38` at C12 | `36;38` | 1370 / 1390 | 1700 / 1700 | exact MVA; closest available MW |
| `30;32;33` at C12 | `30;33;35` | 1532 / 1532 | 2600 / 2600 | exact MW/MVA; closest inertia |

These selections are frozen before new numerical evaluation. Frequency values
for the selected cases and controls are intentionally not copied from an
unmatched scenario; they will be measured in Gate A.

## Prior warnings and open questions

1. The prior tolerance and descriptor checks were stable, but the frozen
   central finite-difference alpha comparison failed for 72/81 cases because
   the rightmost eigenvalue was highly conditioned/non-normal. This campaign
   must report the same warning rather than erase it.
2. The prior high-PLL TDS supported V8 positive growth versus decaying 7/8
   controls, but did not audit the new C12 low-order blockers.
3. The complete holdout census showed H0=42 and H0.05=47 blockers across
   24 conditions; the exact V8 identity did not generalize.
4. A continuous physical SG-to-GFL homotopy was unavailable in the installed
   discrete replacement semantics.
5. Exact PowerDynamics K/D/Q network-closure objects were not exposed or
   documented; no proxy is permitted in this campaign.

## Current state before new numerics

Known: the C12 census contains physically questionable-looking low-order
blockers with very large positive alpha, and aggregate-matched portfolios can
have different alpha values. Unknown: whether the selected blockers survive an
independent descriptor/Jacobian audit and a common nonlinear TDS, whether the
composition effect is strong under the preregistered thresholds, whether an
atlas is justified, and whether exact contextual-return objects exist.

