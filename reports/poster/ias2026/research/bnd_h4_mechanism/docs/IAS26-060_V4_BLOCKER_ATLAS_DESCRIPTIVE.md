# IAS26-060 — V4 blocker atlas (descriptive, read-only)

Date: 2026-09-27  
Source run: `20260927T144629Z_d0fecb32_ias26_060_operating_v1`  
Scope: frozen four-device sublattice `V4 = {30, 33, 35, 37}`; original treatment only.

## Purpose and method

This is a post hoc descriptive reclassification of existing IAS26-060B outputs. It runs no power flow, DAE, or eigenvalue solves; it does not alter the frozen run, experiment, poster figures, or claims ledger.

For each physically valid scenario, the 16 original-treatment `alpha_perp` values were classified using the frozen `tau_dec = 1e-8`: stable when `alpha_perp < -tau_dec`, unstable when `alpha_perp > tau_dec`, and otherwise indeterminate. A collective blocker is an inclusion-minimal unstable mask of cardinality at least two, with every proper subset stable. The frozen bus/bit order is `30, 33, 35, 37`.

The 1,000 frozen IDs were retained. The atlas denominator is 967 valid scenarios; 33 were excluded under the frozen physical-feasibility gate, not replaced. Those 33 IDs have `GENERATOR_Q_LIMIT` failures, with load scale from 1.089882 to 1.099854. No outcome-dependent filtering was applied.

The masks were recomputed directly from `derived/MC_CASES.csv`. All 967 recomputed scenario classifications exactly matched `derived/SCENARIO_METRICS.csv` after filtering on `scenario_valid=True` (zero mismatches). This filter matters: finite spectra/mask labels can be present for physically infeasible IDs, but those IDs are not admitted to the atlas.

## Exact blocker-set atlas

Counts partition the 967 valid scenarios. A semicolon-separated row lists all inclusion-minimal blockers found in that same scenario.

| Exact minimal-blocker set in V4 | Scenarios | Share of valid set |
|---|---:|---:|
| None; all 16 portfolios stable | 412 | 42.61% |
| `{30,33,35}` only | 77 | 7.96% |
| `{30,33,35,37}` only (H4) | 190 | 19.65% |
| `{30,33,35}` and `{30,33,37}` | 151 | 15.62% |
| `{30,33,35}`, `{30,33,37}`, and `{33,35,37}` | 137 | 14.17% |
| **Total** | **967** | **100%** |

Additional exact checks:

- Every singleton and every pair was stable in all 967 valid scenarios; there were no indeterminate portfolio classifications.
- H4 itself was unstable in 555 scenarios (57.39%). Of those, it was minimal in 190 and nonminimal in 365 because at least one unstable triple was already present. No scenario had a minimal H4 and a minimal triple simultaneously.
- The 365 alternative-witness scenarios divide into 77 with `{30,33,35}` alone, 151 with two minimal triples, and 137 with three. Thus 288 scenarios had multiple simultaneous minimal blockers.
- Across those scenarios, the observed minimal triples were `{30,33,35}` (365 memberships), `{30,33,37}` (288), and `{33,35,37}` (137). All include bus 33. The fourth V4 triple, `{30,35,37}`, was never unstable.
- Within this restricted V4 census, the minimum blocker order is 3 in 365 scenarios and 4 in 190; it is undefined where no collective blocker exists. This is **not** a full-V9/network-wide blocker-order result.

The event rates are descriptive frequencies in the declared synthetic ensemble, not real-world probabilities. The frozen event table reports Wilson 95% intervals of 17.27–22.27% for H4 minimality and 54.25–60.48% for any V4 collective blocker.

## Descriptive operating-point context

The following are valid-scenario medians and interquartile ranges (Q1–Q3), grouped by exact blocker set. This is exploratory description, not a causal analysis.

| Exact set | n | Load scale `L` (median [Q1,Q3]) | `epsilon_l2` (median [Q1,Q3]) | Saturated units (median [Q1,Q3]) |
|---|---:|---:|---:|---:|
| None | 412 | 0.9411 [0.9217, 0.9589] | 0.1026 [0.0942, 0.1124] | 0 [0, 0] |
| H4 only | 190 | 1.0000 [0.9881, 1.0103] | 0.1033 [0.0943, 0.1131] | 0 [0, 1] |
| `{30,33,35}` only | 77 | 1.0251 [1.0218, 1.0298] | 0.1009 [0.0946, 0.1100] | 4 [3, 4] |
| `{30,33,35}` + `{30,33,37}` | 151 | 1.0482 [1.0404, 1.0596] | 0.1031 [0.0942, 0.1127] | 6 [5, 7] |
| All three minimal triples | 137 | 1.0763 [1.0669, 1.0834] | 0.1040 [0.0925, 0.1130] | 7 [7, 7] |

The observed category ordering is much clearer in aggregate load scale and generator saturation than in the norm of the spatial perturbation vector. This does not establish that aggregate loading or saturation causes blocker switching; the full spatial pattern, network response, and other operating-point variables have not been isolated.

## Claim boundary

**Supported by this atlas:** in the frozen synthetic envelope and V4 sublattice, the identity/order of the minimal unstable witness changes across valid operating points. In the 365 alternative-witness cases, H4 remains unstable but loses minimality to one or more triples.

**Not tested here:** whether the physical local factors remain regular while `I + Q_B` becomes singular for each changing blocker; whether the collective-closure mechanism is invariant; why the blocker family tracks operating point; whether a gain sweep/eigenlocus explains intervention outcomes; or how the full V9/512-portfolio blocker family behaves. Those claims require new diagnostics/campaign authorization and are not implied by the counts above.

## Provenance

SHA256 of frozen source artifacts:

| Artifact | SHA256 |
|---|---|
| `config/IAS26-060_MC_OPERATING_V1.json` | `75db11e64299757bbe170f5ecadcb8c1a3f2ea969583e3c95bcdebc008bfb4b4` |
| `inputs/IAS26-060_SCENARIOS_V1.csv` | `45c425fa4d146aa3a0a3bf26cfbdda1c99590da67ea6c4318cf788c1006b6eb3` |
| `derived/MC_CASES.csv` | `a8bc6c2d9b9e9dcf4025a3460738775c121819c7c53eb29e5ce9ede410cf0a0c` |
| `derived/SCENARIO_METRICS.csv` | `93daef28f7e804f3e4617f03fa4c3a760ddf6943cf67a01aa01031791bf08e96` |

No source data, frozen configuration, manifest, code, run directory, or poster figure was changed to produce this memo. No solve was run.
