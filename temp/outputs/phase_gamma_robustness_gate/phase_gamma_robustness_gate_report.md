# Gate Gamma Robustness

This gate repeats the equal-cost coordination test using the same networks,
actions, costs, and seeds as Gate 1, but changes the protected-resolvent
input/output weights.

- Commit: `46f48b23d2ca168678e843f32d844eae5f4d93e5`
- Random seed: `20260609`
- Cases: `62`
- Gamma profiles: `12`
- Random diagonal Gamma samples: `8`

## Gamma Robustness Table

| Gamma | kind | ranking agreement S_zeta vs S_g | CV sign agreement | % CV<0 w+D | % CV<0 w+M | % CV<0 D+M |
|---|---:|---:|---:|---:|---:|---:|
| frequency_output | physical | 29.0% | 50.0% | 72.1% | 67.2% | 51.6% |
| identity | physical | 29.0% | 50.0% | 72.1% | 67.2% | 51.6% |
| inertia_output | physical | 21.0% | 49.5% | 73.8% | 73.8% | 48.4% |
| random_diag_00 | random | 19.4% | 52.7% | 67.2% | 78.7% | 43.5% |
| random_diag_01 | random | 27.4% | 49.5% | 63.9% | 75.4% | 37.1% |
| random_diag_02 | random | 27.4% | 56.5% | 67.2% | 68.9% | 25.8% |
| random_diag_03 | random | 32.3% | 57.1% | 62.3% | 68.9% | 22.6% |
| random_diag_04 | random | 16.1% | 50.5% | 67.2% | 73.8% | 35.5% |
| random_diag_05 | random | 29.0% | 52.7% | 70.5% | 73.8% | 32.3% |
| random_diag_06 | random | 14.5% | 58.7% | 65.6% | 78.7% | 24.2% |
| random_diag_07 | random | 27.4% | 50.0% | 70.5% | 68.9% | 35.5% |
| type_input_low_inertia | physical | 22.6% | 48.4% | 73.8% | 70.5% | 48.4% |

## Descriptor Correlations

Correlations use case-level descriptor values and certificate-disagreement scores.

| Gamma | descriptor | target | n | Pearson r | Spearman r |
|---|---|---|---:|---:|---:|
| frequency_output | henrici_departure | ranking_disagreement | 62 | 0.0573 | 0.0695 |
| frequency_output | commutator_chi | ranking_disagreement | 62 | 0.216 | 0.216 |
| identity | henrici_departure | ranking_disagreement | 62 | 0.0573 | 0.0695 |
| identity | commutator_chi | ranking_disagreement | 62 | 0.216 | 0.216 |
| inertia_output | henrici_departure | ranking_disagreement | 62 | -0.164 | -0.192 |
| inertia_output | commutator_chi | ranking_disagreement | 62 | 0.142 | 0.125 |
| random_diag_00 | henrici_departure | ranking_disagreement | 62 | 0.21 | 0.249 |
| random_diag_00 | commutator_chi | ranking_disagreement | 62 | 0.00526 | 0.0183 |
| random_diag_01 | henrici_departure | ranking_disagreement | 62 | -0.237 | -0.221 |
| random_diag_01 | commutator_chi | ranking_disagreement | 62 | 0.204 | 0.217 |
| random_diag_02 | henrici_departure | ranking_disagreement | 62 | -0.031 | 0.00101 |
| random_diag_02 | commutator_chi | ranking_disagreement | 62 | 0.135 | 0.128 |
| random_diag_03 | henrici_departure | ranking_disagreement | 62 | 0.0985 | 0.133 |
| random_diag_03 | commutator_chi | ranking_disagreement | 62 | 0.275 | 0.291 |
| random_diag_04 | henrici_departure | ranking_disagreement | 62 | 0.255 | 0.265 |
| random_diag_04 | commutator_chi | ranking_disagreement | 62 | 0.0593 | 0.0147 |
| random_diag_05 | henrici_departure | ranking_disagreement | 62 | -0.00578 | -0.00596 |
| random_diag_05 | commutator_chi | ranking_disagreement | 62 | 0.213 | 0.214 |
| random_diag_06 | henrici_departure | ranking_disagreement | 62 | 0.0467 | 0.0576 |
| random_diag_06 | commutator_chi | ranking_disagreement | 62 | 0.23 | 0.229 |
| random_diag_07 | henrici_departure | ranking_disagreement | 62 | -0.241 | -0.195 |
| random_diag_07 | commutator_chi | ranking_disagreement | 62 | -0.0485 | -0.0273 |
| type_input_low_inertia | henrici_departure | ranking_disagreement | 62 | -0.128 | -0.14 |
| type_input_low_inertia | commutator_chi | ranking_disagreement | 62 | 0.141 | 0.112 |

## Verdict

- Status: **DISCREPANCY_ROBUST_WEAK_INTERPRETABILITY**
- Implication: The certificate discrepancy survives physical Gamma choices, but descriptor correlations are weak.
- Discrepancy robust over physical Gamma choices: `True`
- Coordination robust over physical Gamma choices for w+D and w+M: `True`
- Max absolute Spearman descriptor/ranking-disagreement correlation: `0.216`

## Scope

- Gamma matrices are fixed at the base case for each network; actions change the system, not the protection standard.
- The frequency-output profile is intentionally equivalent to identity because Gate 1 already measured omega output.
- This is still a reduced dynamic-model gate, not a full ANDES IBR validation.