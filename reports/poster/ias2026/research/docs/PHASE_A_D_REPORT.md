# IEEE-39 SG to PV-GFL replacement: PHASE A-D report

Frozen operating point, matched reactive policy, `rho = 1` for the census.
All numbers regenerate from `experiments/E10..E13` with manifests under
`results/manifests/`.

## PHASE A — the benchmark is sound and eligible

| check | value | reference |
|---|---|---|
| Ybus vs ANDES | `1.14e-13` | andes 2.0.0 |
| bus voltage magnitude vs ANDES | `1.72e-07` | |
| bus angle vs ANDES | `2.44e-07` rad | |
| power-flow mismatch | `2.75e-12` | |
| voltage range | 1.0019 to 1.0730 pu | |
| equilibrium residuals | `7.9e-14` / `2.3e-13` | |
| `cond(gz)` | `5.02e+02` | index-1 confirmed |
| differential / algebraic states | 70 / 78 | 10 machines |
| base spectral abscissa | `-0.126478` | non-reference modes |
| base `zeta_min` (0.1-5 Hz) | `+0.026124` at 0.9101 Hz | 30 oscillatory modes |

Two modes at `|lambda| = 1.8e-05` are excluded as the free common angle and
common frequency, with 100 percent of their participation in machine `delta` and
`omega`. This is a consequence of retaining no governor and is declared, not
discovered after the fact.

Three fidelity decisions, all logged in `results/discovery_log.csv` before any
replacement was evaluated:

1. The **IEEEX1 exciter is not retained**. Six of ten machines carry `KE < 0`,
   and ANDES' own eigenanalysis of the source case returns six unstable real
   modes near `+1.0` with an algebraic Jacobian at `rcond ~ 1e-20`. A
   first-order regulator with the case `KA` and `TE` replaces it.
2. The **TGOV1N governor is not retained**: `T3 = 2.1 s` is outside the 0.2-10 Hz
   band and `Dt = 0` in the case.
3. The preregistered absolute gate `zeta_req = 0.03` was **amended**: the base
   case measures `0.026124`, so an absolute 3 percent gate would fail every
   portfolio vacuously. The primary criterion became threshold-free instability;
   damping is reported relative to base. The amendment was made after the base
   case and the single-replacement atlas, and before any pair or larger portfolio.

## PHASE B — every individual replacement is feasible

All nine candidate machines can be fully replaced. `rho* = 1.0` everywhere under
both reactive policies. Damping change at full replacement:

| bus | SCR | rating MVA | `delta zeta` at `rho=1` |
|---|---|---|---|
| 38 | 0.907 | 1684.1 | **+0.0080** |
| 33 | 2.951 | 1174.8 | +0.0058 |
| 35 | 3.207 | 1085.7 | +0.0035 |
| 36 | 2.115 | 1025.2 | +0.0029 |
| 32 | 3.144 | 843.7 | +0.0007 |
| 31 | 3.691 | 836.0 | +0.0004 |
| 30 | 3.269 | 1040.0 | −0.0007 |
| 37 | 2.760 | 970.2 | −0.0025 |
| 34 | 2.316 | 1080.2 | −0.0038 |

The nodal ranking does not predict the outcome: Spearman between `delta zeta`
and SCR is `-0.233` (`p = 0.55`), and against Thevenin magnitude `+0.033`
(`p = 0.93`). The electrically weakest bus, 38 at SCR 0.907, gives the **largest
damping improvement**. With nine points this is absence of evidence, not
evidence of absence, and is reported as such.

Under the matched policy the power-flow solution is **bit-for-bit identical** to
base at every `rho`, which is asserted in `test_ieee39_baseline.py`. The
frozen-operator campaign is therefore exact, not an approximation.

A separate physical finding: under **unity power factor**, partial replacement is
**infeasible at buses 31, 32 and 35** for some intermediate `rho`, while full
replacement is feasible. The surviving fraction of the machine would have to
carry the entire reactive output on a reduced rating. Feasibility in `rho` is not
monotone, so a single `rho*` is the wrong summary and the feasible set is
reported instead.

## PHASE C — minimum destabilizing order is four

Exhaustive census of the full Boolean lattice, 512 portfolios, 0 infeasible, 51 s.

| size | cases | unstable | below damping threshold | worst `zeta_min` |
|---|---|---|---|---|
| 1 | 9 | 0 | 0 | +0.0223 |
| 2 | 36 | 0 | 0 | +0.0224 |
| 3 | 84 | 0 | 0 | +0.0216 |
| **4** | 126 | **10** | 14 | −0.1139 |
| 5 | 126 | 51 | 59 | −0.4906 |
| 6 | 84 | 78 | 79 | −0.9494 |
| 7 | 36 | 36 | 36 | −0.9717 |
| 8 | 9 | 9 | 9 | −0.2803 |
| 9 | 1 | 1 | 1 | −0.1509 |

**`kappa = 4`.** No single, pair or triple replacement is unstable; ten
quadruples are. 42 portfolios are minimal incompatible cores: every proper subset
stable, the portfolio not.

Mechanism classification, automated and run on the spectral abscissa because its
sign is the criterion:

- **12 `A_GENUINE_HIGH_ORDER`** — no lower-order truncation crosses zero and the
  critical mode keeps the same family across the whole subset lattice.
- **30 `C_MODE_SWITCHING`** — refused. The mode attaining the abscissa changes
  family between subsets, so the Moebius terms mix different physical modes and
  no interaction order may be claimed.

The nine size-4 cores classified `A` all have a synchronous-family critical mode
at 0.42-0.57 Hz, the inter-area mode, and **bus 30 appears in all nine**.

The refused cases are not noise. `33+35+36+38` reaches a spectral abscissa of
`+405` at zero frequency with dominant participation in `theta_pll` of three
converters and an eigenvalue condition number of only 141: a collective loss of
PLL synchronism. Detuning the PLL from `wn = 37.4` to `wn = 5.0 rad/s` reduces it
from `+405` to `+131` but does not remove it, so it is not a bandwidth problem.
No subset of three shows it.

**Two distinct failure mechanisms therefore appear at the same order**: an
inter-area electromechanical crossing and a converter synchronization loss. The
classifier separates them instead of averaging them.

## PHASE D — conventional screens do not identify the failures

ROC AUC for predicting instability, at the orders where both classes exist
(0.5 = no information):

| predictor | size 4 | size 5 | size 6 |
|---|---|---|---|
| replaced inertia fraction | 0.869 | 0.785 | 0.846 |
| replaced MW | 0.813 | 0.773 | 0.861 |
| generalized SCR | 0.639 | 0.626 | 0.720 |
| minimum SCR | 0.605 | 0.653 | 0.678 |
| mean SCR | 0.588 | 0.620 | 0.665 |
| maximum MIIF | 0.471 | 0.552 | 0.750 |
| **additive reconstruction** | **0.144** | 0.146 | 0.090 |
| **pairwise reconstruction** | **0.153** | 0.147 | 0.105 |

Screening at the minimum failing order, flagging the ten worst by each predictor
out of 126 portfolios of which ten are unstable:

| predictor | precision@10 | recall@10 | false negatives |
|---|---|---|---|
| replaced inertia fraction | 0.500 | 0.500 | 5 |
| replaced MW | 0.300 | 0.300 | 7 |
| minimum SCR | 0.200 | 0.200 | 8 |
| generalized SCR | 0.100 | 0.100 | 9 |
| mean SCR | 0.000 | 0.000 | 10 |
| maximum MIIF | 0.000 | 0.000 | 10 |

Two results, of very different character.

**First**, the best conventional screen is simply how much inertia is removed,
and it still misses half the failures. The impedance-based screens do
substantially worse; the multi-infeed interaction factor and mean SCR miss every
one.

**Second, and stronger**, the additive and pairwise reconstructions score
**0.14, far below chance**. They are not uninformative, they are *inverted*. A
planner who ran exhaustive single and pair studies and extrapolated would rank
the dangerous portfolios among the safest. The mechanism is visible in the atlas:
buses 33, 35, 36 and 38 each *improve* damping individually, so any additive
model predicts their combination is unusually safe.

The clearest single pair of cases:

| portfolio | min SCR | gSCR | max MIIF | replaced MW | outcome |
|---|---|---|---|---|---|
| `32+35+36+38` | 0.9073 | 0.8939 | 0.7376 | 4638.7 | stable, `-0.2374` |
| `33+35+36+38` | 0.9073 | 0.8942 | 0.7376 | 4969.8 | **unstable, `+405`** |

Every conventional indicator agrees to three digits. One swaps bus 32 for bus 33
and the system loses synchronism.

## What this does and does not establish

Established at this operating point, this converter tuning and `rho = 1`:

- minimum destabilizing replacement order is four, with 42 minimal incompatible
  cores;
- conventional nodal screens have at best 0.5 precision at the failing order,
  and lower-order extrapolation is anti-predictive;
- two distinct failure mechanisms coexist at the same order, and the classifier
  refuses to label 30 of 42 cores rather than forcing an interaction order.

Not established, and not to be written:

- that any of this survives a different operating point, converter tuning or
  reactive policy. No Monte Carlo has been run.
- that a feedback core has been localized. `K(s)`, `Q(s)`, holonomies and winding
  provenance have not yet been applied to these cores; that is PHASE E onward.
- that a targeted repair exists. No repair has been attempted on IEEE-39.
- that the negative controls behave. **This is the most important gap**: without
  them there is no evidence the method is not an "everything is dangerous"
  detector, and the PHASE D result cannot be promoted to a claim.

## Immediate next step

Negative controls (E14) before anything else, then the interaction calculus on
the twelve `A`-classified cores and on the `33+35+36+38` synchronization case,
which are the two candidate flagships and answer different questions.
