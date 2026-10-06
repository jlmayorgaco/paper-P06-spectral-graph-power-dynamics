# Experiment K — support-aware SG→GFL co-design

## Outcome

**Status: FAIL / NOT CERTIFIED.** The exact analytical fixed-support
model and all 1024 physical SG supports were evaluated. A deterministic
analytical search produced frozen nominal and normalized robust candidates.
The nominal search did **not** recover the better feasible ExpE point, and
the required KKT/global certification is absent. Independent PowerDynamics
and nonlinear TDS results are recorded in Tables K12–K13. Both frozen
candidates fail the independent spectral gate. No
global or robust-optimum claim follows from this experiment.

The work stayed on the existing `research/expH-jordan-mixed-codesign` branch
at the user's direction. No push was made. Existing A–J0 artifacts were not
edited. ExpE and ExpG summary text was seen before the candidate hashes were
created during initial repository inspection; candidate files explicitly
record this **blind-integrity failure**. The candidates were nevertheless
frozen before the later exact-spectrum comparison and any PowerDynamics call.

## Frozen IEEE-39 domain

- Discovered SG buses: 30–39; initial SG power: **5402.761089978776 MW**.
- The static network includes the actual initialized ZIP admittances at
  buses 31 and 39; they differ from CSV nominal load setpoints.
- Every evaluated model removes absent SG state blocks and projects the
  common-angle gauge before computing the full finite spectrum.
- Nominal spectral requirement: `alpha<=-0.05 s^-1`.
- Independent PLL bounds: `Kp,Ki` each from 0.25 to 4 times nominal.

## Architecture sweep and branch evidence

`TABLE_K01_support_catalog.csv` has 1024 distinct masks, dimensions, gauge
residuals, and exact midpoint spectra. 1010 supports are margin feasible at
the midpoint `epsilon_i=0.5` for present SGs. The separate uniform-seed scan
found a feasible point on 1023 supports. **Zero supports were proved
infeasible.** Only 64 supports received local coordinate correction; 959
retain feasible seeds without continuous branch exhaustion. This is an
exhaustive architecture enumeration, not an exhaustive continuous solution.

The branch catalog records active poles at the midpoint using state-energy
participation. The one-sided SG re-entry table independently estimates the
open-loop limits and first slopes. The single-bus nominal-gain scan found two
sampled feasible components for bus 30, confirming that a monotonicity
assumption would be unsafe. Mesh sampling and nearest-pole matching do not
certify all branch components or their remainders.

The analytical derivative formulas are validated at an interior all-bus
point. The maximum complex absolute central-difference error in Table K03 is
`1.46e-7`, above the requested `1e-8` target. The associated test suite
passes 66 tests, including endpoint dimensions, gauge projection, initialized
load integrity, ExpE/ExpG analytical regression, and frozen-hash checks.

## Nominal analytical candidate

The frozen nominal candidate is
[`Z_K_NOMINAL_FINAL.toml`](Z_K_NOMINAL_FINAL.toml), SHA-256
`64989ecfcc93ef471a4f6a5f9fe6be382c634e88a6561faf847dd58871ea65b9`.

| Quantity | Value |
|---|---:|
| Present SG support | 38, 39 |
| Retained SG fractions | `epsilon_38=0.0007208842980`, `epsilon_39=0.01` |
| Converted fractions | `rho_38=0.999279115702`, `rho_39=0.99`; all others 1 |
| Gains | all `Kp=7.853981633974485`, `Ki=986.9604401089358` |
| Retained SG | **3.308756126036 MW** |
| GFL | **5399.452333853 MW**, 99.93876% of initial SG MW |
| Complete-spectrum abscissa | `-0.050001000514425416 s^-1` |
| Active pole | real `-0.050001000514425416 s^-1` |
| Direct normalized radius | `beta_star=1.3464152166764433e-9` |

This point is feasible in the exact analytical quotient. It is **not** a
KKT-certified local minimum. At the reported nonnegative one-mode dual
diagnostic, the interior survival residuals are approximately `+13.53 MW`
at bus 38 and `-218.89 MW` at bus 39. No SOSC/critical-cone certificate is
available.

The post-hash ExpE provisional point retains **1.189211882165 MW** and has
`alpha=-0.05000000218509599 s^-1` under the same ExpK quotient model. It
converts **2.119544244 MW more** than the frozen ExpK nominal candidate.
Thus the answer to “does nominal ExpK recover or beat ExpE?” is **no**, and
ExpE is analytically feasible to the evaluated precision. Its approximately
`2.2e-9 s^-1` margin slack is numerically narrow and would need independent
validation before a physical claim.

## Direct normalized robust frontier

The primary robust metric is the peak singular value of the shifted
full-block resolvent, not ExpG's modal-residue bound. The frozen frontier is
in `TABLE_K09_robust_frontier.csv`. It is the best among 12 explicitly
evaluated candidates and is not a robust Pareto optimum.

| Required normalized beta | Retained SG MW | GFL MW | beta_star | Support |
|---:|---:|---:|---:|---|
| 0 | 3.308756 | 5399.452334 | 1.3464e-9 | 38,39 |
| 1e-6 | 6.690000 | 5396.071090 | 7.4961e-5 | 34,38 |
| ExpG beta=1.6991207e-6 | 6.690000 | 5396.071090 | 7.4961e-5 | 34,38 |
| 3e-6 | 6.690000 | 5396.071090 | 7.4961e-5 | 34,38 |
| 1e-5 | 6.690000 | 5396.071090 | 7.4961e-5 | 34,38 |
| 3e-5 | 6.690000 | 5396.071090 | 7.4961e-5 | 34,38 |
| 1e-4 | 6.950000 | 5395.811090 | 1.1228e-4 | 36,38 |

At ExpG's normalized beta, the frozen comparison candidate
[`Z_K_ROBUST_beta_1p6991206999182038em6.toml`](Z_K_ROBUST_beta_1p6991206999182038em6.toml)
has support 34,38, `epsilon_34=epsilon_38=0.005`, nominal PLL gains,
`alpha=-0.06495339224893061 s^-1`, and a direct radius
`7.496093929784582e-5`. Its analytical GFL share is **3815.081246229 MW**
larger than ExpG's 1580.98984375 MW at the same normalized beta. The
candidate-pool limitation and independent validation result govern any
stronger comparison.

## Transient capacity and BND mechanism

The analytical unit-disturbance scan covers initialized load buses and uses
retained-SG COI frequency. The nominal candidate's worst location is bus 8:
peak RoCoF is `0.0201313 Hz/s/MW` and peak frequency excursion is
`0.2734116 Hz/MW`. With ExpG's declared `0.5 Hz/s` and `0.5 Hz` limits,
the corresponding capacities are **24.837 MW** and **1.829 MW**. The robust
comparison candidate also has worst bus 8, with capacities **94.056 MW**
and **6.876 MW**. These are post-design linear-model metrics. No transient
constraint was inserted into a KKT corrector, and no hard GFL current-limit
claim is made.

At the nominal active pole, the port Schur pivot is bus 38. Its direct and
self-energy sensitivities to `epsilon_38` are approximately `+91268.36` and
`-91325.76 s^-1`, leaving a `-57.40 s^-1` total. This strong cancellation
explains why a simple nodal authority ranking is unreliable. The graph
operator was not substituted for exact closure; no `K=h(L_G)` claim is made.
The pairwise Schur pathway ledger gives cancellation ratio `1.10752`; the
static `L_G/L_B` commutator ratio is `0.0687113`. Both are descriptive
coordinates, not an optimization law.

## Independent PowerDynamics result

The post-hash validator rebuilt both candidates, retained ZIP load components
at buses 31 and 39, applied the frozen independent PLL gains, and initialized
the equilibrium. The ZIP reference parameters at those buses were reconstructed
from the frozen initialized admittances. Residuals were below `2.1e-12`, but
the complete finite spectra disagreed decisively:

| Candidate | Analytical alpha | PowerDynamics alpha | Gate |
|---|---:|---:|---|
| Nominal | `-0.0500010` | `+0.0051060 s^-1` | FAIL |
| Robust at ExpG beta | `-0.0649534` | `+0.0062152 s^-1` | FAIL |

The extra numerical gauge counts (two and three in the PD audit) and the
abscissa gaps require further model reconciliation. A component-level P/Q
sharing audit of the rebuilt mixed buses has not yet closed, so the source
of disagreement is unresolved. The gaps are not treated as small tolerances.
Under the preregistered acceptance rule, candidate
validation **fails**. Nonlinear pulses are a separate falsification check and
cannot overturn the spectral failure.

All eight requested small-pulse runs completed: `0.0005` and `0.001` of
frozen load at buses 16 and 8 for both candidates. Their frequency and
voltage amplitudes scale closely between pulse sizes (largest relative
two-amplitude scaling error `5.1e-5`). This confirms local pulse scaling over
the 12-second window, **not** asymptotic stability. The corresponding frozen
analytical pulse responses differ from the nonlinear COI frequency traces by
up to **38.34%** for the nominal candidate and **31.49%** for the robust
candidate, relative to the nonlinear peaks. Tables K13, K16–K18 and the
linear/nonlinear figure preserve these falsification results. Settling and
local bus frequency/P/Q trajectories were not established; TDS acceptance is
therefore **NO**.

## Required gates still open or failed

- `EXP_K_NOMINAL_PASS = NO` and `EXP_K_ROBUST_PASS(beta) = NO` for the
  frozen comparison candidates, because independent spectral validation fails.
- Full support/branch-complete predictor–corrector, multimode KKT, LICQ,
  SOSC, and perturbation audit: **not completed**.
- Global nominal certificate: **NO**; the feasible ExpE point is a concrete
  better point.
- Robust optimization and direct robust-constraint derivative KKT audit:
  **not completed**; the table is a sampled frontier.
- One-sided branch and single-bus interval *rigorous* certification:
  **not completed**.
- Independent PowerDynamics: **FAIL** for both frozen candidates; see K12.
- Nonlinear TDS: eight finite, well-scaling small pulses, but large
  analytical/nonlinear mismatch and no stable settling certificate; **FAIL**.
- Generic solver falsification benchmark: not run. The ExpE comparison already
  falsifies global optimality of this nominal candidate.

## Reproduction

From the repository root, run in order:

```text
julia --project=. --startup-file=no experiments/bnd_expK/run_catalog_K.jl
julia --project=. --startup-file=no experiments/bnd_expK/run_design_K.jl
julia --project=. --startup-file=no experiments/bnd_expK/finalize_analytic_K.jl
julia --project=. --startup-file=no experiments/bnd_expK/explain_graph_K.jl
julia --project=. --startup-file=no experiments/bnd_expK/validate_powerdynamics_K.jl
julia --project=. --startup-file=no experiments/bnd_expK/compare_tds_K.jl
python experiments/bnd_expK/render_figures_K.py
python experiments/bnd_expK/print_summary_K.py
julia --project=. --startup-file=no test/bnd_expK/runtests.jl
```

`run_design_K.jl` refuses to overwrite frozen candidates. A rerun requires
an isolated output directory or deliberate archival of this run. Do not edit
the frozen candidate values after the hash is written.
