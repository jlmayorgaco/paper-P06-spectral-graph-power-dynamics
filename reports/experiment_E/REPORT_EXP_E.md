# Experiment E — collective IEEE-39 replacement design (in progress)

**Status: INCOMPLETE.** The analytical design has produced a strictly margin-feasible provisional point. It has not passed the mandatory global/branch optimality audit, independent detailed-model validation, or time-domain checks. No final candidate has been frozen for that validation.

## Frozen domain and actual dispatch

The design domain was frozen before any ExpE stability calculation in `experiments/bnd_expE/configs/DESIGN_DOMAIN_FROZEN.toml` (SHA-256 `feb02715bf15a9930d4678c4221535837c2597834fac89e2feed07febe1aae9b`). All ten generator components at buses 30–39 are included. Each has its own `rho`, `Kp`, and `Ki`; the independent gain rectangle is `[0.25,4]` times each nominal gain. No four-bus restriction or beta linkage is used.

Actual SG dispatch was reconstructed from each frozen machine state and its analytical electrical port. It totals **5402.761089978776 MW**. Subtracting nominal load CSV setpoints from the net bus injection would give 6140.811071711622 MW, but that does not describe the initialized detailed-model state at the two generator-plus-load buses. The frozen state and nodal current balance imply 0.107802 MW of load at bus 31 (CSV setpoint 9.2 MW) and 375.042216 MW at bus 39 (CSV setpoint 1104 MW). The resulting actual SG injections are approximately 511.719 and 271.042 MW. Tables E01 and E03 retain both initialized and nominal quantities. This discrepancy must be addressed explicitly in the detailed-model validation; it must not be hidden by treating setpoints as realized dispatch.

## Analytical model checks

The fixed branch network and initialized constant-impedance ZIP loads were reconstructed from frozen voltages and current balance. The resulting all-SG 114-state Schur-reduced spectrum matches `C0_Ared` from ExpC within a maximum bidirectional spectral distance of **2.05e-12**. The independent single-bus replacements at 30, 33, 35, and 37 match the archived ExpC spectra within **1.76e-11**. These checks use the complete finite spectra and distinguish endpoint removal from a zero-rated dormant device.

The independent-gain GFL port has an exact one-scalar Woodbury representation because the `Kp` and `Ki` rows share the same phase-error covector. Across 50 device/frequency checks, the median/p95 relative port errors are **2.64e-16 / 4.81e-16**. The multi-bus determinant-lemma closure has 20 coordinates, against 78 network algebraic coordinates; its maximum sampled log-determinant/phase errors are **1.43e-13 / 8.76e-15**. The design source under `src/bnd_design_e/` has zero imports, calls, or textual references to the detailed simulation package.

## Full replacement and branch structure

At `rho_i=1` for all ten generators, the correctly removed-device model has 90 dynamic states. Its two eigenvalues closest to zero split numerically by about `1e-6` while the analytical uniform-rotation null residual is below `8e-19`. The same near-double-zero structure appears with two independent patterns of `Kp,Ki`. The all-GFL endpoint therefore cannot be certified to satisfy the preregistered `0.05 s^-1` margin after removing only the uniform-angle gauge; the remaining near-zero frequency mode is structurally relevant. An exact all-gain infeasibility proof has not yet been completed.

The graph-metric direction from the all-SG point is essentially uniform and reaches the margin at `rho_i=0.9151377506`. A predictor/corrector with independent gain compensation reached 5017.91 MW before its active-set branch lost correction authority. Explicit one-SG branch enumeration then found a much better nominal-gain branch: nine full replacements and `rho_38=0.9941453795`, or **5397.901755 MW**. Gain compensation on that branch reached the strictly feasible provisional point below. This demonstrates why the earlier continuation stop was not an optimality certificate.

## Provisional analytical point

The best saved point satisfying the strict numerical margin is step 7 of the bus-38 gain-continuation branch:

| Quantity | Value |
| --- | ---: |
| Actual initial SG dispatch | 5402.761089979 MW |
| GFL replacement | **5401.571878097 MW** |
| SG retained | **1.189211882 MW** |
| Replacement fraction | **0.9997798881** |
| Spectral abscissa excluding the angle gauge | **-0.050000000461 s^-1** |
| Maximum nodal KCL residual at frozen operating point | **1.03e-13 pu** |

The full per-generator vector, including `rho_i`, `Kp_i`, `Ki_i`, MW and Mvar shares, is in `tables/TABLE_E14_provisional_per_generator_design.csv`. All buses except 38 are at `rho=1`; its SG component remains physically present. The gain bounds are satisfied. The GFL/SG current split preserves each generator's initialized P/Q injection, and the ZIP loads remain at their reconstructed initialized admittances, so the external voltage/angle operating point is fixed for this analytical model.

The ExpE artifact audit in `test/bnd_expE/runtests.jl` passed **33/33** checks. It covers the frozen-domain hash, actual dispatch and load reconstruction, archived spectral matches, Woodbury/closure residuals, all-GFL neutrality, independent gain bounds, provisional MW split, and nodal KCL.

## Open gates

- The active-set continuation and branch tree are incomplete. A local KKT/SOSC audit, off-branch placement search, and global upper-bound gap remain open. The provisional point is **not** called optimal.
- Closure pole sensitivity validation at the final point, nonnormality, self-energy/pairwise pathway, PLL routing, scaling/cost, and convergence audits remain open.
- The final candidate has not been frozen with its SHA. No detailed-model eigenvalue or time-domain validation has been run for ExpE. The two generator-plus-load buses require special care to preserve the load components and their initialized operating point during that validation.
- `EXP_E_STATUS=INCOMPLETE`; no PASS claim follows from the provisional analytical result.
