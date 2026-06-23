# Phase-0 IEEE39 IBR Bridge Artifacts

This artifact records what is currently validated before any final full-ANDES IBR estimator claim.

## Mix60 replacement map
- GFL generators [9, 7, 5] at buses [34, 36, 38].
- GFM generators [8, 6, 4] at buses [33, 35, 37].
- Retained synchronous generator buses [30, 31, 32, 39].
- Figure: `outputs\ieee39_phase0_bridge\fig_ieee39_mix60_ibr_map.pdf`.

## Phase-0 gates
- PASS: G0.1 baseline ANDES load/PFlow/EIG -- 39 buses, 46 lines, 160 eigenvalues
- PASS-DIAGNOSTIC: G0.2 packaged-case positive-mode diagnostic -- Removing IEEEST/PSS removes all positive-real modes
- PASS: G0.3 SG-to-GFL/GFM case generation -- 5/5 derived cases solved PFlow and EIG
- PARTIAL: G0.4 generic IBR small-signal stability -- 1/5 cases stable under generic gains: ieee39_ibr_mix60
- BLOCKED: G0.5 physical ANDES-to-(L,M,D) IBR bridge -- Strict bridge audit BLOCKED: critical frequency error 12.4% and damping-ratio error 100.0%
- BLOCKED: G0.6 final estimator validation against full ANDES IBR poles -- Requires G0.4 and G0.5 before estimator errors can be claimed.

## Current reduced-model evidence
- IEEE39-derived surrogate estimator errors against the surrogate full QEP: diagonal median 0.004684, second-order median 0.000321, reduced-QEP median 0.000369, adaptive median 0.000306.
- Braess-like surrogate line reversals: 35/46 tested lines.

## Blocking item
The physical bridge from the full calibrated ANDES IBR Jacobian to a scalar `(L,M,D)` model is not yet validated. The manuscript must therefore present these as benchmark-construction and surrogate-validation artifacts, not as final full-ANDES IBR estimator validation.
