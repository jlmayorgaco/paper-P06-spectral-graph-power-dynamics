# J0: one-sided spectrum at the all-GFL endpoint

Run from the repository root:

```text
julia --project=. --startup-file=no experiments/bnd_expJ0/run_j0.jl
```

The runner uses the frozen IEEE-39 network, nominal PLL gains, and one retained-SG fraction at a time on buses 30–39. It removes the exact angular gauge and fits `log10(abs(lambda))` against `log10(epsilon)` for `epsilon = 1e-8, ..., 1e-2`. The local fit uses `1e-8` through `1e-5`. The structural endpoint root is zero; its computed quotient eigenvalue is `2.34e-10` from roundoff. Results are in `J0_ENDPOINT.csv` and `J0_SCALING.csv`.

## Decision

The near-zero physical branch is **linear** in the tested one-sided directions. The ten local exponents range from **0.99923 to 1.00067**; full-range fits range from **0.98120 to 1.02999**. No tested direction supports a square-root law for this branch.

At the all-GFL endpoint, the 90-state matrix has a defective gauge/physical-zero pair (ExpG). After gauge removal, the 89-state quotient has one numerical zero: its smallest singular value is `1.27e-11`, the second is `0.19915`, and the next eigenvalue by distance from zero is `9.10208` away. The gauge residual is `7.06e-19`. Thus the full-matrix Jordan label does not establish a defective physical zero *after* quotienting.

There is a separate spectral obstacle. In `mixed_jacobian`, an SG state block is present for every positive `epsilon` and absent at exactly zero. The dynamic dimension jumps from **90 to 102** for buses 30–38 and from **90 to 96** for bus 39. As `epsilon -> 0+`, an SG block receives the grid voltage but its current feedback vanishes, so its open-loop poles enter the one-sided spectrum. For buses 30, 33, 35, and 37, the SG open-loop rightmost poles are respectively `+0.11319862`, `+0.05992935`, `+0.05077506`, and `+0.05059171` per second. At `epsilon=1e-8`, the full quotient's rightmost poles agree with these values within `2.9e-8` per second. The near-zero linear branch is not the critical branch in those four directions.

A Puiseux expansion of the exact 90-state endpoint cannot predict all poles of the 96- or 102-state one-sided model without a common-dimensional embedding. The exact angular gauge remains present throughout the path; its Jordan coupling does not by itself imply a square-root split of the physical root.

The one-sided family therefore does not justify using a single endpoint derivative as a stability design law. This is a concrete mechanism that the ExpG scalar surrogate did not represent; the present test does not prove it is the sole cause of ExpG's numerical error. The next regular-point test should track the actual rightmost mode and compare its first-order prediction against the full closure, including mode switches and the added SG states, before using an authority/dual LP for design.

## Diagnostic note on ExpH

`BNDDesignH.jordan_audit` currently multiplies its nominal zero tolerance by the largest singular value of the unscaled quotient matrix. At the nominal endpoint that yields `0.357`, above nine nonzero singular values, and its table reports quotient nullity `10` and `UNRESOLVED`. Those fields are invalid as a zero count for this matrix. This J0 test uses the observed spectral gap and reports the explicit absolute threshold `1e-7` per second, which selects one zero. ExpH's scaling values agree with the independent J0 rerun.
