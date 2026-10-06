# D4 exact graph-damping operator gate

**Status: BLOCKED.** After explicitly projecting out the physical global-rotation direction (normalized gauge residual 6.14e-18), the hidden block Δ_cc(0) is singular under the preregistered retention of SG rotor-angle coordinates; its condition estimate is 3.10e18. Its Schur complement and derivative at s=0 are therefore undefined for this partition. The no-gauge calculation is retained only as `TABLE_D04_PRE_GAUGE_DIAGNOSTIC.csv` and is not valid evidence.

No D_G, L_D, or L_τ result is accepted. A future derivation must retain the additional zero-frequency controller/device subspace or choose a nonsingular physically justified partition, then re-run gauge and pole consistency checks. No Laplacian or delay-dressed-damping claim follows from the previous numerical table.

Frozen delay design remains unchanged. All 112 requested coordinates are explicitly marked blocked in `baseline_reproduction/TABLE_D04_GRAPH_DAMPING_OPERATORS.csv`.
