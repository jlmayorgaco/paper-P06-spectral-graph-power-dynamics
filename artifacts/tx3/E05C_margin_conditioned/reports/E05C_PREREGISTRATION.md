# E05C margin-conditioned connected pole externalities preregistration

**Freeze date:** 2026-08-27. **Git SHA before freeze:** `78a34059a2b0641efd1db123e36789cdce7ff6ca`.

E05C preserves C1/C2 supported, C3a rejected, E05B D1/D2/D3 supported, and original cycle/SCC C3
unassessed. It tests a new margin-conditioned claim using all 15 E05B-reproducible coalitions, one
development-locked mode per coalition, uniform constant-power-factor load stress, and eight new Sobol
seeds. A1--A8 and the old `|Delta_S zeta| >= 0.0025` threshold are immutable.

Stress endpoints are calibrated with the EMPTY coalition only. No action vertex may be evaluated until
`E05C_BASELINE_STRESS_LIMITS.parquet` is frozen and committed. The nine tau levels, redispatch factors,
tracking/refinement rules, failure statuses, C3b-A amplification gate, and C3b-B two routes are fixed in
the accompanying JSON files. E05C uses fully re-equilibrated ANDES only; no ParaEMT or E06 execution is
permitted. Even a positive C3b-B result stops for independent review.
