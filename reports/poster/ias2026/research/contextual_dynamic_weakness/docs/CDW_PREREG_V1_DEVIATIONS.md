# CDW preregistration V1: deviations and implementation clarifications

The rule (prereg header): every entry is timestamped, and nothing is applied silently.

## Implementation clarifications fixed BEFORE the campaign launched (2026-09-12)

These entries were written after the code pilots (the first task of each phase, plus
E09 with one iteration) and before the orchestrated campaign produced any result
table. No gate value had been computed.

1. **Critical mode.**
   - The rightmost mode is always identified on the TX4 transverse operator. The
     full-A eigenvector used for derivatives is the one nearest to it.
   - Filtering the full-A eigenvalues by |λ| > 1e-6 is not safe: the structural
     Jordan pair splits numerically to about 1e-5. A pilot showed this, and the
     code was fixed before launch.
2. **Stored mode shapes (E1).** The mode-shape store keeps the transverse modes in
   0.1–2.0 Hz with Re ≥ −1 s⁻¹, plus the critical mode (a storage bound). Tracking
   uses these.
3. **E9 sequential QP.** Implementation details not written in the prereg:
   - per-step trust region |d_k| ≤ 0.25 × range_k;
   - active set capped at the 25 largest-α subsets;
   - stopping at α(T) ≤ −ε/2 (single-boundary tuning) or Φ_T ≤ −ε/2 (plan-level
     tuning), with ε = 0.02;
   - if the QP is infeasible, a least-violation step (min t) is taken after the
     Gordan LP.
   - `ka` and `line` gradients are converted from relative to absolute
     coordinates.
4. **E14 certificate.**
   - It is implemented for the TS (device time-scale) reductions only, in the
     port form `T_S(s)` against the reduced `T̃_S(s)`, which have the same
     dimension (78).
   - For GM/POD (Galerkin network reductions) the port matrices have different
     dimensions. The Gohberg–Sigal comparison does not apply as written, so the
     result is NOT APPLICABLE and GM/POD count with zero abstention in GOLD-C.
5. **E13 reduced α.**
   - GM/POD keep the device state space, so the full-model R_x and w are used in
     the transverse operator.
   - TS uses R_x and w restricted to the retained states. This is exact for TSa
     and TSb; for TSc it is approximate, because PLL angle states are
     eliminated.
6. **E12 Q3 control points.** The midpoint between g = 0 and E1, plus the 8
   midpoints between consecutive frozen events, gives 9 controls.
7. **E16 limiting mode.** This is the rightmost oscillatory transverse mode in
   0.1–2.0 Hz (the same set in which modal energy is compared). It is not the
   global α when α is real or out of band.
8. **E4 Dport.** This is an independent T-form cross-check:
   `−p^H ∂_a T q / p^H ∂_s T q`, with ∂_a T from re-solved Jacobians at fixed s*.
9. **E2 TX4 draws.** They are used only if the regenerated factors match the frozen
   PCV05 `factors` column exactly (checked in code).

## Robustness fix applied AFTER a task crash (2026-09-12, during the campaign)

10. **E09 Gordan-LP solver robustness.** Task (T2, topology, plan) crashed with
    `cvxpy.error.SolverError: Solver 'CLARABEL' failed` inside the Gordan
    conflict-witness LP (`_analysis.gordan`) and the least-violation QP
    (`E09_design.qp_step`), both triggered only when the primary SQP step is
    infeasible. This is an ill-conditioning/numerical-robustness bug in the
    orchestration code, not a change to any hypothesis, threshold, definition,
    baseline or the design problem itself (objective, constraints and stopping
    rule in E09 are unchanged).
    - Fix: both solves now try CLARABEL, then SCS, then OSQP in order; if every
      solver fails, a dependency-free FISTA (accelerated projected gradient)
      solve of the identical convex problem is used, and a NaN/Inf guard was
      added to `gordan()`. If even the least-violation QP fails on every
      solver, the SQP step records `qp: "solver_failed"` and takes no step for
      that iteration (rather than crashing the task); `any_solver_failed` is
      reported per case.
    - Added `test_gordan_projected_gradient_fallback_matches_solver` (8 tests
      now pass).
    - The single failed checkpoint (`raw/E09/911586952a5fc2f1.json`, task
      T2/topology/plan) was deleted and recomputed with the fixed code; every
      other completed E09 checkpoint was left untouched.
