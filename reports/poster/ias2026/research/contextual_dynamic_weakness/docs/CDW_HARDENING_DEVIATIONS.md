# CDW hardening — deviation log

Every deviation from `docs/CDW_HARDENING_PREREG_V1.md` is recorded here with a
timestamp and a reason **before** the affected result is read. Pure execution
fixes (crash, path, resume) are recorded too, with an explicit statement that no
definition, threshold, seed or test set changed.

(no entries at preregistration)

## 2026-09-12T21:40 — execution fix (no definition changed)

- **Relaunch.** The first background launch of `CDWH_MASTER_RUN.py` was terminated
  externally at 325/384 H01 tasks, with no traceback and exit code 127 from the
  shell wrapper.
  - The run was relaunched as a detached process.
  - Completed task checkpoints were reused (resume-safe, deterministic per task).
- **Import order.** Analysis modules now import `_hinfra` first, so that the CDW
  module path is set. This was an import-order crash only.
- **Bound-check tolerance.** H03's numerical verification of the Proposition-D bound
  used an absolute tolerance of 1e-12.
  - It flagged 54 m = 1 pairs at fast-mode magnitudes (abs(D) ≈ 2e4 s⁻¹), where
    witness − D/m = −3.6e-12 is floating-point roundoff.
  - The check now uses a relative tolerance of 1e-9·max(1, abs(D)).
  - The witness definition, thresholds and every reported quantity are unchanged.

## 2026-09-12T21:50 — H17 Q3/Q4 refinements, recorded BEFORE any H18 matrix case was run

**Context.** Qualification of ALT-WECC found the following:

- **Q1 PASS.** The ANDES all-SG base reproduces the internal α_⊥ at all 8 H18 policies,
  with max abs(Δα) = 1.7e-7 s⁻¹ and identical statuses.
- **Q3 case.** The Q3 case was run by mistake at H02, not at P4 as the prereg says.
  H02/H4 is also an H18 matrix case, so that one ALT datum has been seen:
  - under the literal Q4 rule it is BOUNDARY_OR_UNRESOLVED;
  - its rightmost non-structural eigenvalue is +0.257 ± j2.96 (0.472 Hz);
  - two limiter flags are raised.

  No other H18 case was run.

**Diagnosis.** The literal Q4 rule removes the zero eigenvalues of identically zero state
rows, then the structural pair. That leaves 8 extra near-zero eigenvalues in every
converted case:
- the REECB1 Q-PI integrators (PIQ_xi, abs(λ) ≈ 2e-10);
- the V-PI integrators (PIV_xi, abs(λ) ≈ 1e-6).

Under the frozen TX3-GFL-0.1 flags (QFLAG = 0), the PI outputs of these integrators feed
no equation. Their columns in the *reduced* state matrix A_s = T_f⁻¹(f_x − f_y g_y⁻¹ g_x)
are therefore structurally zero. The 4 REPCA1 s2_xi states (Kp = Ki = 0) are the
zero-row dead states already covered by the prereg.

**Refinement 1 (Q4, exact linear algebra; not a tuning).**

- **Rule.** A state is *structurally decoupled* if its row, or its column, of A_s is zero
  within 1e-12·‖A_s‖, diagonal included. The test is applied iteratively.
- **Why this is exact.** Each such state contributes exactly one zero eigenvalue. The
  remaining eigenvalues are exactly those of the submatrix with the state deleted
  (block-triangular permutation).
- **Order of removal.** Remove those zeros first, then the structural pair, then apply
  the unchanged rule: all others ≥ 1e-2, else BOUNDARY_OR_UNRESOLVED.
- **What it changes.** This extends the prereg's "zero row" dead-state rule to its
  column dual. It is fixed by the frozen TX3 parameter set and is identical for every
  case. Nothing in the model or its parameters changes.

**Refinement 2 (limit flags).**

- **Exempt flags** (`LIMIT_ACTIVE` does not count them). The TX3 freeze has no limiter
  check for ANDES, so the following flags are now defined as structurally inactive:
  - REGCP1.HVG: high-voltage reactive-current gain with Khv = 0. Its output is 0 on
    both sides of the limit; zl = 1 simply means v < Volim.
  - REPCA1.feHL and REPCA1.s5: the plant active-power/frequency path, disabled by
    Fflag = 0 and PLflag = 0.
  - REPCA1 Q-path limiters: the plant Q PI is disabled by Kp = Ki = 0.
  - REECB1 PIQ and PIV anti-windup limiters: their outputs are unused when QFLAG = 0.
- **Every other flag counts.** Any zl/zu flag on REGCP1, REECB1 or REPCA1 outside that
  list keeps a case out of the linear verdicts under `LIMIT_ACTIVE`.

**Q3.** Q3 is now also run at P4 (D01) on H4, as the prereg specifies. The H02 run is
kept and reported.

## 2026-09-12T22:05 — H17 refinement 1b, recorded BEFORE any H18 matrix case was run

After refinement 1, the REECB1 V-PI integrators (PIV_xi, one per converter) remain. In A_s
their off-diagonal row **and** column are exactly zero, and their diagonal is +1.2e-6 s⁻¹.
This is the anti-windup tracking term of a PI whose input is multiplied by SWQ_s1 = 0 under
QFLAG = 0.

- **What such a state is.** A fully isolated state is its own 1×1 diagonal block. Its
  eigenvalue is its diagonal, and it neither influences nor is influenced by any physical
  variable.
- **Refinement 1b.** A state whose off-diagonal row and column of A_s are both zero, and whose
  diagonal satisfies abs(a_kk) ≤ 1e-3 s⁻¹, is removed as structurally decoupled. Its
  eigenvalue is recorded in the case output (`isolated_eigs`).
- **What is not removed.** Isolated states with abs(a_kk) > 1e-3 are kept in the spectrum.

Everything else is unchanged:
- the Q1–Q3 results;
- the structural-pair rule;
- the verdict band.

No model or parameter is changed.

**Q3 (for the record).**
- Descriptor QZ vs reduced-A_s spectrum: 1.0e-13 relative at P4.
- ANDES EIG.mu: 3.6e-15.
- Initialization residual: 9.4e-14.

## 2026-09-12T22:10 — execution fix: H06 worker-pool crash (no definition changed)

- **What failed.** 260 of the 440 H06 tasks failed with "A process in the process pool was
  terminated abruptly" (BrokenProcessPool). The failure coincided with a memory-heavy
  analysis running in parallel, so the workers were probably killed for lack of memory.
- **Fix.** The failed-task checkpoints (records with ok = false and no results) were deleted,
  and the 260 tasks were rerun. The code and inputs are unchanged, and each task is
  deterministic.
- **Prevention.** Memory-heavy analyses are no longer run while a large pool is active.

## 2026-09-12T22:20 — H6 eligibility applied literally (interpretation note, no rule changed)

- **The rule.** The prereg declares a GOLD-B condition eligible if "the target solves
  and the derivative engine returns (R0 ≤ 1e-8)", where R0 = max abs(SPR residual) at
  the TX4 solve_case equilibrium.
- **What the data show.** R0 lies between 1.2e-12 and 2.2e-7, with a median of
  4.8e-8. That is the ordinary tolerance of the frozen equilibrium solver; the old E34
  never aggregated R0. The threshold was therefore stricter than the solver the
  campaign uses.
- **Literal application.** Applied as written, only 8 of the 64 H4 primary conditions
  are eligible: 6 draws and 2 policies, in 6 clusters.
- **Primary result.** GOLD-B hardening is decided on these 8 conditions (`GB_prereg_*`).
- **Sensitivity analysis.** The result over all 64 solved conditions (`GB_primary_*` in
  the JSON) is reported as a labelled sensitivity analysis.
- **Why the analysis is still valid.** The implementation-validity check (IFT
  derivative against the small-step finite difference) passes for 100 % of the 4048
  link pairs, so a residual of order 1e-7 does not measurably affect the derivatives.
- **Holm family.** F_conf uses the T2 p-value of the literal (8-condition) set.

## 2026-09-12T23:08:13-0500 — CORRECTION of earlier timestamps, plus post-hoc disclosures (system clock; anchored to commits)

### Earlier timestamps

The clock times in the five earlier entries were written by hand and **not** taken from the system clock. Internal reviewer 2 found that they contradict the git and file times, and it is correct. The times reconstructed from evidence are:

| entry as labelled | true time | evidence |
|---|---|---|
| "21:40" execution fix (relaunch, import order, bound tolerance) | about 21:12–21:20 | H03 run output; the entry is in commit `4de412d2` (21:32:13) |
| "21:50" H17 refinements 1 and 2 | before 21:32:13 | the entry is in commit `4de412d2`, and the ALT matrix case files are dated 21:36:01–21:46:09 |
| "22:05" H17 refinement 1b | about 21:33–21:35 | every ALT matrix case file (21:36–21:46) already contains the `isolated_eigs` field that 1b introduced, so 1b was implemented before the matrix ran; the log text was committed in `eb57d907` (22:06:09) |
| "22:10" H06 worker-pool crash | about 21:36 | `logs/hardening/master_run_H06b.out` starts at 21:36:53 |
| "22:20" H6 literal eligibility | about 22:01 | written after the H07 gate had been computed |

**What the evidence supports.** The ordering claims about the ALT matrix ("recorded before any H18 matrix case was run") are supported by the commit time of `4de412d2` (refinements 1 and 2) and by the `isolated_eigs` field (refinement 1b). The clock times in those entries are wrong.

### H6 eligibility ordering

The literal R0 <= 1e-8 filter was applied **after** the 64-condition GOLD-B result had been computed and read.
- The H07 code committed before the results computed the gate on all 64 conditions.
- The change goes in the conservative direction (a smaller primary set).
- Both sets pass.

### Post-hoc analyses (exploratory, reviewer-requested; no preregistered verdict changes)

- **`regret_anatomy_new`** (H04). The share of material regrets whose ranked choice leaves the EM band. It was added after results commit `eb57d907`.
- **`H31_revision.py` / `H31_revision2.py`.** These cover:
  - the level-D statistics;
  - the tau curve;
  - the regret-optimal and stability-screened rankings;
  - portfolio-conditioned gSCR/SCR sequencing screens;
  - EM-clean one-step squares;
  - the leave-one-cluster-out fixed branch list;
  - gap2 sensitivity;
  - direction counts;
  - the mixing variance terms;
  - within-type topology correlations;
  - the distance between optimal orders;
  - D_tot against doublings and outages;
  - the ALT EM-tracked reversals.
- **`H31_explore.py`.** Participation factors of the fast real modes; fractional-replacement sweeps; eigenvalue continuation of tracked modes; timing.
- **ALT-WECC variant with REECB1 QFLAG = 1.** Library voltage control switched on, every gain unchanged. It lies outside the TX3 validated freeze; results are in `results/hardening/alt/cases_qflag1`.

### Implementation defect found by reviewers (disclosed; the preregistered verdicts are reported as computed)

1. **H18 "tested policies" rule.** `H18_compare.py` applied the rule (ALT base STABLE) only to tests A and B, although the prereg defines it for all tests.
   - Recomputed on the 7 tested policies: C = 0.06 (MODEL-SPECIFIC), C2 = 0.203 (TRANSFERS; 3 of 7 below 0.20), D = 0.71 (PARTIAL, not TRANSFERS), E = 1.00 on 29 pairs.
   - The corrected values are the ones reported in the paper.
2. **H17 refinement 1b and the pinned pole.** Refinement 1b removed isolated states only when abs(a_kk) <= 1e-3.
   - The isolated anti-windup tracking state of REPCA1's disabled plant active-power integrator (s5_xi, eigenvalue −0.100 s⁻¹) was therefore kept.
   - It pins the global α⊥ of most stable ALT portfolios at three policies (13/15 portfolios at H02 and H07), which makes the preregistered global-α reversal test A uninformative there.
   - Test A is reported as computed (0/7), together with this defect and the exploratory EM-tracked result: stabilizing EM marginals at 3/7 policies and nested reversals at 1/7.
