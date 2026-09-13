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
