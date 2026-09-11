# Phase I — complete report and gates (BC00–BC02, NL00–NL02, unit correction)

Branches: `ias2026/binary-certification-v1` (BC) and
`ias2026/nonlinear-portfolio-v2` (NL, created from BC). Nothing was pushed.
Nothing frozen (F1–F12, G1–G6, v2C, poster) was rewritten.

Nothing beyond Phase I has been launched:

- no BC03 or later block;
- no NL03 or later block;
- no threshold sweep, certificate or planning run.

Detailed sources:

- `docs/binary_certification_v1/BC00_DATA_AND_ORIGIN_AUDIT.md` (Spanish);
- `docs/nonlinear_portfolio_v2/NL_PHASE_I_REPORT.md` (Spanish);
- `docs/PHASE_I_UNIT_CORRECTION.md`.

Tests: 31/31. That is `test_bc_quotient` 10, `test_bc_adapter` 5,
`test_nl01_model` 12 and `test_unit_contract` 4.

Gate legend:

- **PASS**: the purpose is met.
- **CONDITIONAL**: met with a declared restriction that the next campaign must
  carry.
- **FAIL**: the tested hypothesis or representation does not hold; it is
  reported, not repaired silently.

## Gate summary

| block | purpose | verdict | blocks next phase? |
|---|---|---|---|
| BC00-A units | MW vs MVA in every quantity | **PASS after correction**: the old wording failed; data contract installed | no; it would have blocked any MW-based planning (BC07) |
| BC00-B origin / zero modes | exact symmetry, zero ledger, four-state classifier | **PASS**, with a **CONDITIONAL** meaning of "stable" | no, but the campaign must declare the stability notion (§B) |
| BC00-C normalization | invariant norms, closure is not a radius | **PASS** | no |
| BC01 physical H | H and kappa without the `abs(s) <= 1e-3` disc | **CONDITIONAL**: Kundur and IEEE-68 PASS; 212 IEEE-39 fast-path labels UNKNOWN | no, unless those low-`g` points enter the campaign (then evaluate them directly) |
| BC01 Kundur crossings | treat the 194 aperiodic crossings at `s = 0` | original port **FAIL** (24/194); relocated port **PASS** 194/194, post hoc | no; the relocated port must be preregistered for BC03+ |
| BC02 representation | which binary representation exists | **CONDITIONAL**: descriptor affine PASS; reduced A affine FAIL; port locality PASS; C1 padding FAIL on IEEE-68 | shapes BC03–BC05 (§C) |
| BC02b Boolean degree | exact Boolean dependence and degree | **PASS** (structure established) | no |
| NL00 audit | provenance, units, rotation, slack, ANDES scope | **PASS** | no |
| NL01 L1 = L0 | matrix rewrite equal to the running code | **PASS** | no |
| NL02 derivatives | chart derivatives, `H_r`, `B`, `C`, `D` | **PASS** on IEEE-39/Kundur; **CONDITIONAL** on IEEE-68 (about 1e-3) | IEEE-68 NL certificates, until rescaled |
| UC unit correction | contract, audit, reruns | **PASS**; several claims changed (§A) | no |

---

## BC00 — data, origin and units

**Purpose.** Make every reported quantity physically labelled. Establish the
structure of the origin: symmetry, zeros, dead states. Freeze a
four-state classifier with error bounds instead of a magnitude cutoff.

**Result and evidence.**

- **Units: see §A.**
- **Symmetry.** `R_x = 1` on `delta` and `theta_pll` with `R_z = jV` gives
  `A R_x = 0`. The residual is ≤ 2.5e-10 on 269 560 linearizations, and on the
  nonlinear residuals under a finite rotation (0.37 rad) the error is 3e-13.
- **Jordan partner.** When every machine has `D = 0`, which holds in all three
  benchmarks, `A w = omega_B R_x` to 3.9e-12. The origin carries a size-2
  Jordan block: T1 removes the rotation, and the second zero (the neutral
  frequency mode) is **physical**.
- **Dead states.** Rows that are identically zero are removed exactly and
  ledgered.
- **Classifier.**
  - Balanced coordinates, with the error matrix transformed by the same
    balance.
  - `eps` = ‖A(h) − A(2h)‖ + backward error.
  - `d_axis = min_w sigma_min(A − iwI) > 10 eps` means the RHP count is exact.
    This is a numerical criterion, not an interval certificate.
- **Synthetic checks.** The 404 `theory_checks` rerun unchanged (synthetic
  only). `test_bc_quotient`: a physical +1e-7 pole stays UNSTABLE; an exact
  zero is BOUNDARY.
- **Normalization.** The metric is transformed by congruence under a base
  change (NL00: 1e-15). The −1 distance is not a non-normal robustness radius
  (`Q_L` example).

**Claims changed.** Every "stable" in F, G and BC means *stable modulo the
declared neutral frequency mode*. Unit claims: §A.

**Blocks next phase?** No. The stability notion is a campaign-level decision
(§B).

## BC01 — physical H and the Kundur crossings

**Purpose.** Recompute H and kappa without the disc. Test the 194 G1 aperiodic
crossings.

**Result and evidence.**

| benchmark | points | subset cases | H_phys = H_G1 | unresolved |
|---|---|---|---|---|
| Kundur (K12A + K12B) | 5 917 | 47 336 | 5 917 / 5 917 | 1 subset |
| IEEE-68 | 121 | 1 936 | 121 / 121 | 0 |
| IEEE-39 (F8, lengua, F7 with `g <= 0.005`, 600 random) | 13 735 | 220 288 | 13 513 / 13 725 | 428 points |

- **IEEE-39.** All 212 label changes are at non-exact fast-path points with
  `g <= 0.005`. Each is exactly one `BOUNDARY_OR_UNRESOLVED` subset: an
  **oscillatory** pair at 2.4–4.4 rad/s whose `d_axis` lies below 10× the
  error of the affine assembly. It is **not** the disc: no eigenvalue of
  `A_qq` has `abs(lambda) < 1e-2`. There are 0 changes at exact points.
- **Kundur crossings.**
  - Mechanism: every crossing is the leaky Q/V integrator, going from −0.05 at
    `g = 0` through the origin at `g` of order 3e-4.
  - Disc offset: the physical boundary sits 3.6e-6 to 2.0e-5 below G1.
  - Original port: 24/194 by the preregistered criterion; it failed because of
    the structural double zero.
  - Relocated port: `T##(0)` flips sign in 194/194 while the device block
    never does, and the result is independent of `beta`.
  - Raw pencil count: it gives 1–2 spurious "RHP" roots in 166/194 on the
    stable side.

**Claims changed.** 212 G1 IEEE-39 labels → UNKNOWN. G1 Kundur boundaries
shift by ≤ 2e-5 in `g`, which is no change at the precision quoted.

**Blocks next phase?** No. Two requirements carry over:

- the relocated port was added after preregistration and must be preregistered
  before BC03 uses it;
- the unresolved IEEE-39 points need direct-path Jacobians if they enter the
  campaign.

## BC02 and BC02b — binary representation, exact Boolean dependence

**Purpose.** Decide which representation of the 0/1 family exists:

- **A**: affinity in control parameters for fixed `S`;
- **B**: port locality;
- **C**: a common realization affine in `delta`.

Then establish the exact polynomial degree of each object.

**Result and evidence.** 245 vertices in BC02, and three families × 4 values of
`s` in BC02b.

| object | dependence on `delta` (exact statement) | evidence |
|---|---|---|
| descriptor Jacobian `J(delta)` (C1/C2) | **degree 1, exactly affine** | Möbius coefficients of order ≥ 2: ≤ 3.7e-16 relative (BC02b); affinity 1e-19 on 245 vertices; common equilibrium `f = g = 0` to 1.5e-12 |
| device factor `h_S(s) = det(sI − f_x)` | **a product of per-device factors**: `log h_S` is additive, with no interaction | order ≥ 2 Möbius coefficients of `log h`: ≤ 3.4e-13 (mod 2πi) |
| port ratio `r_S(s) = det T_S / det T_0 = det(I + M_SS)`, `M = blkdiag(dY_i) E^T T_0^{-1} E` | **multilinear of degree m**. Its Möbius coefficient at `T` is **exactly** the sum of the principal minors of `M` that touch every block of `T` | vertex identity 1.3e-13; Möbius = minor sums 1.3e-13. Coefficients at the full order are non-zero in all families: IEEE-39 order 4 is 0.012–0.062 of `r_0`; Kundur order 3 is 0.18–1.31; IEEE-68 order 4 is 0.028–0.21 |
| reduced `A(delta) = f_x − f_z g_z^{-1} g_x` | **rational, not affine**; full-degree Möbius spectrum | order-k over order-1 norm: IEEE-39 0.24 / 0.14 / 0.32; Kundur 1.03 / 1.17; IEEE-68 1.7e-6 / 1.9e-7 / 2.9e-8 (a norm dominated by stiff entries, not evidence of affinity) |
| B: port locality in the frozen model | `T_S = T_0 + sum(T_i − T_0)` to the power-flow tolerance | ≤ 8e-8. Outside the candidate's two channels the increment is 8e-8 to 1.5e-6 of the inside, located at the slack bus and matching `abs(z_S − z_0) <= 1.5e-7` |
| `det P = h det T` | exact | ≤ 2e-12 |
| ledger | `N_h = 0`: no device block has an RHP pole | IEEE-68 device blocks carry marginal zeros (`Re` ≤ 2e-13) |
| C1 ghost padding | Hurwitz in IEEE-39/Kundur (−0.05); **marginal in IEEE-68** | ghost `Re` up to +7e-14; only the C2 filler is admissible there |

**Consequence.**

    det P_S(s) = h_0(s) · prod_{i in S} kappa_i(s) · det T_0(s) · det(I + M_SS(s))

All Boolean interaction of the characteristic function lives in
`det(I + M_SS)`, equivalently `det(I + Q_SS)` of (20). The order-k interaction
is a sum of principal minors of `M` spanning k actions. No lower-degree
representation is exact on these benchmarks.

**Claims changed.** None of the earlier claims. The theorem applicability
changes: T4 as stated (affine `A(delta)`) is NOT_APPLICABLE; certificates must
be posed on the affine descriptor or on `det(I + M_SS)`.

**Blocks next phase?** It shapes BC03–BC05 but does not block them.

## NL00 — provenance and physics audit

**Purpose.** Identify the running model and its hashes. Audit:

- units and `gamma = Sn/Sbase`;
- the 4270.7 figure;
- the near-zero poles;
- slack versus infinite source;
- the ANDES scope.

**Result and evidence.**

- **Hashes.** Source hashes are in `results/NL/NL00/NL00_hashes.json`, and the
  snapshot manifest in `configs/nonlinear_portfolio_v2/NL_SNAPSHOT_MANIFEST.json`
  (74 files).
- **Base conversion.** A system-base re-expression of every machine gives the
  same current and derivatives to ≤ 2.7e-15.
- **4270.7.** It is the MVA sum 1040 + 1174.8 + 1085.7 + 970.2.
- **Rotation.** The finite rotation holds on the nonlinear residuals to 3e-13.
- **Slack.** No fixed sources exist; every slack bus is a dynamic machine.
- **ANDES.** Static reconciliation only, with the PSS off.

**PASS.** **Blocks?** No.

## NL01 — L1 matrix model

**Purpose.** A matrix rewrite L1 that equals L0, plus the v2 software
contract.

**Result and evidence.** 12/12 tests:

| check | result |
|---|---|
| Ybus of L1 vs L0, three networks | < 1e-12 relative |
| KCL off equilibrium | < 1e-12 |
| per-branch losses with taps | < 1e-10 |
| synthetic phase shifters (`Yft != Ytf`) | 40 cases pass |
| stacked vs interleaved realification | one explicit permutation |
| physical spectrum after bus reordering | < 1e-6 |
| inputs | alter only their declared rows |

The `PhasorModel` contract is implemented, with `active_set`, `guard_values`
and `reset_map` declared ABSENT or NOT_APPLICABLE in L0.

**PASS.** **Blocks?** No.

## NL02 — manifold derivatives (item E)

**Purpose.** Compute the first and second derivatives of the algebraic chart
and the reduced field, and verify them against direct derivatives with `psi`
re-solved at every evaluation.

**Result and evidence.** 11 points, including off-equilibrium points; 9 step
sizes. First-derivative errors:

| quantity | max relative error |
|---|---|
| `J_psi` | 2.8e-8 |
| `A` | 1.4e-8 |
| `B` | 1.0e-7 |
| `C` | 2.8e-8 |
| `D` | 8.9e-8 |

`H_r` (Hessian-vector products, no complex step), best independent agreement:

| benchmark | agreement |
|---|---|
| IEEE-39 | 1.3e-6 to 5.7e-6 |
| Kundur | 1.9e-7 to 5.0e-7 |
| IEEE-68 | 1.1e-3 and 1.7e-4 (stiff scaling, `abs(f)` about 2e5; step sensitivity up to 0.10) |

Chart failures occur only at steps 0.1 and 0.03 (outside the chart) and are
reported. A rerun after reformatting reproduced every output byte for byte.

**Verdict.** PASS on IEEE-39/Kundur; CONDITIONAL on IEEE-68. Everything is
**sampled curvature**: no bound on a domain.

**Blocks?** Certified NL03/NL04 bounds need analytic or interval derivatives
anyway. IEEE-68 NL work needs rescaled coordinates first.

---

## A. Unit audit

See `docs/PHASE_I_UNIT_CORRECTION.md`. In short:

- **Contract.** `replaced_mw` = Σ Sn = **MVA**. New canonical quantities:
  `replaced_sn_mva`, measured `replaced_pg_mw`, `replaced_q_mvar`, and
  `replaced_pmax_mw` only where documented (IEEE-39). There is no PV nameplate.
- **Flagship.**

  | quantity | value |
  |---|---|
  | active dispatch displaced | 2096.6 MW |
  | reactive dispatch displaced | 420.6 Mvar |
  | machine rating | 4270.7 MVA |
  | documented Pmax | 2983 MW |
  | converter rating | 4270.7 MVA |

- **Restorations.** Restoring SG30 or SG37 gives up 436.1 or 321.5 MW of
  converter dispatch (1040 or 970.2 MVA of machine rating).
- **Condensers.** 1067.7 MVA (RD, 0.25 Sn at each retired bus) and 166–270 MVA
  (bus 30 alone).
- **Baselines.**
  - The old "replaced MW" predictor was the retired **rating**.
  - True removed dispatch: AUC 0.57 (p = 0.46) / 0.64 / 0.73, against 0.81 /
    0.77 / 0.86.
  - True MW is **not** a useful baseline at the minimum failing order.
  - The conclusion about nodal/simple metrics is qualitatively unchanged
    except the "megawatts are informative" clause, which is invalidated.
- **Matched analyses** redone on Pg:
  - E14 N4: 107 stable in the window; all top 25 stable.
  - E18: closure 0.156–0.863 vs 1.6e-7; cycle ratio 0.72.
  - E38: OR 0.053, p = 0.18; still NOT SUPPORTED.
- **Planning.** Every "PV MW" axis is INVALID as MW. The re-scored E34 front has
  the **same** 10 points, with M5 costs of 65.4 and 109.0 MW. The E23
  restoration choice (bus 37) is unchanged.
- **Core mathematics.** UNAFFECTED, verified by code path: the quantity is
  read-only and post-solve.

## B. Symmetry and zero-mode audit

- Exact rotation symmetry holds, both in the tangent and in the nonlinear
  model.
- An exact Jordan partner (the neutral frequency mode) exists wherever `D = 0`
  and there is no governor, which is the source data of all three benchmarks.
- Consequence: under the note's definition (8) (`alpha_q < 0`), no portfolio
  is strictly stable. The family would be empty, and T6 trivial.
- Every stability statement is therefore made on `A_qq`: **stability relative
  to the centre-of-inertia frequency**, with the neutral mode declared and
  verified in 100 % of cases.

**Decision needed for the campaign.** Keep this relative notion, or add a
governor or damping as a versioned L2 model. The frozen model must not be
"fixed". The zero ledger is exact: rotation, neutral mode, dead states, and
the physical integrator crossings in Kundur.

## C. Exact Boolean dependence / polynomial degree

See BC02b above:

| object | dependence |
|---|---|
| descriptor `J` | exactly degree 1 |
| `h_S` | exactly multiplicative |
| `det(I + M_SS)` | multilinear, with Möbius coefficients equal to the principal-minor sums of `M`; degree m, non-zero at full order in all three benchmarks |
| reduced `A` | rational, full-degree |

Implication: a "low-degree" Boolean SOS certificate (T5 hierarchy) is a
**hypothesis to test**, not implied by the structure. The exact degree is m.
The decay of the minor sums with order in IEEE-39 (0.53 → 0.20 → 0.05 → 0.012
at `s = 0.5+2j`) is what a low-degree certificate would exploit.

## D. Nonlinear model-to-code traceability

- `docs/nonlinear_portfolio_v2/EQUATION_CODE_TRACEABILITY.csv` has 33 rows
  (equation, code location, units, dimensions, tests, status). By status:

  | status | rows |
  |---|---|
  | verified (V, including V with limits absent) | 27 |
  | T | 1 |
  | T+V | 1 |
  | MODEL_CHANGED (IEEE-39 AVR) | 1 |
  | ABSENT (limits, DC link, governor) | 3 |

- The full equations are in `NONLINEAR_POWER_SYSTEM_MODEL.md`.
- The v2 code package (625 checks) was **not supplied**, so the v2 checks
  could not be reproduced.

## E. First- and second-derivative checks

See NL02 above.

## F. Discrepancies with the frozen IEEE-39 / ANDES model

1. **Excitation.** The IEEE-39 AVR merges KA and TE into a first-order block
   (harmonized IEEEX1): `MODEL_CHANGED`, declared since F1.
2. **PSS.** IEEE-39 uses a power-input washout plus lag (the source's non-zero
   blocks), without limits.
3. **ANDES evidence.**
   - F1 is a static reconciliation only, with the PSS off.
   - E31 reproduces the base inter-area mode (frequency within 4.1 %).
   - E31 does **not** reproduce the machine-removal ordering (C4 FAILED,
     partial).
   - Unchanged by Phase I.
4. **Source data.**
   - `D = 0` on every machine (the origin of the neutral mode).
   - The slack machine (bus 39) has `M = 100` on its own base and represents
     the interconnection; it is excluded from replacement.
   - Pmax comes from the same ANDES workbook (sha `9c2048dc…`).
   - The Kundur Pmax values are placeholders and are not used.
5. **Absent from L0.** Governors, limits, anti-windup, modulation delay, DC
   link and ZIP loads. None was added to the frozen model.

## G. New limitations found in Phase I

- **Strict stability.** It is impossible without a governor or damping: all
  results are relative stability (§B).
- **IEEE-39 fast path.** Below `g = 0.005`, 428 points are unresolved
  (oscillatory pairs within the assembly error). The fast path is **not**
  certifiable there.
- **`d_axis`.** It is a numerical criterion, not an interval certificate, and
  `operator_error_bound` is UNKNOWN.
- **Original port.** It is unusable at `s = 0`. The relocated port is a
  post-hoc method.
- **Reduced A.** It is not affine in `delta`, so T4 as stated is NOT_APPLICABLE.
- **IEEE-68.**
  - C1 ghosts and some device blocks are marginal.
  - The derivative accuracy is only about 1e-3.
  - Relative-norm measures there are dominated by stiff entries.
- **Units.**
  - No PV nameplate and no available PV power are modelled: `P_ref` is the
    dispatch, with unlimited DC.
  - Active-power planning can only use dispatch or documented SG Pmax.
  - Kundur and IEEE-68 have no usable Pmax.
- **E14 N6.** The branch-MAC check has not been recomputed for the Pg-matched
  controls (open).
- **NL02.** Only sampled curvature; there is no automatic differentiation in
  the environment.

## Implications for the final one-run campaign (for review, not executed)

1. Every planning objective uses `replaced_pg_mw` (or documented Pmax) for MW
   and `Sn` for MVA, on separate axes.
2. The stability notion (relative, on `A_qq`, or an L2 governor variant) must
   be fixed before any certificate is run.
3. Certificates are posed on the affine descriptor or on `det(I + M_SS)`,
   never on an assumed-affine reduced `A`. IEEE-68 uses the C2 filler and
   scaled coordinates.
4. The relocated port and the classifier (`d_axis`, `known_zero`) are
   preregistered as used here.
5. IEEE-39 points with `g <= 0.005` enter only through direct-path Jacobians.
