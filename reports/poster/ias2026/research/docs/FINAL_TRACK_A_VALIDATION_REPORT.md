# Final Track-A validation report

Overnight run `outputs/ias2026/final_validation_overnight_20260910T003225/`.
All **12 experiments completed**. Twelve manifests, each carrying its config, git
commit, seed, environment, worker count, accepted and rejected counts with
reasons, raw and summary tables, exact figure source data and a verdict.

Nothing frozen was modified. The flagship (30+33+35+37), the v2C modal-family
observable, the repair definitions, the controller models and the sampling
distributions are as they were. Three code corrections were made and each is
recorded in its experiment's note: an E30 classifier that decided "mechanism
disappears" on a damping comparison rather than on the absence of an instability;
an E34 target test at `1e-9` that reported an optimizer converging *onto* its
constraint as a failure; and an E40 perturbation ensemble that read the closure
eigenvalue off the unperturbed operator.

## Executive table

| # | Claim | Status | Evidence | Independent? | Robust? | Poster-safe? | TPWRS-safe? |
|---|---|---|---|---|---|---|---|
| **C1** | All proper subsets stable, exact portfolio unstable | **SUPPORTED** | E30 (Q-matched), E35 P1 0.947 [0.932, 0.960], E36 1500/1500 | **no** — E31 G2 MIXED | operating point yes, machine data **no** (E37 0.721) | **yes, with dispatch and machine-data caveats** | yes, stated conditionally |
| **C2** | Irreducible order-4 modal-family effect | **SUPPORTED** | E31-audit `μ₄ = +0.433`, orders 1–3 stable; E35 0.871 given unstable; E36 1500/1500 | no | **conditional**: 100 % to controller, 42 % to machine data, 21 % across operating points | **yes, with the conditionality stated** | yes, as a conditional result |
| **C3** | Survives re-equilibrated engineering dispatch | **MIXED** | E30: matched → order 4; unity PF → order **3**, worse; voltage regulation → **stable** | n/a | no — the order is dispatch-dependent | **only with the dispatch named** | yes, as a dispatch-dependence result |
| **C4** | Independent ANDES reproduction | **FAILED (partial)** | base frequency 4.1 %, damping sign agrees; **ordering not reproduced**, shift-sign agreement 73 % < 80 % | yes, by construction | n/a | **base mode only** | **no** — this is the top gap |
| **C5** | Nonlinear TDS consistency | **VALIDATED** | E32: pencil frequency error median 0.00066 Hz, worst 0.0029 Hz; damping sign 25/25 | no (same code base) | 4 disturbances × 7 configurations | **yes** | yes |
| **C6** | Port model reproduces the full model | **VALIDATED** | E41: eigenvalue error ≤ 7.3e-10, `σ_min` ≤ 2.7e-13, RHP counts identical, 18 cases | no | whole subset lattice | **yes** | yes |
| **C7** | Closure margin tracks the full-system boundary | **VALIDATED** | E33: boundary located to 0.0025 in load (median 0 refined), Spearman +0.861, p = 1.5e-86; E35 p = 7.3e-38 | no | 290 continuation points + 1004 MC samples | **yes** | yes |
| **C8** | `f_port` tracks `f_IA` | **VALIDATED** | E33: median 0.0049 Hz near the boundary, 95th pct 0.0093 Hz | no | as above | **yes** | yes |
| **C9** | Controller repair retains all PV | **MIXED** | E34: converter-only reaches **all four margins**, 0 MW lost, 0 MVA. E35: the *frozen* retune works in **117/240** held out | no | **no** as a fixed retune; untested as an adaptive one | **only as "a retune exists at this operating point"** | no, not as a general claim |
| **C10** | Synchronous condenser mitigation | **VALIDATED** | E35 240/240; E34 cheapest is 166–270 MVA at bus 30 alone; E37 686/710 | no | yes across operating points, 96.6 % across machine data | **yes** | yes |
| **C11** | Repair robustness | **MIXED** | condenser robust (E35 100 %, E37 96.6 %); fixed converter retune **not** (E35 48.8 %) | no | split | **yes, if the split is stated** | yes, as a split result |
| **C12** | Bus-30 collective-enabler hypothesis | **NOT SUPPORTED** | unadjusted: 12/12 inter-area cores, Fisher OR 13.2 p = 0.005, matched pairs 41–0 p = 9.1e-13. Adjusted OR 4.09 [0.14, 124], **p = 0.418** | no | underpowered | **descriptive only** | **no** |
| **C13** | Conventional baseline comparison | **SUPPORTED, restated** | E39: lower-order reconstruction AUC 0.09–0.15, and at the declared threshold **0 true positives at every size**; removed inertia AUC 0.78–0.87 | no | 3 sizes, permutation null 0.50002 | **yes, in the corrected form** | yes |

Statuses use only VALIDATED, SUPPORTED, MIXED, FAILED and NON-EXECUTABLE.
NON-EXECUTABLE remains attached to **H5B only**, from v2B, and nothing tonight
changes that.

## Fail-fast gates

| gate | outcome | what actually happened |
|---|---|---|
| **G1** re-equilibrated policies destroy the mechanism | **QUALIFIED, not failed** | unity power factor leaves the portfolio unstable and *worse*, so the mechanism is not destroyed; but the interaction order collapses to 3, and under voltage regulation there is no instability at all |
| **G2** independent implementation cannot reproduce the ordering | **MIXED, effectively failed on the ordering** | ANDES reproduces the base inter-area mode (4.1 % in frequency, damping sign) and **not** the machine-removal ordering; the two disagree on the *sign* of the effect |
| **G3** corrected order-4 effect fails in held-out Monte Carlo | **MIXED** | 0.208 [0.183, 0.235] unconditionally, 0.871 [0.822, 0.911] given the portfolio is unstable |
| **G4** closure does not separate boundary from controls on fresh data | **PASS** | p = 7.3e-38, median `m₄` 0.021 near against 0.187 far |
| **G5** repair works only on discovery operating points | **MIXED** | the fixed retune works on 49 % of held-out unstable samples, all mild ones and none of the severe; the condenser works on all |
| **G6** port reduction fails to reproduce the family | **PASS** | 7.3e-10 eigenvalue agreement, zero invisible band modes |

## The three findings that most change the story

**1. The interaction order depends on the reactive dispatch (E30).** Under
reactive matching the effect is irreducibly fourth order. Under unity power
factor the triple 30+33+35 is already unstable, so it is third order and the
portfolio is worse. Under voltage regulation with either frozen reactive
capability the flagship is stable and there is no failure to explain. The
headline claim must always carry its dispatch.

**2. Uncertainty separates cleanly into three axes (E35, E36, E37).** Converter
tuning: 1500 of 1500 draws reproduce the whole pattern. Machine parameters: 418
of 1000. Operating point: 209 of 1004. The effect is not controller-caused and
that is now measured, not asserted; but it is a property of this machine data,
which is a real limit on generality.

**3. The fixed converter retune is not a robust repair (E35).** It restores
stability in 117 of 240 held-out unstable samples — all 80 mild cases, none of
the 80 severe ones. The 25 % synchronous condenser restores all 240. E34 shows a
converter-only solution exists at every declared margin *at the nominal point*;
whether one exists at every operating point was not tested.

## What was corroborated from a second direction

- The cheapest synchronous mitigation is a condenser at **bus 30 alone** (E34),
  arrived at independently of E38's finding that bus 30 sits in 12 of 12 genuine
  inter-area cores.
- The closure minimum locates the stability boundary to one continuation step
  (E33) and the port operator reproduces the full spectrum to 1e-10 (E41): the
  same conclusion from an 8×8 operator and a 110-state eigenvalue problem.
- The linear prediction is confirmed by nonlinear integration to 0.003 Hz (E32).

## Companion documents

- `POSTER_FINAL_SAFE_CLAIMS.md` — what may go on the poster, with wording
- `TPWRS_REMAINING_GAPS.md` — what a Transactions submission still needs
- `FAILED_FINAL_VALIDATION.md` — the genuine falsifications and downgrades
- `CLAIMS.md` — O48–O53 and N13–N17 added tonight

**The poster has not been edited.**
