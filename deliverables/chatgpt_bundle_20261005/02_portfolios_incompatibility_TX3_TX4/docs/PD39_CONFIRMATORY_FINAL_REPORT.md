# PD39 confirmatory campaign

## 1. Executive verdict

The frozen discovery audit confirms an exact 255+1 pattern, but its meaning is
more precise than “V8 is incompatible.” There are 256 portfolios and nine
discovery conditions. All 255 proper portfolios pass the 0.05 s^-1 requirement
in all nine conditions. V8 is nominally stable, with
`alpha=-0.0014447780 s^-1`, but it is truly unstable (`alpha>0`) in the four
high-PLL conditions and therefore fails all nine robustness-margin tests.

Thus, for the frozen discovery benchmark:

* `H_0 = {V8}`: minimal true-stability blocker.
* `H_0.05 = {V8}`: minimal blocker of the engineering requirement.

The structured-radius result is promising but not confirmatory of a general
hidden-fragility law. The preregistered 128-direction search was reduced to 16
directions for runtime, and the `missing_38` predecessor stalled before its
first direction. The smallest V8 boundary found was `rho_0=0.00190598`; this
is a protocol-limited boundary find, not a global radius. `rho_0.05(V8)=0`
because V8 is already below the engineering target at nominal conditions.

The controller-only repair failed within the frozen bounds. Four single-line
repairs pass the nine discovery conditions at their minimum reinforcement, but
none passes all 24 holdout conditions. The nonlinear TDS rerun is valid and
supports the relative ordering best 7/8 → worst 7/8 → V8, while the fresh
holdout planner comparison does not show a reliable advantage for P3 over P1
or P2. The defensible conclusion is therefore: an auditable PD39 warning
case, not yet a validated universal co-design method.

## 2. Discovery audit

The discovery files were read without rerunning or changing them:

* `portfolio_scenario_results.csv`: 256 × 9 = 2304 rows.
* `link_outage_results.csv`: 46 line outages.
* `summary.toml`, `PD39_RESULTS_REPORT.md`, scripts, and manifests.

The eight replaceable SG buses are 30, 32, 33, 34, 35, 36, 37, and 38. The
nine scenarios are `nominal` plus the eight corners formed by PLL scale,
filter/reactance scale, and current-control scale at 0.8 or 1.2:

`pll0.8_xf0.8_cc0.8`, `pll0.8_xf0.8_cc1.2`,
`pll0.8_xf1.2_cc0.8`, `pll0.8_xf1.2_cc1.2`,
`pll1.2_xf0.8_cc0.8`, `pll1.2_xf0.8_cc1.2`,
`pll1.2_xf1.2_cc0.8`, and `pll1.2_xf1.2_cc1.2`.

The discovery “robust margin” is `dynamic_margin=-max_real`, where
`max_real` is the rightmost real part of the computed spectrum for one
portfolio and scenario. It is not a structured perturbation radius. The
criterion is `dynamic_margin >= 0.05`, equivalently `alpha <= -0.05 s^-1`.

The nine and only nine failures are V8 under the nine scenarios. All 255
proper portfolios pass all nine. The four high-PLL cases are true instability;
the five remaining V8 cases are stable but insufficiently damped.

## 3. Exact V8 audit

| Scenario | alpha (s^-1) | m_alpha | f_critical (Hz) | damping ratio | classification |
|---|---:|---:|---:|---:|---|
| nominal | -0.0014447780 | 0.0014447780 | 0.344642 | 0.000667 | stable, low margin |
| pll0.8_xf0.8_cc0.8 | -0.0464474856 | 0.0464474856 | 0.369479 | 0.020003 | stable, below target |
| pll0.8_xf0.8_cc1.2 | -0.0446749762 | 0.0446750 | 0.369624 | 0.019233 | stable, below target |
| pll0.8_xf1.2_cc0.8 | -0.0470012624 | 0.0470013 | 0.369414 | 0.020245 | stable, below target |
| pll0.8_xf1.2_cc1.2 | -0.0454001542 | 0.0454002 | 0.369548 | 0.019549 | stable, below target |
| pll1.2_xf0.8_cc0.8 | 0.1114144062 | -0.1114144 | 0.318614 | 0.055568 | TRUE instability |
| pll1.2_xf0.8_cc1.2 | 0.1140950119 | -0.1140950 | 0.319047 | 0.056824 | TRUE instability |
| pll1.2_xf1.2_cc0.8 | 0.1104433509 | -0.1104434 | 0.318453 | 0.055113 | TRUE instability |
| pll1.2_xf1.2_cc1.2 | 0.1129068947 | -0.1129069 | 0.318849 | 0.056269 | TRUE instability |

The V8 critical mode is a complex electromechanical/control mode. Its largest
participations are machine speed/angle at bus 39 and AVR/machine states at
bus 31; the nominal low-PLL cases also show a smaller bus-38 GFL PLL
participation. The mode is not a clean converter-only pole. The same mode
family is tracked in the rerun, but the reduced Jacobian is poorly conditioned
(approximately 10^18–10^22) and the interpretation must remain qualified.

## 4. Eight 7/8 predecessors

| Missing SG | worst alpha (s^-1) | worst margin | scenario | critical frequency (Hz) | damping ratio |
|---:|---:|---:|---|---:|---:|
| 30 | -0.126228984 | 0.126228984 | pll0.8_xf0.8_cc0.8 | 0.028083 | 0.581825 |
| 32 | -0.138022315 | 0.138022315 | pll0.8_xf0.8_cc1.2 | 0 | 1 |
| 33 | -0.138295567 | 0.138295567 | pll1.2_xf0.8_cc1.2 | 0 | 1 |
| 34 | -0.085834702 | 0.085834702 | pll0.8_xf0.8_cc1.2 | 0.446512 | 0.030581 |
| 35 | -0.136356176 | 0.136356176 | pll0.8_xf0.8_cc0.8 | 0.026096 | 0.639410 |
| 36 | -0.138057659 | 0.138057659 | pll0.8_xf1.2_cc1.2 | 0 | 1 |
| 37 | -0.138325419 | 0.138325419 | pll1.2_xf0.8_cc1.2 | 0 | 1 |
| 38 | -0.138019716 | 0.138019716 | pll0.8_xf1.2_cc0.8 | 0 | 1 |

The modal assurance calculation used common 39-bus voltage observables. The
worst predecessor-to-V8 MAC values are: missing 30, 0.6938; 32, 0.0088; 33,
0.4481; 34, 0.3903; 35, 0.6938; 36, 0.0183; 37, 0.4999; and 38, 0.0400.
This mixed MAC pattern argues against assuming a single universally tracked
mode. The V8 worst case itself remains physically interpretable as a
machine/control oscillatory mode, but the very high conditioning numbers mean
that a reduced, gauge-aware numerical audit is still needed before claiming a
mechanism independent of implementation details.

## 5. Outage and numerical pathology audit

The 11 discovery outage failures are all `graph_disconnected`, not PF
divergence, dynamic initialization failure, voltage infeasibility, or spectrum
failure. They are lines 5 (2–30), 14 (6–31), 20 (10–32), 27 (16–19), 32
(19–20), 33 (19–33), 34 (20–34), 37 (22–35), 39 (23–36), 41 (25–37), and 46
(29–38). The other 35 outages re-solve with finite PF/state, fixed-point
residual, voltage within 0.90–1.10 pu, and stable spectrum.

The V8/predecessor audit found finite equilibria and spectra for all 81 modal
cases after the classifier fix. It found no evidence of a disabled-state or
duplicate-pole explanation in the recorded modal rows. Nevertheless, the
smallest singular values reach approximately 4.7e-16 and eigenvector
condition numbers are about 4.1e4–7.8e4. The responsible classification is
“physical electromechanical/control mode, numerically ill-conditioned and in
need of gauge-aware replication,” not “cleanly proven converter mechanism.”

## 6. Conventional baselines

The baseline table includes remaining SG MW/MVA/inertia, min static strength,
Thevenin/effective-resistance proxies, voltage, electrical-distance-derived
quantities where meaningful, and nominal/discovery margins. SCR and gSCR are
marked unavailable because this model/campaign has no frozen fault-current
model and does not establish the assumptions needed for a valid generalized
SCR calculation. No fake gSCR was used.

For V8, converted candidate dispatch is 4620 MW, IBR MW penetration is
0.9189255, remaining SG dispatch is about 407.611 MW, remaining SG MVA is
10700, and remaining inertia is about 53030.3 in the model units. Static
strength alone does not uniquely identify the V8 cliff. The discovery static
versus dynamic bus ranking had Spearman rho 0.1190; the link ranking had rho
-0.0506. These are descriptive finite-benchmark comparisons, not population
inferences. Bus 37’s singleton result is consistent with favorable dispatch
and rating as well as location; the data do not isolate those causes.

## 7. Structured dynamic radius

The frozen primary box used controller ±10%, load P/Q ±5%, IBR P ±5%, and line
coordinates ±10%. The executed search used seed 39025 and 16 maximin-LHS
directions, followed by bisection and exact equilibrium/spectrum evaluation.
The following are smallest clean boundary finds; jump crossings are excluded
from the clean count.

| Selection | rho_0 | rho_0.05 | status |
|---|---:|---:|---|
| V8 | 0.00190598 | 0 | clean V8 rho_0 find; already below target for rho_0.05 |
| missing 30 | — | — | jump crossing only; unresolved |
| missing 32 | 0.790983 | 0.754769 | clean direction found |
| missing 33 | — | — | no boundary in tested directions |
| missing 34 | 0.200128 | 0.098158 | clean direction found |
| missing 35 | — | — | jump crossing only; unresolved |
| missing 36 | 0.783359 | 0.686153 | clean direction found |
| missing 37 | 0.711267 | 0.730944 | clean direction found |
| missing 38 | — | — | stalled before first direction; not estimable |

At every clean boundary, one side is feasible and one side violates the target;
the exact side alphas and residuals are recorded in
`PD39_STRUCTURED_RADIUS_SEARCH.csv`. The radius search does not prove global
optimality. No matched-alpha representative portfolios were radius-evaluated,
so H3 was not passed or failed. A radius correlation with nominal alpha is
therefore not promoted as a scientific claim.

## 8. Holdout and hidden-fragility selection

The 24 conditions are deterministic maximin-LHS conditions C01–C24, seed
39024, generated inside the frozen primary box without stability screening.
The core holdout had three numerical exceptions (`InternalLinearSolveFailed`)
in missing 30, missing 37, and missing 38; these were retained as failed cases,
not converted into stability scores. Voltage was not instrumented in the fast
core holdout, so voltage feasibility is not claimed from that table.

The key results are:

| Selection | robust cases | worst alpha (s^-1) |
|---|---:|---:|
| V8 | 12/24 | +0.870750 |
| missing 30 | 23/24 | -0.101954 |
| missing 32 | 22/24 | +21.757133 |
| missing 33 | 23/24 | +16.131602 |
| missing 34 | 13/24 | +19.928510 |
| missing 35 | 23/24 | +61.827526 |
| missing 36 | 20/24 | -0.003174 |
| missing 37 | 19/24 | +26.550743 |
| missing 38 | 16/24 | +0.093897 |
| representative 25% | 24/24 | -0.093767 |
| representative 50% | 24/24 | -0.137772 |
| representative 75% | 22/24 | +50.197711 |
| bus-37 singleton | 24/24 | -0.137793 |
| all-SG | 24/24 | -0.093525 |

Because the required matched-alpha radius pairs were not completed, no claim
of a ≥3× or ≥5× hidden-fragility ratio is justified. H3 remains untested.

## 9. Corrected transition planning

The candidate active dispatch is `P_cand=4620 MW`; targets are 1155, 2310,
3465, and 4620 MW. The discovery planner selections were frozen before the
holdout. The best executed selections, using holdout robust coverage as the
descriptive criterion, are:

* 25%: static P2 selected buses 33;35 (1282 MW), 23/24 robust. The P3
  selection retained 4080 MW but achieved only 19/24; no claim of P3
  superiority is made.
* 50%: P2 selected 33;34;35;38 (2620 MW), 24/24 robust.
* 75%: P2 selected 30;32;33;34;35;38 (3520 MW), 24/24 robust.
* 100%: V8 is fixed; no tested no-intervention design meets the 24-condition
  robust requirement.

P0 penetration-only always selected V8 and obtained 12/24 robust holdout
coverage. P1 nominal selected missing 30 at 25/50/75% and obtained 23/24.
P3 robust-dynamic selected missing 37 at 25/50/75% and obtained 19/24.
This directly triggers the practical-utility warning: in this campaign P3 did
not generalize materially better than P1/P2.

## 10. V8 repairs

### Controller only

The frozen controller coordinates are PLL, filter/reactance, and current
control scale. The best tested grid point was
`delta=(-0.1,+0.1,-0.1)`, with worst discovery alpha `+0.0347404 s^-1`.
No controller-only point inside the frozen ±10% box satisfies all nine
discovery conditions; controller-only repair is FAIL.

### Line only

The minimum single-line reinforcements that pass all nine discovery conditions
are line 1 `gamma=1.0915527`, line 2 `1.1374512`, line 27 `1.1936035`, and
line 46 `1.2114258`. At gamma 1.25 their discovery worst margins are 0.144936,
0.126516, 0.082631, and 0.070223 s^-1 respectively. On the 24 holdout
conditions, the minimum-reinforcement designs obtain 13/24, 13/24, 13/24,
and 12/24 robust cases respectively. None is an accepted holdout repair.

### Topology only and joint action

No one-switch topology action was accepted as a robust V8 repair. The outage
audit shows that 11 candidate openings disconnect the graph and 35 connected
openings are stable in the baseline outage diagnostic; this does not create a
V8 repair result. No joint controller+line Pareto solution was accepted, so
H6 is not tested. No statement that mixed intervention dominates single-family
repair is justified.

## 11. Nonlinear TDS

The preregistered common disturbance was implemented at the deterministic
largest active-load bus: +1% P/Q for 100 ms from 1.0 to 1.1 s, restored, and
simulated to 20 s. The secondary +5% pulse used the same bus and power factor.
All ten traces (T1, T2, T3, reference all-SG, and T5 line repair, each at 1%
and 5%) have finite observables.

For the 1% pulse, estimated rates were approximately: best 7/8 −0.3260 s^-1,
worst 7/8 −0.0975 s^-1, V8 −0.00554 s^-1, all-SG −0.1970 s^-1, and minimum
line-1 repair −0.0923 s^-1. The 5% pulse preserves the same qualitative
ordering. This supports sign/relative-damping consistency with linear alpha,
not exact equality. TDS therefore supports H8 in this bounded case, while it
does not rescue the failed holdout repair claim.

## 12. Figures and tables

The generated vector/PDF/PNG figure package contains F1–F6 and F9–F13. F7/F8
are not generated as evidence because the dynamic weak-node/link radii were
not executed. Figures requiring unexecuted dynamic weak-node/link radii or a
completed matched-pair radius campaign are explicitly data-limited and must
not be presented as evidence of those claims.
Tables T2, T3, and T6 are generated; T1 is documented in the model source and
preregistration. The main artifacts are:

* `results/PD39_PORTFOLIO_BLOCKER_STRUCTURE.csv`
* `results/PD39_7OF8_TO_8OF8_MODAL_ANALYSIS.csv`
* `results/PD39_CLASSICAL_BASELINES.csv`
* `results/PD39_STRUCTURED_RADIUS_SEARCH.csv`
* `results/PD39_HOLDOUT_CORE_SUMMARY.csv`
* `results/PD39_V8_REPAIR_SUMMARY.csv`
* `results/PD39_PLANNER_HOLDOUT_SUMMARY.csv`
* `results/PD39_TDS_SUMMARY.csv`

## 13. Negative results and kill criteria

The campaign does not pass the kill criteria for a strong general method:

* H3 hidden-fragility ratios were not tested with completed matched pairs.
* Dynamic weak-node/link ranking was not completed.
* P3 did not outperform P1/P2 on the 24-condition holdout.
* No repair passed all 24 holdout conditions.
* The modal mechanism is physical-looking but numerically ill-conditioned.

These are preserved as negative results. The appropriate next step is narrow
continuation: complete the 128-direction radius protocol, matched pairs,
weak-node/link campaigns, joint repair, and cross-model check before claiming
an IAS/TPWRS-grade planning method. Redesigning the IAS poster around a strong
“hidden margin/co-design validated” story is not justified yet; a narrower
“255+1 audit, true instability versus engineering margin, and failed
generalization of discovery repairs” poster is defensible.

## 14. Answers to the 32 final questions

1. Yes, the nine failures are exactly V8 × the nine discovery scenarios.
2. Both: V8 is nominally stable, but truly unstable in four high-PLL scenarios
   and below the 0.05 target in the other five.
3. Yes, `H_0={V8}` for the frozen nine-condition benchmark.
4. Yes, `H_0.05={V8}` exactly for the same finite benchmark.
5. The eight margins are 0.126229, 0.138022, 0.138296, 0.085835,
   0.136356, 0.138058, 0.138325, and 0.138020 s^-1 for missing 30, 32, 33,
   34, 35, 36, 37, and 38 respectively.
6. V8 exposes a complex electromechanical/control mode; some 7/8 worst cases
   instead have real/aperiodic critical modes.
7. V8 is physically interpretable but not numerically cleanly isolated;
   conditioning requires a gauge-aware replication.
8. Not reliably. Remaining SG MW/MVA/inertia and static rankings do not uniquely
   identify the cliff; SCR/gSCR were unavailable.
9. `rho_0(V8)=0.00190598` is the smallest boundary found under the 16-direction
   protocol, not a global radius.
10. `rho_0.05(V8)=0` at nominal because V8 already violates the target.
11. Not established: matched-alpha pairs were not radius-evaluated.
12. No strongest hidden-fragility pair can be claimed.
13. No confirmatory alpha–radius correlation is claimed because the matched set
    is incomplete and several crossings are unresolved jumps.
14. Weakest dynamic node: not estimable; node-restricted radius campaign was
    not completed.
15. Weakest dynamic link: not estimable; continuous weakening campaign was not
    completed.
16. Discovery static rankings do not identify the V8 cliff; node/link dynamic
    rank agreement remains untested.
17. The executed 25% P2 choice is buses 33;35, 1282 MW, 23/24 robust.
18. The executed 50% P2 choice is buses 33;34;35;38, 2620 MW, 24/24 robust.
19. The executed 75% P2 choice is buses 30;32;33;34;35;38, 3520 MW, 24/24
    robust.
20. Controller tuning only: no, within the frozen bounds.
21. Line reinforcement only: no tested minimum single-line repair passes all 24
    holdout conditions.
22. Topology only: no accepted one-switch robust repair.
23. Joint action: not solved to an accepted design; no positive claim.
24. The smallest discovery-only intervention found is line 1 at
    `gamma=1.0915527`; it is not a validated holdout solution.
25. The repair holdout cases retain requested load in convergent cases; no
    accepted full robust design exists.
26. No proposed 100% repair passes all 24; the best un-repaired V8 is 12/24.
27. Yes, TDS supports the relative 7/8 → V8 → repair ordering under the frozen
    common disturbance, with qualified rate interpretation.
28. Genuinely new here: the auditable SG→GFL 255+1 blocker distinction and the
    explicit separation of true instability, engineering margin, and bounded
    structured boundary search in this PD39 transition case.
29. Classical: small-signal eigenvalues, modal participation, finite scenario
    screens, static proxies, and conventional repair grids.
30. Useful as a warning/audit workflow; not yet useful as a validated universal
    planner because P3 did not generalize and no 100% repair was validated.
31. Do not redesign the IAS poster around a strong validated-method claim yet;
    a narrower negative-result/audit poster is justified.
32. Research status: **narrow continuation**, not strong completion and not a
    stop. The primary blocker result is solid; the general hidden-radius and
    co-design claims need the unfinished experiments.
