# Final unit and rating audit

This document freezes the unit interpretation for the final campaign and both
papers. The full per-file audit, reruns and claim classes are in
`docs/PHASE_I_UNIT_CORRECTION.md`; the data contract is `src/ibr_cycles/units/`
(tests `tests/test_unit_contract.py`).

## 1. Quantities kept separate

| quantity | unit | canonical name | provenance | documented in the benchmark? |
|---|---|---|---|---|
| dispatched active power displaced | MW | `replaced_pg_mw` | measured at the solved equilibrium | yes (power flow) |
| reactive power displaced | Mvar | `replaced_q_mvar` | measured | yes |
| apparent-power rating of the retired machines | MVA | `retired_sg_sn_mva` | machine `Sn` | yes (IEEE-39, Kundur); for IEEE-68, `Sn` is a per-unit base, not a rating |
| converter rating | MVA | `replaced_sn_mva` | rating rule: `Sn` (IEEE-39, Kundur), `abs(S_gen)/0.8` (IEEE-68) | declared rule |
| synchronous active limit | MW | `replaced_pmax_mw` | ANDES `PV/Slack.pmax` | IEEE-39 only; Kundur has placeholders; IEEE-68 none |
| **PV nameplate MW** | MW | **none** | the converter holds `P_ref` equal to the displaced dispatch, with unlimited DC | **no** |
| condenser rating | MVA | `condenser_sn_mva` | F8/G1 condenser rating fraction × retired `Sn` | declared |

The historical value 4270.7 is the sum of the flagship machines' `Sn`
(1040 + 1174.8 + 1085.7 + 970.2). It is an **MVA rating**. It is never called MW.

## 2. Frozen interpretation (Decision J)

| case | converter active dispatch [MW] | converter rating [MVA] | synchronous MVA |
|---|---|---|---|
| flagship 30+33+35+37 | 2096.6 | 4270.7 | 0 |
| converter retune (RB, RC, E34 M1) | 2096.6 retained | 4270.7 | 0 |
| restore SG30 | gives up 436.1 | 3230.7 | 1040 (machine kept) |
| restore SG37 | gives up 321.5 | 3300.5 | 970.2 (machine kept) |
| condenser 0.25 Sn at each retired bus (RD) | 2096.6 | 4270.7 | 1067.7 condenser |
| condenser at bus 30 only (E34) | 2096.6 | 4270.7 | 166–270 condenser |

Never write "4.27 GW PV retained" or "X MW of PV capacity": the benchmark has no
PV nameplate.

## 3. Figures and tables with ambiguous units

| item | status |
|---|---|
| `outputs/.../E39_baseline_audit/E39_ROC_PR.png` (curve `replaced_mw`) | **prohibited**; replaced by `results/UC/UC02/UC02_ROC_PR_unit_corrected.png` |
| E21/E23/E34 tables with `pv_mw_*` columns | prohibited as MW; the corrected tables are `results/UC/UC03/*` |
| E11 `mw_replaceable`, E12–E14/E18 `replaced_mw` columns | columns are MVA; the corrected quantities are in `results/UC/UC01/UC01_census_quantities.csv` |
| E34 Pareto figures (sync MVA vs margin; controller change vs margin) | unaffected (axes are MVA and controller norm) |
| every figure of this campaign | uses `replaced_pg_mw` [MW] and `Sn` [MVA] on separate axes |

## 4. Consequence for the planning objectives (FC15)

- The objective "MW replaced" is `replaced_pg_mw`, the active dispatch.
- Synchronous support is in MVA.
- The controller change is a norm.
- The three axes are never merged without a stated weight sensitivity.
