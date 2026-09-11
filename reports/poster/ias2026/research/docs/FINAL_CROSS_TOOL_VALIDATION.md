# IEEE-39 final cross-tool validation (Phases A–G)

Date: 2026-09-11. The campaign ran on the repository state at `81ddb963`. Nothing was retuned: no controller parameter, network value, dispatch or gate tolerance was changed. The canonical case was used throughout; stock `pandapower.networks.case39` appears only as an informative note.

- **Inputs:** the handoff zip was extracted to `validation_inputs/cross_tool_handoff/`, a new directory; no frozen output was overwritten.
- **Runs:** all runs live under `validation_runs/cross_tool_final/`.
- **Pack scripts:** the supplied pack scripts were run from byte-identical copies (SHA256 checked).
- **Manifest:** versions, hashes and tolerances are in `results/CROSS_TOOL_MANIFEST.json`.

## Bottom line

| Layer | Verdict |
|---|---|
| Network and equilibrium (Python, pandapower, ANDES) | **PASS on the physical case.** The *as-supplied* pandapower pack **FAILS**, from a translation-convention defect. Two tool-side discrepancies were both root-caused to machine precision. |
| Equation-equivalent synchronous dynamics (Python vs ANDES) | **PASS**: every dynamic gate, 12/12 subsets |
| Branch sensitivity, equation-equivalent same cases (E1, the frozen protocol test) | **PASS**: DAE vs ANDES 12/12 signs, ρ = 1.000; port vs ANDES 11/12, ρ = 1.000 |
| Frozen GFL predictions vs ANDES (E2, cross-model) | **FAIL**: 7/12 signs, ρ = 0.51. The same failure occurs against the project's own static-injection DAE (E2c), so it measures converter-model dependence, not a solver artifact. |
| SG-to-custom-GFL portfolio results | Not testable in ANDES (no equation-equivalent GFL). They remain validated only by the documented project DAE. |

**Hard-stop disclosure.** Under a literal reading of the README hard-stop rule, the campaign stops at Phase B, because the as-supplied pandapower parity fails (Vm 9.3e-2 pu). The failure was diagnosed to its first non-equivalent element: `from_ppc` puts the off-nominal tap on the wrong winding of three transformers. Only that translation was corrected (§B); no data or tolerance changed, and the campaign then continued. **If the correction is not accepted, treat Phases C–E below as not executed** and cite only the Python–ANDES static agreement. Note that Phases C and D do not depend on pandapower.

## Environments

| venv | role | versions |
|---|---|---|
| `.venv/tx3-analysis` (frozen, unchanged) | A, D-supp, E internal, assembly | Python 3.13.14, numpy 2.5.2, scipy 1.18.1, pandas 3.0.5 |
| `.venv/tx3-andes` (frozen, unchanged) | C, D, E ANDES | ANDES 2.0.0 (vendor tag v2.0.0, `eda5163c`, clean), numpy 2.5.2, scipy 1.18.1, pandas 3.0.5, kvxopt 1.3.3.1 |
| `.venv/xtool-pandapower` (new, authorized) | B | pandapower 3.4.0, numpy 2.3.5, scipy 1.16.3, pandas 2.3.3 |

The package lists of `bc-cert`, `tx3-analysis`, `tx3-andes` and `tx3-paraemt` are identical before and after the campaign (`supp/venvs/*.before|after.txt`). The xtool venv is isolated (`include-system-site-packages = false`) and shares nothing with the frozen venvs.

## Phase A: pure-Python canonical reference

This phase used `run_python_reference.py` and `supp_static_checks.py`.

| Check | Measured | Gate | Verdict |
|---|---|---|---|
| PF mismatch | 2.75e-12 | < 1e-10 | PASS |
| pytest (`test_ieee39_network`, `test_ieee39_baseline`) | 9 passed | all | PASS |
| Ybus vs frozen ANDES (stale npy + restored bus-4/5 shunts) | 1.14e-13 pu | ≤ 1e-9 | PASS |
| Vm / Va vs frozen ANDES PF | 1.72e-7 pu / 1.71e-7 rad | ≤ 1e-5 | PASS |
| PV generator P | 1.14e-12 pu | ≤ 1e-8 | PASS |
| Slack generator P | 5.69e-6 pu (0.57 kW) | ≤ 1e-8 | **FAIL (literal), EXPLAINED** (see C) |
| Generator Q (PV / slack) | 2.05e-6 / 3.57e-6 pu | ≤ 1e-5 | PASS |
| Stored F1 R3 (α, f, RHP) | 1.50e-6 1/s, 1.32e-7 Hz, all agree | — | as frozen |

## Phase B: pandapower static parity (canonical JSON → `canonical_ppc` → `from_ppc`)

### Three defects prevent the supplied pack from running on its own recommended stack

`run_pack_pandapower_compat.py` fixes the three defects below in memory. It changes nothing else.

1. **Import path.** `from pandapower.converter import from_ppc` fails in 3.4.0, because `converter/__init__.py` is empty. The wrapper binds `pandapower.converter.pypower.from_ppc` instead.
2. **Key names.** `canonical_ppc` reads the pv keys `p/v/sn`, but the canonical JSON has `p0/v0/Sn`. The wrapper adds the aliases with identical values. The project loader does exactly this mapping (`ieee39_network.py:126-137`).
3. **Upstream bug in pandapower 3.4.0** (`from_ppc.py:303`). For an impedance branch with RATE_A = 0, the code indexes `sn` instead of `sn_mva` and raises IndexError. The pack writes RATE_A = 0 everywhere. The wrapper sets RATE_A to pandapower's own fallback, MAX_VAL = 99999. RATE_A is only a rating and per-unit base, so the electrical parameters are unchanged; the Ybus check below confirms this.

### As supplied (with the three fixes): **FAIL**

| Check | Measured | Gate |
|---|---|---|
| Ybus | 9.383 pu at (6,6) and (31,31), 0.55 at (12,12) | ≤ 1e-9 |
| Vm / Va | 9.26e-2 pu / 1.75e-2 rad | ≤ 1e-5 |
| Generator P (slack) / Q | 3.01e-2 / 5.43 pu | 1e-8 / 1e-5 |
| Branch P/Q: pack formula / pandapower native | 420.3 / 561.3 MVA | ≤ 1e-3 |

The first non-equivalent element is found by `diag_pandapower_ybus.py`. 43 of 46 branch admittances match to 1e-14. The three that don't are the transformers whose canonical *from* bus is the LV bus: 35 (31→6, t = 0.9), 37 (12→11, t = 1.006) and 38 (12→13, t = 1.006).

The cause is a convention mismatch:

- `from_ppc` swaps hv/lv and ends up with the tap on the HV winding.
- The canonical JSON, the project model, ANDES (`Line.build_ybus`, tap at `bus1`) and MATPOWER all put the tap at the from bus.
- Passing `tap_side="lv"` changes the representation but not the Ybus, because `from_ppc` re-refers the impedance.

This is a pack translation defect, not a network error.

### Translation-corrected (`run_pack_pandapower_tapside.py`): **PASS**

The three rows are rewritten into the exactly equivalent two-port with the HV bus as the from bus: tap′ = 1/t, r′ = r·t², x′ = x·t². The derivation is

[[y/t², −y/t], [−y/t, y]] ≡ tap 1/t at the other end, with y′ = y/t².

The script asserts the two-port is unchanged (≤ 1.4e-14) before use. All three rows have b = g = 0 and φ = 0.

| Check | Measured | Gate | Verdict |
|---|---|---|---|
| Ybus (pandapower internal vs project) | 1.14e-13 pu | ≤ 1e-9 | PASS |
| Vm / Va | 2.86e-14 pu / 3.22e-14 rad | ≤ 1e-5 | PASS |
| Generator P: PV / slack | 1.14e-12 / 1.31e-12 pu | ≤ 1e-8 | PASS |
| Generator Q: PV / slack | 1.96e-12 / 1.05e-12 pu | ≤ 1e-5 | PASS |
| Branch P/Q: pack formula / pandapower native | 2.07e-10 / 2.03e-10 MVA | ≤ 1e-3 | PASS |

The pack's own `compare_results.py` exits 1 on the as-supplied variant and 0 on the translation-corrected one, with all 10 rows PASS.

Informative only: stock `case39()` gives bus-1 Vm 1.0394, against 1.0487 in the canonical case. It is a different solved case and is not used as evidence.

## Phase C: live ANDES static reproduction

This phase used `run_andes_dynamic.py` and `supp_andes_static.py`.

- **ANDES version:** 2.0.0.
- **Case file:** `vendor/andes/andes/cases/ieee39/ieee39_full.xlsx`.
- **Case SHA256:** `9c2048dc94201ee48ffe816de65d53831fa1f4b0e7fdf32f24c50f1017edcfe5`. This is identical to the source hash recorded in the canonical JSON.
- **Operating-point SHA256** (bus v, a and generator P, Q rounded to 8 dp): `98373aaf56105e599e38dccc7e55ee83a94687cdbc7ddc21be519fd475807a51`.
- **Repeatability:** the live PF is bit-identical to the frozen reference and to `results/F1/andes_live_powerflow.json`; the error is 0.0 for Vm, Va, P and Q.
- **No stale Ybus:** the live Ybus comes from ANDES itself (`System.build_ybus`: Line + Shunt) and matches Python to **1.14e-13 pu**. `ieee39_ybus_andes.npy` differs from the live Ybus at buses 4 and 5 by exactly the shunts (b = 1.0, 2.0). It is not used for any verdict.
- **Semantics:** shunts, line order, bus1/bus2 and taps are identical to the canonical JSON, and Line b1/b2/g1/g2 are all 0.

| ANDES live vs Python | Measured | Gate | Verdict |
|---|---|---|---|
| Vm / Va | 1.72e-7 pu / 1.71e-7 rad | ≤ 1e-5 | PASS |
| PV generator P | 1.14e-12 pu | ≤ 1e-8 | PASS |
| Slack generator P | 5.69e-6 pu | ≤ 1e-8 | **FAIL (literal), EXPLAINED** |
| Generator Q (PV / slack) | 2.05e-6 / 3.57e-6 pu | ≤ 1e-5 | PASS |

**Discrepancy root cause.** The residual is unchanged at ANDES PFlow tol 1e-12, so it is not a convergence effect. ANDES 2.0.0's Line power-flow/DAE equations use `yhk = u/((r+1e-8) + 1j(x+1e-8))` (`andes/models/line/line.py:200`). Its own `build_ybus` and the project model use u/(r + jx).

Diagnostic (`diag_andes_regularization.py`): a throw-away copy of the network with the same +1e-8 is compared against ANDES at tol 1e-12. The agreement is Vm 2.9e-14, Va 3.1e-14, PV P/Q ≤ 1.9e-12, slack P 1.5e-12 and slack Q 9.3e-13. So the whole Python–ANDES static difference comes from this regularization inside ANDES, and the network and dispatch are identical.

## Phase D: equation-equivalent ANDES dynamic reconciliation (F1 R2/R3, 12 cases)

This phase used `run_andes_dynamic.py` → `F1_andes_equivalent_worker.py` (fresh, 17 s) and `supp_dynamic_checks.py` (internal model re-solved live).

| Check | R3 first-order AVR | R2 manual excitation | Gate | Verdict |
|---|---|---|---|---|
| Critical-band α max error | 1.50e-6 1/s | 1.02e-7 | ≤ 1e-4 | PASS |
| Critical-band frequency max error | 1.32e-7 Hz | 1.58e-7 | ≤ 1e-3 | PASS |
| RHP count agrees, every subset | 6/6 | 6/6 | all | PASS |
| Band 0.3–1.5 Hz nearest-mode mismatch, both directions | 1.61e-6 | 1.62e-6 | ≤ 1e-4 | PASS |
| Fresh vs frozen ANDES spectra | 0.0 | 0.0 | repeatable | PASS |
| Live vs frozen internal α | 8e-17 | 6e-17 | repeatable | PASS |

The flagship R3 critical mode is α = +0.328731 (internal) against +0.328732 (ANDES), at 0.5825 Hz.

**Discrepancy, outside the band.** The whole-dynamic-spectrum worst match is 0.42 1/s (R3) and 0.13 (R2). Every internal eigenvalue that lacks an ANDES partner belongs to the 12 zero-gain stabilizer filter states `pss_w` (−1/4.2) and `pss_l` (−1/0.75). Their column coupling into the rest of A is exactly 0.0; ANDES omits IEEEST in this configuration. Modes present only in ANDES are the fast degraded-GENROU subtransient states (τ = 1e-4). Neither affects the electromechanical band.

## Phase E: frozen 12-line branch-sensitivity holdout

`holdout_lines.json` was loaded first, before any sensitivity result: seed 20260911, canonical line indices [3, 9, 15, 17, 22, 26, 31, 32, 35, 39, 42, 44]. It was not regenerated or replaced.

- **γ_e definition:** γ_e multiplies the complete branch two-port (series and charging, tap unchanged). γ = 1 ± 0.002, central difference.
- **Internal scaling:** `ybus_scaled` is copied verbatim from the frozen F2c script, as are `qs`, `mu` and `trans_eig`.
- **ANDES scaling:** r, x ÷ γ and b, g × γ on the same Line row. The Line sheet order was verified identical to the canonical JSON.

**Why two tests.** The frozen predictions in `F2c_holdout_reequilibrated_port_line_sensitivity.csv` belong to the custom-GFL flagship: Q/V gain 0.20768, k = 1.425, t = 1.5, 0.706 Hz boundary. ANDES cannot express that model. The frozen protocol (`CLAUDE_FINAL_ANDES_NETWORK_VALIDATION.md`) prescribes the F1 equation-equivalent configuration on "exactly the same cases" in both tools, which is **E1**. The comparison against the frozen GFL predictions requested in the handoff is **E2**, labelled cross-model.

**E1 case.** The case is the R3 flagship {30, 33, 35, 37}, replaced by constant-power static injections (`q_policy="matched"`), with first-order AVR, stabilizer gain 0, no governor and constant-power loads.

To solve the same case in ANDES, the replaced units' (P, Q) is re-read at each γ from the all-PV ANDES power flow of the scaled full case. This is the internal "matched" semantics; at γ = 1 it reduces exactly to F1. The port derivative is ds*/dγ = −(∂μ/∂γ)/(∂μ/∂s) at s* = λ₀. This point is admissible because the port operator is built on the base resolvent.

| Comparison (predicted vs reference) | Signs | Spearman (Re) | Top-5 stabilizing | Median / max rel. complex error | Strong branch opposite | Frozen strong gate |
|---|---|---|---|---|---|---|
| E1a internal DAE vs **ANDES** (same cases) | **12/12** | **1.000** | **5/5** | 7.0e-6 / 2.6e-4 (L3) | 0 | **PASS** (preferred met) |
| E1b internal **port** vs **ANDES** (same cases) | 11/12 | 1.000 | 5/5 | 3.0e-3 / 6.6e-2 (L3) | 0 | **PASS** |
| E1c internal port vs internal DAE | 11/12 | 1.000 | 5/5 | 3.1e-3 / 6.5e-2 (L3) | 0 | PASS |
| E2a frozen F2c **GFL** port vs ANDES R3 (cross-model) | 7/12 | 0.510 | 4/5 | 1.02 / 7.43 (L3) | 1 (L44) | **FAIL** |
| E2b frozen F2c GFL DAE vs ANDES R3 (cross-model) | 7/12 | 0.510 | 4/5 | 1.02 / 7.44 (L3) | 1 (L44) | FAIL |
| E2c frozen F2c GFL DAE vs **internal** R3 DAE (same tool) | 7/12 | 0.510 | 4/5 | 1.02 / 7.44 (L3) | 1 (L44) | FAIL |
| E0 frozen F2c GFL port vs frozen F2c GFL DAE (frozen record) | 12/12 | 1.000 | 5/5 | 4.2e-3 / 1.3e-2 (L42) | 0 | PASS |

Observations:

- The only E1 sign split is **L3**. There the reference sensitivity is essentially zero (ANDES +4.4e-5; port −1.4e-4), and L3 is below the median, so it is not a strong branch.
- **Top-5 stabilizing branches.** ANDES R3 gives L42, L35, L26, L9, L44. The GFL port gives L35, L32, L26, L42, L9. Four branches (L42, L35, L26, L9) are among the five most stabilizing in both converter models.
- **Magnitudes and weak branches.** Magnitudes differ by roughly 5–10×. Five branches change sign between models (L3, L17, L31, L39, L44); L44 is strong in R3.
- **Converter dependence, not solver.** E2c reproduces E2a/E2b to three digits. Replacing the GFL by static injections changes the ranking in the project's own DAE exactly as it does in ANDES.

**Port identities in the E1 configuration.** The measured values are |μ(λ₀)+1| = 3.6e-7 and a det(T_S)/det(T_0) = det(I+M) relative error of 1.2e-7, against the 1e-8 theory gate. This is FAIL (literal), EXPLAINED.

`diag_port_identity.py` explains it. In static-injection mode, the base and flagship equilibria differ by 2.5e-7 in bus voltage, independent of solver tolerance (1e-9 and 1e-12 give the same result). The residual T_S − T₀ − U dY Uᵀ therefore sits entirely off the port rows. The frozen GFL boundaries, where the matched equilibria coincide, hold the identities to ≤ 1.4e-13 over 160 checks, with the −1 closure at ≤ 8.5e-14 (FC18).

Per the handoff rule: A–D pass and the equation-equivalent holdout passes, but the holdout against the frozen GFL predictions fails. So **the paper stays and the network-design claim is downgraded**:

- **Supported:** the port derivative reproduces independently computed physical branch sensitivities.
- **Not transferable:** a branch ranking is not a transferable conclusion, because it depends on the converter model.

The TX4 manuscript makes no branch-ranking claim, so its content is unaffected.

## Phase F: theory-consistency audit

| Statement | Verdict | Evidence |
|---|---|---|
| T1 Network and equilibrium are implementation-independent | **CONDITIONAL** (PASS in substance) | Ybus identical in Python, pandapower and ANDES to 1.1e-13. Equilibrium: pandapower 3e-14; ANDES 1.7e-7, exactly explained by ANDES's +1e-8 regularization (→ 3e-14 when mirrored). Condition: acceptance of the pandapower translation correction, since the as-supplied pack fails. |
| T2 Equation-equivalent synchronous dynamics are implementation-independent | **PASS** | α 1.5e-6, f 1.3e-7 Hz, RHP 12/12, band nearest-mode 1.6e-6 both ways, fresh = frozen |
| T3 The custom SG-to-GFL hypergraph remains model-specific | **PASS** (statement confirmed) | ANDES has no equation-equivalent GFL. E2c: replacing the GFL by static injections changes the branch ranking (7/12 signs). |
| T4 Network-closure and line-sensitivity conclusions are not solver artifacts | **CONDITIONAL** | PASS for the equation-equivalent configuration (E1a 12/12, ρ = 1.000; E1b 11/12, ρ = 1.000). The closure identity is algebraic (≤ 1.4e-13 frozen). The GFL-specific sensitivities are reproduced only inside the project DAE (E0 12/12). |
| T5 The framework holds even if the tuning changes the specific coalition | **PASS** (framework level) | Port derivative = DAE derivative under both converter models (E0, E1c), although the rankings differ. Identities hold at all 10 frozen boundaries where κ changes 4→3→2→3→4 (FC18, frozen). This campaign adds device-model robustness only; it did not test new coalitions. |

## Phase G: no overclaiming

**Allowed sentence for TPWRS** (exact):

> The IEEE 39-bus network and operating point were reproduced independently in pandapower 3.4.0 and ANDES 2.0.0 (identical bus-admittance matrices to $10^{-13}$ p.u.; bus voltages within $2\times10^{-7}$ p.u., the residual being ANDES's internal $10^{-8}$ p.u. branch-impedance regularization), and on an equation-equivalent synchronous configuration (first-order AVR, no stabilizer or governor, replaced units as constant-power injections) ANDES reproduces the electromechanical eigenvalues within $2\times10^{-6}$ s$^{-1}$ and the branch-reinforcement sensitivities of the critical mode on a preregistered 12-branch holdout (12/12 signs, Spearman 1.00), which the reduced network-port derivative also reproduces (Spearman 1.00; the one sign difference is on a branch with near-zero sensitivity); the SG-to-GFL portfolio conclusions rest on the documented project DAE, since ANDES has no equation-equivalent model of the converter.

If the pandapower translation correction is not accepted, replace "in pandapower 3.4.0 and ANDES 2.0.0" with "in ANDES 2.0.0".

**Prohibited:**

- "ANDES independently validates the SG-to-GFL portfolio result."
- "ANDES validates the GFL portfolio hypergraph."
- "The line diagnosis / branch ranking of the GFL portfolio is independently validated."
- "pandapower reproduces the canonical case" without the translation note.
- Any claim of independent reproduction of κ, H(θ), P4, or the GFL boundaries.

## Answers to the seven questions

1. **Physical IEEE-39 implementation correct?** Yes. It is identical to the ANDES source file (hash) and to three independent Ybus constructions (1e-13). The dispatch is reproduced exactly.
2. **Does pandapower reproduce the same operating point?** Yes, to 3e-14 pu, once the from-bus tap convention is translated correctly. The supplied translation does not (Vm 9e-2).
3. **Does ANDES reproduce the equation-equivalent synchronous dynamics?** Yes, every gate on all 12 cases.
4. **What remains specific to the custom GFL model?** The GFL converter dynamics and everything built on them: H(θ), κ(θ), the P4 witness, the policy maps, the GFL boundaries and the GFL branch ranking.
5. **Does independent ANDES evidence support the line-diagnosis/design ranking?** It supports the *method*: port derivative = physical sensitivity, confirmed in an independent solver on the same equation-equivalent case. It does *not* support the GFL-specific ranking, which changes with the converter model (7/12 signs). The four most stabilizing branches (L42, L35, L26, L9) recur in both models.
6. **Exact TPWRS sentence:** see Phase G.
7. **Scientific reason to delay submission?** None for the spectral-composability paper. Two author actions remain:
   - accept or reject the pandapower translation correction;
   - optionally replace the Limitations sentence "ANDES reproduces the power flow, the base inter-area mode (within 4.1%) …" (`main.tex:959`) with the Phase G sentence. The manuscript was not edited in this campaign.

## Files

- **Results:**
  - `results/CROSS_TOOL_EVIDENCE_TABLE.csv`
  - `results/PANDAPOWER_STATIC_PARITY.csv`
  - `results/ANDES_DYNAMIC_PARITY.csv`
  - `results/ANDES_LINE_SENSITIVITY_HOLDOUT.csv`
  - `results/CROSS_TOOL_MANIFEST.json`
- **Figure:** `results/CROSS_TOOL_VALIDATION_FIGURE.pdf|.png`, with `_source.csv`. The panels are:
  - (A) Python vs pandapower branch-flow parity;
  - (B) Python vs ANDES electromechanical eigenvalues;
  - (C) port and frozen-GFL predictions vs ANDES line sensitivity.
- **Scripts:** `validation_runs/cross_tool_final/*.py`, with raw outputs in `pack/results`, `pack_tapside/results`, `supp/` and `phaseE/`.

The cross-tool validation is frozen. No further benchmark search is launched.
