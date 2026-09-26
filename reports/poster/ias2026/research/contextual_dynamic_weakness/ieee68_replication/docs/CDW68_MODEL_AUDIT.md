# CDW68 model audit (R1)

Branch `research/cdw-ieee68-replication`, from `e03a5312`. It was written before any portfolio containing a converter was evaluated. The audit computed only the base case (no replacement, documented excitation k = 1); results are in `results/CDW68_R01_audit.json`.

## 1. Network and machine data

| item | value | source |
|---|---|---|
| network | IEEE 68-bus / NETS–NYPS, 68 buses, 83 branch rows (the header of the source says 86) | Singh & Pal (2013), IEEE PES TF report v3.3, Appendix B; transcribed in `configs/ieee68/ieee68_network.json` by `experiments/G3_import_ieee68.py` (Gate 3, commit `d9fa097e`) |
| transformer taps | tap ratio per row (e.g. 1.025 on step-up rows), no phase shift | same |
| base / frequency | 100 MVA / 60 Hz | same |
| machines | 16; sub-transient model with four rotor coils plus the dummy coil E'dc, equations (1)–(6) of the report | `src/ibr_cycles/models/ieee68_devices.py` (read-only) |
| machine per-unit base | the `mac_con` MVA column: 100 MVA (200 MVA for G13 and G16). **This is a per-unit base, not a rating.** | report; Gate 3 note |
| excitation | DC4B on G1–G8 and G10–G12; ST1A on G9; manual (constant Efd) on G13–G16 | report eqs. (7)–(8) |
| PSS | speed-input washout plus three lead-lags on G1–G12 | report eq. (9) |
| loads | constant impedance, converted from the power-flow P/Q at the solved voltage | Gate 3 convention |
| slack | G16 (bus 16, an area equivalent) | report |
| controller limits | not represented; the equilibrium is checked to lie inside every documented limit (smallest margin 3.00) | Gate 3 |

Reproduction (frozen Gate 3 code, rerun here):
- power flow against Table 1: max |dV| = 5.0e-5 pu, max |d angle| = 5.0e-5 deg, **REPRODUCED**;
- the four inter-area modes: **REPRODUCED**;
- all 15 electromechanical modes within 5.0e-4 Hz and 4.5e-4 damping-ratio percentage points.

The Ybus rebuilt from the JSON by `code/_r68.py` equals the frozen network's Ybus to 0.0.

## 2. Damping and governors: provenance

- **Singh & Pal v3.3 (the benchmark as published).** Governors are ignored and the mechanical torque is constant; D = 0 for every machine. Canizares et al. (2017), *IEEE Trans. Power Syst.* 32(1), Benchmark 6, describes the same system without governors.
- **Physical damping still present.** The sixth-order machines keep their damper-winding circuits (T''do, T''qo), which the IEEE-39 two-axis model of the earlier campaigns lacked.
- **PST `data16m.m`** (Power System Toolbox, Chow & Rogers). Copy used: PNNL GridSTAGE, BSD-style licence, sha256 `c16e78f4…`. It is the same system with the same scheduled generation (P = 2.50, 5.45, …, 13.50 pu; slack 40.00), with generator buses numbered 53–68. It documents:
  - machine MVA ratings (300–1900 MVA for G1–G12; 10–12 GVA for the area equivalents G13–G16);
  - rotor damping d_o, zero on G1–G12 and equal to H on the four area equivalents;
  - a turbine-governor block for all 16 machines: PST tg model 1, 1/R = 25 on the machine rating, Ts = 0.1, Tc = 0.5, T3 = 0, T4 = 1.25, T5 = 5.0 s. **In the file this block is commented out** (`tg_con = [];` follows it).
- **PST `tg.m`** (Sandia PSTess, MIT-style licence, sha256 `c1584021…`) supplies the governor equations. They are listed in `inputs/pst_primary_frequency_v1.json`.

No other documented governor or damping data for this benchmark were found (searches listed in the final report).

## 3. Model variants (frozen in `inputs/pst_primary_frequency_v1.json`)

| variant | governors | rotor damping | role |
|---|---|---|---|
| REAL | PST tg model 1 on every surviving machine | PST d_o on every surviving machine (non-zero only on G13–G16) | confirmatory |
| NOGOV | none | as REAL | ablation B (governor off only) |
| SP33 | none | 0 | ablation C: the published Singh & Pal benchmark, identical to Gate 3 |

- **Base conversion.** Damping torque and governor response keep the same MW per pu speed deviation. A coefficient given in pu on the PST rating is multiplied by rating / Sn_SP.
  - Governor gain on the device base: 75 (G1), 200 (G2–G4, G7, G8), 175 (G5), 225 (G6), 250 (G9), 300 (G10), 400 (G11), 475 (G12), 1500 (G13), 2500 (G14, G15), 1375 (G16).
  - Damping on the device base: 244.7 (G13), 300 (G14), 300 (G15), 244.8 (G16).
- **Replaced units.** A replaced machine disappears together with its governor and damping. The converter has no frequency response.

## 4. Base-case audit (S = ∅, k = 1)

| variant | states | equilibrium residual | zeros of full A below 1e-6 | transverse centre | status | α⊥ (s⁻¹) | rightmost EM-band mode |
|---|---|---|---|---|---|---|---|
| REAL | 274 | 7.3e-13 | 1 (rotation) | span{R_x} | STABLE | −0.0669 (real; PSS washout, 1/T_W = 0.0667) | −0.367 s⁻¹ at 0.555 Hz |
| NOGOV | 226 | 7.3e-13 | 1 (rotation) | span{R_x} | STABLE | −0.0461 (slow common-frequency mode, 0.002 Hz) | −0.371 s⁻¹ at 0.518 Hz |
| SP33 | 226 | 7.3e-13 | 0 below 1e-6: the reference Jordan pair splits to about 1e-6 | span{R_x, w} (w derived) | STABLE | −0.0556 (Gate 3 value reproduced) | −0.118 s⁻¹ at 0.520 Hz (Table 4: 0.520 Hz, 3.62 %) |

- **Operating point.** The equilibrium network voltages are identical across the three variants (max difference 1.1e-16). Governors and damping do not move the operating point.
- **What this means for the tests.** With governors or equivalent damping, the rightmost transverse eigenvalue of the base case is a slow real mode, not an electromechanical one.
  - The equivalents' documented damping moves the 0.52 Hz inter-area mode from −0.118 to −0.37 s⁻¹.
  - This is a property of the documented model, known before any portfolio was evaluated. It is why the preregistration adds an EM-tracked secondary test next to the unchanged primary Level-D test.

## 5. Converter models (R2)

**Model A: the frozen TX4 grid-following converter.** `ibr_cycles.models.ieee39_devices.GridFollowingConverter`, `ConverterParameters` defaults.
- Controls: SRF-PLL (kp 53, ki 1400), filtered power measurement (τ 0.03 s), outer PI (0.2/8), inner current PI (0.25/6) behind x_f = 0.15, r_f = 0.01.
- Voltage control: the leaky Q/V regulator q = q_ref + g·2·e + x_v, x_v' = g·20·e − 0.05·(x_v − q_ref).
- Rating |S_gen|/0.8 (Gate 3 rule). Per unit on its own rating.
- Omissions: no current limits, DC link or delays.
- Initialization: matched dispatch; it carries the replaced machine's P and Q at the documented power flow.
- Provenance: TX4 freeze `69f200df`; the same device as on IEEE-39.

**Model B: the frozen WECC library GFL chain of the hardening cross-model study.** TX3-GFL-0.1 parameters (sha256 `f07a6a40…`) in ANDES 2.0.0.
- Blocks: PLL2 + BusFreq + REGCP1 + REECB1 + REPCA1, constant-Q (QFLAG = 0), attached exactly as in `H17_andes_alt.py::add_wecc`.
- Equilibrium: the TX3 forced-PQ formulation.
- Structural treatment: H17 refinements 1, 1b and 2 (the one-dimensional centre applies because REAL is governed).
- The machines must be transcribed into ANDES for the 68-bus case. Qualification criteria are in the preregistration, §5.
- No PLL or voltage-control gain of B may change.

**Grid-forming model.** None is validated in the repository, so the optional third holdout is not run.
