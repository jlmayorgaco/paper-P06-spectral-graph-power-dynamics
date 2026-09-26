# ParaEMT network and operating-point equivalence (EMT01)

- Preregistration: `docs/20260911_PAREMT_EMT_PREREG_V1.md` §2, commit c2947bd8.
- Script: `experiments/paremt_emt/EMT01_network_equivalence.py`.
- Results: `results/EMT01/`.

**Gate G1: PASS.**

## 1. The official ParaEMT IEEE-39 is not the frozen TX4 IEEE-39

`results/EMT01/network_diff.csv` compares the two cases row by row: 195 rows,
of which 171 are identical. `operating_point_diff.csv` compares generator
P/Q, ratings, |V| and θ.

**Identical:**
- base MVA (100) and frequency (60 Hz);
- the 39-bus numbering;
- all 34 line endpoints and R/X/B;
- the 12 transformer endpoints and R/X;
- the two capacitor shunts (buses 4 and 5, 1.0 and 2.0 pu).

**Different:**
- **Transformer taps.** The official case carries only one off-nominal ratio
  (6–31, k = 0.9714, applied on the to-bus). The TX4 case has 11 off-nominal
  taps on bus1, e.g. 31→6 t = 0.9 and 19→33 t = 1.07. Upstream ParaEMT also
  never applies `xfmr_k` (EMT00).
- **Loads.** They differ at buses 3, 4, 7, 8, 12, 16, 31 and 39. For example,
  bus 39 is 11.04 pu officially and 4.0 pu in TX4. The official load 7 Q is
  8.4 pu against 0.84 pu in TX4. The TX4 case is the ANDES `ieee39_full`
  data set, sha256 9c2048dc….
- **Generators.** Ratings are 1000 MVA each officially (2000 at bus 39)
  against the TX4 Sn of 836–1684 MVA. The dispatch, and therefore |V| and θ,
  differ.
- **Units.** The official case is all per unit (`bus_basekV = 1`); so is the
  TX4 case built here.

The official case is therefore **not** used for any TX4 result. The frozen
benchmark was not changed. Instead, a TX4-specific ParaEMT case is built from
the frozen data.

## 2. TX4-specific ParaEMT network

**Construction.** `tx4_emt.build_network` calls ParaEMT's own `numba_InitNet`
(TX4 patch, `experiments/paremt_emt/tx4_patch_paremt.py`):
- 34 lines: series R–L plus C = b/2ω0 per end;
- 12 transformers: series R–L with an ideal from-side tap on bus1;
- 2 capacitor shunts.

**Numerical damping is switched off.** The upstream damping resistors change
the 60-Hz admittance by about 1.4e-3.

**Documented convention translations:**
- **Transformer taps** are applied exactly: y/t² (bus1), y (bus2), −y/t
  (off-diagonals), in the stamps and the history update.
- **Shunts:** `shnt_gb = 100·(g + jb)`.
- **Loads** are not ParaEMT branches (`loadmodel_option = 2`). They are the
  preregistered constant-power realization (prereg §3.3).

**60-Hz positive-sequence admittance** assembled from the element values
stored in ParaEMT's branch table (R, X or L, C, tap):

| Δt | max\|ΔY\|/max\|Y\| vs canonical | trapezoidal warp (discrete vs continuous, not gated) | companion = pure trapezoidal |
|---|---|---|---|
| 25 µs | 1.7e-16 | 7.4e-6 | 2.3e-16 |
| 50 µs | 1.7e-16 | 3.0e-5 | 2.1e-16 |
| 100 µs | 1.7e-16 | 1.2e-4 | 2.2e-16 |

The warp matches (ω0Δt)²/12 exactly. It is a property of the trapezoidal rule,
not a data difference.

## 3. Equilibrium

Test: the base portfolio at P4, a 2 s no-event run at 50 µs, positive-sequence
phasors averaged over the last 0.5 s.

| quantity | EMT vs canonical power flow | gate |
|---|---|---|
| max \|ΔV\| | 1.5e-6 pu | ≤ 1e-4 |
| max \|Δθ\| (relative to bus 39) | 5.9e-6 rad | ≤ 1e-3 |
| max \|ΔP\| (devices) | 1.7e-5 pu | reported |
| max \|ΔQ\| (devices) | 5.8e-5 pu | reported |
| max speed deviation | 8.3e-8 pu | reported |

Per-bus and per-device tables: `equilibrium_bus_voltages.csv` and
`equilibrium_device_pq.csv`.

## 4. Conclusion

The frozen TX4 IEEE-39 network and its canonical operating point are
reproduced in ParaEMT to the preregistered tolerance, after documented
convention translations only:
- the transformer tap, which is a patch to an upstream omission;
- the removal of numerical damping;
- the per-unit shunt scaling.

No physical parameter was changed or tuned.
