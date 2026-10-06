# Gate 3 — IEEE 68-bus / NETS–NYPS replication

- Protocol: `configs/ieee68/G3_preregistration.yaml`, committed in `d9fa097e`
  before any 68-bus code existed.
- Data: `configs/ieee68/ieee68_network.json`, transcribed by
  `experiments/G3_import_ieee68.py` from Singh & Pal (2013), the IEEE PES TF
  benchmark report v3.3, Appendix B. The PDF is not redistributed.
- Model: `src/ibr_cycles/models/ieee68_devices.py`, the report's equations
  (1)–(10) as printed.
- Run: `experiments/G3_ieee68.py`. Post-hoc descriptive scripts:
  `G3b_margin_posthoc.py`, `G3c_family_anatomy_posthoc.py`.
- Data: `results/G3/`. Tests: `tests/test_ieee68_model.py` (7).

## 1. Eligibility (declared checks)

| check | result | verdict |
|---|---|---|
| power flow vs Table 1 (68 buses) | max `|dV| = 5.0e-5` pu, max `|d angle| = 5.0e-5` deg | **REPRODUCED** |
| four inter-area modes vs Table 4 | 0.314 Hz / 33.54 %, 0.520 / 3.62 %, 0.591 / 9.63 %, 0.779 / 3.38 % | **REPRODUCED** |
| all 15 electromechanical modes vs Table 4 | max `|df| = 5.0e-4` Hz, max `|d zeta| = 4.5e-4` percentage points | reproduced |
| base stability | abscissa `-0.0556`; exactly 2 reference zeros | stable |
| distance to documented limits | smallest margin 3.0 (regulator output) | inside every limit |

Our implementation is the benchmark as published. No parameter was changed,
added or tuned. The data file lists 83 branch rows, although its header says
"86-lines".

**Exact assembly.** `A(g, k)` is affine in the policy gain `g` and the
regulator-gain scale `k`, and is assembled from three direct solves per subset.
On 60 random off-grid points the assembled and direct-path hypergraphs agree
60 of 60 times. The matrices agree to a relative `3e-16`.

## 2. Preregistered map (candidates G9, G6, G3, G4; 16 portfolios; 31 x 31 over g in [0,1], k in [0.5, 2])

| quantity | value |
|---|---|
| `H_RHP` | **EMPTY at 961 of 961 nodes** |
| `H_IA` | EMPTY at 961 of 961 nodes |
| RHP boundaries on pure-policy lines | 0 |
| base-unstable nodes | 0 |
| spot checks, direct path vs nearest node | 60 / 60 |

**Preregistered decision: NOT REPRODUCED.** With the four preregistered
candidates, no reactive policy or excitation gain in the window makes any
portfolio small-signal unstable, so there is no policy dependence to observe.

*Post-hoc, descriptive (`G3b`).* At every node, the worst of the 16 portfolios
has spectral abscissa between `-0.0558` and `-0.0500`. The upper value is the
converter regulator's own leak pole at `g = 0`. Otherwise the worst abscissa is
essentially the base's `-0.0556`. The four-candidate map is not barely
composable: its margin is the base system's margin.

## 3. Preregistered secondary check: the 12 physical plants G1–G12 (4 096 portfolios), k = 1

| `g` | `kappa_RHP` | hyperedges | hyperedges inside the 4 candidates | RHP mechanism (post-hoc `G3c`) |
|---|---|---|---|---|
| 0 (fixed Q) | **6** | 259 (sizes 6–9) | 0 | all in band: the 0.58–0.62 Hz inter-area mode dominated by G13's rotor (`delta_sg13`, `omega_sg13`), surviving-machine participation at least 0.89. This is the 0.591 Hz NYPS mode of Table 4. |
| 0.25 | **9** | 1: `1+2+...+9` (every NETS plant) | 0 | oscillatory, outside the band: 2.34 Hz, `Re = +2.28`, converter PLL (`theta_pll`), converter participation 0.92 |
| 1 | **11** | 2 (sizes 11) | 0 | oscillatory, outside the band: 3.50–3.55 Hz, converter PLL, converter participation at least 0.97 |

**What this shows.**

- On the 68-bus system the incompatibility structure exists only at high
  replacement levels: at least 6 of the 12 physical plants.
- Its hypergraph **does** change with the reactive policy (`kappa_RHP`
  6 → 9 → 11), and it changes mechanism along the way. At fixed Q the
  inter-area mode fails, as on IEEE-39. At `g = 0.25` and `g = 1` the failure
  is a converter-synchronization mode at 2–3.5 Hz, which a 0.3–1.5 Hz band
  would miss entirely.
- The four preregistered candidates are a composable sub-fleet at every point
  tested. Restricting to them hides every hyperedge, since the smallest has 6
  members.
- This is three points on one line, not a map. No boundary was localized and
  no port-closure test was run on the 68-bus system.

## 4. Answer to the gate

| question | answer |
|---|---|
| documented data obtainable without inventing parameters? | **yes**; the model reproduces the report's power flow and modes to within `5e-4` |
| preregistered replication (4 candidates, g x k map) | **NOT REPRODUCED**: every portfolio is stable everywhere |
| does `H_RHP` depend on policy on this system at all? | yes, but only in the 12-plant family at `kappa_RHP >= 6`: 259 → 1 → 2 hyperedges between `g = 0, 0.25, 1` (secondary check, three points) |
| band vs RHP | identical on the map (both empty); in the 12-plant family the band would miss the `g = 0.25` and `g = 1` hyperedges (PLL modes at 2.3–3.5 Hz) |
| limitation | one excitation model (the benchmark's), one dispatch, one converter tuning; the converter rating rule `|S|/0.8` was declared, not measured |

The structure is **system-dependent**. On IEEE-39 it appears at `kappa` 1–4
among four plants. On Kundur no policy removes it. On the documented 68-bus
system it needs at least six of twelve plants replaced.
