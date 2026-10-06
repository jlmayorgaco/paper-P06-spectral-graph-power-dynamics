# Experiment F0 — Physical IEEE-39 generator-port graph

**F0_STATUS: PASS**

## Preregistered primary operator

Before any controller analysis, the primary graph was fixed as the real Hermitian part of the passive complex branch admittance after Kron reduction onto generator buses: `Lc = Re(Yport)`. The branch-only network excludes generator, load, and converter device admittances. Generator-port order is `30, 31, 32, 33, 34, 35, 36, 37, 38, 39`; all other buses are eliminated. Values are per-unit conductance on the archived system base. The source-side transformer ratio is retained exactly. No artificial ground, clipping, or absolute-eigenvalue operation is used.

## Physical and numerical checks

- `Ybus` and its Kron port matrix are complex-symmetric to relative residuals `0.0` and `2.057229442163321e-16`.

- `Lc` eigenvalue range: `0.001127537540072563` to `5.0133384196549695`; zero modes: `0`; condition number: `4446.271845931033`.

- The complex Kron map reproduces port currents for seven deterministic port-voltage trials with relative residual `4.0750477862428105e-16` (maximum absolute residual `2.5635335727383778e-14`).

- Minimum tested real power over the full passive network: `205.64926878309018`; minimum port quadratic form: `1.9485493444618582`.

- The candidate table classifies reactive susceptance and ExpC C0/C33 synchronizing backbones separately. The C33 backbone is tested directly; no indefinite operator is modified to make it PSD.

- ExpC critical modes are compared with the graph basis only as a diagnostic: `computed from frozen ExpC C0 critical poles with M-weighted q projection`. The M-weighted mode-energy and principal-angle tables are retained separately.

## Gate

**F0_PASS.** A network-derived PSD operator exists and the port-current relation is validated. This establishes a mathematically admissible graph spectrum; it does not establish that the detailed GFL model is modal in this basis.

Inputs are the frozen CSV copies under `reports/experiment_D/inputs`; SHA-256 values are in `RESULTS_EXP_F0.toml`. The experiment code imports no PowerDynamics module.


## Synchronizing-backbone classification

The frozen ExpC all-SG C0 synchronizing backbone has spectrum [-1.18515e-17, 0.302968] and is PSD within tolerance. The bus-33 GFL C33 backbone has spectrum [-19.1495, 1439.58] and is indefinite. The C33 operator is not sent through a real square root.
