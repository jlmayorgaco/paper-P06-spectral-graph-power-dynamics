# H4 mechanism campaign run

Status: **PASS**

Run root: `C:\Users\walla\Documents\Github\paper-P06-spectral-graph-power-dynamics\reports\poster\ias2026\research\bnd_h4_mechanism\results\20260925T203454_265f773a_h4_crossmode_v1`

The primary TX4 contextual-return campaign was executed through the frozen implementation with output redirection set before evaluation.

## Gate status

- G0 isolated immutable run root: **PASS** — observed `C:\Users\walla\Documents\Github\paper-P06-spectral-graph-power-dynamics\reports\poster\ias2026\research\bnd_h4_mechanism\results\20260925T203454_265f773a_h4_crossmode_v1`
- G1 frozen source hashes: **PASS** — observed `{'models/ieee39_case.py': {'observed': '1f38a8e541780ab7205809dc1ba58bdafb6c923933f5772d46407064edc8c722', 'expected': '1f38a8e541780ab7205809dc1ba58bdafb6c923933f5772d46407064edc8c722', 'match': True}, 'models/ieee39_devices.py': {'observed': 'bcf0c76be65b171b580684065cde4eb88d0433c6b4d7beefec966f8711036d90', 'expected': 'bcf0c76be65b171b580684065cde4eb88d0433c6b4d7beefec966f8711036d90', 'match': True}, 'models/ieee39_network.py': {'observed': 'f42144cd02a6dc779ed690fd866ad026f33eb7844584c7f10281fc7d98e3fe13', 'expected': 'f42144cd02a6dc779ed690fd866ad026f33eb7844584c7f10281fc7d98e3fe13', 'match': True}, 'models/port_admittance.py': {'observed': 'b0fa694dd0111be158133a9cf8c4e36be670e50ef562dcd63f89d7ff87b8d4ef', 'expected': 'b0fa694dd0111be158133a9cf8c4e36be670e50ef562dcd63f89d7ff87b8d4ef', 'match': True}, 'models/port_core.py': {'observed': '38c54b1a70f87fbbe2c301a717e6dd6b4112b3a9cf8deb663426ce08b729bb9b', 'expected': '38c54b1a70f87fbbe2c301a717e6dd6b4112b3a9cf8deb663426ce08b729bb9b', 'match': True}, 'experiments/tx4_contextual_return.py': {'observed': '77d4ed0f76f8c4a6ec369bc30b306361cff9b789c6f432e287d5f180d9210b74', 'expected': '77d4ed0f76f8c4a6ec369bc30b306361cff9b789c6f432e287d5f180d9210b74', 'match': True}}`
- G2 equilibrium residuals: **PASS** — observed `{'max_f': 5.580589596813811e-13, 'max_g': 9.379164112033322e-13}`
- G3 H4 dynamic-state dimension: **PASS** — observed `86`
- G4 DAE/index-1 conditioning: **PASS** — observed `2793.5943065470565`
- G5 transverse target mode: **PASS** — observed `{'frequency_hz': 0.6222796695036534, 'numpy_scipy_error': 0.0}`
- G6 H4 P4 instability: **PASS** — observed `{'alpha_s-1': 0.12700646782830968, 'frequency_hz': 0.6222796695036534}`
- G7 proper-subset minimality: **PASS** — observed `{'rows': 15, 'all_stable': True}`
- G8 eigenvalue boundary: **PASS** — observed `{'g_root': 0.20768140519037842, 'nearest_grid_alpha': 6.472444463139398e-11}`
- G9 port/Schur identities: **PASS** — observed `{'max_schur': 1.343782689223772e-15, 'max_port': 1.4304896245381993e-16}`
- G10 collective not local: **PASS** — observed `{'min_local_sigma': 0.9999999999999997, 'min_collective_sigma': 2.6081328337715837e-08}`
- G11 derivative and numerical audit: **PASS** — observed `{'audit': {'n_portfolios': 9, 'tolerance_rows': 27, 'max_numpy_scipy_eig_error': 0.0, 'pass': True, 'descriptor': 'NOT_AVAILABLE'}, 'derivative_rows': 4, 'max_derivative_error': np.float64(2.662806050493975e-07)}`

## Core summary

```json
{
  "audit": {
    "n_portfolios": 9,
    "tolerance_rows": 27,
    "max_numpy_scipy_eig_error": 0.0,
    "pass": true,
    "descriptor": "NOT_AVAILABLE"
  },
  "core": {
    "g_eigen_boundary": 0.20768140519037842,
    "g_return_minimum": 0.20768140481809155,
    "g_return_grid_minimum": 0.2076814045,
    "g_boundary_difference": 3.7228686800006017e-10,
    "root_iterations": 26,
    "alpha_P4": 0.12700646782830968,
    "alpha_after": -0.017384676061556643,
    "frequency_root_hz": 0.706424788006417,
    "local_sigma_min_root": 1.0,
    "collective_sigma_min_root": 4.3936703686787105e-08,
    "eta_H": 0.12083120539406678,
    "tau_H": 1.3069396516158314,
    "nearest_Q_eigenvalue_to_minus_one": [
      -1.0000000690223476,
      8.670802298382796e-09
    ],
    "max_boundary_schur_residual": 2.2959283690665653e-16,
    "max_return_derivative_error": 2.662806050493975e-07,
    "minimality_all_proper_stable": true,
    "newton_or_bisection_iterations": 26
  }
}
```
