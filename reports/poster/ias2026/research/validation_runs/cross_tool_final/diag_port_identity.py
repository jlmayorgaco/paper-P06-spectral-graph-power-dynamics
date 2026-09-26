# ruff: noqa: E501  -- evidence labels, docstrings and verbatim source quotes kept on one line
"""DIAGNOSTIC: what limits the port identities in the E1 (R3 static_power) configuration?

T_S = T_0 + U dY U^T holds exactly only if the two operators share every non-port
block. T_0 and T_S are built at separately solved equilibria (solve_case tol), so the
load linearizations at non-port buses can differ at the equilibrium-tolerance level.
This script reports the residual R = T_S - T_0 - U dY U^T and both identities at the
default equilibrium tolerance (1e-9) and at 1e-12. The model is unchanged; only the
equilibrium solver tolerance differs.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

repo = Path(sys.argv[1]).resolve()
sys.path[:0] = [str(repo / "src"), str(repo / "experiments")]
from _v2c_common import BAND_HZ  # noqa: E402
from ibr_cycles.dynamics.modal_family import band_candidates  # noqa: E402
from ibr_cycles.dynamics.modes import eigen_analysis  # noqa: E402
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case  # noqa: E402
from ibr_cycles.models.port_admittance import build_action_space  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
CORE = (30, 33, 35, 37)
S = {"pss": 0.0}
res = {}
for tol in (1e-9, 1e-12):
    b = solve_case(ReplacementPlan.of({}), machine_services=S, tol=tol)
    f = solve_case(
        ReplacementPlan.of({k: 1.0 for k in CORE}, device="static_power"),
        machine_services=S,
        tol=tol,
    )
    sp = build_action_space(b, f, CORE)
    s = 1j * 2 * np.pi * 0.7
    u = sp.selector()
    Rp = sp.ts.evaluate(s) - sp.t0.evaluate(s) - u @ sp.update(s) @ u.T
    Rm = sp.ts.evaluate(s) - sp.t0.evaluate(s) + u @ sp.update(s) @ u.T
    R = Rp if np.max(np.abs(Rp)) < np.max(np.abs(Rm)) else Rm
    port_rows = np.zeros(R.shape[0], bool)
    for bus in CORE:
        p = sp.t0.bus_index[bus]
        port_rows[2 * p : 2 * p + 2] = True
    res.setdefault("residual_off_port_rows", {})[f"{tol:g}"] = float(
        np.max(np.abs(R[~port_rows][:, ~port_rows]))
    )
    lhs = np.linalg.det(sp.ts.evaluate(s)) / np.linalg.det(sp.t0.evaluate(s))
    rhs = np.linalg.det(np.eye(sp.dimension) + sp.m(s))
    band = band_candidates(eigen_analysis(f.system.A), BAND_HZ)
    lam = complex(max(band, key=lambda m: m.real).value)
    q = sp.split(lam)["closest_to_minus_one"]
    zmax = float(np.max(np.abs(b.equilibrium.z - f.equilibrium.z)))
    res[f"tol_{tol:g}"] = {
        "max_abs_T_S_minus_T0_minus_UdYU": float(np.max(np.abs(R))),
        "det_identity_rel_error": float(abs(lhs - rhs) / abs(rhs)),
        "abs_q_eig_plus_1_at_lambda": float(abs(q + 1)),
        "max_abs_voltage_diff_base_vs_flagship_equilibrium": zmax,
    }
print(json.dumps(res, indent=2))
(HERE / "phaseE" / "port_identity_diag.json").write_text(json.dumps(res, indent=2))
