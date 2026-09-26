# ruff: noqa: E501  -- evidence labels, docstrings and verbatim source quotes kept on one line
"""Phase E1, internal side: equation-equivalent branch sensitivity on the frozen holdout.

Configuration (the F1 R3 equation-equivalent case, as prescribed by the frozen
protocol CLAUDE_FINAL_ANDES_NETWORK_VALIDATION.md): first-order AVR, stabilizer gain
0, no governor, constant-power loads, and the flagship {30,33,35,37} replaced by
static constant-power injections ("static_power", q_policy "matched"). No parameter
is tuned.

For each frozen holdout line (loaded from holdout_lines.json; never regenerated):
  - gamma_e multiplies the complete branch two-port admittance (series and charging,
    tap unchanged), gamma = 1 -/+ 0.002, exactly as in the frozen F2c script
    (ybus_scaled is copied verbatim);
  - the case is re-solved (PF, equilibrium, linearization) at each gamma;
  - DAE:  d lambda / d gamma_e by central difference, tracking the transverse
    eigenvalue nearest the gamma = 1 critical band mode (trans_eig, verbatim);
  - PORT: the action space is rebuilt at each perturbed equilibrium and
    ds*/d gamma_e = -(d mu / d gamma) / (d mu / d s) at s* = lambda_0 (qs and mu
    verbatim from F2c). Evaluating at the replaced-system eigenvalue is admissible
    because the port operator is built on the BASE resolvent.

Also records the port identities used by Phase F:
  - |mu(lambda_0) + 1|  (the characteristic zero of the port closure);
  - det(T_S)/det(T_0) vs det(I + M) at s = j 2 pi 0.7.

Writes phaseE/internal_E1.csv and phaseE/internal_E1.json.
Usage (tx3-analysis venv): python phaseE_internal.py <research_root> <holdout_lines.json>
"""

from __future__ import annotations

import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

repo = Path(sys.argv[1]).resolve()
holdout = json.loads(Path(sys.argv[2]).read_text())
SEL = [int(x) for x in holdout["lines"]]
assert holdout["seed"] == 20260911 and SEL == [
    3,
    9,
    15,
    17,
    22,
    26,
    31,
    32,
    35,
    39,
    42,
    44,
]
sys.path[:0] = [str(repo / "src"), str(repo / "experiments")]
from _v2c_common import BAND_HZ  # noqa: E402
from ibr_cycles.certification.symmetry import (  # noqa: E402
    frequency_partner,
    rotation_generator,
)
from ibr_cycles.certification.transverse import transverse_operator  # noqa: E402
from ibr_cycles.dynamics.modal_family import band_candidates  # noqa: E402
from ibr_cycles.dynamics.modes import eigen_analysis  # noqa: E402
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case  # noqa: E402
from ibr_cycles.models.ieee39_network import load_network  # noqa: E402
from ibr_cycles.models.port_admittance import build_action_space  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "phaseE"
OUT.mkdir(exist_ok=True)
CORE = (30, 33, 35, 37)
SERVICES = {"pss": 0.0}  # F1 stage "R3 first-order AVR"
EPS = 0.002
payload = json.loads((repo / "configs/ias2026/ieee39_network.json").read_text())
net = load_network(repo / "configs/ias2026/ieee39_network.json")


# ---- verbatim from dynamic_forest_line_port_holdout.py (frozen F2c) ----
def ybus_scaled(line_idx, gamma):
    order = [int(b["idx"]) for b in payload["buses"]]
    index = {b: i for i, b in enumerate(order)}
    n = len(order)
    y = np.zeros((n, n), complex)
    for ell, line in enumerate(payload["lines"]):
        if float(line["u"]) == 0:
            continue
        f, t = index[int(line["bus1"])], index[int(line["bus2"])]
        scale = gamma if ell == line_idx else 1.0
        series = scale / complex(line["r"], line["x"])
        charging = scale * complex(line["g"], line["b"]) / 2
        m = float(line["tap"]) * np.exp(1j * float(line["phi"]))
        m2 = abs(m) ** 2
        y[f, f] += (series + charging) / m2
        y[t, t] += series + charging
        y[f, t] += -series / np.conj(m)
        y[t, f] += -series / m
    for sh in payload["shunts"]:
        y[index[int(sh["bus"])], index[int(sh["bus"])]] += complex(sh["g"], sh["b"])
    return y


def qs(space, s):
    m = space.m(s)
    eye = np.eye(m.shape[0], dtype=complex)
    total = eye + m
    sb = np.zeros_like(total)
    for k in range(space.order):
        sl = slice(2 * k, 2 * k + 2)
        sb[sl, sl] = total[sl, sl]
    return np.linalg.solve(sb, total) - eye


def mu(space, s, target=-1):
    vals = np.linalg.eigvals(qs(space, s))
    return vals[np.argmin(np.abs(vals - target))]


def trans_eig(case, target):
    rx, _ = rotation_generator(case.dae, case.equilibrium.z)
    w = frequency_partner(case.dae).w
    ev = np.linalg.eigvals(transverse_operator(case.system.A, rx, w).a_perp)
    c = ev[ev.imag >= 0]
    return c[np.argmin(np.abs(c - target))]


# ---- end verbatim ----


def solve(members, network):
    plan = (
        ReplacementPlan.of({b: 1.0 for b in members}, device="static_power")
        if members
        else ReplacementPlan.of({})
    )
    return solve_case(plan, network=network, machine_services=SERVICES)


assert float(np.max(np.abs(ybus_scaled(-1, 1.0) - net.ybus))) < 1e-12
base0, full0 = solve((), net), solve(CORE, net)
band = band_candidates(eigen_analysis(full0.system.A), BAND_HZ)
crit = max(band, key=lambda m: m.real)
lam0 = trans_eig(full0, complex(crit.value))
sp0 = build_action_space(base0, full0, CORE)
hs = 1e-5
mu0 = mu(sp0, lam0)
mu_s = (mu(sp0, lam0 + hs) - mu(sp0, lam0 - hs)) / (2 * hs)
# determinant identity off both spectra
s_test = 1j * 2 * np.pi * 0.7
lhs = np.linalg.det(sp0.ts.evaluate(s_test)) / np.linalg.det(sp0.t0.evaluate(s_test))
rhs = np.linalg.det(np.eye(sp0.dimension) + sp0.m(s_test))
info = {
    "configuration": "F1 R3: first-order AVR, pss gain 0, no governor, constant-power loads, "
    "flagship 30+33+35+37 as static_power (q matched)",
    "holdout_seed": holdout["seed"],
    "holdout_lines": SEL,
    "eps": EPS,
    "lambda0_real": lam0.real,
    "lambda0_imag": lam0.imag,
    "lambda0_freq_hz": lam0.imag / (2 * np.pi),
    "band_critical_alpha": float(crit.real),
    "port_mu_at_lambda0_plus_1_abs": float(abs(mu0 + 1)),
    "det_identity_rel_error_at_j2pi0p7": float(abs(lhs - rhs) / abs(rhs)),
    "mu_s_real": mu_s.real,
    "mu_s_imag": mu_s.imag,
}
print(json.dumps(info, indent=2), flush=True)

rows = []
for li in SEL:
    spaces, lams, lams_full = [], [], []
    for gam in (1 - EPS, 1 + EPS):
        nn = replace(net, ybus=ybus_scaled(li, gam))
        b, f = solve((), nn), solve(CORE, nn)
        spaces.append(build_action_space(b, f, CORE))
        lams.append(trans_eig(f, lam0))
        ev = np.linalg.eigvals(f.system.A)
        lams_full.append(ev[np.argmin(np.abs(ev - lam0))])
    mug = (mu(spaces[1], lam0) - mu(spaces[0], lam0)) / (2 * EPS)
    ds = -mug / mu_s
    dl = (lams[1] - lams[0]) / (2 * EPS)
    dl_full = (lams_full[1] - lams_full[0]) / (2 * EPS)
    line = payload["lines"][li]
    rows.append(
        {
            "line_index": li,
            "line": f"L{li:02d}:{int(line['bus1'])}-{int(line['bus2'])}",
            "internal_dae_dlambda_real": dl.real,
            "internal_dae_dlambda_imag": dl.imag,
            "internal_fullA_dlambda_real": dl_full.real,
            "internal_fullA_dlambda_imag": dl_full.imag,
            "internal_port_ds_real": ds.real,
            "internal_port_ds_imag": ds.imag,
            "port_vs_dae_abs": abs(ds - dl),
        }
    )
    print(rows[-1], flush=True)

pd.DataFrame(rows).to_csv(OUT / "internal_E1.csv", index=False)
(OUT / "internal_E1.json").write_text(json.dumps(info, indent=2))
