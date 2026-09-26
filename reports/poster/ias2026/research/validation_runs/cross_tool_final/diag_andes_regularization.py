"""DIAGNOSTIC ONLY: is the Python-vs-ANDES static residual ANDES's +1e-8 impedance
regularization?

ANDES 2.0.0 Line equations use yhk = u/((r+1e-8) + 1j*(x+1e-8)) (andes/models/line/
line.py:200), while Line.build_ybus (line.py:268) and the project model use u/(r+jx).
Here the same +1e-8 is applied to r and x of every branch in a TEMPORARY copy of the
canonical JSON (written to supp/, never to configs/). Python is re-solved and compared
with ANDES solved at tol 1e-12. The frozen canonical model is not changed.

Writes supp/andes_regularization_diag.json. Usage (tx3-andes venv):
    python diag_andes_regularization.py <research_root>
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import andes
import numpy as np

repo = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(repo / "src"))
from ibr_cycles.models.ieee39_network import (  # noqa: E402
    load_network,
    solve_power_flow,
)

HERE = Path(__file__).resolve().parent
out = HERE / "supp"
out.mkdir(exist_ok=True)
cfg = repo / "configs/ias2026/ieee39_network.json"
payload = json.loads(cfg.read_text())
for line in payload["lines"]:
    line["r"] = float(line["r"]) + 1e-8
    line["x"] = float(line["x"]) + 1e-8
tmp = out / "ieee39_network_andes_regularized_DIAGNOSTIC.json"
tmp.write_text(json.dumps(payload))

andes.config_logger(stream_level=40)
case = Path(andes.__file__).resolve().parent / "cases/ieee39/ieee39_full.xlsx"
ss = andes.load(str(case), setup=True, no_output=True)
ss.PFlow.config.tol = 1e-12
ss.PFlow.config.max_iter = 50
ss.PFlow.run()
bus_ids = [int(b) for b in ss.Bus.idx.v]
v = np.asarray(ss.Bus.v.v)
a = np.asarray(ss.Bus.a.v)
g = {
    int(b): complex(p, q)
    for b, p, q in zip(ss.PV.bus.v, ss.PV.p.v, ss.PV.q.v, strict=False)
}
g.update(
    {
        int(b): complex(p, q)
        for b, p, q in zip(ss.Slack.bus.v, ss.Slack.p.v, ss.Slack.q.v, strict=False)
    }
)


def compare(path: Path) -> dict:
    net = load_network(path)
    pf = solve_power_flow(net)
    pos = [net.position(b) for b in bus_ids]
    V = pf.voltages[pos]
    gen = {
        b: pf.injection(b, net.ybus) + net.loads.get(b, 0j)
        for b in list(net.pv) + [net.slack_bus]
    }
    off = np.median(a - np.angle(V))
    return {
        "python_pf_mismatch": float(pf.max_mismatch),
        "vm_max_abs_error_pu": float(np.max(np.abs(v - np.abs(V)))),
        "va_max_gauge_aligned_error_rad": float(np.max(np.abs(a - off - np.angle(V)))),
        "pv_p_max_abs_error_pu": float(
            max(abs(g[b].real - gen[b].real) for b in net.pv)
        ),
        "pv_q_max_abs_error_pu": float(
            max(abs(g[b].imag - gen[b].imag) for b in net.pv)
        ),
        "slack_p_abs_error_pu": float(
            abs(g[net.slack_bus].real - gen[net.slack_bus].real)
        ),
        "slack_q_abs_error_pu": float(
            abs(g[net.slack_bus].imag - gen[net.slack_bus].imag)
        ),
    }


res = {
    "andes_tol": 1e-12,
    "andes_final_mismatch": float(ss.PFlow.mis[-1]),
    "python_canonical_vs_andes": compare(cfg),
    "python_with_andes_1e-8_regularization_vs_andes": compare(tmp),
    "source": "andes/models/line/line.py:200 yhk = u/((r+1e-8) + 1j*(x+1e-8))",
}
(out / "andes_regularization_diag.json").write_text(json.dumps(res, indent=2))
print(json.dumps(res, indent=2))
