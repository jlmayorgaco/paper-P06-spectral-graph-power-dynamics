"""Supplementary Phase C checks: live ANDES static reproduction.

The pack's run_andes_dynamic.py writes live bus voltages but compares nothing, hashes
nothing and does not look at shunts or taps. This script:

- records the ANDES version, the exact case file and its SHA256, and compares that
  hash with the source hash stored in the canonical JSON;
- runs ANDES PFlow with the case's own configuration (tol 1e-6, as in F1 and in the
  frozen reference) and compares Vm, Va and generator P/Q with (i) the frozen ANDES
  reference JSON, (ii) the frozen F1 live power flow and (iii) the project Python PF;
- builds the Ybus LIVE from ANDES (System.build_ybus: Line + Shunt) and compares it
  with the project Ybus; the stale configs/ias2026/ieee39_ybus_andes.npy is only
  reported, never used for a verdict;
- checks shunt and tap semantics: every ANDES bus power balance recomputed from the
  live Ybus must equal ANDES's own device injections, and Line b1/b2/g1/g2 must be 0;
- hashes the operating point;
- DIAGNOSTIC ONLY: re-solves the same case with PFlow tol 1e-12 to see whether the
  ANDES-vs-Python slack-P residual is set by ANDES's convergence tolerance.

Writes supp/andes_static_supp.json. Usage (tx3-andes venv):
    python supp_andes_static.py <research_root>
"""

from __future__ import annotations

import hashlib
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
cfg = repo / "configs/ias2026"
payload = json.loads((cfg / "ieee39_network.json").read_text())
ref = json.loads((cfg / "ieee39_andes_powerflow_reference.json").read_text())
f1_live = json.loads((repo / "results/F1/andes_live_powerflow.json").read_text())

andes.config_logger(stream_level=40)
case = Path(andes.__file__).resolve().parent / "cases/ieee39/ieee39_full.xlsx"
case_sha = hashlib.sha256(case.read_bytes()).hexdigest()


def solve(tol: float | None):
    ss = andes.load(str(case), setup=True, no_output=True)
    if tol is not None:
        ss.PFlow.config.tol = tol
        ss.PFlow.config.max_iter = 50
    ss.PFlow.run()
    return ss


def gens(ss) -> dict[int, complex]:
    g = {}
    for b, p, q in zip(ss.PV.bus.v, ss.PV.p.v, ss.PV.q.v, strict=False):
        g[int(b)] = complex(p, q)
    for b, p, q in zip(ss.Slack.bus.v, ss.Slack.p.v, ss.Slack.q.v, strict=False):
        g[int(b)] = complex(p, q)
    return g


def gauge_err(a, b):
    off = np.median(a - b)
    return float(np.max(np.abs(a - off - b)))


ss = solve(None)
res: dict = {
    "andes_version": andes.__version__,
    "case_file": str(case),
    "case_sha256": case_sha,
    "canonical_source_sha256": payload["source"]["sha256"],
    "case_hash_matches_canonical_source": case_sha == payload["source"]["sha256"],
    "pflow_tol": float(ss.PFlow.config.tol),
    "pflow_converged": bool(ss.PFlow.converged),
    "pflow_iterations": int(ss.PFlow.niter),
    "pflow_final_mismatch": float(ss.PFlow.mis[-1]),
}
bus_ids = [int(b) for b in ss.Bus.idx.v]
v = np.asarray(ss.Bus.v.v, float)
a = np.asarray(ss.Bus.a.v, float)

# project Python PF in the same bus order
pynet = load_network(cfg / "ieee39_network.json")
pf = solve_power_flow(pynet)
pos = [pynet.position(b) for b in bus_ids]
Vpy = pf.voltages[pos]
py_gen = {
    b: pf.injection(b, pynet.ybus) + pynet.loads.get(b, 0j)
    for b in list(pynet.pv) + [pynet.slack_bus]
}
g_live = gens(ss)

assert [int(b) for b in ref["bus_idx"]] == bus_ids
res["vs_frozen_reference"] = {
    "vm_max_abs_error_pu": float(np.max(np.abs(v - np.asarray(ref["v"])))),
    "va_max_abs_error_rad": float(np.max(np.abs(a - np.asarray(ref["a"])))),
    "pv_p_max_abs_error_pu": float(
        max(
            abs(g_live[int(b)].real - p)
            for b, p in zip(ref["pv_bus"], ref["pv_p"], strict=False)
        )
    ),
    "pv_q_max_abs_error_pu": float(
        max(
            abs(g_live[int(b)].imag - q)
            for b, q in zip(ref["pv_bus"], ref["pv_q"], strict=False)
        )
    ),
    "slack_p_abs_error_pu": float(
        abs(g_live[int(ref["slack_bus"][0])].real - ref["slack_p"][0])
    ),
    "slack_q_abs_error_pu": float(
        abs(g_live[int(ref["slack_bus"][0])].imag - ref["slack_q"][0])
    ),
}
fl_v = np.asarray(f1_live.get("v", f1_live.get("bus_v", [])), float)
fl_a = np.asarray(f1_live.get("a", f1_live.get("bus_a", [])), float)
if fl_v.size == v.size:
    res["vs_frozen_F1_live_powerflow"] = {
        "vm_max_abs_error_pu": float(np.max(np.abs(v - fl_v))),
        "va_max_abs_error_rad": float(np.max(np.abs(a - fl_a))),
    }
else:
    res["vs_frozen_F1_live_powerflow"] = {
        "note": f"keys {sorted(f1_live)}; not compared"
    }


def vs_python(v, a, g):
    return {
        "vm_max_abs_error_pu": float(np.max(np.abs(v - np.abs(Vpy)))),
        "va_max_gauge_aligned_error_rad": gauge_err(a, np.angle(Vpy)),
        "pv_p_max_abs_error_pu": float(
            max(abs(g[b].real - py_gen[b].real) for b in pynet.pv)
        ),
        "pv_q_max_abs_error_pu": float(
            max(abs(g[b].imag - py_gen[b].imag) for b in pynet.pv)
        ),
        "slack_p_abs_error_pu": float(
            abs(g[pynet.slack_bus].real - py_gen[pynet.slack_bus].real)
        ),
        "slack_q_abs_error_pu": float(
            abs(g[pynet.slack_bus].imag - py_gen[pynet.slack_bus].imag)
        ),
    }


res["vs_python_pf"] = vs_python(v, a, g_live)

# live Ybus from ANDES itself (Line + Shunt)
Y_live = np.array(andes.shared.matrix(ss.build_ybus()), dtype=complex)
Y_py = pynet.ybus[np.ix_(pos, pos)]
res["ybus_live_andes_vs_python_max_abs_error"] = float(np.max(np.abs(Y_live - Y_py)))
ybus_models = sorted(ss.exist.ybus)
res["ybus_contributing_models"] = ybus_models
stale = np.load(cfg / "ieee39_ybus_andes.npy")
d = np.abs(stale - Y_live)
res["stale_npy_vs_live_max_abs_diff"] = float(d.max())
res["stale_npy_differs_at_buses"] = [
    bus_ids[i] for i in range(len(bus_ids)) if d[i, i] > 1e-9
]
res["stale_npy_used_for_any_verdict"] = False

# shunt / tap semantics
res["andes_shunts"] = [
    {"bus": int(b), "g": float(g), "b": float(bb)}
    for b, g, bb in zip(ss.Shunt.bus.v, ss.Shunt.g.v, ss.Shunt.b.v, strict=False)
]
res["canonical_shunts"] = payload["shunts"]
res["line_b1_b2_g1_g2_all_zero"] = bool(
    all(
        np.all(np.asarray(getattr(ss.Line, k).v) == 0.0)
        for k in ("b1", "b2", "g1", "g2")
    )
)
tap_rows = [
    (int(b1), int(b2), float(t))
    for b1, b2, t in zip(ss.Line.bus1.v, ss.Line.bus2.v, ss.Line.tap.v, strict=False)
]
canon_rows = [
    (int(x["bus1"]), int(x["bus2"]), float(x["tap"])) for x in payload["lines"]
]
res["line_order_bus1_bus2_tap_identical_to_canonical"] = tap_rows == canon_rows
V = v * np.exp(1j * a)
S_net = V * np.conj(Y_live @ V)  # injection computed from the live Ybus
S_dev = np.zeros(len(bus_ids), complex)
for b, s in g_live.items():
    S_dev[bus_ids.index(b)] += s
for b, p, q in zip(ss.PQ.bus.v, ss.PQ.p0.v, ss.PQ.q0.v, strict=False):
    S_dev[bus_ids.index(int(b))] -= complex(p, q) * 1.0  # constant-power at PF
res["bus_balance_live_ybus_vs_andes_devices_max_abs_pu"] = float(
    np.max(np.abs(S_net - S_dev))
)

# operating-point hash
op = {
    "bus_idx": bus_ids,
    "v": [round(x, 8) for x in v],
    "a": [round(x, 8) for x in a],
    "gen": {
        str(b): [round(s.real, 8), round(s.imag, 8)] for b, s in sorted(g_live.items())
    },
}
res["operating_point_sha256_8dp"] = hashlib.sha256(
    json.dumps(op, sort_keys=True).encode()
).hexdigest()
op_ref = {
    "bus_idx": bus_ids,
    "v": [round(x, 8) for x in ref["v"]],
    "a": [round(x, 8) for x in ref["a"]],
}
op_live = {"bus_idx": bus_ids, "v": op["v"], "a": op["a"]}
res["voltage_hash_live_equals_frozen_reference_8dp"] = (
    hashlib.sha256(json.dumps(op_live).encode()).hexdigest()
    == hashlib.sha256(json.dumps(op_ref).encode()).hexdigest()
)

# diagnostic: tight PFlow tolerance
ss2 = solve(1e-12)
v2 = np.asarray(ss2.Bus.v.v, float)
a2 = np.asarray(ss2.Bus.a.v, float)
res["DIAGNOSTIC_tol_1e-12"] = {
    "converged": bool(ss2.PFlow.converged),
    "iterations": int(ss2.PFlow.niter),
    "final_mismatch": float(ss2.PFlow.mis[-1]),
    **vs_python(v2, a2, gens(ss2)),
}
(out / "andes_static_supp.json").write_text(json.dumps(res, indent=2))
print(json.dumps(res, indent=2))
