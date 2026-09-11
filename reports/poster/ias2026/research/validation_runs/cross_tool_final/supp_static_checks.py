"""Supplementary static checks the supplied pack does not perform (Phases A and B).

A-supp  Python generator P/Q vs the frozen ANDES 2.0.0 power-flow reference.
B-supp  For both pandapower variants (as supplied + compat fixes; B-corrected):
        - pandapower internal Ybus vs project Ybus
        - generator P/Q from res_gen / res_ext_grid vs Python
        - pandapower NATIVE branch terminal flows (res_line/res_trafo/res_impedance)
          vs the Python canonical pi-model flows. The pack's own branch metric applies
          one formula to both voltage vectors, so it only tests voltages; this one
          uses pandapower's own branch models.

Writes supp/static_supp.json and supp/pandapower_static_parity_long.csv.
Usage (xtool-pandapower venv): python supp_static_checks.py <research_root>
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pandapower as pp  # noqa: E402
import run_pack_pandapower_tapside as tapside  # noqa: E402

compat = tapside.compat
cp = compat.canonical_ppc
repo = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(repo / "src"))
from ibr_cycles.models.ieee39_network import (  # noqa: E402
    load_network,
    solve_power_flow,
)

BASE = 100.0
out = HERE / "supp"
out.mkdir(exist_ok=True)
cfg = repo / "configs/ias2026"
pynet = load_network(cfg / "ieee39_network.json")
pf = solve_power_flow(pynet)
V = pf.voltages
payload = json.loads((cfg / "ieee39_network.json").read_text())

gen_buses = sorted(pynet.pv) + [pynet.slack_bus]
py_gen = {b: pf.injection(b, pynet.ybus) + pynet.loads.get(b, 0j) for b in gen_buses}
rows: list[dict] = []
summary: dict = {}

# ---------------- A-supp: Python vs frozen ANDES reference ----------------
ref = json.loads((cfg / "ieee39_andes_powerflow_reference.json").read_text())
a_err = {"pv_p": [], "pv_q": [], "slack_p": [], "slack_q": []}
for b, p, q in zip(ref["pv_bus"], ref["pv_p"], ref["pv_q"], strict=False):
    s = py_gen[int(b)]
    a_err["pv_p"].append(abs(s.real - p))
    a_err["pv_q"].append(abs(s.imag - q))
    rows += [
        {
            "comparison": "A_python_vs_frozen_andes",
            "quantity": "gen_P_pu",
            "element": int(b),
            "reference": p,
            "other": s.real,
            "abs_error": abs(s.real - p),
        },
        {
            "comparison": "A_python_vs_frozen_andes",
            "quantity": "gen_Q_pu",
            "element": int(b),
            "reference": q,
            "other": s.imag,
            "abs_error": abs(s.imag - q),
        },
    ]
for b, p, q in zip(ref["slack_bus"], ref["slack_p"], ref["slack_q"], strict=False):
    s = py_gen[int(b)]
    a_err["slack_p"].append(abs(s.real - p))
    a_err["slack_q"].append(abs(s.imag - q))
    rows += [
        {
            "comparison": "A_python_vs_frozen_andes",
            "quantity": "slack_P_pu",
            "element": int(b),
            "reference": p,
            "other": s.real,
            "abs_error": abs(s.real - p),
        },
        {
            "comparison": "A_python_vs_frozen_andes",
            "quantity": "slack_Q_pu",
            "element": int(b),
            "reference": q,
            "other": s.imag,
            "abs_error": abs(s.imag - q),
        },
    ]
summary["A_supp"] = {k: float(max(v)) for k, v in a_err.items()}

# ---------------- B-supp ----------------
py_flows = cp.branch_flows_from_voltages(repo, V)  # idx, f, t, Pf, Qf, Pt, Qt (MW/Mvar)


def run_variant(name: str, builder) -> dict:
    ppc = builder(repo)
    net = compat.from_ppc(ppc, f_hz=60, validate_conversion=False)
    pp.runpp(
        net,
        algorithm="nr",
        calculate_voltage_angles=True,
        init="flat",
        enforce_q_lims=False,
        tolerance_mva=1e-10,
        max_iteration=50,
        numba=False,
    )
    res: dict = {"converged": bool(net.converged)}
    # Ybus
    yb = net._ppc["internal"]["Ybus"].toarray()
    lk = net._pd2ppc_lookups["bus"]
    order = np.asarray([lk[b] for b in pynet.bus_idx])
    res["ybus_max_abs_error"] = float(
        np.max(np.abs(yb[np.ix_(order, order)] - pynet.ybus))
    )
    # generators
    pg: dict[int, complex] = {}
    for tbl, rtbl in (
        ("gen", "res_gen"),
        ("ext_grid", "res_ext_grid"),
        ("sgen", "res_sgen"),
    ):
        t, r = getattr(net, tbl), getattr(net, rtbl)
        for i in t.index:
            b = int(t.at[i, "bus"])
            pg[b] = pg.get(b, 0j) + complex(r.at[i, "p_mw"], r.at[i, "q_mvar"]) / BASE
    ep, eq, es_p, es_q = [], [], [], []
    for b in gen_buses:
        s_pp, s_py = pg[b], py_gen[b]
        kind = "slack" if b == pynet.slack_bus else "gen"
        (es_p if kind == "slack" else ep).append(abs(s_pp.real - s_py.real))
        (es_q if kind == "slack" else eq).append(abs(s_pp.imag - s_py.imag))
        rows.extend(
            {
                "comparison": f"B_{name}",
                "quantity": f"{kind}_{c}_pu",
                "element": b,
                "reference": getattr(s_py, part),
                "other": getattr(s_pp, part),
                "abs_error": abs(getattr(s_pp, part) - getattr(s_py, part)),
            }
            for c, part in (("P", "real"), ("Q", "imag"))
        )
    res.update(
        gen_P_max_abs_error_pu=float(max(ep)),
        gen_Q_max_abs_error_pu=float(max(eq)),
        slack_P_abs_error_pu=float(max(es_p)),
        slack_Q_abs_error_pu=float(max(es_q)),
    )
    # voltages
    ids = np.asarray(pynet.bus_idx)
    rb = net.res_bus.loc[ids]
    Vpp = rb.vm_pu.to_numpy() * np.exp(1j * np.radians(rb.va_degree.to_numpy()))
    off = np.median(np.angle(Vpp) - np.angle(V))
    res["vm_max_abs_error_pu"] = float(np.max(np.abs(np.abs(Vpp) - np.abs(V))))
    res["va_max_gauge_aligned_error_rad"] = float(
        np.max(np.abs(np.angle(Vpp) - off - np.angle(V)))
    )
    for k, b in enumerate(ids):
        rows.append(
            {
                "comparison": f"B_{name}",
                "quantity": "Vm_pu",
                "element": int(b),
                "reference": abs(V[k]),
                "other": abs(Vpp[k]),
                "abs_error": abs(abs(Vpp[k]) - abs(V[k])),
            }
        )
    # native branch flows, mapped to canonical from/to terminals
    lookup = net._from_ppc_lookups["branch"]
    worst = 0.0
    for k, x in enumerate(payload["lines"]):
        f, t = int(x["bus1"]), int(x["bus2"])
        et, e = lookup.loc[k, "element_type"], int(lookup.loc[k, "element"])
        if et == "line":
            r = net.res_line.loc[e]
            ends = {
                int(net.line.at[e, "from_bus"]): (r.p_from_mw, r.q_from_mvar),
                int(net.line.at[e, "to_bus"]): (r.p_to_mw, r.q_to_mvar),
            }
        elif et == "trafo":
            r = net.res_trafo.loc[e]
            ends = {
                int(net.trafo.at[e, "hv_bus"]): (r.p_hv_mw, r.q_hv_mvar),
                int(net.trafo.at[e, "lv_bus"]): (r.p_lv_mw, r.q_lv_mvar),
            }
        else:
            r = net.res_impedance.loc[e]
            ends = {
                int(net.impedance.at[e, "from_bus"]): (r.p_from_mw, r.q_from_mvar),
                int(net.impedance.at[e, "to_bus"]): (r.p_to_mw, r.q_to_mvar),
            }
        ref4 = py_flows[k, 3:7]
        got4 = np.array([*ends[f], *ends[t]])
        err = float(np.max(np.abs(got4 - ref4)))
        worst = max(worst, err)
        rows.append(
            {
                "comparison": f"B_{name}",
                "quantity": "branch_terminal_PQ_MVA_maxabs",
                "element": k,
                "reference": float(np.max(np.abs(ref4))),
                "other": float(np.max(np.abs(got4))),
                "abs_error": err,
                "note": f"{f}-{t} {et}",
            }
        )
    res["native_branch_terminal_flow_max_error_MVA"] = worst
    return res


summary["B_as_supplied"] = run_variant("as_supplied", tapside._rated_canonical_to_ppc)
summary["B_corrected"] = run_variant("corrected", tapside.canonical_to_ppc_hv_from)
summary["B_corrected_reoriented_branches"] = [
    {k: (float(v) if isinstance(v, (float, np.floating)) else v) for k, v in r.items()}
    for r in tapside.REORIENTED[:3]
]
pd.DataFrame(rows).to_csv(out / "pandapower_static_parity_long.csv", index=False)
(out / "static_supp.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
