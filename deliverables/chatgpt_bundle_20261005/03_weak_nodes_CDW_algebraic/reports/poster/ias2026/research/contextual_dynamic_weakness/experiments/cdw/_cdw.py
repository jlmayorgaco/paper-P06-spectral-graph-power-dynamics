# ruff: noqa: E501
"""CDW model library: frozen inputs, portfolio evaluation (TX4 machinery), network tools.

Reuses TX4 code READ-ONLY: ibr_cycles (solve_case, devices, classifier, transverse
operator) and experiments/_f7_common.LEAK. Definitions follow docs/CDW_PREREG_V1.md.
"""

from __future__ import annotations

import json
import math
import sys
from dataclasses import replace
from functools import lru_cache
from itertools import combinations

from _infra import PROJECT, RESEARCH

for p in (RESEARCH / "src", RESEARCH / "experiments"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import numpy as np  # noqa: E402

from ibr_cycles.certification.classify import SAFETY, classify_spectrum  # noqa: E402
from ibr_cycles.certification.physical import physical_matrices  # noqa: E402
from ibr_cycles.certification.symmetry import frequency_partner, rotation_generator  # noqa: E402
from ibr_cycles.certification.transverse import transverse_operator  # noqa: E402
from ibr_cycles.models.ieee39_case import InfeasibleReplacement, ReplacementPlan, _controller_payload, solve_case  # noqa: E402
from ibr_cycles.models.ieee39_devices import ConverterParameters  # noqa: E402
from ibr_cycles.models.ieee39_network import load_network, solve_power_flow  # noqa: E402

LEAK = 0.05  # TX4 _f7_common.LEAK (frozen)
V9 = (30, 31, 32, 33, 34, 35, 36, 37, 38)
V4 = (30, 33, 35, 37)
SG_BUSES = (30, 31, 32, 33, 34, 35, 36, 37, 38, 39)
TAU_MAT = 0.01
TAU_RES = 1e-3
MAC_MIN = 0.80
DF_MAX = 0.15
MODE_BAND = (0.1, 2.0)
MODE_RE_MIN = -1.0
STRUCT_ZERO = 1e-6

DISCOVERY = {
    "D01": (0.03625, 1.425, 1.5, 1.0),
    "D02": (1.0, 0.5, 1.5, 1.0),
    "D03": (0.25, 1.425, 1.5, 1.0),
    "D04": (0.0, 1.425, 1.5, 1.0),
    "D05": (0.1, 1.425, 1.5, 1.0),
    "D06": (0.5, 1.425, 1.5, 1.0),
    "D07": (1.0, 1.425, 1.5, 1.0),
    "D08": (0.1, 0.75, 1.5, 1.0),
    "D09": (0.5, 1.75, 1.5, 1.0),
    "D10": (0.9, 2.2, 1.5, 1.0),
    "D11": (0.2, 1.0, 0.8515625, 1.0),
    "D12": (0.6, 1.0, 0.8515625, 1.0),
    "D13": (0.3, 1.0, 2.5, 1.0),
    "D14": (0.2, 1.0, 1.0, 0.5),
    "D15": (0.7, 1.0, 1.0, 1.5),
}
BOUNDARY_THETA = (0.20768, 1.425, 1.5, 1.0)  # TX4 frozen H4 boundary on the P4 line


@lru_cache(maxsize=1)
def holdout() -> dict:
    rows = json.loads((PROJECT / "results/prereg_inputs/holdout_policies.json").read_text())
    return {r["id"]: (r["g"], r["k"], r["t"], r["h"]) for r in rows}


def policy(pid: str) -> tuple:
    if pid in DISCOVERY:
        return DISCOVERY[pid]
    if pid == "BSTAR":
        return BOUNDARY_THETA
    return holdout()[pid]


@lru_cache(maxsize=1)
def cdw_draws() -> list:
    return json.loads((PROJECT / "results/prereg_inputs/cdw_envelope_draws.json").read_text())


MACHINE_GROUPS = {
    "M": ("m",), "XP": ("xd1", "xq1"), "KA": ("ka",), "TE": ("ta",),
    "PSS_K": ("pss_gain",), "PSS_T": ("pss_washout", "pss_wash_lag", "pss_lag"),
}
CONV_GROUPS = {
    "PLL": ("kp_pll", "ki_pll"), "OUTER": ("kp_p", "ki_p", "kp_q", "ki_q"),
    "CURRENT": ("kp_i", "ki_i"), "TAU_P": ("tau_p",), "XF": ("xf",),
}
DEFAULT_CONV = ConverterParameters()


def tx4_draws() -> list:
    """The TX4 PCV05 draws, regenerated with the frozen generator (seed 20260917)."""

    import importlib

    sys.path.insert(0, str(RESEARCH / "experiments/post_cumulant_validation"))
    # PCV05 creates results/PCV/PCV05 on import (out_dir); that directory exists already
    # (frozen) and mkdir(exist_ok) writes nothing.
    mod = importlib.import_module("PCV05_robustness")
    out = []
    for d in mod.make_draws():
        if d["envelope"] in ("EM-f", "EM-u", "EC", "EMC"):
            out.append({
                "envelope": d["envelope"], "draw": d["draw"], "fleet": d["fleet"],
                "unit": {str(b): g for b, g in d["unit"].items()},
                "conv": {str(b): g for b, g in d["conv"].items()},
            })
    return out


def te_native() -> dict:
    return {int(b): float(r["TE"]) for b, r in _controller_payload()["avr_by_bus"].items()}


def te_gmean() -> float:
    v = list(te_native().values())
    return float(math.exp(sum(math.log(x) for x in v) / len(v)))


# ------------------------------------------------------------------ kwargs --
def case_kwargs(theta, *, draw=None, conv_extra=None, mbs_extra=None, network=None, device="gfl"):
    """solve_case keyword arguments for policy theta, an envelope draw and extra
    per-bus perturbations (conv_extra {bus: {field: factor or ('abs', value)}},
    mbs_extra {bus: {field: factor}})."""

    g, k, t, h = theta
    ms = {"ka": k, "ta": t}
    mbs: dict[int, dict] = {}
    if h != 1.0:
        mean = te_gmean()
        for b, te in te_native().items():
            mbs.setdefault(b, {})["ta"] = (mean / te) ** (1.0 - h)
    draw = draw or {}
    for grp, f in (draw.get("fleet") or {}).items():
        for field in MACHINE_GROUPS[grp]:
            ms[field] = ms.get(field, 1.0) * f
    for b, groups in (draw.get("unit") or {}).items():
        for grp, f in groups.items():
            for field in MACHINE_GROUPS[grp]:
                d = mbs.setdefault(int(b), {})
                d[field] = d.get(field, 1.0) * f
    for b, fields in (mbs_extra or {}).items():
        for field, f in fields.items():
            d = mbs.setdefault(int(b), {})
            d[field] = d.get(field, 1.0) * f
    base = ConverterParameters(voltage_control=True, voltage_gain=g, voltage_leak=LEAK)
    convs = {}
    conv_draw = draw.get("conv") or {}
    for b in set(map(int, conv_draw)) | set(map(int, conv_extra or {})):
        kw = {}
        for grp, f in conv_draw.get(str(b), {}).items():
            for field in CONV_GROUPS[grp]:
                kw[field] = getattr(DEFAULT_CONV, field) * f
        for field, spec in (conv_extra or {}).get(b, {}).items():
            if isinstance(spec, tuple) and spec[0] == "abs":
                kw[field] = spec[1]
            else:
                kw[field] = kw.get(field, getattr(base, field)) * spec
        convs[b] = replace(base, **kw)
    out = {"converter": base, "machine_scaling": ms}
    if convs:
        out["converters"] = convs
    if mbs:
        out["machine_bus_scaling"] = mbs
    if network is not None:
        out["network"] = network
    return out


def solve(members, theta, *, device="gfl", **kw):
    plan = ReplacementPlan.of({b: 1.0 for b in members}, device=device)
    return solve_case(plan, **case_kwargs(theta, device=device, **kw))


# -------------------------------------------------------------- evaluation --
def transverse_eval(case, *, modes=True) -> dict:
    """TX4 direct path (FC01/PCV B0): status, alpha; plus rightmost mode and band mode shapes."""

    a, d, jac = physical_matrices(case)
    r_x, _ = rotation_generator(case.dae, case.equilibrium.z)
    w = frequency_partner(case.dae).w
    tr = transverse_operator(a, r_x, w)
    ev_t = np.linalg.eigvals(tr.a_perp)
    status = classify_spectrum(tr.a_perp, tr.z.T @ d @ tr.z, SAFETY).status
    i_c = int(np.argmax(ev_t.real))
    lam_c = complex(ev_t[i_c])
    if lam_c.imag < 0:
        lam_c = lam_c.conjugate()
    out = {
        "status": status,
        "alpha": float(ev_t.real.max()),
        "lam_re": lam_c.real,
        "lam_hz": abs(lam_c.imag) / (2 * np.pi),
        "rhp": int((ev_t.real > 0).sum()),
        "n_x": int(a.shape[0]),
    }
    srt = np.sort(ev_t.real)[::-1]
    # second rightmost distinct mode (skip the conjugate partner)
    others = ev_t[np.abs(ev_t - lam_c) > 1e-9]
    others = others[np.abs(others - lam_c.conjugate()) > 1e-9]
    out["gap2"] = float(lam_c.real - others.real.max()) if others.size else float("inf")
    del srt
    if modes:
        vals, vecs = np.linalg.eig(a)
        phi_all = -np.linalg.solve(jac.gz, jac.gx @ vecs)  # bus-voltage mode shapes (78 x n)
        f = np.abs(vals.imag) / (2 * np.pi)
        keep = []
        for j in range(vals.size):
            lam = vals[j]
            if abs(lam) < STRUCT_ZERO or lam.imag < -1e-12:
                continue
            in_band = MODE_BAND[0] <= f[j] <= MODE_BAND[1] and lam.real >= MODE_RE_MIN
            if in_band:
                keep.append(j)
        j_c = int(np.argmin(np.abs(vals - lam_c)))
        if j_c not in keep:
            keep.append(j_c)
        modes_out = []
        for j in keep:
            ph = phi_all[:, j]
            n = np.linalg.norm(ph)
            ph = ph / n if n > 0 else ph
            # fix the arbitrary complex phase: largest entry real positive
            k = int(np.argmax(np.abs(ph)))
            ph = ph * np.exp(-1j * np.angle(ph[k]))
            modes_out.append({
                "re": float(vals[j].real), "hz": float(abs(vals[j].imag) / (2 * np.pi)),
                "crit": j == j_c,
                "phi": np.round(np.concatenate([ph.real, ph.imag]), 5).astype(float).tolist(),
            })
        out["modes"] = modes_out
    return out


def eval_portfolio(members, theta, *, modes=True, device="gfl", **kw) -> dict:
    try:
        case = solve(tuple(members), theta, device=device, **kw)
    except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError) as e:
        return {"status": "INFEASIBLE", "alpha": float("nan"), "error": str(e)[:200]}
    rec = transverse_eval(case, modes=modes)
    rec["gz_cond"] = float(case.gz_condition)
    return rec


def mac(phi1, phi2) -> float:
    a = np.asarray(phi1[: len(phi1) // 2]) + 1j * np.asarray(phi1[len(phi1) // 2:])
    b = np.asarray(phi2[: len(phi2) // 2]) + 1j * np.asarray(phi2[len(phi2) // 2:])
    num = abs(np.vdot(a, b)) ** 2
    den = (np.vdot(a, a).real * np.vdot(b, b).real) or 1e-300
    return float(num / den)


def label(members) -> str:
    return "+".join(map(str, sorted(members))) or "BASE"


def subsets(cands):
    return [tuple(sorted(s)) for r in range(len(cands) + 1) for s in combinations(cands, r)]


# ---------------------------------------------------------------- network --
@lru_cache(maxsize=1)
def network_payload() -> dict:
    return json.loads((RESEARCH / "configs/ias2026/ieee39_network.json").read_text(encoding="utf-8"))


def build_ybus(scale: dict | None = None, outage: set | None = None) -> np.ndarray:
    """Exact Ybus (TX4 _build_ybus) with whole-two-port branch scaling rho_e.

    ``scale`` {branch index: rho}; ``outage`` removes branches (rho = 0)."""

    payload = network_payload()
    order = [int(b["idx"]) for b in payload["buses"]]
    index = {bus: i for i, bus in enumerate(order)}
    y = np.zeros((len(order), len(order)), dtype=np.complex128)
    for e, line in enumerate(payload["lines"]):
        if float(line["u"]) == 0.0:
            continue
        rho = (scale or {}).get(e, 1.0)
        if outage and e in outage:
            rho = 0.0
        if rho == 0.0:
            continue
        f, t = index[int(line["bus1"])], index[int(line["bus2"])]
        series = 1.0 / complex(line["r"], line["x"])
        charging = complex(line["g"], line["b"]) / 2.0
        m = float(line["tap"]) * np.exp(1j * float(line["phi"]))
        m2 = float(abs(m) ** 2)
        y[f, f] += rho * (series + charging) / m2
        y[t, t] += rho * (series + charging)
        y[f, t] += rho * (-series / np.conj(m))
        y[t, f] += rho * (-series / m)
    for sh in payload["shunts"]:
        y[index[int(sh["bus"])], index[int(sh["bus"])]] += complex(sh["g"], sh["b"])
    return y


def network_with(scale=None, outage=None, loads_scale: dict | None = None, vset: dict | None = None):
    net = load_network()
    kw = {}
    if scale or outage:
        kw["ybus"] = build_ybus(scale, outage)
    if loads_scale:
        kw["loads"] = {b: v * loads_scale.get(b, 1.0) for b, v in net.loads.items()}
    if vset:
        pv = {b: dict(s) for b, s in net.pv.items()}
        for b, dv in vset.items():
            if b in pv:
                pv[b]["v"] = pv[b]["v"] + dv
        kw["pv"] = pv
        if net.slack_bus in vset:
            kw["slack_voltage"] = net.slack_voltage + vset[net.slack_bus]
    return replace(net, **kw) if kw else net


def branches() -> list[dict]:
    payload = network_payload()
    out = []
    for e, line in enumerate(payload["lines"]):
        out.append({"e": e, "f": int(line["bus1"]), "t": int(line["bus2"]), "r": float(line["r"]),
                    "x": float(line["x"]), "b": float(line["b"]), "tap": float(line["tap"]),
                    "trans": int(float(line["trans"]))})
    return out


def laplacian_B(scale=None, outage=None) -> tuple[np.ndarray, list]:
    """Exact coupling Laplacian L_B (weights b^s/t), see theory C8."""

    payload = network_payload()
    order = [int(b["idx"]) for b in payload["buses"]]
    index = {bus: i for i, bus in enumerate(order)}
    n = len(order)
    L = np.zeros((n, n))
    for e, br in enumerate(branches()):
        rho = (scale or {}).get(e, 1.0)
        if outage and e in outage:
            rho = 0.0
        y = 1.0 / complex(br["r"], br["x"])
        w = rho * (-y.imag) / br["tap"]
        i, j = index[br["f"]], index[br["t"]]
        L[i, i] += w
        L[j, j] += w
        L[i, j] -= w
        L[j, i] -= w
    return L, order


def connected(outage: set) -> bool:
    import networkx as nx

    G = nx.MultiGraph()
    G.add_nodes_from(int(b["idx"]) for b in network_payload()["buses"])
    for br in branches():
        if br["e"] not in outage:
            G.add_edge(br["f"], br["t"])
    return nx.is_connected(G)


def pf(net):
    return solve_power_flow(net)
