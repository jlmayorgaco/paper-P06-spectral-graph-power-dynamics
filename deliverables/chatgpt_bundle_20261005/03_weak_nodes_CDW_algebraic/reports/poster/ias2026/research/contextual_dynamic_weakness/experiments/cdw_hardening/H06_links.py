# ruff: noqa: E501
"""H06 link sensitivity on the new holdout (+ H08 node coordinates): SPR frozen / total derivatives,
small-step validation and finite re-equilibrated truth at gamma in {1.10, 1.25, 1.50} with PF
convergence and bus-voltage recording (prereg H6, H8)."""

from __future__ import annotations

import _hinfra as HI
import numpy as np

import _cdw as C
import _sens as S
import E34_sens as E34

PHASE = "H_H06"
PHASE_N = "H_H08n"
GAMMAS = (1.10, 1.25, 1.50)
CHUNK = 12
TARGETS = {"H4": C.V4, "V9": C.V9}
VGUARD = (0.90, 1.10)


def conditions():
    out = [{"pid": p, "target": t, "source": None, "env": None, "draw": None} for p in HI.HPOL for t in ("H4", "V9")]
    out += [{"pid": "D01", "target": "H4", "source": "FRESH", "env": d["envelope"], "draw": d["draw"]} for d in HI.fresh_draws()]
    return out


def tasks():
    specs = E34.link_specs()
    out = []
    for c in conditions():
        for k in range(0, len(specs), CHUNK):
            out.append({"phase": PHASE, "family": "link", **c, "specs": specs[k: k + CHUNK]})
    seen = set()
    for c in conditions():
        key = (c["pid"], c["env"], c["draw"])
        if key in seen:
            continue
        seen.add(key)
        out.append({"phase": PHASE, "family": "conv", **{**c, "target": "BASE"}, "specs": specs})
    for p in HI.HPOL:
        out.append({"phase": PHASE, "family": "port", "pid": p, "target": "H4", "source": None, "env": None, "draw": None, "specs": specs})
    return out


def finite_rec(eng, sp, value, lam_ref):
    case = S.spr_case(eng, sp, value)
    if case is None:
        return {"ok": False}
    A = C.physical_matrices(case)[0]
    ev = S.transverse_spectrum(case, A)
    j = int(np.argmin(np.abs(ev - lam_ref)))
    vm = np.abs(case.dae.voltages(case.equilibrium.z))
    ic = int(np.argmax(ev.real))
    return {"ok": True, "alpha": float(ev.real.max()), "alpha_hz": float(abs(ev[ic].imag) / (2 * np.pi)),
            "tracked_re": float(ev[j].real), "tracked_hz": float(abs(ev[j].imag) / (2 * np.pi)),
            "vmin": float(vm.min()), "vmax": float(vm.max()),
            "admissible": bool(vm.min() >= VGUARD[0] and vm.max() <= VGUARD[1])}


def run_task(task):
    theta = HI.theta_of(task["pid"])
    draw = HI.draw_of(task["source"], task["env"], task["draw"])
    members = () if task["target"] == "BASE" else TARGETS[task["target"]]
    eng = S.Engine(members, theta, draw=draw)
    if task["family"] == "port":
        return {"items": E34.port_items(eng, task["specs"])}
    d = eng.derivatives(task["specs"], "SPR")
    lam = complex(*d["lam"])
    A00, _ = eng.A_of(eng.w0("SPR"), None, "SPR")
    alpha0 = float(eng.transverse(A00).real.max())
    items = []
    for it in d["items"]:
        rec = dict(it)
        if task["family"] in ("link", "node"):
            small, _ = E34.steps(it, "SPR")
            fp, fm = eng.finite(it, small[0], "SPR", lam), eng.finite(it, small[1], "SPR", lam)
            rec["fd_small"] = (fp["tracked_re"] - fm["tracked_re"]) / (small[0] - small[1]) if (fp and fm) else None
            rec["alpha0"] = alpha0
        if task["family"] == "link":
            for gma in GAMMAS:
                f = finite_rec(eng, it, gma, lam)
                for k, v in f.items():
                    rec[f"g{gma:.2f}_{k}"] = v
        items.append(rec)
    return {"lam": d["lam"], "gap2": d["gap2"], "R0": d["R0"], "match_err": d.get("match_err"), "items": items,
            "alpha0_transverse": alpha0}


# ------------------------------------------------------------------ H08 node specs --
def tasks_node():
    out = []
    for p in HI.HPOL:
        specs = E34.node_specs(C.V4, HI.theta_of(p))
        for k in range(0, len(specs), CHUNK):
            out.append({"phase": PHASE_N, "family": "node", "pid": p, "target": "H4", "source": None, "env": None,
                        "draw": None, "specs": specs[k: k + CHUNK]})
    return out
