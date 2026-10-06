# ruff: noqa: E501
"""E8 local stability-equivalent exchange rates (control <-> topology), prereg E8."""

from __future__ import annotations

import json

import _infra as I
import numpy as np
import pandas as pd

import _cdw as C
import _sens as S

PHASE_A = "E08a"
PHASE_B = "E08b"
POINTS = ("BSTAR", "D01")
SURV = tuple(b for b in C.SG_BUSES if b not in C.V4)


def specs(theta):
    sp = [{"kind": "g", "idx": i, "value": theta[0]} for i in C.V4]
    sp += [{"kind": "ka", "idx": i, "value": 1.0} for i in SURV]
    sp += [{"kind": "line", "idx": e, "value": 1.0} for e in range(46)]
    return sp


def tasks_a():
    return [{"phase": PHASE_A, "pid": p} for p in POINTS]


def run_task_a(task):
    theta = C.policy(task["pid"])
    eng = S.Engine(C.V4, theta)
    d = eng.derivatives(specs(theta), "SPR")
    return {"lam": d["lam"], "gap2": d["gap2"], "items": d["items"]}


def derivs():
    out = {}
    for r in I.Store(PHASE_A).all():
        if r.get("ok"):
            out[r["task"]["pid"]] = r
    return out


def tasks_b():
    der = derivs()
    b = der["BSTAR"]
    lines = sorted((it for it in b["items"] if it["kind"] == "line"), key=lambda it: -abs(it["d_total"][0]))[:5]
    top = [it["idx"] for it in lines]
    I.atomic_write_json(I.RESULTS / "CDW_E8_selected_lines.json", {"rule": "top-5 |d alpha/d gamma| total at BSTAR", "lines": top})
    out = []
    for pid in POINTS:
        items = {(it["kind"], it["idx"]): it["d_total"][0] for it in der[pid]["items"]}
        theta = C.policy(pid)
        for ck, ci in [("g", i) for i in C.V4] + [("ka", i) for i in SURV]:
            for e in top:
                dth, dga = items[(ck, ci)], items[("line", e)]
                rate = -dth / dga if dga != 0 else np.nan
                for size, delta in (("small", 0.02), ("large", 0.1)):
                    out.append({"phase": PHASE_B, "pid": pid, "control": ck, "ci": ci, "line": e, "rate": rate,
                                "d_theta": dth, "d_gamma": dga, "size": size, "delta": delta,
                                "g0": theta[0]})
    return out


def run_task_b(task):
    theta = C.policy(task["pid"])
    base = C.eval_portfolio(C.V4, theta, modes=False)["alpha"]
    ck, ci, d = task["control"], task["ci"], task["delta"]
    if ck == "g":
        design_c = {"g": {ci: task["g0"] + d}}
    else:
        design_c = {"ka": {ci: 1.0 + d}}
    kw_c = S.design_kwargs(design_c)
    a_c = C.eval_portfolio(C.V4, theta, modes=False, **kw_c)["alpha"]
    gam = 1.0 + task["rate"] * d
    if not np.isfinite(gam) or gam <= 0:
        return {"alpha0": base, "alpha_control": a_c, "alpha_pair": None, "gamma": gam}
    design_p = {**design_c, "line": {task["line"]: gam}}
    a_p = C.eval_portfolio(C.V4, theta, modes=False, **S.design_kwargs(design_p))["alpha"]
    return {"alpha0": base, "alpha_control": a_c, "alpha_pair": a_p, "gamma": gam}


def aggregate():
    rows = []
    for r in I.Store(PHASE_B).all():
        if not r.get("ok"):
            continue
        t = r["task"]
        dc = r["alpha_control"] - r["alpha0"]
        dp = (r["alpha_pair"] - r["alpha0"]) if r["alpha_pair"] is not None else np.nan
        rows.append({**{k: t[k] for k in ("pid", "control", "ci", "line", "rate", "d_theta", "d_gamma", "size", "delta")},
                     "gamma": r["gamma"], "d_alpha_control": dc, "d_alpha_pair": dp,
                     "compensation_ratio": abs(dp) / abs(dc) if dc != 0 and np.isfinite(dp) else np.nan})
    df = pd.DataFrame(rows)
    df.to_csv(I.RESULTS / "CDW_E8_exchange.csv", index=False)
    g = {}
    for (pid, size), gg in df.groupby(["pid", "size"]):
        g[f"{pid}_{size}_frac_ratio_le_0.2"] = float((gg.compensation_ratio <= 0.2).mean())
        g[f"{pid}_{size}_median_ratio"] = float(gg.compensation_ratio.median())
        g[f"{pid}_{size}_n"] = int(len(gg))
    g["E8_pass_BSTAR_small"] = bool(g.get("BSTAR_small_frac_ratio_le_0.2", 0) >= 0.8)
    I.atomic_write_json(I.RESULTS / "CDW_E8_gates.json", g)
    return g


if __name__ == "__main__":
    print(json.dumps(aggregate(), indent=1))
