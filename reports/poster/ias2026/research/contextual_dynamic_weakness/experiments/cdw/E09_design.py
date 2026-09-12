# ruff: noqa: E501
"""E9 plan-level multi-constraint design vs single-boundary tuning (prereg E9, GOLD-D)."""

from __future__ import annotations

import json
import time

import _infra as I
import numpy as np
import pandas as pd

import _analysis as AN
import _cdw as C
import _sens as S

PHASE = "E09"
TARGETS = {"T1": ("H4", "D01"), "T2": ("V9", "D01"), "T3": ("H4", "D11"), "T4": ("V9", "D03")}
SETS = {"H4": C.V4, "V9": C.V9}
FAMILIES = ("control", "topology", "joint")
METHODS = ("single", "plan")
EPS = 0.02
ITERS = 6
BOUNDS = {"g": (0.0, 1.0), "ka": (0.7, 1.5), "line": (1.0, 2.0)}
RANGE = {"g": 1.0, "ka": 0.8, "line": 1.0}
TRUST = 0.25  # |d| <= TRUST * range per step (implementation detail, fixed before E9 ran)
MAX_ACTIVE = 25


def params(T, family):
    ctrl = [("g", i) for i in T] + [("ka", i) for i in C.SG_BUSES if i not in T]
    topo = [("line", e) for e in range(46)]
    return {"control": ctrl, "topology": topo, "joint": ctrl + topo}[family]


def tasks():
    return [{"phase": PHASE, "target": tid, "family": f, "method": m} for tid in TARGETS for f in FAMILIES for m in METHODS]


def design_of(vals):
    d = {"g": {}, "ka": {}, "line": {}}
    for (k, i), v in vals.items():
        d[k][i] = v
    return d


def evaluate_all(T, theta, vals):
    kw = S.design_kwargs(design_of(vals))
    out = {}
    for s in C.subsets(T):
        r = C.eval_portfolio(s, theta, modes=False, **kw)
        out[s] = (r["alpha"], r["status"])
    return out


def gradient(Ssub, theta, vals, plist):
    rel = []
    for k, i in plist:
        if k == "g" and i not in Ssub:
            continue
        if k == "ka" and i in Ssub:
            continue
        rel.append((k, i))
    eng = S.Engine(Ssub, theta, design=design_of(vals))
    specs = [{"kind": k, "idx": i, "value": (vals[(k, i)] if k == "g" else 1.0)} for k, i in rel]
    d = eng.derivatives(specs, "SPR") if specs else {"items": []}
    g = np.zeros(len(plist))
    pos = {p: j for j, p in enumerate(plist)}
    for sp, it in zip(rel, d["items"], strict=True):
        v = it["d_total"][0]
        if sp[0] in ("ka", "line"):
            v = v / vals[sp]  # relative -> absolute coordinate
        g[pos[sp]] = v
    return g


def qp_step(active, alphas, grads, vals, plist):
    import cvxpy as cp

    n = len(plist)
    a = np.array([vals[p] for p in plist])
    lo = np.array([BOUNDS[k][0] for k, _ in plist]) - a
    hi = np.array([BOUNDS[k][1] for k, _ in plist]) - a
    rng = np.array([RANGE[k] for k, _ in plist])
    tr = TRUST * rng
    lo, hi = np.maximum(lo, -tr), np.minimum(hi, tr)
    d = cp.Variable(n)
    cons = [d >= lo, d <= hi]
    G = np.array([grads[s] for s in active])
    al = np.array([alphas[s] for s in active])
    cons.append(al + G @ d <= -EPS)
    prob = cp.Problem(cp.Minimize(0.5 * cp.sum_squares(cp.multiply(d, 1 / rng))), cons)
    try:
        prob.solve(solver=cp.CLARABEL)
    except Exception:  # noqa: BLE001
        prob = None
    if prob is not None and prob.status in ("optimal", "optimal_inaccurate"):
        return np.asarray(d.value).ravel(), "feasible", None
    # infeasible: Gordan witness on the active gradients, then the least-violation step
    gv, lam = AN.gordan(G) if len(active) > 1 else (float(np.sum(G**2)), np.ones(1))
    t = cp.Variable()
    d2 = cp.Variable(n)
    prob2 = cp.Problem(cp.Minimize(t + 1e-3 * cp.sum_squares(cp.multiply(d2, 1 / rng))),
                       [d2 >= lo, d2 <= hi, al + G @ d2 <= t])
    prob2.solve(solver=cp.CLARABEL)
    return np.asarray(d2.value).ravel(), "infeasible", {"gordan_value": gv, "lambda": lam.tolist()}


def run_task(task):
    tname, pid = TARGETS[task["target"]]
    T = SETS[tname]
    theta = C.policy(pid)
    plist = params(T, task["family"])
    vals = {p: (theta[0] if p[0] == "g" else 1.0) for p in plist}
    hist = []
    ev = evaluate_all(T, theta, vals)
    phi0 = max(a for a, _ in ev.values())
    rec = {"phi_nominal": phi0, "alpha_T_nominal": ev[T][0], "n_unstable_nominal": sum(st == "UNSTABLE" for _, st in ev.values())}
    if phi0 < 0:
        rec["outcome"] = "ALREADY_SAFE"
        return rec
    active = set()
    t0 = time.time()
    for it in range(ITERS):
        alphas = {s: a for s, (a, _) in ev.items()}
        phi = max(alphas.values())
        aT = alphas[T]
        done = (aT <= -EPS / 2) if task["method"] == "single" else (phi <= -EPS / 2)
        hist.append({"iter": it, "phi": phi, "alpha_T": aT, "n_unstable": sum(st == "UNSTABLE" for _, st in ev.values()),
                     "n_active": len(active)})
        if done:
            break
        if task["method"] == "single":
            active = {T}
        else:
            cand = sorted((s for s in alphas if alphas[s] >= max(phi - 0.05, -0.05)), key=lambda s: -alphas[s])[:MAX_ACTIVE]
            active |= set(cand)
            active = set(sorted(active, key=lambda s: -alphas[s])[:MAX_ACTIVE])
        hist[-1]["n_active"] = len(active)
        grads = {s: gradient(s, theta, vals, plist) for s in active}
        d, status, conflict = qp_step(sorted(active), alphas, grads, vals, plist)
        hist[-1].update(qp=status, conflict=conflict, step_norm=float(np.abs(d).max()))
        for j, p in enumerate(plist):
            lo, hi = BOUNDS[p[0]]
            vals[p] = float(min(max(vals[p] + d[j], lo), hi))
        ev = evaluate_all(T, theta, vals)
    alphas = {s: a for s, (a, _) in ev.items()}
    statuses = {s: st for s, (_, st) in ev.items()}
    phi = max(alphas.values())
    unstable = [C.label(s) for s, st in statuses.items() if st == "UNSTABLE"]
    rec.update(outcome="DESIGNED", history=hist, phi_final=phi, alpha_T_final=alphas[T], status_T_final=statuses[T],
               unstable_subsets_final=unstable, n_unresolved_final=sum(st == "BOUNDARY_OR_UNRESOLVED" for st in statuses.values()),
               final_values={f"{k}{i}": v for (k, i), v in vals.items() if abs(v - (theta[0] if k == "g" else 1.0)) > 1e-9},
               wall_design_s=time.time() - t0)
    return rec


def aggregate():
    rows = []
    for r in I.Store(PHASE).all():
        t = r["task"]
        if not r.get("ok"):
            rows.append({**{k: t[k] for k in ("target", "family", "method")}, "outcome": "ERROR", "error": r.get("error")})
            continue
        row = {**{k: t[k] for k in ("target", "family", "method")}, **{k: v for k, v in r.items() if k not in ("task", "ok", "history", "final_values", "unstable_subsets_final")}}
        row["final_values"] = json.dumps(r.get("final_values", {}))
        row["unstable_subsets_final"] = "|".join(r.get("unstable_subsets_final", []))
        row["n_iters"] = len(r.get("history", []))
        row["any_infeasible_qp"] = any(h.get("qp") == "infeasible" for h in r.get("history", []))
        row["min_gordan"] = min((h["conflict"]["gordan_value"] for h in r.get("history", []) if h.get("conflict")), default=np.nan)
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(I.RESULTS / "CDW_E9_design.csv", index=False)
    gold, cases = False, []
    for (tg, fam), g in df[df.outcome == "DESIGNED"].groupby(["target", "family"]):
        s = g[g.method == "single"]
        p = g[g.method == "plan"]
        if len(s) and len(p):
            s, p = s.iloc[0], p.iloc[0]
            single_final_only = (s.status_T_final == "STABLE") and (s.phi_final >= 0)
            plan_safe = (p.phi_final < 0) and (p.n_unresolved_final == 0) and (not p.unstable_subsets_final)
            cases.append({"target": tg, "family": fam, "single_T_stable": s.status_T_final == "STABLE", "single_phi": s.phi_final,
                          "single_leaves_unsafe_subset": bool(single_final_only), "plan_phi": p.phi_final, "plan_safe": bool(plan_safe)})
            gold |= bool(single_final_only and plan_safe)
    gate = {"cases": cases, "GOLD_D_pass": gold,
            "single_boundary_suffices": [c for c in cases if c["single_T_stable"] and c["single_phi"] < 0],
            "already_safe": df[df.outcome == "ALREADY_SAFE"][["target", "family", "method"]].to_dict("records")}
    I.atomic_write_json(I.RESULTS / "CDW_E9_gates.json", gate)
    return gate


if __name__ == "__main__":
    print(json.dumps(aggregate(), indent=1, default=str))
