# ruff: noqa: E501
"""H09 equal-budget corridors and H10 size-matched corridor nulls (prereg H9, H10, H11).

Truth is the finite re-equilibrated alpha(H4) under SPR (TX4 solve_case) with the branch
scaling applied; budgets: L1 0.50 (primary), L1 1.00, L2 0.50 (uniform within a group)."""

from __future__ import annotations

import json
import math

import _hinfra as HI

import _cdw as C

PHASE9 = "H_H09"
PHASE10 = "H_H10"
CHUNK10 = 50
BUDGETS = {"L1_0.50": ("L1", 0.50), "L1_1.00": ("L1", 1.00), "L2_0.50": ("L2", 0.50)}


def corridor_defs() -> dict:
    return json.loads((C.RESEARCH / "contextual_dynamic_weakness/results/CDW_E6_corridor_definitions.json").read_text())


def scale_of(group, budget) -> dict:
    kind, b = BUDGETS[budget]
    k = len(group)
    dg = b / k if kind == "L1" else b / math.sqrt(k)
    return {int(e): 1.0 + dg for e in group}


def cond_list(sets=("discovery", "old", "new")):
    out = []
    if "discovery" in sets:
        out += [{"pid": p, "set": "discovery", "source": None, "env": None, "draw": None} for p in ("D01", "D03", "D08", "D11", "D14")]
    if "old" in sets:
        out += [{"pid": f"H{i:02d}", "set": "old", "source": None, "env": None, "draw": None} for i in range(1, 13)]
        out += [{"pid": "D01", "set": "old_draws", "source": "CDW", "env": e, "draw": k} for e in HI.ENVS for k in range(5)]
    if "new" in sets:
        out += [{"pid": p, "set": "new", "source": None, "env": None, "draw": None} for p in HI.HPOL]
        out += [{"pid": "D01", "set": "fresh_draws", "source": "FRESH", "env": d["envelope"], "draw": d["draw"]} for d in HI.fresh_draws()]
    return out


def cond_key(c) -> str:
    return f"{c['pid']}|{c['env'] or '-'}|{c['draw'] if c['draw'] is not None else -1}"


def tasks9():
    cor = corridor_defs()
    return [{"phase": PHASE9, **c, "corridors": cor} for c in cond_list()]


def _alpha(theta, draw, scale):
    net = C.network_with(scale=scale) if scale else None
    r = C.eval_portfolio(C.V4, theta, modes=False, draw=draw, network=net)
    return r


def run_task9(task):
    theta = HI.theta_of(task["pid"])
    draw = HI.draw_of(task["source"], task["env"], task["draw"])
    base = _alpha(theta, draw, None)
    done, items = {}, []
    for name, cut in task["corridors"].items():
        for bud in BUDGETS:
            sc = scale_of(cut, bud)
            key = json.dumps(sorted(sc.items()))
            if key not in done:
                done[key] = _alpha(theta, draw, sc)
            r = done[key]
            items.append({"corridor": name, "size": len(cut), "budget": bud, "alpha": r["alpha"], "status": r["status"],
                          "lam_hz": r.get("lam_hz")})
    return {"alpha0": base["alpha"], "status0": base["status"], "lam_hz0": base.get("lam_hz"), "items": items}


# ------------------------------------------------------------------------ nulls --
def unique_groups() -> list:
    ng = HI.null_groups()
    seen, out = set(), []
    for fam in ("A", "B"):
        for k, lst in ng[fam].items():
            for g in lst:
                t = tuple(sorted(g))
                if t not in seen:
                    seen.add(t)
                    out.append(list(t))
    return out


def tasks10():
    groups = unique_groups()
    out = []
    for c in cond_list(("old", "new")):
        for i in range(0, len(groups), CHUNK10):
            out.append({"phase": PHASE10, **c, "g0": i, "groups": groups[i: i + CHUNK10]})
    return out


def run_task10(task):
    theta = HI.theta_of(task["pid"])
    draw = HI.draw_of(task["source"], task["env"], task["draw"])
    base = _alpha(theta, draw, None)
    res = []
    for g in task["groups"]:
        r = _alpha(theta, draw, scale_of(g, "L1_0.50"))
        res.append([r["alpha"], r["status"], r.get("lam_hz")])
    return {"alpha0": base["alpha"], "status0": base["status"], "lam_hz0": base.get("lam_hz"), "res": res}
