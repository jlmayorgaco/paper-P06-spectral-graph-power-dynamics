# ruff: noqa: E501
"""H01: V9 census (512 portfolios, mode shapes) at HARDENING_H01-H24; H01e: core lattice at P4
for the 40 fresh envelope draws (prereg H1, L7). Same recording as the old E1 run_task."""

from __future__ import annotations

import _hinfra as HI

import _cdw as C

PHASE = "H_H01"
PHASE_E = "H_H01e"
CHUNK = 32


def tasks():
    subs = C.subsets(C.V9)
    return [{"phase": PHASE, "pid": pid, "subsets": [list(s) for s in subs[c: c + CHUNK]]}
            for pid in HI.HPOL for c in range(0, len(subs), CHUNK)]


def run_task(task):
    theta = HI.theta_of(task["pid"])
    recs = []
    for s in task["subsets"]:
        r = C.eval_portfolio(tuple(s), theta, modes=True)
        r["S"] = C.label(s)
        recs.append(r)
    return {"records": recs}


def tasks_e():
    return [{"phase": PHASE_E, "pid": "D01", "source": "FRESH", "env": d["envelope"], "draw": d["draw"]}
            for d in HI.fresh_draws()]


def run_task_e(task):
    theta = HI.theta_of(task["pid"])
    draw = HI.draw_of(task["source"], task["env"], task["draw"])
    recs = []
    for s in C.subsets(C.V4):
        r = C.eval_portfolio(s, theta, modes=True, draw=draw)
        r["S"] = C.label(s)
        recs.append(r)
    return {"records": recs}
