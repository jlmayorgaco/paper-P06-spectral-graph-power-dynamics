# ruff: noqa: E501
"""H14 topology with electromechanical mode tracking (prereg H14-H16): core lattice with mode
shapes at the old E7 policies plus HARDENING_H01-H06, nominal topology and every old-admissible
action; H18c: custom-GFL H4 evaluations for the six H18 topology actions at the 8 H18 policies."""

from __future__ import annotations

import pandas as pd

import _cdw as C
import _hinfra as HI
import _infra as I
import E07_topology as E7

PHASE = "H_H14"
PHASE_C = "H_H18c"
POLICIES = ("D01", "D03", "D11", "D14", "H01", "H02", "H03", "H04", "H05", "H06") + tuple(HI.HPOL[:6])
H18_POLICIES = ("H02", "H04", "H05", "H07") + tuple(HI.HPOL[:4])
H18_ACTIONS = ("dbl45", "dbl41", "out2", "out33", "out26", "out20")


def admissible_actions():
    adm = pd.read_csv(I.RESULTS / "CDW_E7_admissibility.csv")
    ok = set(adm[adm.admissible].action)
    return [a for a in E7.actions() if a["action"] in ok]


def tasks():
    acts = [{"action": "NOMINAL", "outage": [], "scale": {}}] + admissible_actions()
    return [{"phase": PHASE, "pid": p, "action": a["action"], "act": a} for p in POLICIES for a in acts]


def run_task(task):
    theta = HI.theta_of(task["pid"])
    act = task["act"]
    net = E7.network_of(act) if (act["outage"] or act["scale"]) else None
    recs = []
    for s in C.subsets(C.V4):
        r = C.eval_portfolio(s, theta, modes=True, network=net)
        r["S"] = C.label(s)
        recs.append(r)
    return {"records": recs}


def tasks_c():
    acts = {a["action"]: a for a in E7.actions()}
    acts["NOMINAL"] = {"action": "NOMINAL", "outage": [], "scale": {}}
    return [{"phase": PHASE_C, "pid": p, "action": a, "act": acts[a]} for p in H18_POLICIES for a in ("NOMINAL",) + H18_ACTIONS]


def run_task_c(task):
    theta = HI.theta_of(task["pid"])
    act = task["act"]
    net = E7.network_of(act) if (act["outage"] or act["scale"]) else None
    r = C.eval_portfolio(C.V4, theta, modes=True, network=net)
    return {"record": r}
