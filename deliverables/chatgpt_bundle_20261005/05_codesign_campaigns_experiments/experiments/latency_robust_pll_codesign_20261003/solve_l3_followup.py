"""One further frozen-trust multimode LP from the accepted event-active design.

The event Jacobian is reused from the preceding measured probes. This is only a
predictor; the full exact spectrum and all frozen events decide acceptance.
"""

from __future__ import annotations

import csv
import json
import math
import sys
import tomllib
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


HERE = Path(__file__).resolve().parent
BASE = sys.argv[1] if len(sys.argv)>1 else "N_step1_multimode_lp_event_active_lp"
NAME = BASE + "_followup_lp"


def rows(path: Path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def main() -> None:
    original = "N_step1_multimode_lp"
    d = tomllib.loads((HERE / "designs" / f"{BASE}.toml").read_text())
    events = tomllib.loads((HERE / "event_validations" / BASE / "Q0_RESULT.toml").read_text())
    assert events["all_five_events_pass"]
    slack = min(float(e["min_SG_actuator_fraction_slack"]) for e in events["events"])
    probes = rows(HERE / "L3_FROZEN_EVENT_GRADIENT_PROBES.csv")
    measured = {r["candidate_id"]: r for r in rows(HERE / "L3_EVENT_ACTUATOR_SENSITIVITY.csv")}
    old_slack = float(measured[original]["actuator_slack"])
    E = np.zeros(20)
    allowed = []
    for p in probes[1:]:
        j = int(p["bus"]) - 30 + (0 if p["gain_kind"] == "Kp" else 10)
        E[j] = (float(measured[p["candidate_id"]]["actuator_slack"]) - old_slack) / float(p["delta_log_gain"])
        allowed.append(j)
    modes = rows(HERE / "evaluations" / BASE / "ROOTS.csv")[:3]
    grads = rows(HERE / "evaluations" / BASE / "GRADIENTS.csv")
    tau = np.array([float(r["local_crossing_ms"]) for r in modes])
    G = np.zeros((len(modes), 20))
    for g in grads:
        m = int(g["mode_rank"]) - 1
        if m >= len(modes):
            continue
        b = int(g["bus"]) - 30
        G[m,b] = float(g["d_margin_ms_d_logKp"])
        G[m,10+b] = float(g["d_margin_ms_d_logKi"])
    base = np.r_[d["Kp"],d["Ki"]]
    lim = json.loads((HERE / "EXPERIMENT_MANIFEST.json").read_text())["gain_bounds"]
    low = np.r_[np.full(10,lim["Kp_min"]),np.full(10,lim["Ki_min"])]
    high = np.r_[np.full(10,lim["Kp_max"]),np.full(10,lim["Ki_max"])]
    lb = np.maximum(-0.01,np.log(low/base))
    ub = np.minimum(0.01,np.log(high/base))
    for j in range(20):
        if j not in allowed:
            lb[j]=ub[j]=0
    A=[];b=[]
    for m in range(len(modes)):
        a=np.zeros(41);a[:20]=-G[m];a[-1]=1
        A.append(a);b.append(tau[m])
    a=np.zeros(41);a[:20]=-E
    A.append(a);b.append(slack-0.002-3e-6)
    for j in range(20):
        a=np.zeros(41);a[j]=1;a[20+j]=-1
        c=np.zeros(41);c[j]=-1;c[20+j]=-1
        A.extend([a,c]);b.extend([0,0])
    a=np.zeros(41);a[20:40]=1
    A.append(a);b.append(0.05)
    result=linprog(np.r_[np.zeros(40),-1.0],A_ub=np.asarray(A),b_ub=np.asarray(b),
                   bounds=[(float(lb[j]),float(ub[j])) for j in range(20)]+[(0,0.05)]*20+[(None,None)],
                   method="highs")
    assert result.success,result.message
    du=result.x[:20]
    new=base*np.exp(du)
    (HERE/"designs"/f"{NAME}.toml").write_text(
        "rho = ["+", ".join(repr(float(x)) for x in d["rho"])+"]\n"+
        "Kp = ["+", ".join(repr(float(x)) for x in new[:10])+"]\n"+
        "Ki = ["+", ".join(repr(float(x)) for x in new[10:])+"]\n"+
        f'parent = "{BASE}"\npredictor = "reuse_measured_event_jacobian_multimode_lp"\n',
        encoding="utf-8")
    with (HERE/f"L3_FOLLOWUP_PREDICTION_{BASE}.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=["candidate_id","parent","base_tau_ms","predicted_tau_ms",
            "base_slack","predicted_slack","guard_above_limit","log_step_l1","status"])
        w.writeheader();w.writerow(dict(candidate_id=NAME,parent=BASE,base_tau_ms=float(min(tau)),
            predicted_tau_ms=float(min(tau+G@du)),base_slack=slack,
            predicted_slack=float(slack+E@du),guard_above_limit=3e-6,
            log_step_l1=float(np.sum(np.abs(du))),
            status="REUSED_EVENT_JACOBIAN_PREDICTOR_REQUIRES_EXACT_CORRECTOR"))
    print(NAME,"predicted",float(min(tau+G@du)),"event slack",float(slack+E@du))


if __name__ == "__main__":
    main()
