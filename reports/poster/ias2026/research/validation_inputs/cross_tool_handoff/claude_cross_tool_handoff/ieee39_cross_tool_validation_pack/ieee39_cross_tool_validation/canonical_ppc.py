
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

BASE_MVA = 100.0

# PYPOWER/MATPOWER indices (minimal subset)
BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, VA, BASE_KV, ZONE, VMAX, VMIN = range(13)
GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN = range(10)
F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, TAP, SHIFT, BR_STATUS, ANGMIN, ANGMAX = range(13)

def load_payload(repo: Path) -> dict:
    return json.loads((repo / "configs/ias2026/ieee39_network.json").read_text())

def canonical_to_ppc(repo: Path) -> dict:
    p = load_payload(repo)
    buses = sorted(p["buses"], key=lambda x: int(x["idx"]))
    bus_ids = [int(x["idx"]) for x in buses]
    pos = {b:i for i,b in enumerate(bus_ids)}
    pv = {int(x["bus"]):x for x in p["pv"]}
    slack = p["slack"][0]
    slack_bus = int(slack["bus"])
    loads = {}
    for x in p["loads"]:
        b = int(x["bus"])
        loads[b] = loads.get(b, 0j) + complex(float(x["p0"]), float(x["q0"]))
    shunts = {}
    for x in p["shunts"]:
        b = int(x["bus"])
        shunts[b] = shunts.get(b, 0j) + complex(float(x["g"]), float(x["b"]))

    bus = np.zeros((len(buses), 13), float)
    for row, x in enumerate(buses):
        b = int(x["idx"])
        typ = 3 if b == slack_bus else (2 if b in pv else 1)
        ld = loads.get(b, 0j)
        sh = shunts.get(b, 0j)
        bus[row] = [
            b, typ,
            ld.real * BASE_MVA, ld.imag * BASE_MVA,
            sh.real * BASE_MVA, sh.imag * BASE_MVA,
            float(x.get("area", 1)), float(x["v0"]),
            np.degrees(float(x["a0"])), float(x["Vn"]),
            1, 1.2, 0.8
        ]

    gens = []
    for b in sorted(pv):
        x = pv[b]
        gens.append([
            b, x["p"]*BASE_MVA, x["qmax"]*0.0, # QG is an initial guess only
            x["qmax"]*BASE_MVA, x["qmin"]*BASE_MVA,
            x["v"], x["sn"], 1,
            x["pmax"]*BASE_MVA, x["pmin"]*BASE_MVA,
        ] + [0.0]*11)
    gens.append([
        slack_bus, float(slack["p0"])*BASE_MVA, float(slack["q0"])*BASE_MVA,
        9999.0, -9999.0, float(slack["v0"]), float(slack["Sn"]), 1,
        float(slack["pmax"])*BASE_MVA, float(slack["pmin"])*BASE_MVA,
    ] + [0.0]*11)
    gen = np.asarray(gens, float)

    branches = []
    for x in p["lines"]:
        branches.append([
            int(x["bus1"]), int(x["bus2"]),
            float(x["r"]), float(x["x"]), float(x["b"]),
            0.0, 0.0, 0.0,
            float(x["tap"]), np.degrees(float(x["phi"])),
            float(x["u"]), -360.0, 360.0,
        ] + [0.0]*8)  # room for PF/QF/PT/QT etc.
    branch = np.asarray(branches, float)

    return {"version":"2", "baseMVA":BASE_MVA, "bus":bus, "gen":gen, "branch":branch}

def branch_flows_from_voltages(repo: Path, V: np.ndarray) -> np.ndarray:
    p = load_payload(repo)
    buses = sorted(p["buses"], key=lambda x: int(x["idx"]))
    pos = {int(x["idx"]): i for i,x in enumerate(buses)}
    rows = []
    for k, x in enumerate(p["lines"]):
        i, j = pos[int(x["bus1"])], pos[int(x["bus2"])]
        y = 1.0 / complex(float(x["r"]), float(x["x"]))
        yc = complex(float(x["g"]), float(x["b"])) / 2.0
        tap = float(x["tap"]) * np.exp(1j*float(x["phi"]))
        Yff = (y + yc)/(abs(tap)**2)
        Yft = -y/np.conj(tap)
        Ytf = -y/tap
        Ytt = y + yc
        If = Yff*V[i] + Yft*V[j]
        It = Ytf*V[i] + Ytt*V[j]
        Sf = BASE_MVA * V[i] * np.conj(If)
        St = BASE_MVA * V[j] * np.conj(It)
        rows.append((k, int(x["bus1"]), int(x["bus2"]), Sf.real, Sf.imag, St.real, St.imag))
    return np.asarray(rows, float)
