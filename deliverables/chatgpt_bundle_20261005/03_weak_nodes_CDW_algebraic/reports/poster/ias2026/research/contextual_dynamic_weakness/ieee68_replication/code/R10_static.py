# ruff: noqa: E501
"""R10 static baselines on IEEE-68 (base network and power flow; policy-independent), mirroring
experiments/cdw/E34_sens.static_link_baselines: |P_e|, |S_e|, |z_e|, electrical distance (sub-transient short-circuit
impedance), effective resistance and Fiedler edge score (coupling Laplacian), weighted edge betweenness, endpoint dV/dQ,
Delta gSCR of the target V68 at gamma = 1.5 (normalized by the converter rating |S_gen|/0.8, the 68-bus rule).
Writes results/CDW68_R10_static.csv.
"""

from __future__ import annotations

import _r68 as R  # noqa: I001

import networkx as nx
import numpy as np
import pandas as pd

from ibr_cycles.diagnosis.baselines import short_circuit_ybus
from ibr_cycles.models.ieee39_network import solve_power_flow


def laplacian_B(scale=None):
    net = R.base_network()
    order = list(net.bus_idx)
    idx = {b: i for i, b in enumerate(order)}
    L = np.zeros((len(order), len(order)))
    for b in R.branches():
        rho = (scale or {}).get(b["e"], 1.0)
        tap = b["tap"] if b["tap"] != 0.0 else 1.0
        w = rho * (-(1.0 / complex(b["r"], b["x"])).imag) / tap
        i, j = idx[b["f"]], idx[b["t"]]
        L[i, i] += w
        L[j, j] += w
        L[i, j] -= w
        L[j, i] -= w
    return L, order


def gscr(net, buses, flow, ratings):
    y = short_circuit_ybus(net, exclude=tuple(buses))
    keep = [net.position(b) for b in buses]
    drop = [p for p in range(net.n_bus) if p not in keep]
    red = y[np.ix_(keep, keep)] - y[np.ix_(keep, drop)] @ np.linalg.solve(y[np.ix_(drop, drop)], y[np.ix_(drop, keep)])
    sus = np.abs(np.imag(red))
    volt = np.array([abs(flow.at(b)) ** 2 for b in buses])
    rat = np.array([ratings[b] for b in buses])
    return float(np.min(np.abs(np.linalg.eigvals(np.diag(volt / rat) @ sus))))


def qv_sensitivity(net, flow):
    order = list(net.bus_idx)
    n = len(order)
    V = flow.voltages
    slack = order.index(net.slack_bus)
    pvset = {order.index(b) for b in net.pv}
    ang = [i for i in range(n) if i != slack]
    pq = [i for i in range(n) if i != slack and i not in pvset]
    Y = net.ybus

    def mism(th, vm):
        v = vm * np.exp(1j * th)
        s = v * np.conj(Y @ v)
        return s.real, s.imag

    th0, vm0 = np.angle(V), np.abs(V)
    h = 1e-7
    cols = [("t", i) for i in ang] + [("v", i) for i in pq]
    J = np.zeros((len(cols), len(cols)))
    for c, (kind, i) in enumerate(cols):
        thp, vmp, thm, vmm = th0.copy(), vm0.copy(), th0.copy(), vm0.copy()
        if kind == "t":
            thp[i] += h
            thm[i] -= h
        else:
            vmp[i] += h
            vmm[i] -= h
        pp, qp = mism(thp, vmp)
        pm, qm = mism(thm, vmm)
        J[:, c] = np.concatenate([(pp - pm)[ang], (qp - qm)[pq]]) / (2 * h)
    inv = np.linalg.inv(J)
    nq = len(pq)
    sv = np.diag(inv[len(ang):, len(ang):]) if nq else []
    return {order[i]: float(abs(s)) for i, s in zip(pq, sv, strict=True)}


def main():
    net = R.base_network()
    flow = solve_power_flow(net)
    V = flow.voltages
    order = list(net.bus_idx)
    idx = {b: i for i, b in enumerate(order)}
    Z = np.linalg.inv(short_circuit_ybus(net))
    L, _ = laplacian_B()
    Lp = np.linalg.pinv(L)
    _, U = np.linalg.eigh(L)
    u2 = U[:, 1]
    G = nx.Graph()
    for b in R.branches():
        if G.has_edge(b["f"], b["t"]):
            G[b["f"]][b["t"]]["weight"] = min(G[b["f"]][b["t"]]["weight"], abs(b["x"]))
        else:
            G.add_edge(b["f"], b["t"], weight=abs(b["x"]))
    eb = nx.edge_betweenness_centrality(G, weight="weight")
    dvdq = qv_sensitivity(net, flow)
    ratings = {b: abs(flow.injection(b, net.ybus) + net.loads.get(b, 0j)) / net.converter_loading for b in R.V68}
    g0 = gscr(net, R.V68, flow, ratings)
    rows = []
    for b in R.branches():
        f, t = idx[b["f"]], idx[b["t"]]
        y = 1 / complex(b["r"], b["x"])
        m = b["tap"] if b["tap"] != 0.0 else 1.0
        bsh = float(R.network_payload()["lines"][b["e"]]["b"])
        i_f = ((y + 1j * bsh / 2) / m**2) * V[f] + (-y / m) * V[t]
        s_f = V[f] * np.conj(i_f)
        ev = np.zeros(len(order))
        ev[f], ev[t] = 1, -1
        net2 = R.network_with(scale={b["e"]: 1.5})
        fl2 = solve_power_flow(net2)
        rows.append({"e": b["e"], "f": b["f"], "t": b["t"], "transformer": m != 1.0,
                     "S1_absP": abs(s_f.real), "S2_absS": abs(s_f), "S3_absz": abs(complex(b["r"], b["x"])),
                     "S4_elecdist": abs(Z[f, f] + Z[t, t] - Z[f, t] - Z[t, f]), "S5_reff": float(ev @ Lp @ ev),
                     "S6_fiedler": float((u2[f] - u2[t]) ** 2), "S7_betweenness": eb.get((b["f"], b["t"]), eb.get((b["t"], b["f"]), np.nan)),
                     "S8_dvdq": max(dvdq.get(b["f"], 0.0), dvdq.get(b["t"], 0.0)),
                     "S9_dgscr": gscr(net2, R.V68, fl2, ratings) - g0})
    df = pd.DataFrame(rows)
    df.to_csv(R.RESULTS / "CDW68_R10_static.csv", index=False)
    print(df.describe().T[["mean", "min", "max"]])


if __name__ == "__main__":
    main()
