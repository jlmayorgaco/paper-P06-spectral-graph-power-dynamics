# ruff: noqa: E501
"""R14 (POST HOC, exploratory): nonlinear phasor TDS of EM-tracked nested pairs on IEEE-68, Model A REAL.

The preregistered R14 selection (3 strongest level-D reversals, 2 near-threshold, 2 level-D negative controls) yields
no case: no level-D or level-T reversal and no level-D-eligible nested pair exists. This exploratory check takes the
two EM-tracked (level-T conditions) nested pairs with the largest SAME-SIGN marginals at the R15 reference policy and
asks whether the time response reproduces the linear sign of the intervention effect in both contexts.

Integrator: the validated tool's scheme (scipy BDF on the reduced state, network solved at every right-hand side by
ibr_cycles.nonlinear.tds._Network). The IEEE-39-specific observer/guards of tds.simulate are not used.
Disturbance (preregistered): 0.5 pu shunt reactor at bus 3 for 10 s, then free response to t = 30 s.
Decay estimate: least-squares slope of log|q(t)|, q = v^H (x - x*), v the left eigenvector of the tracked EM mode,
over 11-29 s. Output: results/CDW68_R14_tds_posthoc.json
"""

from __future__ import annotations

import _r68 as R  # noqa: I001

import json
from dataclasses import replace

import numpy as np
import scipy.linalg as sla
from scipy.integrate import BDF

import _rev68 as V
from ibr_cycles.dynamics.linearize import central_difference_jacobians
from ibr_cycles.nonlinear.tds import _Network


def em_top_left(case):
    A = R.matrices(case)[0]
    vals, vl, vr = sla.eig(A, left=True, right=True)
    f = np.abs(vals.imag) / (2 * np.pi)
    band = [j for j in range(vals.size) if R.MODE_BAND[0] <= f[j] <= R.MODE_BAND[1] and vals[j].imag > 0]
    j = max(band, key=lambda i: vals[i].real)
    return complex(vals[j]), vl[:, j], A


def simulate(case, t_pulse=10.0, t_end=30.0, max_step=0.02):
    dae0 = case.dae
    net = dae0.network
    pos = net.position(3)
    yb = net.ybus.copy()
    yb[pos, pos] += -1j * 0.5
    dae_p = replace(dae0, network=replace(net, ybus=yb))
    x0, z0 = case.equilibrium.x.copy(), case.equilibrium.z.copy()
    nw = _Network(dae0, x0, z0)
    ts, xs = [], []
    x = x0
    for t0, t1, dae in ((0.0, t_pulse, dae_p), (t_pulse, t_end, dae0)):
        def rhs(t, xx, dae=dae):
            z, _ = nw.solve(dae, xx)
            return dae.f(xx, z, {}) if z is not None else np.full(xx.shape, np.nan)

        def jac(t, xx, dae=dae):
            z, _ = nw.solve(dae, xx)
            j = central_difference_jacobians(dae, xx, z, {})
            return j.fx - j.fz @ np.linalg.solve(j.gz, j.gx)

        sol = BDF(rhs, t0, x, t1, jac=jac, rtol=1e-7, atol=1e-9, max_step=max_step)
        while sol.status == "running":
            sol.step()
            if not np.all(np.isfinite(sol.y)):
                raise RuntimeError("non-finite state")
            ts.append(sol.t)
            xs.append(sol.y.copy())
        x = sol.y.copy()
    return np.array(ts), np.array(xs), x0


def decay(ts, xs, x0, vl, win=(11.0, 29.0)):
    q = (xs - x0) @ vl.conj()
    m = (ts >= win[0]) & (ts <= win[1])
    y = np.log(np.maximum(np.abs(q[m]), 1e-300))
    from scipy.signal import find_peaks

    pk, _ = find_peaks(np.abs(q[m]))
    if pk.size >= 4:
        return float(np.polyfit(ts[m][pk], y[pk], 1)[0])
    return float(np.polyfit(ts[m], y, 1)[0])


def select_pairs(n=2):
    import R06_census as C6

    ref = C6.reference_policy()
    recs = [(r["task"]["pid"], r["task"]["S"], r["record"]) for r in R.Store("R06A").all() if r["task"]["variant"] == "REAL" and r["task"]["pid"] == ref["id"]]
    cen = V.census_from_records(recs)
    d = cen[ref["id"]]
    mg = V.marginals(ref["id"], d)
    e = mg[mg.lvT & np.isfinite(mg.t_delta) & (mg.cls_T != 0)]
    cands = []
    for i in V.V:
        ei = e[e.i == i]
        for c in (1, -1):
            ms = ei[ei.cls_T == c]
            dl = dict(zip(ms["mask"], ms.t_delta, strict=True))
            for S1, S2 in V._nested(ms["mask"].to_numpy(), ms["mask"].to_numpy()):
                c1, c2 = V.pick(d["modes"][int(S1)], "em_top"), V.pick(d["modes"][int(S2)], "em_top")
                if c1 is None or c2 is None or abs(c1["hz"] - c2["hz"]) > R.DF_MAX or V.mac_c(c1["phi"], c2["phi"]) < R.MAC_MIN:
                    continue
                cands.append((min(abs(dl[S1]), abs(dl[S2])), int(i), int(S1), int(S2), float(dl[S1]), float(dl[S2])))
    cands.sort(reverse=True)
    out, seen = [], set()
    for c in cands:
        if c[1] in seen:
            continue
        seen.add(c[1])
        out.append(c)
        if len(out) == n:
            break
    return ref, out


def main():
    R.env_threads()
    ref, pairs = select_pairs()
    res = {"label": "POST HOC exploratory (preregistered R14 selection yields no case)", "reference_policy": ref["id"], "pairs": []}
    for mag, i, S1, S2, d1, d2 in pairs:
        rec = {"i": i, "S1": R.label(R.members_of(S1)), "S2": R.label(R.members_of(S2)), "lin_d1": d1, "lin_d2": d2, "portfolios": {}}
        for tag, m in (("S1", S1), ("S1i", S1 | (1 << V.V.index(i))), ("S2", S2), ("S2i", S2 | (1 << V.V.index(i)))):
            case = R.build68(R.members_of(m), variant="REAL", g=ref["g"], k=ref["k"])
            lam, vl, _ = em_top_left(case)
            ts, xs, x0 = simulate(case)
            q = np.abs((xs - x0) @ vl.conj())
            grid = np.linspace(0.0, 30.0, 601)
            rec["portfolios"][tag] = {"S": R.label(R.members_of(m)), "lin_re": lam.real, "lin_hz": abs(lam.imag) / (2 * np.pi), "tds_decay": decay(ts, xs, x0, vl),
                                      "n_steps": int(ts.size), "trace_t": grid.tolist(), "trace_absq": np.interp(grid, ts, q).tolist()}
        p = rec["portfolios"]
        rec["tds_d1"] = p["S1i"]["tds_decay"] - p["S1"]["tds_decay"]
        rec["tds_d2"] = p["S2i"]["tds_decay"] - p["S2"]["tds_decay"]
        rec["sign_agree_1"] = bool(np.sign(rec["tds_d1"]) == np.sign(d1))
        rec["sign_agree_2"] = bool(np.sign(rec["tds_d2"]) == np.sign(d2))
        res["pairs"].append(rec)
        print(json.dumps({k: v for k, v in rec.items() if k != "portfolios"}), flush=True)
    R.write_json(R.RESULTS / "CDW68_R14_tds_posthoc.json", res)


if __name__ == "__main__":
    main()
