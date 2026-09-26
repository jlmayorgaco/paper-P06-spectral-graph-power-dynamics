# ruff: noqa: E501  -- library tables kept on one line
"""SQ1-1: scientific measurement domain from the frozen phasor models (.venv/tx3-analysis).

For every frozen scientific phasor model that the frozen G2 builder can construct, EXCEPT the nine
T1-T9 instrument-holdout models (kept independent), compute the linear modal response of the v1
section-5 estimator channels to the preregistered pulse (+2 % bus-20 P for 0.2 s):
    channel c after the pulse (t >= 1.2 s):  y_c(t) = sum_m r_cm exp(lambda_m (t - 1.2)),
    r_cm = (C_c v_m)(w_m . B) dP (exp(0.2 lambda_m) - 1)/lambda_m,
with C = d(channels)/dx and B = d(xdot)/d(P_load20) by central differences through the algebraic
equations, (lambda, v, w) the eigen-decomposition of the frozen A. Structural modes (|lambda| < 1e-3)
and modes decaying faster than 5 s^-1 are dropped (gone before the 2.2-s window start).
Also validates the linear responses against the non-holdout development TDS trace DEV_P4_BASE.
Writes results/EMTSQ1/domain/modal_library.npz and domain_summary.json. No estimator is run.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.linalg import lu_factor, lu_solve  # noqa: E402

import G2_tds as G2  # noqa: E402

RESEARCH = HERE.parents[2]
REPO = HERE.parents[6]
OUT = RESEARCH / "results" / "EMTSQ1" / "domain"
H4 = (30, 33, 35, 37)
GEN_BUSES = (30, 31, 32, 33, 34, 35, 36, 37, 38, 39)
W_B = 2 * np.pi * 60
HOLDOUT_IDS = {"EMT05_P4_30+33+35+37", "EMT05_P4_30+33+35", "EMT05_P4_33", "EMT05_P4_37", "EMT08_g0.180", "EMT08_g0.250", "EMT08_g0.205",
               "EMT13_cond2.0pct", "EMT13_cond3.0pct", "EMT04_G_S_H4"}
T_PULSE, DP_FRAC = 0.2, 0.02


def spec_of(row, cases):
    c = cases[row.case_id]
    g, k, t, h = c["theta"]
    var = c.get("variant", {"kind": "none"})
    if var.get("kind") not in ("none", "condenser"):
        return None
    cond = ("R_0010000", var["fraction"]) if var.get("kind") == "condenser" else None
    return tuple(c["members"]), {"g": g, "k": k}, cond


def solve_z(dae, x, z0):
    z = z0.copy()
    for _ in range(40):
        r = dae.g(x, z, {})
        if np.abs(r).max() < 1e-12:
            return z
        cols = []
        for kk in range(z.size):
            dz = np.zeros_like(z)
            dz[kk] = 1e-7
            cols.append((dae.g(x, z + dz, {}) - dae.g(x, z - dz, {})) / 2e-7)
        z = z - lu_solve(lu_factor(np.column_stack(cols)), r)
    raise RuntimeError("algebraic solve failed")


def channel_fn(case):
    dae = case.dae
    labels = list(case.system.labels)
    net = dae.network
    sg = sorted(int(n[len("omega_sg"):]) for n in labels if n.startswith("omega_sg"))
    gf = sorted(int(n[len("theta_pll_gfl"):]) for n in labels if n.startswith("theta_pll_gfl"))
    i_w = {b: labels.index(f"omega_sg{b}") for b in sg}
    i_th = {b: labels.index(f"theta_pll_gfl{b}") for b in gf}
    pos = {b: net.position(b) for b in GEN_BUSES}
    names = [f"relspeed_{b}" for b in sg if b != 39] + [f"gflfreq_{b}" for b in gf] + [f"angdiff_{b}" for b in GEN_BUSES if b != 39] + [f"vmag_{b}" for b in GEN_BUSES]

    def h(x, z):
        v = dae.voltages(z)
        dx = dae.f(x, z, {})
        w39 = x[i_w[39]]
        out = [x[i_w[b]] - w39 for b in sg if b != 39]
        out += [dx[i_th[b]] / W_B - (w39 - 1.0) for b in gf]
        out += [np.angle(v[pos[b]] / v[pos[39]]) for b in GEN_BUSES if b != 39]
        out += [abs(v[pos[b]]) for b in GEN_BUSES]
        return np.array(out)

    return h, names


def modal(case):
    dae = case.dae
    x0, z0 = case.equilibrium.x.copy(), case.equilibrium.z.copy()
    h, names = channel_fn(case)
    n = x0.size
    C = np.zeros((len(names), n))
    for j in range(n):
        e = 1e-6 * max(1.0, abs(x0[j]))
        xp, xm = x0.copy(), x0.copy()
        xp[j] += e
        xm[j] -= e
        C[:, j] = (h(xp, solve_z(dae, xp, z0)) - h(xm, solve_z(dae, xm, z0))) / (2 * e)
    net = dae.network
    s20 = net.loads[20]
    d = 1e-6 * abs(s20.real)
    net.loads[20] = complex(s20.real + d, s20.imag)
    fp = dae.f(x0, solve_z(dae, x0, z0), {})
    net.loads[20] = complex(s20.real - d, s20.imag)
    fm = dae.f(x0, solve_z(dae, x0, z0), {})
    net.loads[20] = s20
    B = (fp - fm) / (2 * d)
    dP = DP_FRAC * s20.real
    lam, V = np.linalg.eig(case.system.A)
    W = np.linalg.inv(V)
    keep = (np.abs(lam) > 1e-3) & (lam.real > -5.0) & (lam.imag >= 0)
    lam_k = lam[keep]
    fac = (np.exp(T_PULSE * lam_k) - 1.0) / lam_k * dP
    R = (C @ V[:, keep]) * (W[keep] @ B)[None, :] * fac[None, :]
    return lam_k, R, names, x0, z0, h, C, B, dP


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    pred = pd.read_csv(RESEARCH / "results" / "EMT_PRED" / "phasor_predictions.csv")
    cases = {c["case_id"]: c for c in json.loads((RESEARCH / "results" / "EMT_PRED" / "emt_cases.json").read_text())["cases"]}
    lib, summary = {}, {"cases": {}}
    for row in pred.itertuples():
        if row.case_id in HOLDOUT_IDS or row.experiment not in ("EMT05", "EMT08", "EMT09", "EMT10", "EMT11", "EMT13"):
            continue
        spec = spec_of(row, cases)
        if spec is None:
            continue
        members, point, cond = spec
        case = G2.build("IEEE-39", members, point, cond)
        lam, R, names, *_ = modal(case)
        target = int(np.argmin(np.abs(lam - complex(row.band_re, 2 * np.pi * row.band_hz))))
        lib[row.case_id] = (lam, R, names, target)
        summary["cases"][row.case_id] = {"n_ch": len(names), "n_modes": int(lam.size), "target_lambda": [float(lam[target].real), float(lam[target].imag / (2 * np.pi))],
                                         "band_re": row.band_re, "band_hz": row.band_hz, "target_match_err": float(abs(lam[target] - complex(row.band_re, 2 * np.pi * row.band_hz)))}
        print(row.case_id, len(names), lam.size, flush=True)
    # ---- channel observability and nuisance structure over the window [2.2, 30] s ----
    tau = np.linspace(1.0, 28.8, 2781)  # t - 1.2 for t in [2.2, 30]
    obs, nuis = [], []
    for cid, (lam, R, names, tgt) in lib.items():
        en = np.zeros(R.shape)
        for m in range(lam.size):
            mult = 2.0 if abs(lam[m].imag) > 1e-9 else 1.0
            sig = mult * np.abs(R[:, m])[:, None] * np.exp(lam[m].real * tau)[None, :]
            en[:, m] = np.sum(sig**2, axis=1) / 2.0 if mult == 2.0 else np.sum(sig**2, axis=1)
        tot = en.sum(axis=1)
        share = en[:, tgt] / np.maximum(tot, 1e-300)
        obs += share.tolist()
        f = lam.imag / (2 * np.pi)
        for m in range(lam.size):
            if m == tgt:
                continue
            ms = float(np.median(en[:, m] / np.maximum(tot, 1e-300)))
            if ms >= 0.01:
                nuis.append({"case": cid, "d_alpha": float(lam[m].real - lam[tgt].real), "d_f": float(f[m] - f[tgt]), "f": float(f[m]), "alpha": float(lam[m].real),
                             "median_energy_share": ms, "real_mode": bool(abs(lam[m].imag) < 1e-9)})
    obs = np.array(obs)
    nd = pd.DataFrame(nuis)
    summary["n_library_cases"] = len(lib)
    summary["channels"] = {"n_min": int(min(v[2].__len__() for v in lib.values())), "n_max": int(max(len(v[2]) for v in lib.values()))}
    summary["observability_target_energy_share"] = {q: float(np.quantile(obs, p)) for q, p in (("p05", 0.05), ("p25", 0.25), ("median", 0.5), ("p75", 0.75), ("p95", 0.95))}
    summary["observability_fraction_below_0.1"] = float(np.mean(obs < 0.1))
    summary["nuisance_modes_ge_1pct"] = {"count": len(nd), "per_case_mean": float(len(nd) / max(len(lib), 1)),
                                         "d_f_range": [float(nd.d_f.min()), float(nd.d_f.max())] if len(nd) else None,
                                         "d_alpha_range": [float(nd.d_alpha.min()), float(nd.d_alpha.max())] if len(nd) else None,
                                         "real_modes": int(nd.real_mode.sum()) if len(nd) else 0}
    nd.to_csv(OUT / "nuisance_modes.csv", index=False)
    # ---- validation against the non-holdout development TDS trace DEV_P4_BASE ----
    dev = np.load(REPO / "external" / "paremt_runs" / "v3_tds" / "DEV_P4_BASE.npz")
    lam, R, names, tgt = lib["EMT05_P4_BASE"]
    t = dev["t"]
    m = (t >= 2.2) & (t <= 30.0)
    ylin = np.zeros((len(names), m.sum()))
    for k in range(lam.size):
        term = R[:, k][:, None] * np.exp(lam[k] * (t[m] - 1.2))[None, :]
        ylin += 2.0 * term.real if abs(lam[k].imag) > 1e-9 else term.real
    ytds = dev["y"][:, m] - dev["y"][:, :1]
    ytds = ytds - ytds[:, -1:] * 0.0
    # compare oscillatory parts (remove each channel's affine fit on both)
    tt = t[m]
    A = np.column_stack([np.ones_like(tt), tt])
    def detr(y):
        c, *_ = np.linalg.lstsq(A, y.T, rcond=None)
        return y - (A @ c).T
    dl, dt_ = detr(ylin), detr(ytds)
    rel = np.linalg.norm(dl - dt_, axis=1) / np.maximum(np.linalg.norm(dt_, axis=1), 1e-300)
    summary["validation_DEV_P4_BASE"] = {"median_rel_err_oscillatory": float(np.median(rel)), "max_rel_err_oscillatory": float(rel.max()),
                                         "channel_names_equal": list(dev["names"]) == names}
    np.savez_compressed(OUT / "modal_library.npz", **{f"{cid}__lam": v[0] for cid, v in lib.items()}, **{f"{cid}__R": v[1] for cid, v in lib.items()},
                        **{f"{cid}__names": np.array(v[2]) for cid, v in lib.items()}, **{f"{cid}__target": np.array(v[3]) for cid, v in lib.items()})
    (OUT / "domain_summary.json").write_text(json.dumps(summary, indent=1, default=float), encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k != "cases"}, indent=1, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
