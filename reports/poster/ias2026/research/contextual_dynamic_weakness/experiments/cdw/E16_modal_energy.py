# ruff: noqa: E501
"""E16 modal energy vs stability-limiting mode (phasor domain, prereg E16)."""

from __future__ import annotations

from dataclasses import replace

import _infra as I
import numpy as np
import pandas as pd

import _cdw as C
from ibr_cycles.models.ieee39_case import Ieee39Dae

PHASE = "E16"
POLICIES = list(C.DISCOVERY) + [f"H{i:02d}" for i in range(1, 25)]
PULSE_BUSES = (16, 20)
AMP, DUR, T0, T1 = 0.1, 0.2, 2.0, 30.0


def tasks():
    return [{"phase": PHASE, "pid": p, "S": list(s)} for p in POLICIES for s in C.subsets(C.V4)]


def energy(a, lam, t0, t1):
    """integral_{t0}^{t1} (2 Re(a e^{lam t}))^2 dt."""

    def prim(t):
        s = 2 * lam.real
        v1 = 2 * abs(a) ** 2 * (np.exp(s * t) / s if s != 0 else t)
        l2 = 2 * lam
        v2 = 2 * (a * a * np.exp(l2 * t) / l2).real
        return v1 + v2

    return float(prim(t1) - prim(t0))


def run_task(task):
    theta = C.policy(task["pid"])
    try:
        case = C.solve(tuple(task["S"]), theta)
    except (C.InfeasibleReplacement, ValueError, np.linalg.LinAlgError):
        return {"status": "INFEASIBLE"}
    ev0 = C.transverse_eval(case, modes=False)
    if ev0["status"] != "STABLE":
        return {"status": ev0["status"], "skipped": True}
    a, _, jac = C.physical_matrices(case)
    dae = case.dae
    x0, z0 = case.equilibrium.x, case.equilibrium.z
    vals, vr = np.linalg.eig(a)
    vl = np.linalg.inv(vr)  # rows: left eigenvectors with vl @ vr = I
    ev_t = C.transverse_eval(case, modes=False)
    del ev_t
    # modal set: oscillatory transverse modes in 0.1-2.0 Hz (upper half plane, non-structural)
    f = vals.imag / (2 * np.pi)
    keep = [k for k in range(vals.size) if 0.1 <= f[k] <= 2.0 and abs(vals[k]) > 1e-4]
    labels = dae.labels
    speed_rows = [k for k, lab in enumerate(labels) if lab.startswith("omega_sg")]
    V0 = z0[0::2] + 1j * z0[1::2]
    Gz_inv_gx = np.linalg.solve(jac.gz, jac.gx)
    out = {"status": "STABLE", "cases": []}
    for bus in PULSE_BUSES:
        loads = dict(dae.network.loads)
        loads[bus] = loads.get(bus, 0j) + 1e-4
        net2 = replace(dae.network, loads=loads)
        dae2 = Ieee39Dae(network=net2, power_flow=dae.power_flow, slots=dae.slots, plan=dae.plan)
        gP = (dae2.g(x0, z0, {}) - dae.g(x0, z0, {})) / 1e-4
        b = -jac.fz @ np.linalg.solve(jac.gz, gP)
        E = np.zeros((len(keep),))
        shares = []
        for c_kind in ("speed", "vmag"):
            if c_kind == "speed":
                Cm = np.zeros((len(speed_rows), a.shape[0]))
                for r, k in enumerate(speed_rows):
                    Cm[r, k] = 1.0
            else:
                dv = -Gz_inv_gx  # 78 x n: bus-voltage perturbation per state
                dvc = dv[0::2] + 1j * dv[1::2]
                Cm = np.real(np.conj(V0)[:, None] * dvc) / np.abs(V0)[:, None]
            Ek = np.zeros((Cm.shape[0], len(keep)))
            for j, k in enumerate(keep):
                lam = vals[k]
                ck = AMP * (vl[k] @ b) * (np.exp(DUR * lam) - 1) / lam
                amp = (Cm @ vr[:, k]) * ck
                for r in range(Cm.shape[0]):
                    Ek[r, j] = energy(amp[r], lam, T0 - DUR, T1 - DUR)
            tot = Ek.sum(1, keepdims=True)
            sh = np.where(tot > 0, Ek / np.where(tot > 0, tot, 1), 0)
            shares.append(sh)
            E += Ek.sum(0)
        sh_all = np.vstack(shares).mean(0)
        re = vals[keep].real
        k_lim = int(np.argmax(re))
        k_dom = int(np.argmax(sh_all))
        k_dom_raw = int(np.argmax(E))
        out["cases"].append({"bus": bus, "n_modes": len(keep), "limit_hz": float(f[keep][k_lim]), "limit_re": float(re[k_lim]),
                             "dom_hz": float(f[keep][k_dom]), "dom_re": float(re[k_dom]), "dom_share": float(sh_all[k_dom]),
                             "limit_share": float(sh_all[k_lim]), "differ": bool(k_dom != k_lim), "differ_raw": bool(k_dom_raw != k_lim),
                             "re_gap": float(re[k_lim] - re[k_dom])})
    return out


def aggregate():
    rows = []
    for r in I.Store(PHASE).all():
        if not r.get("ok") or r.get("status") != "STABLE":
            continue
        t = r["task"]
        for c in r["cases"]:
            rows.append({"pid": t["pid"], "S": C.label(t["S"]), **c})
    df = pd.DataFrame(rows)
    df.to_csv(I.RESULTS / "CDW_E16_modal_energy.csv", index=False)
    per = df.groupby("pid").differ.mean()
    hold = per[per.index.str.startswith("H")]
    gate = {"frac_differ_all": float(df.differ.mean()), "frac_differ_raw_all": float(df.differ_raw.mean()),
            "holdout_frac_policies_ge_0.25": float((hold >= 0.25).mean()), "systematic": bool((hold >= 0.25).mean() >= 0.5),
            "median_re_gap_when_differ": float(df[df.differ].re_gap.median()) if df.differ.any() else np.nan,
            "n_cases": int(len(df))}
    I.atomic_write_json(I.RESULTS / "CDW_E16_gates.json", gate)
    return gate


if __name__ == "__main__":
    print(aggregate())
