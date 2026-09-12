# ruff: noqa: E501  -- case tables kept on one line
"""Gate E0 phasor-TDS instrument holdout traces (prereg V3 section 4.2; .venv/tx3-analysis).

Imports the frozen experiments/G2_tds.py unchanged (build, critical, observable, simulate, D2 pulse)
and records the v1 section-5 estimator channels from the phasor states: omega_sg,b - omega_sg39;
GFL theta'_pll/w_B - (omega_39 - 1) (theta' from the DAE right-hand side); angle(V_b) - angle(V_39)
and |V_b| at the ten generator buses. Sampling 5 ms; time shifted by +1.0 s (pulse at [1.0, 1.2)).
Usage: EMTV3_tds_traces.py dev      -> development case DEV_P4_BASE only (not a holdout)
       EMTV3_tds_traces.py holdout  -> T1-T9 (only after the implementation commit)
Writes external/paremt_runs/v3_tds/<id>.npz and results/EMTV3/E0/tds_manifest_<mode>.json.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))

import G2_tds as G2  # noqa: E402
import numpy as np  # noqa: E402
from scipy.linalg import lu_factor, lu_solve  # noqa: E402

from ibr_cycles.dynamics.linearize import central_difference_jacobians  # noqa: E402

RESEARCH = HERE.parents[2]
REPO = HERE.parents[6]
BIG = REPO / "external" / "paremt_runs" / "v3_tds"
OUT = RESEARCH / "results" / "EMTV3" / "E0"
H4 = (30, 33, 35, 37)
GEN_BUSES = (30, 31, 32, 33, 34, 35, 36, 37, 38, 39)
W_B = 2 * np.pi * 60
DT_S = 0.005
HORIZON = 89.0
SHIFT = 1.0
HOLDOUT = {
    "T1": (H4, {"g": 0.03625, "k": 1.425}, None, "EMT05_P4_30+33+35+37"),
    "T2": ((30, 33, 35), {"g": 0.03625, "k": 1.425}, None, "EMT05_P4_30+33+35"),
    "T3": ((33,), {"g": 0.03625, "k": 1.425}, None, "EMT05_P4_33"),
    "T4": ((37,), {"g": 0.03625, "k": 1.425}, None, "EMT05_P4_37"),
    "T5": (H4, {"g": 0.18, "k": 1.425}, None, "EMT08_g0.180"),
    "T6": (H4, {"g": 0.25, "k": 1.425}, None, "EMT08_g0.250"),
    "T7": (H4, {"g": 0.205, "k": 1.425}, None, "EMT08_g0.205"),
    "T8": (H4, {"g": 0.03625, "k": 1.425}, ("R_0010000", 0.020), "EMT13_cond2.0pct"),
    "T9": (H4, {"g": 0.03625, "k": 1.425}, ("R_0010000", 0.030), "EMT13_cond3.0pct"),
}
DEV = {"DEV_P4_BASE": ((), {"g": 0.03625, "k": 1.425}, None, "EMT05_P4_BASE")}


def solve_z(dae, x, z, lu, gz_full):
    for _ in range(20):
        r = dae.g(x, z, {})
        if np.abs(r).max() < 1e-11:
            return z, lu
        z = z - lu_solve(lu, r)
    for _ in range(12):
        r = dae.g(x, z, {})
        if np.abs(r).max() < 1e-11:
            return z, lu
        lu = lu_factor(gz_full(x, z))
        z = z - lu_solve(lu, r)
    raise RuntimeError("algebraic solve failed")


def channels(case, sol, t_end):
    dae = case.dae
    labels = list(case.system.labels)
    net = dae.network
    ts = np.arange(0.0, t_end - 1e-9, DT_S)
    xs = sol(ts)
    z = case.equilibrium.z.copy()
    lu = lu_factor(central_difference_jacobians(dae, case.equilibrium.x, z, {}).gz)

    def gz_full(x, zz, h=1e-7):
        cols = []
        for k in range(zz.size):
            dz = np.zeros_like(zz)
            dz[k] = h
            cols.append((dae.g(x, zz + dz, {}) - dae.g(x, zz - dz, {})) / (2 * h))
        return np.column_stack(cols)

    sg = sorted(int(n[len("omega_sg"):]) for n in labels if n.startswith("omega_sg"))
    gf = sorted(int(n[len("theta_pll_gfl"):]) for n in labels if n.startswith("theta_pll_gfl"))
    i_w = {b: labels.index(f"omega_sg{b}") for b in sg}
    i_th = {b: labels.index(f"theta_pll_gfl{b}") for b in gf}
    pos = {b: net.position(b) for b in GEN_BUSES}
    out = np.zeros((len(sg) - 1 + len(gf) + 2 * len(GEN_BUSES) - 1, ts.size))
    names = [f"relspeed_{b}" for b in sg if b != 39] + [f"gflfreq_{b}" for b in gf] + [f"angdiff_{b}" for b in GEN_BUSES if b != 39] + [f"vmag_{b}" for b in GEN_BUSES]
    ang_prev = None
    for k in range(ts.size):
        x = xs[:, k]
        z, lu = solve_z(dae, x, z, lu, gz_full)
        v = dae.voltages(z)
        dx = dae.f(x, z, {})
        w39 = x[i_w[39]]
        row = [x[i_w[b]] - w39 for b in sg if b != 39]
        row += [dx[i_th[b]] / W_B - (w39 - 1.0) for b in gf]
        ang = np.array([np.angle(v[pos[b]]) - np.angle(v[pos[39]]) for b in GEN_BUSES if b != 39])
        if ang_prev is not None:
            ang = ang_prev + np.angle(np.exp(1j * (ang - ang_prev)))
        ang_prev = ang
        row += list(ang)
        row += [abs(v[pos[b]]) for b in GEN_BUSES]
        out[:, k] = row
    # shift by +1.0 s: equilibrium samples on [0, 1.0)
    npre = int(round(SHIFT / DT_S))
    full = np.concatenate([np.repeat(out[:, :1], npre, axis=1), out], axis=1)
    t_full = np.arange(full.shape[1]) * DT_S
    return t_full, full, names


def main() -> int:
    mode = sys.argv[1] if len(sys.argv) > 1 else "dev"
    cases = DEV if mode == "dev" else HOLDOUT
    BIG.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    man = {}
    for cid, (members, point, cond, pred_id) in cases.items():
        t0 = time.time()
        case = G2.build("IEEE-39", members, point, cond)
        lam, vec, rhp = G2.critical(case)
        (i, j), speeds = G2.observable(case, vec)
        sol, t_end, status = G2.simulate(case, G2.DISTURBANCES["IEEE-39"]["D2"], HORIZON, speeds)
        t, y, names = channels(case, sol, t_end)
        p = BIG / f"{cid}.npz"
        np.savez_compressed(p, t=t, y=y, names=np.array(names))
        man[cid] = {"members": list(members), "point": point, "condenser": list(cond) if cond else None, "pred_case_id": pred_id,
                    "g2_critical_re": float(lam.real), "g2_critical_hz": float(abs(lam.imag) / (2 * np.pi)), "rhp": rhp,
                    "run_status": status, "t_end_after_pulse_start_s": t_end, "record_end_s": float(t[-1]), "n_channels": len(names),
                    "file": str(p.relative_to(REPO)).replace("\\", "/"), "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
                    "wall_s": round(time.time() - t0, 1)}
        print(cid, status, round(t_end, 2), len(names), man[cid]["wall_s"], flush=True)
    (OUT / f"tds_manifest_{mode}.json").write_text(json.dumps(man, indent=1), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
