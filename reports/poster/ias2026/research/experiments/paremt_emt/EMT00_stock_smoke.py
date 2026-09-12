# ruff: noqa: E501  -- long diagnostic strings and code-patch literals kept on one line
"""EMT00b - stock ParaEMT IEEE-39 smoke test (installation / tool check only; no TX4 meaning).

Runs the UNMODIFIED official ParaEMT (upstream d79d735a) IEEE-39 case (systemN = 3) from a
disposable copy of the upstream tree (external/paremt_runs/EMT00_stock), so that the pristine
clone never receives output files. Serial LU network solver. Stock loads (const RLC) and
stock models, nothing modified.

Checks: initialization, no-event stationarity, a small governor-reference step (the public
exciter step is a no-op: StepChange edits ini.Init_mac_vref, which the solver never reads), pickle
serialization / reload, snapshot creation and a run resumed from the snapshot.

Usage (xtool-paremt):  python EMT00_stock_smoke.py
Writes results/EMT00/stock_smoke_summary.json and stock_smoke_traces.npz.
"""

from __future__ import annotations

import filecmp
import json
import os
import pickle
import shutil
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
REPO = HERE.parents[5]
UP = REPO / "external" / "ParaEMT_upstream"
RUN = REPO / "external" / "paremt_runs" / "EMT00_stock"
OUT = RESEARCH / "results" / "EMT00"
TS = 50e-6
DSRATE = 20  # 1 ms stored interval
T_NOEVENT = 10.0
T_STEP = 6.0
T_RESUME = 1.0


def prepare_run_dir() -> None:
    if RUN.exists():
        shutil.rmtree(RUN)
    shutil.copytree(UP, RUN, ignore=shutil.ignore_patterns(".git"))
    # provenance: every copied file is byte-identical to the pristine upstream
    for p in RUN.rglob("*"):
        if p.is_file():
            assert filecmp.cmp(p, UP / p.relative_to(RUN), shallow=False), p


def simulate(emt, pfd, dyd, ini, tlen, net_mod="lu"):
    """The main_step1_simulation.py time loop, verbatim in its calls, without prints."""

    tn = 0
    tsave = 0
    t0 = time.time()
    while tn * TS < tlen:
        tn += 1
        emt.StepChange(dyd, ini, tn)
        emt.GenTrip(pfd, dyd, ini, tn, net_mod)
        emt.predictX(pfd, dyd, emt.ts)
        emt.Igs = emt.Igs * 0
        emt.updateIg(pfd, dyd, ini)
        emt.Igi = emt.Igi * 0
        emt.Iibr = emt.Iibr * 0
        emt.updateIibr(pfd, dyd, ini)
        if emt.loadmodel_option != 1:
            emt.Il = emt.Il * 0
            emt.updateIl(pfd, dyd, tn)
        emt.solveV(ini)
        emt.BusMea(pfd, dyd, tn)
        emt.updateX(pfd, dyd, ini, tn)
        emt.updateXibr(pfd, dyd, ini, TS)
        if emt.loadmodel_option != 1:
            emt.updateXl(pfd, dyd, tn)
        emt.x_pred = {0: emt.x_pred[1], 1: emt.x_pred[2], 2: emt.x_pv_1}
        if np.mod(tn, DSRATE) == 0:
            tsave += 1
            emt.t.append(tn * TS)
            emt.x[tsave] = emt.x_pv_1.copy()
            if len(pfd.ibr_bus) > 0:
                emt.x_ibr[tsave] = emt.x_ibr_pv_1.copy()
            if len(pfd.bus_num) > 0:
                emt.x_bus[tsave] = emt.x_bus_pv_1.copy()
            if len(pfd.load_bus) > 0:
                emt.x_load[tsave] = emt.x_load_pv_1.copy()
            emt.v[tsave] = emt.Vsol.copy()
        if (emt.flag_gentrip == 0) & (emt.flag_reinit == 1):
            emt.Re_Init(pfd, dyd, ini)
        else:
            emt.updateIhis(ini)
    return time.time() - t0, tn


def stack(d):
    return np.array([np.asarray(d[k]) for k in sorted(d)])


def vmag(v, nbus):
    va, vb, vc = v[:, :nbus], v[:, nbus : 2 * nbus], v[:, 2 * nbus :]
    return np.sqrt((2.0 / 3.0) * (va**2 + vb**2 + vc**2))


def run_case(label, tlen, step=None, snapshot=None):
    from psutils import initialize_emt, initialize_from_snp

    t_init = time.time()
    if snapshot is None:
        pfd, ini, dyd, emt = initialize_emt(".", 3, 1, 1, TS, tlen, mode="lu", nparts=2)
    else:
        pfd, ini, dyd, emt = initialize_from_snp(snapshot, "lu", 2)
    t_init = time.time() - t_init
    emt.t_gentrip = 0  # falsy: GenTrip disabled (no event)
    emt.t_sc = 1e9
    emt.flag_sc = 0
    emt.t_release_f = 0.0
    emt.loadmodel_option = 1
    if step is not None:
        # governor-reference step (the public exciter step is a no-op: StepChange edits ini.Init_mac_vref)
        emt.t_sc, emt.i_gen_sc, emt.flag_exc_gov, emt.dsp, emt.flag_sc = (
            step["t"],
            step["gen"],
            1,
            step["dgref"],
            1,
        )
    wall, nsteps = simulate(emt, pfd, dyd, ini, tlen)
    return pfd, ini, dyd, emt, wall, nsteps, t_init


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    prepare_run_dir()
    os.chdir(RUN)
    sys.path.insert(0, str(RUN))
    summary = {
        "paremt_commit": "d79d735a4a587d56c5b88187d1a499195b6b2b84",
        "systemN": 3,
        "ts_s": TS,
        "stored_interval_s": TS * DSRATE,
        "net_solver": "lu (serial)",
        "load_model": "stock const RLC (option 1)",
    }

    # --- A: no-event run, then snapshot ---
    pfd, ini, dyd, emt, wall, nsteps, t_init = run_case("noevent", T_NOEVENT)
    nbus, ngen = len(pfd.bus_num), len(pfd.gen_bus)
    x = stack(emt.x)
    v = stack(emt.v)
    t = np.asarray(emt.t)
    vm = vmag(v, nbus)
    speed = x[:, 1 : 18 * ngen : 18] / pfd.ws
    delta = x[:, 0 : 18 * ngen : 18]
    rel_delta = delta - delta[:, [0]]
    summary["A_noevent"] = {
        "tlen_s": T_NOEVENT,
        "steps": nsteps,
        "wall_s": round(wall, 2),
        "init_s": round(t_init, 2),
        "us_per_step": round(1e6 * wall / nsteps, 1),
        "network_dimension": int(ini.Init_net_N),
        "bus_count": nbus,
        "line_count": int(len(pfd.line_from)),
        "transformer_count": int(len(pfd.xfmr_from)),
        "generator_count": ngen,
        "load_count": int(len(pfd.load_bus)),
        "shunt_count": int(len(pfd.shnt_bus)),
        "ibr_count": int(len(pfd.ibr_bus)),
        "max_speed_drift_pu": float(np.abs(speed - speed[0]).max()),
        "max_relative_angle_drift_rad": float(np.abs(rel_delta - rel_delta[0]).max()),
        "max_voltage_magnitude_drift_pu": float(np.abs(vm - vm[0]).max()),
        "max_state_drift_any": float(np.abs(x - x[0]).max()),
        # drift after an initial settling window (t >= 5 s), relative to the value at t = 5 s
        "settled_max_speed_drift_pu": float(
            np.abs(speed[t >= 5.0] - speed[t >= 5.0][0]).max()
        ),
        "settled_max_voltage_magnitude_drift_pu": float(
            np.abs(vm[t >= 5.0] - vm[t >= 5.0][0]).max()
        ),
        "settled_max_relative_angle_drift_rad": float(
            np.abs(rel_delta[t >= 5.0] - rel_delta[t >= 5.0][0]).max()
        ),
        "initial_voltage_step_first_1ms_pu": float(np.abs(vm[1] - vm[0]).max()),
    }
    np.savez_compressed(
        OUT / "stock_smoke_traces.npz",
        t_noevent=t,
        speed_noevent=speed,
        vm_noevent=vm[:, :nbus],
    )
    emt.dump_res(pfd, dyd, ini, 0, "snp_ful.pkl", "snp_1pt.pkl", "res.pkl")

    # --- B: serialization reload ---
    with open("snp_ful.pkl", "rb") as f:
        pfd_r, dyd_r, ini_r, emt_r = pickle.load(f)
    summary["B_serialization"] = {
        "reloaded": True,
        "x_shape": list(np.asarray(emt_r.x).shape),
        "v_shape": list(np.asarray(emt_r.v).shape),
        "x_identical": bool(np.array_equal(np.asarray(emt_r.x).T, x)),
    }

    # --- C: resume from the snapshot, no event ---
    last_state = x[-1].copy()
    pfd2, ini2, dyd2, emt2, wall2, n2, _ = run_case(
        "resume", T_RESUME, snapshot="snp_1pt.pkl"
    )
    x2 = stack(emt2.x)
    summary["C_snapshot_resume"] = {
        "first_state_equals_snapshot": bool(
            np.allclose(x2[0], last_state, atol=0, rtol=0)
        ),
        "max_state_drift_over_1s": float(np.abs(x2 - x2[0]).max()),
        "max_speed_drift_pu": float(
            np.abs(
                (x2[:, 1 : 18 * ngen : 18] - x2[0, 1 : 18 * ngen : 18]) / pfd2.ws
            ).max()
        ),
        "wall_s": round(wall2, 2),
    }

    # --- D: small disturbance: +0.01 pu governor-reference step on generator index 0 at t = 1 s ---
    step = {"t": 1.0, "gen": 0, "dgref": 0.01}
    pfd3, ini3, dyd3, emt3, wall3, n3, _ = run_case("step", T_STEP, step=step)
    x3 = stack(emt3.x)
    t3 = np.asarray(emt3.t)
    v3 = stack(emt3.v)
    vm3 = vmag(v3, nbus)
    gbus = int(np.where(pfd3.bus_num == pfd3.gen_bus[0])[0][0])
    sp3 = x3[:, 1 : 18 * ngen : 18] / pfd3.ws
    pe3 = x3[:, 16]
    summary["D_governor_step"] = {
        "step": step,
        "gen_bus": int(pfd3.gen_bus[0]),
        "pe_gen0_before_machine_pu": float(pe3[t3 < 1.0][-1]),
        "pe_gen0_final_machine_pu": float(pe3[-1]),
        "pe_gen0_change_machine_pu": float(pe3[-1] - pe3[t3 < 1.0][-1]),
        "terminal_vm_change_pu": float(vm3[-1, gbus] - vm3[t3 < 1.0, gbus][-1]),
        "max_relative_speed_excursion_pu": float(
            np.abs(sp3 - sp3[:, [0]] - (sp3[0] - sp3[0, 0])).max()
        ),
        "finite": bool(np.isfinite(x3).all()),
        "wall_s": round(wall3, 2),
    }
    np.savez_compressed(OUT / "stock_smoke_step.npz", t=t3, speed=sp3, vm=vm3[:, :nbus])
    gate = (
        summary["A_noevent"]["max_voltage_magnitude_drift_pu"] < 1e-3
        and summary["A_noevent"]["max_speed_drift_pu"] < 1e-4
        and summary["B_serialization"]["x_identical"]
        and summary["C_snapshot_resume"]["first_state_equals_snapshot"]
        and summary["D_governor_step"]["finite"]
        and summary["D_governor_step"]["pe_gen0_change_machine_pu"] > 0
    )
    summary["GATE_EMT00"] = "PASS" if gate else "FAIL"
    (OUT / "stock_smoke_summary.json").write_text(
        json.dumps(summary, indent=1), encoding="utf-8"
    )
    print(json.dumps(summary, indent=1))
    return 0 if gate else 2


if __name__ == "__main__":
    raise SystemExit(main())
