# ruff: noqa: E501
"""CDW68 Model B (xtool-andes-gfl ONLY): IEEE-68 with SG68 machines and the frozen WECC library GFL chain.

Composition frozen in docs/CDW68_PREREG_V1.md section 5:
- network rebuilt from the handoff (Bus, Line with taps, PV, Slack, PQ -> constant impedance at the power-flow
  voltage for TDS/linearization; no PF-time conversion);
- SG68D/S/M (runtime-registered, private pycode) at every non-converted generator bus;
- converted buses: TX3-GFL-0.1 chain exactly as H17 add_wecc, REGCP1 Sn = |S_gen|/0.8 * 100 MVA, forced-PQ
  equilibrium holding the unit's base P and Q;
- spectrum: H17 descriptor-free reduced solve, H17 refinements 1/1b (decoupled states), one-dimensional
  structural centre (the governed model has exactly one rotation zero), TX4 verdict band 2e-3.

    python run68.py setup                 private pycode with SG68 (single process)
    python run68.py ybus                  Ybus check against the internal network
    python run68.py qual                  Q1-Q4 -> results/b68/B68_qualification.json
    python run68.py census                16 k x 64 portfolios -> raw/B68C/*.json
    python run68.py branches              branch actions on V68 -> raw/B68L/*.json
    python run68.py draws                 20 envelope draws x 64 -> raw/B68D/*.json
Launch with USERPROFILE/HOME set to .venv/xtool-andes-gfl/home.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import shutil
import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")
assert importlib.util.find_spec("ibr_cycles") is None, "internal package must not be importable"

import andes  # noqa: E402
import andes.models as AM  # noqa: E402
from kvxopt import matrix  # noqa: E402

HERE = Path(__file__).resolve().parent
P68 = HERE.parents[1]
REPO = P68.parents[5]
# the isolated ANDES venv lives in the main checkout (worktrees share it read-only; only home/pycode_cdw68 is written)
VENV = Path(os.environ.get("B68_VENV", r"C:\Users\walla\Documents\Github\paper-P06-spectral-graph-power-dynamics\.venv\xtool-andes-gfl"))
HOME = VENV / "home"
PYCODE_OLD = HOME / "pycode"
PYCODE = HOME / "pycode_cdw68"
B68 = P68 / "results" / "b68"
RAW = P68 / "raw"
EM = (0.1, 2.0)
TAU_MAT = 0.01
MODELS = ["SG68D", "SG68S", "SG68M"]
H: dict = {}
WORKERS = int(os.environ.get("B68_WORKERS", "12"))


def register():
    if "andes.models.cdw68_models" not in sys.modules:
        spec = importlib.util.spec_from_file_location("andes.models.cdw68_models", HERE / "cdw68_models.py")
        mod = importlib.util.module_from_spec(spec)
        sys.modules["andes.models.cdw68_models"] = mod
        spec.loader.exec_module(mod)
    if not any(f == "cdw68_models" for f, _ in AM.file_classes):
        k = next(i for i, (f, _) in enumerate(AM.file_classes) if f == "static")
        AM.file_classes.insert(k + 1, ("cdw68_models", list(MODELS)))


def new_system():
    andes.config_logger(stream_level=50)
    return andes.System(default_config=True, no_output=True, pycode_path=str(PYCODE), autogen_stale=False)


def setup():
    register()
    if not PYCODE.exists():
        shutil.copytree(PYCODE_OLD, PYCODE)
    andes.config_logger(stream_level=20)
    ss = andes.System(default_config=True, no_output=True, pycode_path=str(PYCODE), no_undill=True)
    ss.prepare(quick=True, incremental=True, models=list(MODELS), nomp=True)
    ss2 = new_system()
    files = {m: (PYCODE / f"{m}.py").exists() for m in MODELS}
    print("pycode files:", files, "calls loaded:", {m: ss2.models[m].calls.f is not None for m in MODELS}, PYCODE)


# ----------------------------------------------------------------- building --
def _incident_line(ss, bus):
    for idx, b1, b2 in zip(ss.Line.idx.v, ss.Line.bus1.v, ss.Line.bus2.v, strict=True):
        if int(b1) == bus or int(b2) == bus:
            return idx
    raise ValueError(bus)


def add_wecc(ss, bus, gen_idx, rating, cfg):
    """H17 add_wecc (TX3 _add_gfl), verbatim parameter mapping; QFLAG as frozen (0)."""
    regca, reecb, repca = cfg["regca"], dict(cfg["reecb"]), cfg["repca"]
    pll_idx, reg_idx, ree_idx, freq_idx = f"PLL2_{bus}", f"REGCP1_{bus}", f"REECB1_{bus}", f"BusFreq_{bus}"
    ss.add("PLL2", {"idx": pll_idx, "bus": bus, "Kp": cfg["pll"]["andes"]["Kp"], "Ki": cfg["pll"]["andes"]["Ki"]})
    ss.add("BusFreq", {"idx": freq_idx, "bus": bus, "Tf": 0.02, "Tw": 0.1})
    ss.add("REGCP1", {"idx": reg_idx, "bus": bus, "gen": gen_idx, "Sn": rating, "pll": pll_idx, "Tg": regca["Tg_s"],
                      "Rrpwr": regca["Rrpwr_pu_per_s"], "Brkpt": regca["Brkpt_pu"], "Zerox": regca["Zerox_pu"],
                      "Lvplsw": regca["Lvplsw"], "Lvpl1": regca["Lvpl1_pu"], "Volim": regca["Volim_pu"],
                      "Lvpnt1": regca["Lvpnt1_pu"], "Lvpnt0": regca["Lvpnt0_pu"], "Iolim": regca["Iolim_pu"],
                      "Tfltr": regca["Tfltr_s"], "Khv": regca["Khv"], "Iqrmax": regca["Iqrmax_pu_per_s"],
                      "Iqrmin": regca["Iqrmin_pu_per_s"], "Accel": regca["Accel"]})
    ss.add("REECB1", {"idx": ree_idx, "reg": reg_idx, "busr": bus, "PFFLAG": reecb["PFFLAG"], "VFLAG": reecb["VFLAG"],
                      "QFLAG": reecb["QFLAG"], "PQFLAG": reecb["PQFLAG"], "Vdip": reecb["Vdip_pu"], "Vup": reecb["Vup_pu"],
                      "Trv": reecb["Trv_s"], "dbd1": reecb["dbd1_pu"], "dbd2": reecb["dbd2_pu"], "Kqv": reecb["Kqv"],
                      "Iqh1": reecb["Iqhl_pu"], "Iql1": reecb["Iqll_pu"], "Vref0": reecb["Vref0_pu"], "Tp": reecb["Tp_s"],
                      "QMax": reecb["QMax_pu"], "QMin": reecb["QMin_pu"], "VMAX": reecb["VMAX_pu"], "VMIN": reecb["VMIN_pu"],
                      "Kqp": reecb["Kqp"], "Kqi": reecb["Kqi"], "Kvp": reecb["Kvp"], "Kvi": reecb["Kvi"], "Tiq": reecb["Tiq_s"],
                      "dPmax": reecb["dPmax_pu_per_s"], "dPmin": reecb["dPmin_pu_per_s"], "PMAX": reecb["PMAX_pu"],
                      "PMIN": reecb["PMIN_pu"], "Imax": reecb["Imax_pu"], "Tpord": reecb["Tpord_s"]})
    ss.add("REPCA1", {"idx": f"REPCA1_{bus}", "ree": ree_idx, "line": _incident_line(ss, bus), "busr": bus, "busf": freq_idx,
                      "VCFlag": repca["VCFlag"], "RefFlag": repca["RefFlag"], "Fflag": repca["Fflag"], "PLflag": repca["PLflag"],
                      "Tfltr": repca["Tfltr_s"], "Kp": repca["Kp"], "Ki": repca["Ki"], "Tft": repca["Tft_s"], "Tfv": repca["Tfv_s"],
                      "Vfrz": repca["Vfrz_pu"], "Rc": repca["Rc_pu"], "Xc": repca["Xc_pu"], "Kc": repca["Kc"],
                      "emax": repca["emax_pu"], "emin": repca["emin_pu"], "dbd1": repca["dbd1_pu"], "dbd2": repca["dbd2_pu"],
                      "Qmax": repca["Qmax_pu"], "Qmin": repca["Qmin_pu"], "Kpg": repca["Kpg"], "Kig": repca["Kig"],
                      "Tp": repca["Tp_s"], "fdbd1": repca["fdbd1_pu"], "fdbd2": repca["fdbd2_pu"], "femax": repca["femax_pu"],
                      "femin": repca["femin_pu"], "Pmax": repca["Pmax_pu"], "Pmin": repca["Pmin_pu"], "Tg": repca["Tg_s"],
                      "Ddn": repca["Ddn"], "Dup": repca["Dup"]})


def enforce_pq(ss, targets):
    """H17 enforce_pq (dicgrid forced-PQ logic): hold the registered P and Q at converted buses."""
    pv = ss.PV
    pv.qlim.enable = False
    pv.qlim.zi[:] = 1.0
    pv.qlim.zl[:] = 0.0
    pv.qlim.zu[:] = 0.0
    for bus, (tp, tq) in targets.items():
        uid = [u for u, b in enumerate(pv.bus.v) if int(b) == int(bus)]
        assert len(uid) == 1
        uid = uid[0]
        assert np.isclose(float(pv.p0.v[uid]), tp, atol=1e-8), (bus, pv.p0.v[uid], tp)
        pv.q0.v[uid] = tq
        pv.qmin.v[uid] = tq
        pv.qmax.v[uid] = tq
        pv.q.v[uid] = tq
        pv.qlim.zi[uid] = 0.0
        pv.qlim.zu[uid] = 1.0


def build(members, k, *, scale=None, mdraw=None):
    """One IEEE-68 case. scale = {branch e: gamma}; mdraw = machine multipliers (h, xd1q1, ka, pss_k)."""
    register()
    ss = new_system()
    net = H["network"]
    for b in net["buses"]:
        ss.add("Bus", {"idx": int(b["idx"]), "name": f"B{int(b['idx'])}", "Vn": 100.0, "v0": float(b["v0"])})
    for e, ln in enumerate(net["lines"]):
        g = float((scale or {}).get(e, 1.0))
        tap = float(ln["tap"]) if float(ln["tap"]) != 0.0 else 1.0
        ss.add("Line", {"idx": f"L{e}", "bus1": int(ln["bus1"]), "bus2": int(ln["bus2"]), "r": float(ln["r"]) / g, "x": float(ln["x"]) / g,
                        "b": float(ln["b"]) * g, "g": float(ln["g"]) * g, "tap": tap, "phi": float(ln["phi"]), "Vn1": 100.0, "Vn2": 100.0, "u": float(ln["u"])})
    for p in net["pv"]:
        ss.add("PV", {"idx": f"G{int(p['bus'])}", "bus": int(p["bus"]), "p0": float(p["p0"]), "v0": float(p["v0"]), "Sn": 100.0,
                      "qmax": 999.0, "qmin": -999.0, "pmax": 999.0, "pmin": -999.0, "Vn": 100.0})
    for s in net["slack"]:
        ss.add("Slack", {"idx": f"G{int(s['bus'])}", "bus": int(s["bus"]), "v0": float(s["v0"]), "a0": float(s.get("a0", 0.0)), "Sn": 100.0,
                         "p0": 0.0, "qmax": 999.0, "qmin": -999.0, "pmax": 999.0, "pmin": -999.0, "Vn": 100.0})
    for i, ld in enumerate(net["loads"]):
        ss.add("PQ", {"idx": f"PQ{i}", "bus": int(ld["bus"]), "p0": float(ld["p0"]), "q0": float(ld["q0"]), "Vn": 100.0})
    cfg = H["wecc_config"]
    for bstr, p in H["sg"].items():
        bus = int(bstr)
        gidx = f"G{bus}"
        if bus in members:
            add_wecc(ss, bus, gidx, H["targets"][bstr]["rating_mva"], cfg)
            continue
        q = {kk: vv for kk, vv in p.items() if kk != "cls"}
        if "KA" in q:
            q["GS"] = float(k) * float((mdraw or {}).get("ka", 1.0))
        if mdraw:
            q["H"] = q["H"] * mdraw.get("h", 1.0)
            q["xd1"] = q["xd1"] * mdraw.get("xd1q1", 1.0)
            q["xq1"] = q["xq1"] * mdraw.get("xd1q1", 1.0)
            if "K" in q:
                q["K"] = q["K"] * mdraw.get("pss_k", 1.0)
        ss.add(p["cls"], {"idx": f"SG_{bus}", "name": f"SG_{bus}", "u": 1, "bus": bus, "gen": gidx, **q})
    ss.setup()
    ss.PQ.config.pq2z = 0
    ss.PQ.config.p2p, ss.PQ.config.q2q = 0.0, 0.0
    ss.PQ.config.p2i, ss.PQ.config.q2i = 0.0, 0.0
    ss.PQ.config.p2z, ss.PQ.config.q2z = 1.0, 1.0
    ss.PV.config.pv2pq = 0
    ss.PFlow.config.tol = 1e-12
    ss.TDS.config.criteria = 0
    if members:
        enforce_pq(ss, {b: (H["targets"][str(b)]["p"], H["targets"][str(b)]["q"]) for b in members})
    ss.PFlow.run()
    if not ss.PFlow.converged:
        return ss, False
    ss.TDS.init()
    return ss, True


def residual(ss) -> float:
    ss.TDS.fg_update(ss.exist.tds)
    return float(max(np.abs(ss.dae.f).max(), np.abs(ss.dae.g).max()))


def _dense(v):
    return np.array(matrix(v))


def decoupled_states(As, tol=1e-12):
    """H17 refinements 1 and 1b, verbatim."""
    keep = list(range(As.shape[0]))
    thr = tol * max(np.abs(As).max(), 1.0)
    removed, isolated_eigs = [], []
    changed = True
    while changed:
        changed = False
        sub = As[np.ix_(keep, keep)]
        for pos in range(len(keep)):
            row, col = sub[pos].copy(), sub[:, pos].copy()
            d = row[pos]
            row[pos] = col[pos] = 0.0
            zero_row, zero_col = np.abs(row).max() <= thr, np.abs(col).max() <= thr
            dead = (abs(d) <= thr) and (zero_row or zero_col)
            isolated = zero_row and zero_col and abs(d) <= 1e-3
            if dead or isolated:
                removed.append(keep[pos])
                if isolated and not dead:
                    isolated_eigs.append(float(d))
                keep.pop(pos)
                changed = True
                break
    return keep, removed, isolated_eigs


def spectrum(ss, want_modes=True):
    dae = ss.dae
    fx, fy, gx, gy = _dense(dae.fx), _dense(dae.fy), _dense(dae.gx), _dense(dae.gy)
    Tf = np.asarray(dae.Tf, float)
    assert np.all(Tf > 0)
    gyinv_gx = np.linalg.solve(gy, gx)
    As = (fx - fy @ gyinv_gx) / Tf[:, None]
    keep, removed, isolated_eigs = decoupled_states(As)
    Ak = As[np.ix_(keep, keep)]
    vals, vecs = np.linalg.eig(Ak)
    order = np.argsort(np.abs(vals), kind="stable")
    small = np.abs(vals[order[:1]])
    rest_idx = order[1:]
    rest = vals[rest_idx]
    centre_ok = bool(small.size == 1 and small.max() < 1e-3 and np.abs(rest).min() >= 1e-2)
    ic = int(np.argmax(rest.real))
    lam = complex(rest[ic])
    if lam.imag < 0:
        lam = lam.conjugate()
    min_abs_re = float(np.abs(rest.real).min())
    resolved = centre_ok and min_abs_re >= 2e-3
    rhp = int((rest.real > 0).sum())
    status = ("UNSTABLE" if rhp else "STABLE") if resolved else "BOUNDARY_OR_UNRESOLVED"
    names = list(dae.x_name)
    others = rest[(np.abs(rest - lam) > 1e-9) & (np.abs(rest - lam.conjugate()) > 1e-9)]
    out = {"status": status, "alpha": float(rest.real.max()), "lam_re": lam.real, "lam_hz": abs(lam.imag) / (2 * np.pi), "rhp": rhp,
           "n_dead": len(removed), "removed_states": [names[i] for i in removed], "isolated_eigs": isolated_eigs, "centre": small.tolist(),
           "centre_ok": centre_ok, "min_abs_re": min_abs_re, "n_x": int(dae.n), "n_y": int(dae.m),
           "gap2": float(lam.real - others.real.max()) if others.size else float("inf")}
    if want_modes:
        a_idx = np.asarray(ss.Bus.a.a)
        v_idx = np.asarray(ss.Bus.v.a)
        a0, v0 = np.asarray(ss.Bus.a.v, float), np.asarray(ss.Bus.v.v, float)
        full = np.zeros((As.shape[0], vecs.shape[1]), complex)
        full[keep, :] = vecs
        f = np.abs(vals.imag) / (2 * np.pi)
        band = [j for j in rest_idx if EM[0] <= f[j] <= EM[1] and vals[j].imag >= -1e-12]
        j_em = max(band, key=lambda j: vals[j].real) if band else None
        modes = []
        for j in rest_idx:
            z = vals[j]
            if z.imag < -1e-12:
                continue
            hz = f[j]
            crit = abs(z - lam) < 1e-9 or abs(z - lam.conjugate()) < 1e-9
            if not ((EM[0] <= hz <= EM[1] and z.real >= -1.0) or crit):
                continue
            y = -gyinv_gx @ full[:, j]
            dv = np.exp(1j * a0) * (y[v_idx] + 1j * v0 * y[a_idx])
            n = np.linalg.norm(dv)
            dv = dv / n if n > 0 else dv
            jj = int(np.argmax(np.abs(dv)))
            dv = dv * np.exp(-1j * np.angle(dv[jj]))
            modes.append({"re": float(z.real), "hz": float(hz), "crit": bool(crit), "em_top": bool(j == j_em),
                          "phi": np.round(np.concatenate([dv.real, dv.imag]), 6).tolist()})
        out["modes"] = modes
        out["em_top_re"] = float(vals[j_em].real) if j_em is not None else float("nan")
        out["em_top_hz"] = float(f[j_em]) if j_em is not None else float("nan")
    return out


EXEMPT = {"REGCP1": ("HVG",), "REPCA1": ("dbd", "eHL", "s2", "fdbd", "feHL", "s5"), "REECB1": ("PIQ", "PIV")}


def limits_active(ss):
    act, allf = [], []
    for mname in ("REGCP1", "REECB1", "REPCA1"):
        mdl = getattr(ss, mname, None)
        if mdl is None or mdl.n == 0:
            continue
        for dname, disc in mdl.discrete.items():
            for flag in ("zl", "zu"):
                arr = getattr(disc, flag, None)
                if arr is not None and np.any(np.asarray(arr) > 0.5):
                    tag = f"{mname}.{dname}.{flag}"
                    allf.append(tag)
                    if not dname.startswith(EXEMPT.get(mname, ())):
                        act.append(tag)
    return act, allf


def run_case(members, k, *, scale=None, mdraw=None, modes=True):
    t0 = time.perf_counter()
    rec = {"S": "+".join(map(str, sorted(members))) or "BASE", "k": k}
    try:
        ss, ok = build(tuple(members), k, scale=scale, mdraw=mdraw)
        if not ok:
            rec.update(status="PF_FAIL")
            return rec
        res = residual(ss)
        lim, lim_all = limits_active(ss) if members else ([], [])
        sp = spectrum(ss, want_modes=modes)
        # forced-PQ check (nominal network only): converted-bus voltages equal the internal power flow
        # (ANDES 2.0 does not populate PV.q.v, so the registered targets are checked through the solution)
        v_err = float("nan")
        if members and not scale:
            bidx = [int(b) for b in ss.Bus.idx.v]
            v_err = max(abs(float(ss.Bus.v.v[bidx.index(int(b))]) - float(H["pf_voltage"][str(b)])) for b in members)
        rec.update(sp, init_residual=res, limits_active=lim, limit_flags_all=lim_all, pq_voltage_error=v_err,
                   vmin=float(np.min(ss.Bus.v.v)), vmax=float(np.max(ss.Bus.v.v)))
        if lim:
            rec["status_linear"] = rec["status"]
            rec["status"] = "LIMIT_ACTIVE"
    except Exception as e:  # noqa: BLE001
        rec.update(status="ERROR", error=f"{type(e).__name__}: {e}"[:300])
    rec["wall_s"] = round(time.perf_counter() - t0, 3)
    return rec


# -------------------------------------------------------------------- phases --
def subsets():
    from itertools import combinations
    v = H["V68"]
    return [tuple(sorted(s)) for r in range(len(v) + 1) for s in combinations(v, r)]


def _jdefault(o):
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def key(task):
    return hashlib.sha1(json.dumps(task, sort_keys=True).encode()).hexdigest()[:16]


def _init():
    for kk in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[kk] = "1"
    H.update(json.loads((B68 / "B68_handoff.json").read_text(encoding="utf-8")))
    register()


def _work(task):
    members = () if task["S"] == "BASE" else tuple(int(b) for b in task["S"].split("+"))
    scale = {int(e): g for e, g in (task.get("scale") or {}).items()} or None
    mdraw = task.get("mdraw")
    rec = run_case(members, task["k"], scale=scale, mdraw=mdraw, modes=task.get("modes", True))
    return task, rec


def run_shard(phase, tasks, n_shards, shard):
    """Single-process worker over every n_shards-th pending task (Windows spawn is unusable in this venv)."""
    out = RAW / phase
    out.mkdir(parents=True, exist_ok=True)
    mine = [t for j, t in enumerate(tasks) if j % n_shards == shard]
    todo = [t for t in mine if not (out / f"{key(t)}.json").exists()]
    print(time.strftime("%Y-%m-%dT%H:%M:%S"), phase, f"shard {shard}/{n_shards}", len(mine), "tasks", len(todo), "to run", flush=True)
    t0 = time.time()
    for i, task in enumerate(todo, 1):
        task, rec = _work(task)
        tmp = out / f"{key(task)}.tmp"
        tmp.write_text(json.dumps({"task": task, "record": rec}, default=_jdefault))
        tmp.replace(out / f"{key(task)}.json")
        if i % 25 == 0:
            print(time.strftime("%Y-%m-%dT%H:%M:%S"), phase, f"shard {shard}", i, "/", len(todo), f"{time.time() - t0:.0f}s", flush=True)
    print(time.strftime("%Y-%m-%dT%H:%M:%S"), phase, f"shard {shard}", "COMPLETE", f"{time.time() - t0:.0f}s", flush=True)


def run_phase(phase, tasks):
    out = RAW / phase
    out.mkdir(parents=True, exist_ok=True)
    todo = [t for t in tasks if not (out / f"{key(t)}.json").exists()]
    print(time.strftime("%Y-%m-%dT%H:%M:%S"), phase, len(tasks), "tasks", len(todo), "to run", flush=True)
    t0 = time.time()
    with Pool(WORKERS, initializer=_init) as pool:
        for i, (task, rec) in enumerate(pool.imap_unordered(_work, todo, chunksize=1), 1):
            (out / f"{key(task)}.json").write_text(json.dumps({"task": task, "record": rec}, default=_jdefault))
            if i % 50 == 0:
                print(time.strftime("%Y-%m-%dT%H:%M:%S"), phase, i, "/", len(todo), f"{time.time() - t0:.0f}s", flush=True)
    print(time.strftime("%Y-%m-%dT%H:%M:%S"), phase, "COMPLETE", f"{time.time() - t0:.0f}s", flush=True)


def tasks_census():
    return [{"phase": "B68C", "pid": p["id"], "k": p["k"], "S": "+".join(map(str, s)) or "BASE"} for p in H["policies_B"] for s in subsets()]


def tasks_draws(ref_pid):
    p = next(q for q in H["policies_B"] if q["id"] == ref_pid)
    return [{"phase": "B68D", "pid": p["id"], "k": p["k"], "S": "+".join(map(str, s)) or "BASE", "draw": d["id"], "mdraw": d["machine"]}
            for d in H["draws_B"] for s in subsets()]


def tasks_branches():
    S = "+".join(map(str, sorted(H["V68"])))
    out = []
    for p in H["policies_B"]:
        out.append({"phase": "B68L", "pid": p["id"], "k": p["k"], "S": S, "scale": {}, "e": -1, "gamma": 1.0, "modes": False})
        for e in H["eligible_branches"]:
            for g in list(H["gammas"]) + [1.001, 0.999]:
                out.append({"phase": "B68L", "pid": p["id"], "k": p["k"], "S": S, "scale": {str(e): g}, "e": e, "gamma": g, "modes": False})
        if p["split"] == "discovery":
            for e in H["r11_branches"]:
                for g in (1.0001, 0.9999):
                    out.append({"phase": "B68L", "pid": p["id"], "k": p["k"], "S": S, "scale": {str(e): g}, "e": e, "gamma": g, "modes": False})
    return out


def main(mode):
    _init()
    if mode == "setup":
        setup()
    elif mode == "ybus":
        ss, ok = build((), 1.0)
        Y = ss.Line.build_ybus() if hasattr(ss.Line, "build_ybus") else ss.Line.build_y()
        Ya = np.array(matrix(Y))
        order = [int(b) for b in ss.Bus.idx.v]
        Yi = np.array(H["ybus_internal_re"]) + 1j * np.array(H["ybus_internal_im"])
        ordi = [int(b["idx"]) for b in H["network"]["buses"]]
        perm = [ordi.index(b) for b in order]
        Yi = Yi[np.ix_(perm, perm)]
        print("ybus max abs diff", float(np.abs(Ya - Yi).max()), "rel", float(np.abs(Ya - Yi).max() / np.abs(Yi).max()), "pf", ok)
    elif mode == "qual":
        out = {"andes_version": andes.__version__, "wecc_config_sha256": H["wecc_config_sha256"], "Q1": {}}
        for p in H["policies_B"]:
            r = run_case((), p["k"])
            r.pop("modes", None)
            out["Q1"][p["id"]] = r
        r = run_case(tuple(H["V68"]), H["policies_B"][0]["k"])
        r.pop("modes", None)
        out["Q2Q4_full_V68"] = r
        # Q3: descriptor QZ of the full Jacobian pencil against the reduced-state spectrum (H17 cross-check)
        from scipy.linalg import eig as _eig
        ss, ok = build(tuple(H["V68"]), H["policies_B"][0]["k"])
        dae = ss.dae
        fx, fy, gx, gy = _dense(dae.fx), _dense(dae.fy), _dense(dae.gx), _dense(dae.gy)
        Tf = np.asarray(dae.Tf, float)
        As = (fx - fy @ np.linalg.solve(gy, gx)) / Tf[:, None]
        red = np.linalg.eigvals(As)
        J = np.block([[fx, fy], [gx, gy]])
        E = np.zeros_like(J)
        E[: dae.n, : dae.n] = np.diag(Tf)
        dv_ = _eig(J, E, right=False)
        dv_ = dv_[np.isfinite(dv_) & (np.abs(dv_) < 1e8)]
        crit = red[np.argsort(-red.real)[:40]]
        out["Q3_descriptor_rel_err"] = float(max(np.min(np.abs(dv_ - z)) / max(1.0, abs(z)) for z in crit))
        (B68 / "B68_qualification.json").write_text(json.dumps(out, indent=1))
        print(json.dumps({k: (v if k != "Q1" else {p: (x["status"], x["alpha"], x.get("init_residual")) for p, x in v.items()}) for k, v in out.items()}, indent=1, default=str)[:5000])
    elif mode == "census":
        run_phase("B68C", tasks_census())
    elif mode == "shard":
        phase, n_sh, sh = sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
        tasks = {"B68C": tasks_census, "B68L": tasks_branches}.get(phase)
        tasks = tasks() if tasks else tasks_draws(sys.argv[5])
        run_shard(phase, tasks, n_sh, sh)
    elif mode == "branches":
        run_phase("B68L", tasks_branches())
    elif mode == "draws":
        run_phase("B68D", tasks_draws(sys.argv[2]))


if __name__ == "__main__":
    assert Path(os.path.expanduser("~")).resolve() == HOME.resolve(), "USERPROFILE must point to the private home"
    main(sys.argv[1])
