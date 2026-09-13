# ruff: noqa: E501
"""H17/H18 step 2 (xtool-andes-gfl ONLY): the ALT-WECC model in ANDES 2.0.0 (prereg H17, H18).

Composition (frozen in docs/CDW_HARDENING_PREREG_V1.md, H17):
- TX4 network (ANDES ieee39_full.xlsx static sheets, identical to the TX4 JSON) and the SG2AX
  transcription of the TX4 machine/AVR/PSS (PCV06, validated against the internal DAE);
- at converted buses the SG is removed and the TX3 WECC library GFL chain is attached with the
  TX3-GFL-0.1 parameters unchanged (PLL2 + BusFreq + REGCP1 + REECB1 + REPCA1, exactly the
  construction of experiments/tx3/E02_converter_validation/build_andes_gfl_cases.py::_add_gfl);
  REGCP1 Sn = the StaticGen Sn of the replaced unit;
- converted-bus equilibrium: TX3 forced-PQ formulation (logic of
  dicgrid.adapters.andes.enforce_gfl_pq_equilibrium, reproduced verbatim below), holding the
  unit's base (all-SG, nominal network) power-flow P and Q;
- constant-power loads; policies enter only through SG2AX KA/TE (exact exported values).

Spectrum: generalized-QZ descriptor solve of the ANDES Jacobians (TX3 readiness approach),
dead-state zero removal, structural pair removal, TX4/PCV06 verdict band (min |Re| >= 2e-3).

    python H17_andes_alt.py qual     -> results/hardening/alt/H17_qualification_andes.json
    python H17_andes_alt.py matrix   -> results/hardening/alt/cases/*.json (resume-safe)
Launch with USERPROFILE/HOME set to .venv/xtool-andes-gfl/home.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
assert importlib.util.find_spec("ibr_cycles") is None, "internal package must not be importable"

import andes  # noqa: E402
from kvxopt import matrix  # noqa: E402
from scipy.linalg import eig  # noqa: E402

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
REPO = PROJECT.parents[4]
VENV = REPO / ".venv" / "xtool-andes-gfl"
HOME = VENV / "home"
PYCODE = HOME / "pycode"
WORK = VENV / "work" / "cdw_h17"
OUT = PROJECT / "results" / "hardening" / "alt"
CASES = OUT / "cases"
KEEP = ("Bus", "PQ", "PV", "Slack", "Shunt", "Line", "Area")
EM = (0.1, 2.0)
INPUTS: dict = {}


# ------------------------------------------------------------------ networks --
def static_variant(name: str, spec: dict) -> Path:
    WORK.mkdir(parents=True, exist_ok=True)
    target = WORK / f"static_{name}.xlsx"
    if target.exists():
        return target
    src = Path(andes.__file__).parent / "cases" / "ieee39" / "ieee39_full.xlsx"
    sheets = pd.read_excel(src, sheet_name=None)
    ln = sheets["Line"].copy()
    for e, g in spec["scale"].items():
        row = int(e)
        for col in ("r", "x"):
            ln.loc[row, col] = float(ln.loc[row, col]) / g
        for col in ("b", "g", "b1", "b2", "g1", "g2"):
            if col in ln.columns:
                ln.loc[row, col] = float(ln.loc[row, col]) * g
    for e in spec["outage"]:
        ln.loc[int(e), "u"] = 0
    sheets["Line"] = ln
    with pd.ExcelWriter(target) as writer:
        for nm in KEEP:
            if nm in sheets and len(sheets[nm]):
                sheets[nm].to_excel(writer, sheet_name=nm, index=False)
    return target


# ------------------------------------------------------------------- devices --
def _incident_line(ss, bus):
    for idx, b1, b2 in zip(ss.Line.idx.v, ss.Line.bus1.v, ss.Line.bus2.v, strict=True):
        if int(b1) == bus or int(b2) == bus:
            return idx
    raise ValueError(bus)


def add_wecc(ss, bus, gen_idx, rating, cfg):
    """TX3 _add_gfl, verbatim parameter mapping (build_andes_gfl_cases.py)."""

    regca, reecb, repca = cfg["regca"], cfg["reecb"], cfg["repca"]
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
    """dicgrid.adapters.andes.enforce_gfl_pq_equilibrium (verbatim logic): hold registered P, Q."""

    pv = ss.PV
    pv.qlim.enable = False
    pv.qlim.zi[:] = 1.0
    pv.qlim.zl[:] = 0.0
    pv.qlim.zu[:] = 0.0
    for bus, (tp, tq) in targets.items():
        uid = [u for u, b in enumerate(pv.bus.v) if int(b) == int(bus)]
        assert len(uid) == 1
        uid = uid[0]
        assert np.isclose(float(pv.p0.v[uid]), tp, atol=1e-10), (bus, pv.p0.v[uid], tp)
        pv.q0.v[uid] = tq
        pv.qmin.v[uid] = tq
        pv.qmax.v[uid] = tq
        pv.q.v[uid] = tq
        pv.qlim.zi[uid] = 0.0
        pv.qlim.zu[uid] = 1.0


def build(pid, members, net_name, q_targets=None):
    andes.config_logger(stream_level=50)
    static = static_variant(net_name, INPUTS["networks"][net_name])
    ss = andes.load(str(static), setup=False, no_output=True, default_config=True, pycode_path=str(PYCODE))
    gen_of = {int(b): (i, float(sn)) for b, i, sn in zip(ss.PV.bus.v, ss.PV.idx.v, ss.PV.Sn.v, strict=True)}
    gen_of.update({int(b): (i, float(sn)) for b, i, sn in zip(ss.Slack.bus.v, ss.Slack.idx.v, ss.Slack.Sn.v, strict=True)})
    cfg = INPUTS["wecc_config"]
    for bus in sorted(gen_of):
        gidx, rating = gen_of[bus]
        if bus in members:
            add_wecc(ss, bus, gidx, rating, cfg)
        else:
            p = INPUTS["sg"][pid][str(bus)]
            ss.add("SG2AX", {"idx": f"SG_{bus}", "name": f"SG_{bus}", "u": 1, "bus": bus, "gen": gidx, "w": p["w"], "ra": p["ra"],
                             "xd": p["xd"], "xq": p["xq"], "xd1": p["xd1"], "xq1": p["xq1"], "Td10": p["Td10"], "Tq10": p["Tq10"],
                             "M": p["M"], "D": p["D"], "KA": p["KA"], "TE": p["TE"], "KS": p["KS"], "T4": p["T4"], "T5": p["T5"],
                             "T6": p["T6"]})
    ss.setup()
    ss.PQ.config.pq2z = 0
    ss.PQ.config.p2p, ss.PQ.config.q2q = 1.0, 1.0
    ss.PQ.config.p2z, ss.PQ.config.q2z = 0.0, 0.0
    ss.PQ.config.p2i, ss.PQ.config.q2i = 0.0, 0.0
    ss.PV.config.pv2pq = 0
    ss.PFlow.config.tol = 1e-12
    ss.TDS.config.criteria = 0
    if members:
        enforce_pq(ss, {b: q_targets[b] for b in members})
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


def spectrum(ss, want_modes=True, cross_check=False):
    dae = ss.dae
    fx, fy, gx, gy = _dense(dae.fx), _dense(dae.fy), _dense(dae.gx), _dense(dae.gy)
    J = np.block([[fx, fy], [gx, gy]])
    E = np.zeros_like(J)
    Tf = np.asarray(dae.Tf, float)
    E[: dae.n, : dae.n] = np.diag(Tf)
    vals, vecs = eig(J, E, right=True)
    fin = np.where(np.isfinite(vals) & (np.abs(vals) < 1e8))[0]
    ev, V = vals[fin], vecs[:, fin]
    dead_rows = [i for i in range(dae.n) if not np.any(fx[i]) and not np.any(fy[i])]
    n_dead = len(dead_rows)
    order = np.argsort(np.abs(ev), kind="stable")
    drop = list(order[: n_dead + 2])
    small = np.abs(ev[order[n_dead: n_dead + 2]])
    rest_idx = np.array([k for k in range(ev.size) if k not in drop])
    rest = ev[rest_idx]
    dead_ok = bool(n_dead == 0 or np.abs(ev[order[:n_dead]]).max() < 1e-8)
    pair_ok = bool(small.size == 2 and small.max() < 1e-3 and np.abs(rest).min() >= 1e-2)
    ic = int(np.argmax(rest.real))
    lam = complex(rest[ic])
    if lam.imag < 0:
        lam = lam.conjugate()
    min_abs_re = float(np.abs(rest.real).min())
    resolved = pair_ok and dead_ok and min_abs_re >= 2e-3
    rhp = int((rest.real > 0).sum())
    status = ("UNSTABLE" if rhp else "STABLE") if resolved else "BOUNDARY_OR_UNRESOLVED"
    out = {"status": status, "alpha": float(rest.real.max()), "lam_re": lam.real, "lam_hz": abs(lam.imag) / (2 * np.pi), "rhp": rhp,
           "n_dead": n_dead, "dead_ok": dead_ok, "pair": small.tolist(), "pair_ok": pair_ok, "min_abs_re": min_abs_re,
           "n_x": int(dae.n), "n_y": int(dae.m)}
    if cross_check:
        ss.EIG.run()
        mu = np.asarray(ss.EIG.mu, complex)
        crit = rest[np.argsort(-rest.real)[:12]]
        out["eig_cross_rel_err"] = float(max(np.min(np.abs(mu - z)) / max(1.0, abs(z)) for z in crit))
    if want_modes:
        a_idx = np.asarray(ss.Bus.a.a) + dae.n
        v_idx = np.asarray(ss.Bus.v.a) + dae.n
        a0, v0 = np.asarray(ss.Bus.a.v, float), np.asarray(ss.Bus.v.v, float)
        modes = []
        for k in rest_idx:
            z = ev[k]
            if z.imag < -1e-12:
                continue
            hz = abs(z.imag) / (2 * np.pi)
            crit = abs(z - lam) < 1e-9 or abs(z - lam.conjugate()) < 1e-9
            if not ((EM[0] <= hz <= EM[1] and z.real >= -1.0) or crit):
                continue
            x = vecs[:, fin[k]]
            dv = np.exp(1j * a0) * (x[v_idx] + 1j * v0 * x[a_idx])
            n = np.linalg.norm(dv)
            dv = dv / n if n > 0 else dv
            j = int(np.argmax(np.abs(dv)))
            dv = dv * np.exp(-1j * np.angle(dv[j]))
            modes.append({"re": float(z.real), "hz": float(hz), "crit": bool(crit),
                          "phi": np.round(np.concatenate([dv.real, dv.imag]), 6).tolist()})
        out["modes"] = modes
    return out


def limits_active(ss) -> list:
    """Report limiter flags that are active at the initial point (never adjusted)."""

    act = []
    for mname in ("REGCP1", "REECB1", "REPCA1"):
        mdl = getattr(ss, mname, None)
        if mdl is None or mdl.n == 0:
            continue
        for dname, disc in mdl.discrete.items():
            for flag in ("zl", "zu"):
                arr = getattr(disc, flag, None)
                if arr is not None and np.any(np.asarray(arr) > 0.5):
                    if dname.lower().startswith(("dbd", "db")):
                        continue
                    act.append(f"{mname}.{dname}.{flag}")
    return act


def q_targets_nominal():
    """Base (all-SG, nominal network) PF P and Q of every generator bus (policy-independent)."""

    pid0 = next(iter(INPUTS["policies"]))
    ss, ok = build(pid0, (), "NOMINAL")
    assert ok
    out = {}
    for b, p, q in zip(ss.PV.bus.v, ss.PV.p.v, ss.PV.q.v, strict=True):
        out[int(b)] = (float(p), float(q))
    return out


def run_case(pid, members, net, qt, cross=False):
    t0 = time.perf_counter()
    rec = {"pid": pid, "S": "+".join(map(str, sorted(members))) or "BASE", "net": net}
    try:
        ss, ok = build(pid, tuple(members), net, qt)
        if not ok:
            rec.update(status="PF_FAIL")
            return rec
        res = residual(ss)
        lim = limits_active(ss) if members else []
        sp = spectrum(ss, want_modes=True, cross_check=cross)
        pq_err = 0.0
        for b in members:
            uid = [u for u, bb in enumerate(ss.PV.bus.v) if int(bb) == int(b)][0]
            pq_err = max(pq_err, abs(float(ss.PV.q.v[uid]) - qt[b][1]), abs(float(ss.PV.p.v[uid]) - qt[b][0]))
        rec.update(sp, init_residual=res, limits_active=lim, pq_error=pq_err,
                   vmin=float(np.min(ss.Bus.v.v)), vmax=float(np.max(ss.Bus.v.v)))
        if lim:
            rec["status_linear"] = rec["status"]
            rec["status"] = "LIMIT_ACTIVE"
    except Exception as e:  # noqa: BLE001
        rec.update(status="ERROR", error=f"{type(e).__name__}: {e}"[:300])
    rec["wall_s"] = round(time.perf_counter() - t0, 3)
    return rec


def key_of(pid, S, net):
    return hashlib.sha1(f"{pid}|{S}|{net}".encode()).hexdigest()[:16]


def main(mode):
    CASES.mkdir(parents=True, exist_ok=True)
    qt = q_targets_nominal()
    (OUT / "H17_q_targets.json").write_text(json.dumps({str(k): v for k, v in qt.items()}, indent=1))
    pols = list(INPUTS["policies"])
    if mode == "qual":
        out = {"andes_version": andes.__version__, "wecc_config_sha256": INPUTS["wecc_config_sha256"], "Q1": {}}
        for p in pols:
            out["Q1"][p] = run_case(p, (), "NOMINAL", qt)
            out["Q1"][p].pop("modes", None)
        out["Q3_P4like"] = run_case(pols[0], (30, 33, 35, 37), "NOMINAL", qt, cross=True)
        out["Q3_P4like"].pop("modes", None)
        (OUT / "H17_qualification_andes.json").write_text(json.dumps(out, indent=1))
        print(json.dumps(out, indent=1)[:4000])
        return
    jobs = []
    for p in pols:
        for lab in INPUTS["subsets"]:
            mem = () if lab == "BASE" else tuple(int(b) for b in lab.split("+"))
            jobs.append((p, mem, "NOMINAL"))
        for net in INPUTS["networks"]:
            if net != "NOMINAL":
                jobs.append((p, (30, 33, 35, 37), net))
    n_new = 0
    for p, mem, net in jobs:
        lab = "+".join(map(str, sorted(mem))) or "BASE"
        f = CASES / f"{key_of(p, lab, net)}.json"
        if f.exists():
            continue
        rec = run_case(p, mem, net, qt)
        f.write_text(json.dumps(rec))
        n_new += 1
        if n_new % 20 == 0:
            print(time.strftime("%H:%M:%S"), n_new, "new cases", flush=True)
    print("done; new", n_new, "total", len(jobs))


if __name__ == "__main__":
    assert Path(os.path.expanduser("~")).resolve() == HOME.resolve(), "USERPROFILE must point to the private home"
    INPUTS.update(json.loads((OUT / "H17_handoff.json").read_text(encoding="utf-8")))
    assert INPUTS["wecc_config_sha256"].startswith("f07a6a40")
    main(sys.argv[1])
