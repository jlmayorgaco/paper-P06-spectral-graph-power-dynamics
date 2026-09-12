# ruff: noqa: E501  -- gate tables kept on one line
"""Preregistration V2 gates G3a -> G3b -> G4a -> G4b (docs/20260911_PAREMT_EMT_PREREG_V2.md section 4).

Requires gate SW PASS (results/EMTV2/SW/sw_tests.json) and the references of EMTV2_refs.py.
Stop rule: the first failing gate stops the sequence (later gates are recorded NOT RUN).
Writes results/EMTV2/<gate>/*.csv|json, results/EMTV2/EMTV2_gates.json; every-step G4b traces go to
external/paremt_runs/v2_raw/ (hashed). Run with .venv/xtool-paremt.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))

import _emt  # noqa: E402
import EMT02_03_diagnostics as D  # noqa: E402
import EMT02_03_unit_tests as U  # noqa: E402
import emtv2_rules as R  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import tx4_case as C  # noqa: E402
import tx4_emt as K  # noqa: E402
import tx4_emt_v2 as K2  # noqa: E402
from emt_estimator import matrix_pencil  # noqa: E402

OUT = _emt.RESULTS / "EMTV2"
REFS = OUT / "refs"
RAW = _emt.REPO / "external" / "paremt_runs" / "v2_raw"
SG_IDS = ("R-SG30", "B-SG36")
GFL_IDS = ("R-GFL30", "B1-GFL35", "B2-GFL37")
GS = (0.0, 0.03625, 0.25)
SG_NAMES = U.SG_NAMES
GFL_NAMES = K.GFL_STATES
EV = {"a": [K.EV_PREF, 0, 1.0, 1e9, 0.01], "b": [K.EV_EMAG, 0, 1.0, 1e9, -0.01], "c": [K.EV_EPHASE, 0, 1.0, 1e9, 0.02]}
META = json.loads((REFS / "refs_manifest.json").read_text())["meta"]
E_I, E_C = np.zeros(0, np.int64), np.zeros(0, np.complex128)


def op_of(meta):
    vt = complex(*meta["V_T"])
    s = complex(*meta["S"])
    vinf = 1.0 + 1e-4j * np.conj(s / vt)
    return {"V_T": [vt.real, vt.imag], "V_INF": [vinf.real, vinf.imag], "S_dev": [s.real, s.imag], "z_s": [0.0, 1e-4], "z_line": [0.0, 0.05],
            "w": meta["w"]}


# ------------------------------------------------------------------ SG ----------
def sg_run(cid, algebraic, dt):
    meta = META[cid]
    op = op_of(meta)
    U.DT, U.DS = dt, int(round(1e-3 / dt))
    vt, s, net, g, ys = D.quasi_static_g(op) if algebraic else U.common(op)
    data = C.Data(_emt.RESEARCH)
    p = C.machine_params(data, meta["bus"], (0.03625, 1.425, 1.5, 1.0), {}, {})
    w = meta["w"]
    ini = C.sg_init(p, w, vt, s)
    row = np.zeros(len(K.SG_COLS))
    vals = {"bus": 0, "w": w, "ra": p["ra"], "xd": p["xd"], "xq": p["xq"], "x1": p["xd1"], "td10": p["td10"], "tq10": p["tq10"], "m": p["m"],
            "d": p["d"], "ka": p["ka"], "ta": p["ta"], "ks": p["pss_gain"], "t4": p["pss_lag"], "t5": p["pss_washout"], "t6": p["pss_wash_lag"],
            "pm": ini["pm"], "vref": ini["vref"], "flux_blend": 1.0, "avr_blend": 1.0, "has_gov": 0.0, "gov_r": 1.0, "gov_t1": 1.0, "gov_t2": 1.0,
            "gov_t3": 1.0, "gov_dt": 0.0, "gov_pref": ini["pm"]}
    for n, v in vals.items():
        row[K.SG[n]] = v
    ysg = w / complex(p["ra"], p["xd1"])
    K.stamp(g, 2, 0, ysg)
    x0 = np.array([[ini["delta"], 1.0, ini["eq1"], ini["ed1"], ini["efd"], ini["pm"], 0.0, ini["pm"], ini["pm"]]])
    t, sg, gf, v, fin = U.run_kernel(net, g, row[None, :], np.array([ysg]), x0, np.zeros((0, 20)), np.zeros((0, 11)), ys,
                                     np.array([[K.EV_PM, 0, 1.0, 1.2, 0.02]]), 10.0)
    return t, sg[:, 0, :7], sg[:, 0, 10] + 1j * sg[:, 0, 11], bool(fin), complex(*meta["S"])


def sg_gate(ref_kind, algebraic, name):
    rows, checks = [], []
    for cid in SG_IDS:
        ref = np.load(REFS / f"{cid}_{ref_kind}.npz")
        t, x, s, fin, s0 = sg_run(cid, algebraic, 50e-6)
        rr = R.v1_rule(t, x, ref["t"], ref["x"], SG_NAMES, 1.0, f"{name}_{cid}")
        rows += rr
        i09 = int(np.argmin(np.abs(t - 0.9)))
        me = matrix_pencil(t, (x[:, 1] - 1.0)[None, :], 1e-3) if fin and len(t) == len(ref["t"]) else {"alpha": np.nan, "freq": np.nan}
        mp = matrix_pencil(ref["t"], (ref["x"][:, 1] - 1.0)[None, :], 1e-3)
        c = {"case": cid, "finite": fin, "dP": float(s[i09].real - s0.real), "dQ": float(s[i09].imag - s0.imag), "alpha_emt": me["alpha"],
             "f_emt": me["freq"], "alpha_ref": mp["alpha"], "f_ref": mp["freq"], "d_alpha": me["alpha"] - mp["alpha"], "d_f": me["freq"] - mp["freq"]}
        c["pass_states"] = all(r["pass"] for r in rr)
        c["pass"] = bool(fin and c["pass_states"] and abs(c["dP"]) <= 1e-4 and abs(c["dQ"]) <= 1e-4 and abs(c["d_alpha"]) <= 0.005 and abs(c["d_f"]) <= 0.005)
        checks.append(c)
    return rows, checks


def sg_reported():
    """G3b reported quantities: EMT(EMT line) vs quasi-static; dyn vs quasi-static; 25 us."""

    out = []
    for cid in SG_IDS:
        qs = np.load(REFS / f"{cid}_qs.npz")
        dyn = np.load(REFS / f"{cid}_dyn.npz")
        t, x, s, fin, s0 = sg_run(cid, False, 50e-6)
        t25, x25, *_ = sg_run(cid, False, 25e-6)
        for lab, (ta, xa, tb, xb) in {"emt_vs_quasi_static": (t, x, qs["t"], qs["x"]), "dyn_vs_quasi_static": (dyn["t"], dyn["x"], qs["t"], qs["x"]),
                                      "emt25_vs_dyn": (t25, x25, dyn["t"], dyn["x"])}.items():
            rr = R.v1_rule(ta, xa, tb, xb, SG_NAMES, 1.0, lab)
            out.append({"case": cid, "comparison": lab, "max_ratio": float(max(r["ratio"] for r in rr)),
                        "worst_state": max(rr, key=lambda r: r["ratio"])["state"]})
    return out


# ------------------------------------------------------------------ GFL (G4a harness) --
def gfl_v1_alg(cid, g, case):
    meta = META[f"{cid}_g{g:g}_{case}"]
    op = op_of(meta)
    U.DT, U.DS = 50e-6, 20
    vt, s, net, gm, ys = D.quasi_static_g(op)
    par = dict(C.GFL_DEFAULTS)
    x, p_ref, q_ref, v_ref = C.gfl_init(par, meta["w"], vt, s)
    row = np.zeros(len(K.GFL_COLS))
    for n, v in dict(par, bus=0, w=meta["w"], g=g, leak=C.LEAK, p_ref=p_ref, q_ref=q_ref, v_ref=v_ref).items():
        row[K.GF[n]] = v
    t, sg, gf, vv, fin = U.run_kernel(net, gm, np.zeros((0, 27)), np.zeros(0, np.complex128), np.zeros((0, 9)), row[None, :], np.array([x]),
                                      ys, np.array([EV[case]]), 5.0)
    return t, gf[:, 0], vv[:, 0, 0] + 1j * vv[:, 0, 1], bool(fin), s, q_ref


def g4a():
    rows, checks = [], []
    for cid in GFL_IDS:
        ref = np.load(REFS / f"{cid}_qs.npz")
        for g in GS:
            for case in "abc":
                key = f"g{g:g}_{case}"
                t, gf, vt, fin, s, q_ref = gfl_v1_alg(cid, g, case)
                rr = R.v1_rule(t, gf[:, :11], ref[f"{key}_t"], ref[f"{key}_x"], GFL_NAMES, 1.0, f"G4a_{cid}_{key}")
                rows += rr
                done = fin and t[-1] >= 5.0 - 1e-9
                i09 = int(np.argmin(np.abs(t - 0.9)))
                pll = float(abs(np.angle(np.exp(1j * (gf[-1, 0] - np.angle(vt[-1])))))) if done else float("inf")
                qf = (float(abs(gf[-1, 3] - q_ref)) if done else float("inf")) if (g == 0.0 and case == "b") else float("nan")
                c = {"case": cid, "run": key, "finite": done, "dP": float(gf[i09, 11] - s.real), "dQ": float(gf[i09, 12] - s.imag), "pll_err": pll,
                     "g0_qf_err": qf, "pass_states": all(r["pass"] for r in rr)}
                c["pass"] = bool(done and c["pass_states"] and abs(c["dP"]) <= 1e-4 and abs(c["dQ"]) <= 1e-4 and pll <= 1e-3
                                 and (np.isnan(qf) or qf <= 1e-3))
                checks.append(c)
    return rows, checks


# ------------------------------------------------------------------ GFL V2 (G4b) ------
def gfl_v2(cid, g, case, dt, tlen=5.0):
    meta = META[f"{cid}_g{g:g}_{case if case else 'a'}"]
    op = op_of(meta)
    U.DT = dt
    vt, s, net, gm, ys = U.common(op)
    par = dict(C.GFL_DEFAULTS)
    w = meta["w"]
    x, p_ref, q_ref, v_ref = C.gfl_init(par, w, vt, s)
    row = np.zeros(len(K.GFL_COLS))
    for n, v in dict(par, bus=0, w=w, g=g, leak=C.LEAK, p_ref=p_ref, q_ref=q_ref, v_ref=v_ref).items():
        row[K.GF[n]] = v
    filt = np.array(K2.filter_coeffs(w, par["rf"], par["xf"], dt))
    K2.stamp_filter(gm, 2, 0, filt[2])
    ibr, ihis, _ = K2.filter_init(row, vt, complex(x[6], x[7]), x[0], filt)
    events = np.array([EV[case]]) if case else np.zeros((0, 5))
    t0 = time.time()
    t, sg, gf, vv, nyq, fin = K2.simulate_v2(np.linalg.inv(gm), net["coe0"], net["vsol0"], net["brch_ihis"], net["node_ihis"], 2,
                                             np.zeros((0, 27)), np.zeros(0, np.complex128), np.zeros((0, 9)), row[None, :], np.array([x]),
                                             filt[None, :], ibr[None, :], ihis[None, :], E_I, E_C, E_C, 1e-3,
                                             np.array([1], np.int64), np.array([ys], np.complex128), np.array([1.0 + 0j]), events, dt,
                                             int(round(tlen / dt)), 1, np.arange(2, dtype=np.int64))
    return {"t": t, "gf": gf[:, 0], "vt": vv[:, 0, 0] + 1j * vv[:, 0, 1], "nyq": nyq[:, 0], "finite": bool(fin) and t[-1] >= tlen - 1e-9,
            "s": s, "q_ref": q_ref, "wall": time.time() - t0}


def save_raw(name, r):
    RAW.mkdir(parents=True, exist_ok=True)
    p = RAW / f"{name}.npz"
    np.savez_compressed(p, t=r["t"], gf=r["gf"], vt=r["vt"], nyq=r["nyq"])
    return str(p.relative_to(_emt.REPO)).replace("\\", "/"), hashlib.sha256(p.read_bytes()).hexdigest()


def lp_ms(r, dt):
    """LP'd quantities (11 + P + Q) on the 1-ms grid."""

    q = np.column_stack([r["gf"][:, :11], r["gf"][:, 11], r["gf"][:, 12]])
    return R.on_ms_grid(R.lp(q, dt), dt, 5.0)


def g4b():
    labels = list(GFL_NAMES) + ["P", "Q"]
    rows, checks, raw = [], [], {}
    t_ms = np.arange(5001) * 1e-3
    for cid in GFL_IDS:
        big = np.load(_emt.REPO / "external" / "paremt_runs" / "v2_refs" / f"{cid}_dyn_50us.npz")
        for g in GS:
            # no-event run (checks 1, 2)
            r0 = gfl_v2(cid, g, None, 50e-6)
            raw[f"{cid}_g{g:g}_noevent_50us"] = save_raw(f"{cid}_g{g:g}_noevent_50us", r0)
            tb, ab = R.nyquist_blocks(r0["nyq"], 50e-6)
            vm = np.abs(r0["vt"])
            dev_x = np.abs(r0["gf"][:, :11] - r0["gf"][0, :11]).max(axis=0)
            c1 = {"max_dVmag": float(np.abs(vm - vm[0]).max()), "max_dP": float(np.abs(r0["gf"][:, 11] - r0["gf"][0, 11]).max()),
                  "max_dQ": float(np.abs(r0["gf"][:, 12] - r0["gf"][0, 12]).max()), "max_dx": float(dev_x.max()),
                  "worst_x": GFL_NAMES[int(np.argmax(dev_x))]}
            c1["pass"] = bool(r0["finite"] and max(c1["max_dVmag"], c1["max_dP"], c1["max_dQ"], c1["max_dx"]) <= 1e-4)
            c2 = R.nyquist_noevent_check(tb, ab) if r0["finite"] else {"pass": False}
            checks.append({"case": cid, "run": f"g{g:g}_noevent", "check": "1_stationarity", **c1})
            checks.append({"case": cid, "run": f"g{g:g}_noevent", "check": "2_nyquist", **c2})
            for case in "abc":
                key = f"g{g:g}_{case}"
                meta = META[f"{cid}_{key}"]
                r50 = gfl_v2(cid, g, case, 50e-6)
                r25 = gfl_v2(cid, g, case, 25e-6)
                raw[f"{cid}_{key}_50us"] = save_raw(f"{cid}_{key}_50us", r50)
                raw[f"{cid}_{key}_25us"] = save_raw(f"{cid}_{key}_25us", r25)
                for dt, rr_ in ((50e-6, r50), (25e-6, r25)):
                    if rr_["finite"]:
                        tb, ab = R.nyquist_blocks(rr_["nyq"], dt)
                        c2 = R.nyquist_event_check(tb, ab)
                    else:
                        c2 = {"pass": False}
                    checks.append({"case": cid, "run": key, "check": f"2_nyquist_{int(dt * 1e6)}us", **c2})
                if not (r50["finite"] and r25["finite"]):
                    checks.append({"case": cid, "run": key, "check": "finite", "pass": False})
                    continue
                qref = np.column_stack([big[f"{key}_x"], big[f"{key}_s"].real, big[f"{key}_s"].imag])
                ref_lp = R.on_ms_grid(R.lp(qref, 50e-6), 50e-6, 5.0)
                x0_ref = qref[0]
                e50, e25 = lp_ms(r50, 50e-6), lp_ms(r25, 25e-6)
                # check 3: 50 vs 25 us
                rr3 = R.lp_compare(t_ms, e50, e25, x0_ref, labels, R.TOL_CONV, f"{cid}_{key}_50vs25")
                # hmm: tolerance uses the reference excursion
                for r_ in rr3:
                    j = labels.index(r_["quantity"])
                    m = R.window_mask(t_ms)
                    exc = float(np.abs(ref_lp[m, j] - x0_ref[j]).max())
                    r_["excursion"], r_["tolerance"] = exc, R.TOL_CONV * exc + R.FLOOR
                    r_["ratio"] = r_["max_abs_err"] / exc if exc > 0 else float("nan")
                    r_["pass"] = bool(r_["max_abs_err"] <= r_["tolerance"])
                rows += [{"check": "3_convergence", **r_} for r_ in rr3]
                checks.append({"case": cid, "run": key, "check": "3_convergence", "pass": all(r_["pass"] for r_ in rr3),
                               "worst": max(rr3, key=lambda q: q["max_abs_err"] / q["tolerance"])["quantity"]})
                # checks 4 and 5 at 50 us
                rr45 = R.lp_compare(t_ms, e50, ref_lp, x0_ref, labels, R.TOL_TRAJ, f"{cid}_{key}")
                rows += [{"check": "4_5_trajectory", **r_} for r_ in rr45]
                t50 = r50["t"]
                m_eq = (t50 >= 0.8 - 1e-12) & (t50 <= 0.9 + 1e-12)
                s0 = r50["s"]
                dp = float(r50["gf"][m_eq, 11].mean() - s0.real)
                dq = float(r50["gf"][m_eq, 12].mean() - s0.imag)
                pq_ok = all(r_["pass"] for r_ in rr45 if r_["quantity"] in ("P", "Q"))
                checks.append({"case": cid, "run": key, "check": "4_PQ", "dP_eq": dp, "dQ_eq": dq, "pass": bool(abs(dp) <= 1e-4 and abs(dq) <= 1e-4 and pq_ok)})
                m_end = (t50 >= 4.9 - 1e-12) & (t50 <= 5.0 + 1e-12)
                pll = float(abs(np.mean(np.angle(np.exp(1j * (r50["gf"][m_end, 0] - np.angle(r50["vt"][m_end])))))))
                qf = float(abs(np.mean(r50["gf"][m_end, 3]) - meta["q_ref"])) if (g == 0.0 and case == "b") else float("nan")
                st_ok = all(r_["pass"] for r_ in rr45 if r_["quantity"] not in ("P", "Q"))
                checks.append({"case": cid, "run": key, "check": "5_states", "pll_err": pll, "g0_qf_err": qf,
                               "worst": max((r_ for r_ in rr45 if r_["quantity"] not in ("P", "Q")), key=lambda q: q["max_abs_err"] / q["tolerance"])["quantity"],
                               "pass": bool(st_ok and pll <= 1e-3 and (np.isnan(qf) or qf <= 1e-3))})
                # check 6: ringdown
                re_ = R.estimator_r(t_ms, e50[:, :11])
                rf_ = R.estimator_r(t_ms, ref_lp[:, :11])
                c6 = R.ringdown_check(re_, rf_)
                checks.append({"case": cid, "run": key, "check": "6_ringdown", "alpha_emt": re_["alpha"], "f_emt": re_["freq"], "resolved_emt": re_["resolved"],
                               "alpha_ref": rf_["alpha"], "f_ref": rf_["freq"], "resolved_ref": rf_["resolved"], **c6})
                print(cid, key, "done", flush=True)
    return rows, checks, raw


def g4b_reported():
    """Reported: unfiltered EMT vs dyn and vs quasi-static, dyn vs quasi-static (1-ms grid)."""

    out = []
    for cid in GFL_IDS:
        qs = np.load(REFS / f"{cid}_qs.npz")
        dyn = np.load(REFS / f"{cid}_dyn_1ms.npz")
        for g in GS:
            for case in "abc":
                key = f"g{g:g}_{case}"
                d = np.load(RAW / f"{cid}_{key}_50us.npz")
                x_ms = d["gf"][::20, :11]
                t_ms = d["t"][::20]
                for lab, (ta, xa, tb, xb) in {"emt_raw_vs_dyn": (t_ms, x_ms, dyn[f"{key}_t"], dyn[f"{key}_x"]),
                                              "emt_raw_vs_quasi_static": (t_ms, x_ms, qs[f"{key}_t"], qs[f"{key}_x"]),
                                              "dyn_vs_quasi_static": (dyn[f"{key}_t"], dyn[f"{key}_x"], qs[f"{key}_t"], qs[f"{key}_x"])}.items():
                    rr = R.v1_rule(ta, xa, tb, xb, GFL_NAMES, 1.0, lab)
                    fin = [r for r in rr if np.isfinite(r["ratio"])]
                    out.append({"case": cid, "run": key, "comparison": lab, "max_ratio": float(max(r["ratio"] for r in fin)) if fin else float("nan"),
                                "worst_state": max(fin, key=lambda r: r["ratio"])["state"] if fin else ""})
    return out


def main() -> int:
    sw = json.loads((OUT / "SW" / "sw_tests.json").read_text())
    assert sw["GATE_SW"] == "PASS", "gate SW must pass first"
    summary = {"GATE_SW": "PASS", "provenance": _emt.provenance()}
    order = ["G3a", "G3b", "G4a", "G4b"]
    stopped = None
    for gate in order:
        if stopped:
            summary[f"GATE_{gate}"] = f"NOT RUN (stopped at {stopped})"
            continue
        d = OUT / gate
        d.mkdir(parents=True, exist_ok=True)
        t0 = time.time()
        if gate == "G3a":
            rows, checks = sg_gate("qs", True, "G3a")
        elif gate == "G3b":
            rows, checks = sg_gate("dyn", False, "G3b")
            pd.DataFrame(sg_reported()).to_csv(d / "G3b_reported.csv", index=False)
        elif gate == "G4a":
            rows, checks = g4a()
        else:
            rows, checks, raw = g4b()
            (d / "G4b_raw_manifest.json").write_text(json.dumps(raw, indent=1), encoding="utf-8")
        pd.DataFrame(rows).to_csv(d / f"{gate}_trajectories.csv", index=False)
        pd.DataFrame(checks).to_csv(d / f"{gate}_checks.csv", index=False)
        ok = bool(all(c["pass"] for c in checks))
        summary[f"GATE_{gate}"] = "PASS" if ok else "FAIL"
        summary[f"{gate}_failed_checks"] = [c for c in checks if not c["pass"]][:20]
        summary[f"{gate}_wall_s"] = round(time.time() - t0, 1)
        print(gate, summary[f"GATE_{gate}"], flush=True)
        if gate == "G4b":
            pd.DataFrame(g4b_reported()).to_csv(d / "G4b_reported.csv", index=False)
        if not ok:
            stopped = gate
    (OUT / "EMTV2_gates.json").write_text(json.dumps(summary, indent=1, default=float), encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k.startswith("GATE")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
