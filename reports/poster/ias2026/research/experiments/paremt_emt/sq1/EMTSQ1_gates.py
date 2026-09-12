# ruff: noqa: E501  -- gate tables kept on one line
# STATUS: PREPARED BUT NEVER EXECUTED - SQ1 stopped at SQ1-2 (blind synthetic qualification FAIL); the ParaEMT line is stopped permanently (prereg SQ1).
"""SQ1-4/5/6 device gates G3a -> G3b -> G4a -> G4b (prereg SQ1 sections 5-6).

Requires SQ1-2 and SQ1-3 PASS. Device equivalence = equilibrium + exact parameters/bases + state
trajectories; the frozen V3 estimator is SECONDARY (modal labels MODAL AGREES / MODAL DISAGREES /
MODAL MEASUREMENT UNRESOLVED; DEVICE PASS / MODAL MEASUREMENT UNRESOLVED is a pass).
SG: R-SG30, R-SG36 (seen) and blind B-SG38. GFL: R-GFL30, blind B1-GFL35, B2-GFL37; V2 harnesses and
the V2 voltage interface unchanged; V2 G4b checks 1-5 unchanged. First failing gate stops.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "v2"))
sys.path.insert(0, str(HERE.parent / "v3"))

import _emt  # noqa: E402
import EMTV2_gates as G2  # noqa: E402
import emtv2_rules as R  # noqa: E402
import emtv3_estimator as V  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import tx4_case as C  # noqa: E402
import tx4_emt as K  # noqa: E402

FROZEN_SHA = "849472d8804c694eee69701dd251aaa69d506878007a3a4fc88fde1fb3b0817f"
OUT = _emt.RESULTS / "EMTSQ1"
REFS2 = _emt.RESULTS / "EMTV2" / "refs"
REFS = OUT / "refs"
G2.RAW = _emt.REPO / "external" / "paremt_runs" / "sq1_raw"
CANON = json.loads((REFS / "canonical_params.json").read_text())
META_SQ1 = json.loads((REFS / "refs_manifest.json").read_text())["meta"]
SG_CASES = {"R-SG30": (G2.META["R-SG30"], REFS2 / "R-SG30"), "R-SG36": (G2.META["B-SG36"], REFS2 / "B-SG36"), "B-SG38": (META_SQ1["B-SG38"], REFS / "B-SG38")}
SG_MAP = {"ra": "ra", "xd": "xd", "xq": "xq", "x1": "xd1", "td10": "td10", "tq10": "tq10", "m": "m", "d": "d", "ka": "ka", "ta": "ta", "ks": "pss_gain",
          "t4": "pss_lag", "t5": "pss_washout", "t6": "pss_wash_lag", "w": "w", "flux_blend": "flux_blend", "avr_blend": "avr_blend"}
GFL_MAP = {k: k for k in ("kp_pll", "ki_pll", "tau_p", "kp_p", "ki_p", "kp_q", "ki_q", "kp_i", "ki_i", "xf", "rf", "kp_v", "ki_v", "w")}


def rel_err(a, b):
    return abs(a - b) / max(abs(b), 1e-300) if b != 0 else abs(a)


def sg_row(meta):
    data = C.Data(_emt.RESEARCH)
    p = C.machine_params(data, meta["bus"], (0.03625, 1.425, 1.5, 1.0), {}, {})
    return p


def sg_run(meta, algebraic, dt):
    op = G2.op_of(meta)
    U = G2.U
    U.DT, U.DS = dt, int(round(1e-3 / dt))
    vt, s, net, g, ys = G2.D.quasi_static_g(op) if algebraic else U.common(op)
    p = sg_row(meta)
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
    return t, sg[:, 0, :7], sg[:, 0, 10] + 1j * sg[:, 0, 11], bool(fin), complex(*meta["S"]), row


def modal_label(re_, rf_, tol_a=0.005, tol_f=0.005, rel=0.0):
    if not (re_.get("resolved") and rf_.get("resolved")):
        return "MODAL MEASUREMENT UNRESOLVED", np.nan, np.nan
    da, dfq = abs(re_["alpha"] - rf_["alpha"]), abs(re_["freq"] - rf_["freq"])
    ok = da <= tol_a + rel * abs(rf_["alpha"]) and dfq <= tol_f + rel * rf_["freq"]
    return ("MODAL AGREES" if ok else "MODAL DISAGREES"), da, dfq


def sg_gate(kind, algebraic, name):
    rows, checks = [], []
    for cid, (meta, base) in SG_CASES.items():
        ref = np.load(f"{base}_{kind}.npz")
        t, x, s, fin, s0, row = sg_run(meta, algebraic, 50e-6)
        canon = CANON["sg"][str(meta["bus"])]
        perr = max(rel_err(row[K.SG[k]], canon[c]) for k, c in SG_MAP.items())
        rr = R.v1_rule(t, x, ref["t"], ref["x"], G2.SG_NAMES, 1.0, f"{name}_{cid}")
        rows += rr
        i09 = int(np.argmin(np.abs(t - 0.9)))
        full = fin and len(t) == len(ref["t"])
        re_ = V.unit_ringdown(t, (x[:, 1] - 1.0)[None, :]) if full else {"resolved": False}
        rf_ = V.unit_ringdown(ref["t"], (ref["x"][:, 1] - 1.0)[None, :])
        lab, da, dfq = modal_label(re_, rf_)
        c = {"case": cid, "finite": fin, "dP": float(s[i09].real - s0.real), "dQ": float(s[i09].imag - s0.imag), "param_max_rel_err": float(perr),
             "worst_state_ratio": float(max(r["ratio"] for r in rr)), "pass_states": all(r["pass"] for r in rr), "modal_label": lab,
             "alpha_emt": re_.get("alpha"), "alpha_ref": rf_.get("alpha"), "f_emt": re_.get("freq"), "f_ref": rf_.get("freq"), "d_alpha": da, "d_f": dfq}
        c["pass"] = bool(fin and c["pass_states"] and abs(c["dP"]) <= 1e-4 and abs(c["dQ"]) <= 1e-4 and perr <= 1e-12)
        c["label"] = ("DEVICE PASS / " + lab) if c["pass"] else "DEVICE FAIL"
        checks.append(c)
    return rows, checks


def g3b_reported():
    out = []
    for cid, (meta, base) in SG_CASES.items():
        qs, dyn = np.load(f"{base}_qs.npz"), np.load(f"{base}_dyn.npz")
        t, x, *_ = sg_run(meta, False, 50e-6)
        t25, x25, *_ = sg_run(meta, False, 25e-6)
        for lab, (ta, xa, tb, xb) in {"emt_vs_dynamic_line (transcription error in EMT embedding)": (t, x, dyn["t"], dyn["x"]),
                                      "dynamic_line_vs_quasi_static (electromagnetic line correction)": (dyn["t"], dyn["x"], qs["t"], qs["x"]),
                                      "emt_vs_quasi_static": (t, x, qs["t"], qs["x"]), "emt25_vs_dynamic_line": (t25, x25, dyn["t"], dyn["x"])}.items():
            rr = R.v1_rule(ta, xa, tb, xb, G2.SG_NAMES, 1.0, lab)
            fin = [r for r in rr if np.isfinite(r["ratio"])]
            out.append({"case": cid, "comparison": lab, "max_ratio": float(max(r["ratio"] for r in fin)), "worst_state": max(fin, key=lambda r: r["ratio"])["state"]})
    return out


def g4a():
    rows, checks = G2.g4a()
    for c in checks:
        meta = G2.META[f"{c['case']}_{c['run']}"]
        canon = CANON["gfl"][str(meta["bus"])]
        par = dict(C.GFL_DEFAULTS, w=meta["w"])
        c["param_max_rel_err"] = float(max(max(rel_err(par[k], canon[k]) for k in GFL_MAP), rel_err(C.LEAK, canon["voltage_leak"])))
        c["pass"] = bool(c["pass"] and c["param_max_rel_err"] <= 1e-12)
    return rows, checks


def g4b():
    """V2 G4b checks 1-5 unchanged (as EMTV2_gates.g4b); ringdown with the frozen V3 estimator reported only."""

    labels = list(G2.GFL_NAMES) + ["P", "Q"]
    rows, checks, raw = [], [], {}
    t_ms = np.arange(5001) * 1e-3
    for cid in G2.GFL_IDS:
        big = np.load(_emt.REPO / "external" / "paremt_runs" / "v2_refs" / f"{cid}_dyn_50us.npz")
        for g in G2.GS:
            r0 = G2.gfl_v2(cid, g, None, 50e-6)
            raw[f"{cid}_g{g:g}_noevent_50us"] = G2.save_raw(f"{cid}_g{g:g}_noevent_50us", r0)
            tb, ab = R.nyquist_blocks(r0["nyq"], 50e-6)
            vm = np.abs(r0["vt"])
            dev_x = np.abs(r0["gf"][:, :11] - r0["gf"][0, :11]).max(axis=0)
            c1 = {"max_dVmag": float(np.abs(vm - vm[0]).max()), "max_dP": float(np.abs(r0["gf"][:, 11] - r0["gf"][0, 11]).max()),
                  "max_dQ": float(np.abs(r0["gf"][:, 12] - r0["gf"][0, 12]).max()), "max_dx": float(dev_x.max()), "worst_x": G2.GFL_NAMES[int(np.argmax(dev_x))]}
            c1["pass"] = bool(r0["finite"] and max(c1["max_dVmag"], c1["max_dP"], c1["max_dQ"], c1["max_dx"]) <= 1e-4)
            c2 = R.nyquist_noevent_check(tb, ab) if r0["finite"] else {"pass": False}
            checks.append({"case": cid, "run": f"g{g:g}_noevent", "check": "1_stationarity", **c1})
            checks.append({"case": cid, "run": f"g{g:g}_noevent", "check": "2_nyquist", **c2})
            for case in "abc":
                key = f"g{g:g}_{case}"
                meta = G2.META[f"{cid}_{key}"]
                r50 = G2.gfl_v2(cid, g, case, 50e-6)
                r25 = G2.gfl_v2(cid, g, case, 25e-6)
                raw[f"{cid}_{key}_50us"] = G2.save_raw(f"{cid}_{key}_50us", r50)
                raw[f"{cid}_{key}_25us"] = G2.save_raw(f"{cid}_{key}_25us", r25)
                for dt, rr_ in ((50e-6, r50), (25e-6, r25)):
                    c2 = R.nyquist_event_check(*R.nyquist_blocks(rr_["nyq"], dt)) if rr_["finite"] else {"pass": False}
                    checks.append({"case": cid, "run": key, "check": f"2_nyquist_{int(dt * 1e6)}us", **c2})
                if not (r50["finite"] and r25["finite"]):
                    checks.append({"case": cid, "run": key, "check": "finite", "pass": False})
                    continue
                qref = np.column_stack([big[f"{key}_x"], big[f"{key}_s"].real, big[f"{key}_s"].imag])
                ref_lp = R.on_ms_grid(R.lp(qref, 50e-6), 50e-6, 5.0)
                x0_ref = qref[0]
                e50, e25 = G2.lp_ms(r50, 50e-6), G2.lp_ms(r25, 25e-6)
                m = R.window_mask(t_ms)
                rr3 = []
                for j, lab in enumerate(labels):
                    err = float(np.abs(e50[m, j] - e25[m, j]).max())
                    exc = float(np.abs(ref_lp[m, j] - x0_ref[j]).max())
                    rr3.append({"test": f"{cid}_{key}_50vs25", "quantity": lab, "max_abs_err": err, "excursion": exc, "tolerance": R.TOL_CONV * exc + R.FLOOR,
                                "pass": bool(err <= R.TOL_CONV * exc + R.FLOOR)})
                rows += [{"check": "3_convergence", **r_} for r_ in rr3]
                checks.append({"case": cid, "run": key, "check": "3_convergence", "pass": all(r_["pass"] for r_ in rr3),
                               "worst": max(rr3, key=lambda q: q["max_abs_err"] / q["tolerance"])["quantity"]})
                rr45 = R.lp_compare(t_ms, e50, ref_lp, x0_ref, labels, R.TOL_TRAJ, f"{cid}_{key}")
                rows += [{"check": "4_5_trajectory", **r_} for r_ in rr45]
                t50 = r50["t"]
                m_eq = (t50 >= 0.8 - 1e-12) & (t50 <= 0.9 + 1e-12)
                dp = float(r50["gf"][m_eq, 11].mean() - r50["s"].real)
                dq = float(r50["gf"][m_eq, 12].mean() - r50["s"].imag)
                pq_ok = all(r_["pass"] for r_ in rr45 if r_["quantity"] in ("P", "Q"))
                checks.append({"case": cid, "run": key, "check": "4_PQ", "dP_eq": dp, "dQ_eq": dq, "pass": bool(abs(dp) <= 1e-4 and abs(dq) <= 1e-4 and pq_ok)})
                m_end = (t50 >= 4.9 - 1e-12) & (t50 <= 5.0 + 1e-12)
                pll = float(abs(np.mean(np.angle(np.exp(1j * (r50["gf"][m_end, 0] - np.angle(r50["vt"][m_end])))))))
                qf = float(abs(np.mean(r50["gf"][m_end, 3]) - meta["q_ref"])) if (g == 0.0 and case == "b") else float("nan")
                st_ok = all(r_["pass"] for r_ in rr45 if r_["quantity"] not in ("P", "Q"))
                checks.append({"case": cid, "run": key, "check": "5_states", "pll_err": pll, "g0_qf_err": qf,
                               "worst": max((r_ for r_ in rr45 if r_["quantity"] not in ("P", "Q")), key=lambda q: q["max_abs_err"] / q["tolerance"])["quantity"],
                               "pass": bool(st_ok and pll <= 1e-3 and (np.isnan(qf) or qf <= 1e-3))})
                re_ = V.unit_ringdown(t_ms, e50[:, :11].T, 1.05, 4.8, (0.2, 4.0))
                rf_ = V.unit_ringdown(t_ms, ref_lp[:, :11].T, 1.05, 4.8, (0.2, 4.0))
                lab, da, dfq = modal_label(re_, rf_, rel=0.02)
                checks.append({"case": cid, "run": key, "check": "6_modal_secondary", "modal_label": lab, "d_alpha": da, "d_f": dfq, "pass": True})
                print(cid, key, "done", flush=True)
    return rows, checks, raw


def main() -> int:
    assert hashlib.sha256((HERE.parent / "v3" / "emtv3_estimator.py").read_bytes()).hexdigest() == FROZEN_SHA, "frozen V3 estimator changed"
    s2 = json.loads((OUT / "SQ1-2" / "sq1_synthetic.json").read_text())
    s3 = json.loads((OUT / "SQ1-3" / "sq1_tds.json").read_text())
    assert s2["gates"]["pass"] and s3["pass"], "SQ1-2 and SQ1-3 must pass first"
    summary = {"SQ1-2": "PASS", "SQ1-3": "PASS", "provenance": _emt.provenance()}
    stopped = None
    for gate in ("G3a", "G3b", "G4a", "G4b"):
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
            pd.DataFrame(g3b_reported()).to_csv(d / "G3b_reported.csv", index=False)
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
        summary[f"{gate}_modal_labels"] = sorted({c.get("label", c.get("modal_label")) for c in checks if c.get("label") or c.get("modal_label")})
        summary[f"{gate}_wall_s"] = round(time.time() - t0, 1)
        print(gate, summary[f"GATE_{gate}"], flush=True)
        if not ok:
            stopped = gate
    (OUT / "EMTSQ1_gates.json").write_text(json.dumps(summary, indent=1, default=float), encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k.startswith("GATE")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
