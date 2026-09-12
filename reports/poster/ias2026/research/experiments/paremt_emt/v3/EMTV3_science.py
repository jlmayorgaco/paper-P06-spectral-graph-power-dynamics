# ruff: noqa: E501  -- result tables kept on one line
# STATUS: PREPARED BUT NEVER EXECUTED - V3 stopped at gate E0 (commit 229613de). Any use needs a new preregistration.
"""V3 scientific blocks EMT04 -> EMT05 -> EMT06 -> EMT07/08 -> EMT10 (prereg V3 section 7), then STOP.

Requires E0, G3a, G3b, G4a, G4b all PASS (results/EMTV3/EMTV3_gates.json). IEEE-39 = v1 realization with
every GFL on the V2 voltage interface (v3_case.build_v3); verdicts from the qualified V3 instrument
(targeted tracker at the committed band_hz, sentinel, adaptive 30/60/90 s). Frozen predictions,
portfolios, policies, pulse, dt, thresholds as in prereg v1.
Usage: EMTV3_science.py run | verify   (verify reruns EMT05 and the EMT08 grid for byte identity)
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
import emtv3_estimator as V  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import tx4_case as C  # noqa: E402
import v3_case as B  # noqa: E402

OUT = _emt.RESULTS / "EMTV3"
RAW = _emt.REPO / "external" / "paremt_runs" / "v3_sci"
PRED = pd.read_csv(_emt.RESULTS / "EMT_PRED" / "phasor_predictions.csv").set_index("case_id")
DATA = C.Data(_emt.RESEARCH)
H4 = [30, 33, 35, 37]
P4 = [0.03625, 1.425, 1.5, 1.0]
_CACHE = {}
MANIFESTS = []


def case_key(case, dt, T, amp, tau_m, native):
    return json.dumps([case["members"], case["theta"], case.get("variant", {"kind": "none"}), dt, T, amp, tau_m, native], sort_keys=True)


def run(case, dt, T, amp=0.02, tau_m=1e-3, native=False, tag=""):
    key = case_key(case, dt, T, amp, tau_m, native)
    if key in _CACHE:
        return _CACHE[key]
    t0 = time.time()
    if native:
        b = B.build_native(DATA, case, dt, amplitude=amp, tau_m=tau_m)
        t, sg, gf, v, nyq, fin = B.run_native(b, dt, T, int(round(1e-3 / dt)))
    else:
        b = B.build_v3(DATA, case, dt, amplitude=amp, tau_m=tau_m)
        t, sg, gf, v, nyq, fin = B.run_v3(b, dt, T, int(round(1e-3 / dt)))
    out = {"t": t, "sg": sg, "gf": gf, "v": v[..., 0] + 1j * v[..., 1], "labels_sg": b["labels_sg"], "labels_gf": b["labels_gf"], "order": b["order"],
           "finite": bool(fin) and t[-1] >= T - 1e-9}
    chans, names = _emt.estimator_channels(out)
    out["chans"], out["names"] = chans, names
    rid = f"{case['case_id']}{tag}_dt{int(round(dt * 1e6))}us_T{int(T)}_a{amp:g}_tm{tau_m * 1e3:g}{'_native' if native else ''}"
    RAW.mkdir(parents=True, exist_ok=True)
    p = RAW / f"{rid}.npz"
    np.savez_compressed(p, t=t, chans=chans.astype(np.float64), names=np.array(names), sg=sg.astype(np.float32), gf=gf.astype(np.float32))
    MANIFESTS.append({"run_id": rid, **_emt.provenance(), **{f"{k}_hash": v_ for k, v_ in b["hashes"].items()}, "dt_s": dt, "duration_s": T,
                      "disturbance": {"bus": C.PULSE["bus"], "fraction": amp, "t0": C.PULSE["t0"], "t1": C.PULSE["t1"]}, "tau_m_s": tau_m, "native_stator": native,
                      "seed": None, "wall_s": round(time.time() - t0, 1), "finite": out["finite"], "output_path": str(p.relative_to(_emt.REPO)).replace("\\", "/"),
                      "output_sha256": hashlib.sha256(p.read_bytes()).hexdigest()})
    _CACHE[key] = out
    return out


def classify(case, f_pred, dt=50e-6, amp=0.02, tau_m=1e-3, native=False, tag=""):
    def rec(T):
        o = run(case, dt, T, amp, tau_m, native, tag)
        if not o["finite"]:
            raise FloatingPointError("non-finite EMT run")
        return o["t"], o["chans"]

    try:
        r = V.classify_adaptive(rec, f_pred)
    except FloatingPointError:
        return {"verdict": "NUMERICAL_FAILURE", "resolved": False, "alpha": np.nan, "freq": np.nan, "ci_alpha": [np.nan, np.nan], "duration": np.nan,
                "reasons": ["non-finite"], "sentinel": {"flag": False}, "history": []}
    return r


def row(label, case, r, extra=None):
    cid = case["case_id"]
    a_ph = float(PRED.loc[cid, "band_re"]) if cid in PRED.index else np.nan
    f_ph = float(PRED.loc[cid, "band_hz"]) if cid in PRED.index else np.nan
    return {"label": label, "case_id": cid, "members": "+".join(map(str, case["members"])) or "BASE", "g": case["theta"][0], "verdict": r["verdict"],
            "resolved": bool(r.get("resolved")), "alpha_emt": r.get("alpha"), "f_emt": r.get("freq"), "ci_lo": r.get("ci_alpha", [np.nan, np.nan])[0],
            "ci_hi": r.get("ci_alpha", [np.nan, np.nan])[1], "duration_s": r.get("duration"), "reasons": ";".join(r.get("reasons", [])),
            "sentinel_flag": bool(r.get("sentinel", {}).get("flag", False)), "alpha_phasor": a_ph, "f_phasor": f_ph,
            "phasor_verdict": "UNSTABLE" if a_ph > 0 else "STABLE", **(extra or {})}


def emt04():
    cases = ["EMT05_P4_BASE", "EMT05_P4_30+33+35", "EMT05_P4_30+33+35+37", "EMT04_G_S_H4"]
    rows = []
    for cid in cases:
        case = dict(_emt.CASE[cid])
        f_pred = float(PRED.loc[cid, "band_hz"])
        for dt in (100e-6, 50e-6, 25e-6):
            rows.append(row(f"dt{int(dt * 1e6)}us", case, classify(case, f_pred, dt), {"dt_us": int(dt * 1e6)}))
        if cid in ("EMT05_P4_BASE", "EMT05_P4_30+33+35+37"):
            rows.append(row("tau_m_0.5ms", case, classify(case, f_pred, 50e-6, tau_m=0.5e-3), {"dt_us": 50}))
            rows.append(row("native_stator", case, classify(case, f_pred, 50e-6, native=True), {"dt_us": 50}))
        print("EMT04", cid, [r["verdict"] for r in rows[-5:]], flush=True)
    df = pd.DataFrame(rows)

    def acc(a, b):
        same = a["verdict"] == b["verdict"]
        both = a["resolved"] and b["resolved"]
        da = abs(a["alpha_emt"] - b["alpha_emt"]) if both else np.nan
        dfreq = abs(a["f_emt"] - b["f_emt"]) if both else np.nan
        return {"verdict_invariant": bool(same), "d_alpha": da, "d_f": dfreq, "pass": bool(same and (not both or (da <= 0.01 and dfreq <= 0.01)))}

    g5, hold = [], []
    for cid in cases:
        s = df[df.case_id == cid].set_index("label")
        g5.append({"case_id": cid, **acc(s.loc["dt50us"], s.loc["dt25us"])})
        for lab in ("tau_m_0.5ms", "native_stator"):
            if lab in s.index:
                hold.append({"case_id": cid, "holdout": lab, **acc(s.loc["dt50us"], s.loc[lab])})
    noevent = []
    for cid in cases:
        case = dict(_emt.CASE[cid])
        o = run(case, 50e-6, 5.0, amp=0.0)
        vm = np.abs(o["v"])
        rel = o["chans"][: len([b for b in o["labels_sg"] if b != 39])]
        noevent.append({"case_id": cid, "finite": o["finite"], "max_dVmag": float(np.abs(vm - vm[0]).max()),
                        "max_rel_speed": float(np.abs(rel).max()) if rel.size else 0.0})
    return df, pd.DataFrame(g5), pd.DataFrame(hold), pd.DataFrame(noevent)


def emt05(df04):
    ids = [c["case_id"] for c in _emt.CASES if c["experiment"] == "EMT05"]
    rows = []
    for cid in ids:
        case = dict(_emt.CASE[cid])
        rows.append(row("P4_2pct", case, classify(case, float(PRED.loc[cid, "band_hz"]))))
    df = pd.DataFrame(rows)
    proper = df[df.case_id != "EMT05_P4_30+33+35+37"]
    h4 = df[df.case_id == "EMT05_P4_30+33+35+37"].iloc[0]
    primary = bool((proper.verdict == "STABLE").all() and h4.verdict == "UNSTABLE")
    res = df[df.resolved]
    sec = {"mae_alpha": float((res.alpha_emt - res.alpha_phasor).abs().mean()) if len(res) else np.nan,
           "max_abs_df": float((res.f_emt - res.f_phasor).abs().max()) if len(res) else np.nan}
    sec["pass"] = bool(sec["mae_alpha"] <= 0.03 and sec["max_abs_df"] <= 0.05)
    amp_rows = []
    for cid in ("EMT05_P4_30+33+35", "EMT05_P4_30+33+35+37"):
        case = dict(_emt.CASE[cid])
        base = df[df.case_id == cid].iloc[0]
        for amp in (0.01, 0.04):
            r = classify(case, float(PRED.loc[cid, "band_hz"]), amp=amp)
            ok = bool(r.get("resolved") and base.resolved and abs(r["alpha"] - base.alpha_emt) <= 0.01)
            amp_rows.append({"case_id": cid, "amplitude": amp, "verdict": r["verdict"], "alpha_emt": r.get("alpha"), "alpha_2pct": base.alpha_emt,
                             "d_alpha": abs(r["alpha"] - base.alpha_emt) if r.get("resolved") else np.nan, "pass": ok})
    return df, {"primary_pass": primary, "kappa_emt": 4 if primary else None, "secondary": sec}, pd.DataFrame(amp_rows)


def emt06(df05):
    a = {tuple(sorted(int(x) for x in m.split("+"))) if m != "BASE" else (): r for m, r in zip(df05.members, df05.to_dict("records"), strict=True)}
    full = tuple(H4)
    rows = []
    for i in H4:
        for kind, s_with, s_without in (("alone", (i,), ()), ("last", full, tuple(b for b in full if b != i))):
            rw, ro = a[s_with], a[s_without]
            usable = rw["resolved"] and ro["resolved"]
            d_emt = rw["alpha_emt"] - ro["alpha_emt"] if usable else np.nan
            d_ph = rw["alpha_phasor"] - ro["alpha_phasor"]
            rows.append({"device": i, "context": kind, "d_alpha_emt": d_emt, "d_alpha_phasor": d_ph, "usable": bool(usable),
                         "sign_agrees": bool(usable and np.sign(d_emt) == np.sign(d_ph))})
    eff = pd.DataFrame(rows)
    lattice = []
    import itertools

    for i in H4:
        others = [b for b in H4 if b != i]
        for k in range(4):
            for s in itertools.combinations(others, k):
                rw, ro = a[tuple(sorted(s + (i,)))], a[tuple(sorted(s))]
                if rw["resolved"] and ro["resolved"]:
                    lattice.append(bool(np.sign(rw["alpha_emt"] - ro["alpha_emt"]) == np.sign(rw["alpha_phasor"] - ro["alpha_phasor"])))
    status = "UNRESOLVED" if not eff.usable.all() else ("PASS" if eff.sign_agrees.all() else "FAIL")
    return eff, {"status": status, "lattice_sign_agreement": float(np.mean(lattice)) if lattice else np.nan, "lattice_usable": len(lattice)}


def h4_case(g, cid):
    return {"case_id": cid, "experiment": "EMT08", "members": H4, "theta": [g, 1.425, 1.5, 1.0], "variant": {"kind": "none"}, "note": "bisection"}


def emt08(df05):
    grid = [("EMT05_P4_30+33+35+37", 0.03625)] + [(c["case_id"], c["theta"][0]) for c in _emt.CASES if c["experiment"] == "EMT08"]
    rows = []
    for cid, g in grid:
        rows.append(row("grid", dict(_emt.CASE[cid]), classify(dict(_emt.CASE[cid]), float(PRED.loc[cid, "band_hz"])), {"g": g}))
    df = pd.DataFrame(rows).sort_values("g").reset_index(drop=True)
    pts = df[df.verdict.isin(["STABLE", "UNSTABLE"])].sort_values("g")
    bracket = None
    for (_, lo), (_, hi) in zip(pts.iloc[:-1].iterrows(), pts.iloc[1:].iterrows(), strict=True):
        if lo.verdict == "UNSTABLE" and hi.verdict == "STABLE":
            bracket = [lo.to_dict(), hi.to_dict()]
            break
    bis = []
    if bracket is not None:
        gpts = df[["g", "f_phasor"]].dropna().sort_values("g")
        for k in range(4):
            lo, hi = bracket
            if hi["g"] - lo["g"] <= 0.0025:
                break
            gm = 0.5 * (lo["g"] + hi["g"])
            f_pred = float(np.interp(gm, gpts.g, gpts.f_phasor))
            case = h4_case(gm, f"EMT08_bisect{k + 1}_g{gm:.6f}")
            r = classify(case, f_pred)
            rr = row(f"bisect{k + 1}", case, r, {"g": gm, "f_pred_interpolated": f_pred})
            bis.append(rr)
            if r["verdict"] == "UNSTABLE":
                bracket = [rr, hi]
            elif r["verdict"] == "STABLE":
                bracket = [lo, rr]
            else:
                break
    out = {"bracket": None, "g_star_emt": None, "primary_pass": False, "target_pass": False}
    if bracket is not None:
        lo, hi = bracket
        if lo["resolved"] and hi["resolved"]:
            gs = lo["g"] + (hi["g"] - lo["g"]) * lo["alpha_emt"] / (lo["alpha_emt"] - hi["alpha_emt"])
            out.update(bracket=[lo["g"], hi["g"]], g_star_emt=float(gs), primary_pass=bool(lo["verdict"] == "UNSTABLE" and hi["verdict"] == "STABLE"),
                       target_pass=bool(abs(gs - 0.20768) <= 0.02))
    return pd.concat([df, pd.DataFrame(bis)], ignore_index=True), out


def emt10(df05, df08, g_star):
    seq = [("EMT05_P4_30+33+35+37", 0.03625), ("EMT10_newton1", 0.08806), ("EMT10_gstar", 0.20768), ("EMT08_g0.250", 0.25)]
    rows = []
    for cid, g in seq:
        case = dict(_emt.CASE[cid])
        rows.append(row("newton_sequence", case, classify(case, float(PRED.loc[cid, "band_hz"])), {"g": g}))
    df = pd.DataFrame(rows)
    mono = bool(df.resolved.all() and (np.diff(df.alpha_emt.to_numpy()) < 0).all())
    last_stable = bool(df.iloc[-1].verdict == "STABLE")
    gs_ok = bool(g_star is not None and abs(g_star - 0.20768) <= 0.02)
    return df, {"monotone_decrease": mono, "g025_stable": last_stable, "g_star_within_0.02": gs_ok, "pass": bool(mono and last_stable and gs_ok)}


def main() -> int:
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    gates = json.loads((OUT / "EMTV3_gates.json").read_text())
    assert all(gates.get(f"GATE_{g}") == "PASS" for g in ("E0", "G3a", "G3b", "G4a", "G4b")), "all V3 model gates must pass first"
    if mode == "verify":
        d = OUT / "verify"
        d.mkdir(parents=True, exist_ok=True)
        df05, s05, amp = emt05(None)
        df05.to_csv(d / "p4_16_subsets.csv", index=False)
        df08, s08 = emt08(df05)
        df08[df08.label == "grid"].to_csv(d / "g_grid.csv", index=False)
        same = {f: (d / f).read_bytes() == (OUT / ref / f2).read_bytes() for f, ref, f2 in (("p4_16_subsets.csv", "EMT05", "p4_16_subsets.csv"), ("g_grid.csv", "EMT08", "g_grid.csv"))}
        (OUT / "EMTV3_science_rerun_identity.json").write_text(json.dumps(same, indent=1), encoding="utf-8")
        print(same)
        return 0
    summary = {"provenance": _emt.provenance()}
    t0 = time.time()
    df04, g5, hold, noev = emt04()
    d = OUT / "EMT04"
    d.mkdir(parents=True, exist_ok=True)
    df04.to_csv(d / "emt04_runs.csv", index=False)
    g5.to_csv(d / "emt04_g5.csv", index=False)
    hold.to_csv(d / "emt04_holdouts.csv", index=False)
    noev.to_csv(d / "emt04_noevent.csv", index=False)
    summary["EMT04"] = {"G5_pass": bool(g5["pass"].all()), "G2_tau_m_pass": bool(hold[hold.holdout == "tau_m_0.5ms"]["pass"].all()),
                        "native_stator_descriptive": hold[hold.holdout == "native_stator"].to_dict("records")}
    if not summary["EMT04"]["G5_pass"]:
        summary["STOPPED"] = "G5 failed (prereg v1 section 6: stop and amend)"
    else:
        df05, s05, amp = emt05(df04)
        d = OUT / "EMT05"
        d.mkdir(parents=True, exist_ok=True)
        df05.to_csv(d / "p4_16_subsets.csv", index=False)
        amp.to_csv(d / "amplitude_checks.csv", index=False)
        s05["amplitude_checks_pass"] = bool(amp["pass"].all())
        summary["EMT05"] = s05
        eff, s06 = emt06(df05)
        d = OUT / "EMT06"
        d.mkdir(parents=True, exist_ok=True)
        eff.to_csv(d / "contextual_effects.csv", index=False)
        summary["EMT06"] = s06
        df08, s08 = emt08(df05)
        d = OUT / "EMT08"
        d.mkdir(parents=True, exist_ok=True)
        df08.to_csv(d / "g_boundary.csv", index=False)
        df08[df08.label == "grid"].to_csv(d / "g_grid.csv", index=False)
        summary["EMT08"] = s08
        df10, s10 = emt10(df05, df08, s08.get("g_star_emt"))
        d = OUT / "EMT10"
        d.mkdir(parents=True, exist_ok=True)
        df10.to_csv(d / "newton_sequence.csv", index=False)
        summary["EMT10"] = s10
        summary["STOPPED"] = "after EMT10 by prereg V3 section 1.5 (report before EMT09, EMT11-EMT17)"
    summary["wall_s"] = round(time.time() - t0, 1)
    (OUT / "EMTV3_science.json").write_text(json.dumps(summary, indent=1, default=float), encoding="utf-8")
    with (OUT / "EMTV3_science_manifests.jsonl").open("w", encoding="utf-8") as f:
        for m in MANIFESTS:
            f.write(json.dumps(m, default=float) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "provenance"}, indent=1, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
