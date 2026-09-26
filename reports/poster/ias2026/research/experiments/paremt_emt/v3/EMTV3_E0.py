# ruff: noqa: E501  -- stratum tables kept on one line
"""Gate E0, synthetic part (prereg V3 section 4.1) and DV2-2 regression.

Usage: EMTV3_E0.py dev    -> development set, seed 1111, 5 cases per stratum (never qualification)
       EMTV3_E0.py blind  -> blind qualification set, seed 20260912, 140 cases (only after the
                             implementation commit); writes results/EMTV3/E0/E0_synthetic.{csv,json}
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import emtv3_estimator as V  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

OUT = HERE.parents[2] / "results" / "EMTV3" / "E0"
T_ON = 1.2
FS = 1000.0
COUNTS = {"S1": 40, "S2": 25, "S3": 20, "S4": 20, "S5": 20, "S6": 15}


def draw_case(rng, stratum):
    c = {"stratum": stratum}
    if stratum == "S6":
        c.update(n_ch=1, alpha=rng.uniform(-0.5, -0.05), f=rng.uniform(1.0, 1.6), beta=rng.uniform(0.02, 0.2), b_ratio=rng.uniform(0.2, 1.0),
                 amp=rng.uniform(0.2, 1.0) * rng.choice([-1, 1]), phase=rng.uniform(0, 2 * np.pi))
        return c
    c["n_ch"] = int(rng.integers(10, 41))
    if stratum == "S3":
        c["alpha"] = rng.uniform(-0.010, 0.010)
    else:
        c["alpha"] = rng.uniform(-0.30, 0.20)
    c["f"] = rng.uniform(0.25, 1.15)
    c["delta"] = rng.uniform(-0.03, 0.03)
    c["amps"] = (rng.uniform(0.2, 1.0, c["n_ch"]) * rng.choice([-1, 1], c["n_ch"])).tolist()
    c["phases"] = rng.uniform(0, 2 * np.pi, c["n_ch"]).tolist()
    if stratum == "S2":
        while True:
            f2 = rng.uniform(0.2, 1.2)
            if abs(f2 - c["f"]) >= 0.20:
                break
        c.update(f2=f2, alpha2=c["alpha"] - rng.uniform(0.15, 0.60), ratio=rng.uniform(0.3, 1.5))
    if stratum == "S5":
        c.update(f2=c["f"] + rng.choice([-1, 1]) * rng.uniform(0.04, 0.12), alpha2=c["alpha"] - rng.uniform(0.05, 0.30), ratio=rng.uniform(0.5, 1.5))
    if stratum in ("S2", "S5"):
        c["amps2"] = (rng.uniform(0.2, 1.0, c["n_ch"]) * rng.choice([-1, 1], c["n_ch"])).tolist()
        c["phases2"] = rng.uniform(0, 2 * np.pi, c["n_ch"]).tolist()
    if stratum == "S4":
        c["trend_c"] = rng.uniform(-1, 1, c["n_ch"]).tolist()
        c["trend_d"] = rng.uniform(-0.05, 0.05, c["n_ch"]).tolist()
        low = rng.choice(c["n_ch"], size=max(1, int(round(0.2 * c["n_ch"]))), replace=False)
        c["low_obs"] = sorted(int(i) for i in low)
    c["noise_seed"] = int(rng.integers(0, 2**31 - 1))
    return c


def signal(c, t_end):
    t = np.arange(0.0, t_end - 1e-9, 1.0 / FS)
    on = t >= T_ON
    dt = np.where(on, t - T_ON, 0.0)
    rng = np.random.default_rng(c.get("noise_seed", 0))
    if c["stratum"] == "S6":
        y = on * (c["amp"] * np.exp(c["alpha"] * dt) * np.cos(2 * np.pi * c["f"] * dt + c["phase"]) + c["b_ratio"] * abs(c["amp"]) * np.exp(-c["beta"] * dt))
        y = y + 1e-4 * abs(c["amp"]) * rng.standard_normal(t.size)
        return t, y[None, :]
    n = c["n_ch"]
    y = np.zeros((n, t.size))
    for i in range(n):
        a = c["amps"][i]
        if c["stratum"] == "S4" and i in c["low_obs"]:
            a = a * 1e-3
        y[i] = on * a * np.exp(c["alpha"] * dt) * np.cos(2 * np.pi * c["f"] * dt + c["phases"][i])
        if "f2" in c:
            y[i] += on * c["ratio"] * c["amps2"][i] * np.exp(c["alpha2"] * dt) * np.cos(2 * np.pi * c["f2"] * dt + c["phases2"][i])
        sigma = (1e-3 if c["stratum"] == "S4" else 1e-4) * abs(c["amps"][i])
        y[i] += sigma * rng.standard_normal(t.size)
        if c["stratum"] == "S4":
            y[i] += abs(c["amps"][i]) * (c["trend_c"][i] + c["trend_d"][i] * t)
    return t, y


def evaluate(c):
    if c["stratum"] == "S6":
        t, y = signal(c, 10.0 + 1e-3)
        r = V.unit_ringdown(t, y)
        r["duration"] = 10.0
        r["sentinel"] = {"flag": False}
    else:
        full_t, full_y = signal(c, 90.0)

        def rec(T):
            n = int(round(T * FS))
            return full_t[:n], full_y[:, :n]

        r = V.classify_adaptive(rec, c["f"] + c["delta"])
    return r


def dv2_2_cases():
    rng = np.random.default_rng(7)
    t = np.arange(0, 30.001, 1e-3)
    out = []
    for alpha_true, f_true in ((0.127, 0.622), (-0.144, 0.637), (-0.0174, 0.71), (0.0036, 0.705)):
        chans = []
        for _ in range(20):
            a1, a2 = rng.uniform(0.2, 1.0), rng.uniform(0.2, 1.0)
            env = np.where(t >= 1.0, 1.0, 0.0)
            chans.append(env * (a1 * np.exp(alpha_true * (t - 1.0)) * np.cos(2 * np.pi * f_true * (t - 1.0) + rng.uniform(0, 6.28))
                                + a2 * np.exp(-0.6 * (t - 1.0)) * np.cos(2 * np.pi * 0.95 * (t - 1.0) + rng.uniform(0, 6.28))))
        out.append((alpha_true, f_true, t, np.array(chans)))
    return out


def row_of(c, r):
    v = r["verdict"]
    a, f = c["alpha"], c["f"]
    wrong = (v == "STABLE" and a > 0) or (v == "UNSTABLE" and a < 0)
    correct = (v == "STABLE" and a < 0) or (v == "UNSTABLE" and a > 0)
    return {"stratum": c["stratum"], "n_ch": c["n_ch"], "alpha_true": a, "f_true": f, "f_pred": f + c.get("delta", 0.0), "verdict": v,
            "resolved": bool(r.get("resolved")), "alpha_hat": r.get("alpha"), "f_hat": r.get("freq"), "ci_lo": r.get("ci_alpha", [np.nan])[0],
            "ci_hi": r.get("ci_alpha", [np.nan, np.nan])[1], "duration": r.get("duration"), "reasons": ";".join(r.get("reasons", [])),
            "resid_median": r.get("resid_median"), "sentinel_flag": bool(r.get("sentinel", {}).get("flag", False)), "wrong_sign": bool(wrong),
            "correct_conclusive": bool(correct), "err_alpha": abs(r["alpha"] - a) if r.get("resolved") else np.nan,
            "err_f": abs(r["freq"] - f) if r.get("resolved") else np.nan}


def criteria(df, dv):
    big = df.alpha_true.abs() >= 0.005
    a2 = int((df.wrong_sign & big).sum())
    s14 = df[df.stratum.isin(["S1", "S4"]) & (df.alpha_true.abs() >= 0.01)]
    s2 = df[(df.stratum == "S2") & (df.alpha_true.abs() >= 0.01)]
    res = df[df.resolved]
    s6 = df[df.stratum == "S6"]
    p95a = float(np.percentile(res.err_alpha, 95)) if len(res) else np.inf
    p95f = float(np.percentile(res.err_f, 95)) if len(res) else np.inf
    dv_ok = all(d["no_crash"] for d in dv) and all(d["ok"] for d in dv)
    cov = res[(res.ci_lo <= res.alpha_true) & (res.alpha_true <= res.ci_hi)]
    crit = {"A1_no_crash": True, "A2_wrong_sign_count": a2, "A3_S1S4_correct_fraction": float(s14.correct_conclusive.mean()) if len(s14) else np.nan,
            "A3b_S2_correct_fraction": float(s2.correct_conclusive.mean()) if len(s2) else np.nan, "A4_p95_err_alpha": p95a, "A4_p95_err_f": p95f,
            "A5_S6_resolved_fraction": float(s6.resolved.mean()) if len(s6) else np.nan, "A6_dv2_2_pass": dv_ok,
            "reported_ci_coverage": float(len(cov) / max(len(res), 1)), "reported_S3_resolved": float(df[df.stratum == "S3"].resolved.mean()) if (df.stratum == "S3").any() else np.nan,
            "reported_S5_resolved": float(df[df.stratum == "S5"].resolved.mean()) if (df.stratum == "S5").any() else np.nan,
            "reported_sentinel_flags": int(df.sentinel_flag.sum()), "n_cases": int(len(df)), "n_resolved": int(len(res))}
    crit["pass"] = bool(a2 == 0 and crit["A3_S1S4_correct_fraction"] >= 0.95 and crit["A3b_S2_correct_fraction"] >= 0.80 and p95a <= 0.005
                        and p95f <= 0.005 and crit["A5_S6_resolved_fraction"] >= 0.95 and dv_ok)
    return crit


def run_dv2_2():
    out = []
    for alpha_true, f_true, t, y in dv2_2_cases():
        d = {"alpha_true": alpha_true, "f_true": f_true}
        try:
            r = V.classify_adaptive(lambda T, t=t, y=y: (t[: int(round(T * FS))], y[:, : int(round(T * FS))]), f_true)
            d.update(no_crash=True, verdict=r["verdict"], alpha_hat=r["alpha"], f_hat=r["freq"], resolved=r["resolved"], ci=r.get("ci_alpha"))
            if abs(alpha_true) >= 0.01:
                corr = (r["verdict"] == "STABLE") == (alpha_true < 0) and r["verdict"] in ("STABLE", "UNSTABLE")
                d["ok"] = bool(corr and abs(r["alpha"] - alpha_true) <= 0.005 and abs(r["freq"] - f_true) <= 0.005)
            else:
                d["ok"] = not ((r["verdict"] == "STABLE" and alpha_true > 0) or (r["verdict"] == "UNSTABLE" and alpha_true < 0))
        except Exception as exc:
            d.update(no_crash=False, ok=False, error=f"{type(exc).__name__}: {exc}")
        out.append(d)
    return out


def main() -> int:
    mode = sys.argv[1] if len(sys.argv) > 1 else "dev"
    seed, scale = (1111, 5) if mode == "dev" else (20260912, None)
    rng = np.random.default_rng(seed)
    rows, t0 = [], time.time()
    crashed = []
    for stratum, n in COUNTS.items():
        for k in range(scale if scale else n):
            c = draw_case(rng, stratum)
            try:
                rows.append(row_of(c, evaluate(c)))
            except Exception as exc:
                crashed.append({"stratum": stratum, "k": k, "error": f"{type(exc).__name__}: {exc}"})
            print(stratum, k, rows[-1]["verdict"] if rows else "", flush=True)
    df = pd.DataFrame(rows)
    dv = run_dv2_2()
    crit = criteria(df, dv)
    crit["A1_no_crash"] = not crashed
    crit["pass"] = bool(crit["pass"] and not crashed)
    out = {"mode": mode, "seed": seed, "criteria": crit, "crashed": crashed, "dv2_2": dv, "wall_s": round(time.time() - t0, 1)}
    if mode == "blind":
        OUT.mkdir(parents=True, exist_ok=True)
        df.to_csv(OUT / "E0_synthetic.csv", index=False)
        (OUT / "E0_synthetic.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    else:
        dev = HERE.parents[2] / "results" / "EMTV3" / "dev"
        dev.mkdir(parents=True, exist_ok=True)
        df.to_csv(dev / "E0_dev.csv", index=False)
        (dev / "E0_dev.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print(json.dumps({k: v for k, v in out.items() if k != "dv2_2"}, indent=1, default=float))
    print(json.dumps(dv, indent=1, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
