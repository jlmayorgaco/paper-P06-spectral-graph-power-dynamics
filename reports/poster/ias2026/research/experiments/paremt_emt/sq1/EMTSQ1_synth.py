# ruff: noqa: E501  -- stratum tables kept on one line
"""SQ1-2: fresh blind synthetic qualification of the FROZEN V3 estimator (prereg SQ1 section 3).

Synthetic multichannel records built from the frozen-model modal library (results/EMTSQ1/domain/
modal_library.npz, 26 non-holdout frozen IEEE-39 phasor models): the channel residue patterns and
nuisance-mode structure of a randomly chosen template are kept (with bounded random perturbation);
the template's target mode is replaced by (alpha, f) drawn from the preregistered scientific domain;
affine trends and white noise are added; 90-s records at 1 kHz; the frozen V3 targeted tracker +
sentinel + adaptive 30/60/90 s rule (emtv3_estimator, sha256 849472d8..., NOT modified).
Usage: EMTSQ1_synth.py check  -> generator self-check with dev seed 2222 (no estimator is run)
       EMTSQ1_synth.py blind  -> blind qualification, seed 20260913 (only after the SQ1 prereg commit)
Writes results/EMTSQ1/SQ1-2/sq1_synthetic.{csv,json}.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "v3"))

import emtv3_estimator as V  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

RESEARCH = HERE.parents[2]
LIBF = RESEARCH / "results" / "EMTSQ1" / "domain" / "modal_library.npz"
OUT = RESEARCH / "results" / "EMTSQ1" / "SQ1-2"
FROZEN_SHA = "849472d8804c694eee69701dd251aaa69d506878007a3a4fc88fde1fb3b0817f"
A_DOM, F_DOM = (-0.23, 0.20), (0.55, 1.05)
T_ON, FS, T_REC = 1.2, 1000.0, 90.0
COUNTS = {"SD1": 130, "SD2": 40, "SD3": 40, "ST1": 20, "ST2": 15}
TAU = np.linspace(1.0, 28.8, 2781)


def load_library():
    z = np.load(LIBF)
    cids = sorted({k.split("__")[0] for k in z.files})
    lib = {}
    for cid in cids:
        lam, R, tgt = z[f"{cid}__lam"], z[f"{cid}__R"], int(z[f"{cid}__target"])
        en = np.zeros(R.shape)
        for m in range(lam.size):
            cpx = abs(lam[m].imag) > 1e-9
            sig = (2.0 if cpx else 1.0) * np.abs(R[:, m])[:, None] * np.exp(lam[m].real * TAU)[None, :]
            en[:, m] = np.sum(sig**2, axis=1) / (2.0 if cpx else 1.0)
        share = np.median(en / np.maximum(en.sum(axis=1, keepdims=True), 1e-300), axis=0)
        lib[cid] = (lam, R, tgt, share)
    return cids, lib


def draw(rng, stratum, cids, lib):
    c = {"stratum": stratum}
    if stratum == "ST1":
        c.update(template="single-channel", alpha=rng.uniform(-0.5, -0.3), f=rng.uniform(1.0, 1.6), beta=rng.uniform(0.02, 0.2), b_ratio=rng.uniform(0.2, 1.0),
                 amp=rng.uniform(0.2, 1.0) * rng.choice([-1, 1]), phase=rng.uniform(0, 2 * np.pi), noise_seed=int(rng.integers(0, 2**31 - 1)))
        return c
    cid = str(rng.choice(cids))
    lam, R, tgt, share = lib[cid]
    if stratum == "SD2":
        alpha = rng.uniform(-0.01, 0.01)
    elif stratum == "ST2":
        alpha = rng.uniform(-0.6, -0.3)
    else:
        alpha = rng.uniform(*A_DOM)
    others = [m for m in range(lam.size) if m != tgt and abs(lam[m].imag) > 1e-9 and share[m] >= 0.05]
    for _ in range(200):
        f = rng.uniform(*F_DOM)
        if all(abs(lam[m].imag / (2 * np.pi) - f) >= 0.03 for m in others):
            break
    c.update(template=cid, alpha=float(alpha), f=float(f), delta=float(rng.uniform(-0.03, 0.03)), mag=rng.uniform(0.8, 1.25, R.shape).tolist(),
             ph=rng.uniform(-0.2, 0.2, R.shape).tolist(), nu=float(10 ** rng.uniform(-4, np.log10(3e-3))), trend_c=rng.uniform(-1, 1, R.shape[0]).tolist(),
             trend_d=rng.uniform(-0.01, 0.01, R.shape[0]).tolist(), noise_seed=int(rng.integers(0, 2**31 - 1)))
    if stratum == "SD3":
        k = max(1, int(round(0.3 * R.shape[0])))
        c["low_obs"] = sorted(int(i) for i in rng.choice(R.shape[0], size=k, replace=False))
        c["low_fac"] = rng.uniform(1e-3, 1e-1, k).tolist()
    c["n_ch"] = int(R.shape[0])
    return c


def signal(c, lib):
    t = np.arange(0.0, T_REC - 1e-9, 1.0 / FS)
    on = t >= T_ON
    tt = np.where(on, t - T_ON, 0.0)
    rng = np.random.default_rng(c["noise_seed"])
    if c["stratum"] == "ST1":
        y = on * (c["amp"] * np.exp(c["alpha"] * tt) * np.cos(2 * np.pi * c["f"] * tt + c["phase"]) + c["b_ratio"] * abs(c["amp"]) * np.exp(-c["beta"] * tt))
        y = y + 1e-4 * abs(c["amp"]) * rng.standard_normal(t.size)
        return t, y[None, :]
    lam, R, tgt, share = lib[c["template"]]
    R = R * np.array(c["mag"]) * np.exp(1j * np.array(c["ph"]))
    lam = lam.copy()
    lam[tgt] = complex(c["alpha"], 2 * np.pi * c["f"])
    if "low_obs" in c:
        for i, fac in zip(c["low_obs"], c["low_fac"], strict=True):
            R[i, tgt] *= fac
    y = np.zeros((R.shape[0], t.size))
    for m in range(lam.size):
        term = R[:, m][:, None] * np.exp(lam[m] * tt[on])[None, :]
        y[:, on] += 2.0 * term.real if abs(lam[m].imag) > 1e-9 else term.real
    w = (t >= 2.2) & (t <= 30.0)
    rms = np.sqrt(np.mean(y[:, w] ** 2, axis=1))
    y += (np.array(c["trend_c"]) * rms)[:, None] + (np.array(c["trend_d"]) * rms)[:, None] * t[None, :]
    y += (c["nu"] * rms)[:, None] * rng.standard_normal(y.shape)
    return t, y


def evaluate(c, lib):
    t, y = signal(c, lib)
    if c["stratum"] == "ST1":
        n = int(round(10.0 * FS)) + 1
        r = V.unit_ringdown(t[:n], y[:, :n])
        r["duration"] = 10.0
        r["sentinel"] = {"flag": False}
        return r

    def rec(T):
        n = int(round(T * FS))
        return t[:n], y[:, :n]

    return V.classify_adaptive(rec, c["f"] + c["delta"])


def row_of(c, r):
    v, a, f = r["verdict"], c["alpha"], c["f"]
    return {"stratum": c["stratum"], "template": c["template"], "n_ch": c.get("n_ch", 1), "alpha_true": a, "f_true": f, "f_pred": f + c.get("delta", 0.0),
            "verdict": v, "resolved": bool(r.get("resolved")), "alpha_hat": r.get("alpha"), "f_hat": r.get("freq"), "ci_lo": r.get("ci_alpha", [np.nan, np.nan])[0],
            "ci_hi": r.get("ci_alpha", [np.nan, np.nan])[1], "duration": r.get("duration"), "reasons": ";".join(r.get("reasons", [])),
            "resid_median": r.get("resid_median"), "sentinel_flag": bool(r.get("sentinel", {}).get("flag", False)),
            "wrong_sign": bool((v == "STABLE" and a > 0) or (v == "UNSTABLE" and a < 0)),
            "correct_conclusive": bool((v == "STABLE" and a < 0) or (v == "UNSTABLE" and a > 0)),
            "err_alpha": abs(r["alpha"] - a) if r.get("resolved") else np.nan, "err_f": abs(r["freq"] - f) if r.get("resolved") else np.nan}


def gates(df):
    dom = df[df.stratum.isin(["SD1", "SD2", "SD3"])]
    big = dom[dom.alpha_true.abs() >= 0.01]
    res = dom[dom.resolved]
    g = {"GA_wrong_sign_count": int((dom.wrong_sign & (dom.alpha_true.abs() >= 0.005)).sum()),
         "GB_resolved_fraction_abs_alpha_ge_0.01": float(big.resolved.mean()) if len(big) else np.nan,
         "reported_correct_conclusive_fraction_abs_alpha_ge_0.01": float(big.correct_conclusive.mean()) if len(big) else np.nan,
         "GC_p95_err_alpha": float(np.percentile(res.err_alpha, 95)) if len(res) else np.inf,
         "GD_p95_err_f": float(np.percentile(res.err_f, 95)) if len(res) else np.inf, "n_domain": int(len(dom)), "n_domain_abs_alpha_ge_0.01": int(len(big)),
         "n_domain_resolved": int(len(res))}
    g["pass"] = bool(g["GA_wrong_sign_count"] == 0 and g["GB_resolved_fraction_abs_alpha_ge_0.01"] >= 0.95 and g["GC_p95_err_alpha"] <= 0.005 and g["GD_p95_err_f"] <= 0.005)
    st = df[df.stratum.isin(["ST1", "ST2"])]
    g["stress_reported"] = {s: {"n": int((st.stratum == s).sum()), "resolved": float(st[st.stratum == s].resolved.mean()) if (st.stratum == s).any() else np.nan,
                                "wrong_sign": int(st[st.stratum == s].wrong_sign.sum())} for s in ("ST1", "ST2")}
    return g


def main() -> int:
    mode = sys.argv[1] if len(sys.argv) > 1 else "check"
    est = HERE.parent / "v3" / "emtv3_estimator.py"
    sha = hashlib.sha256(est.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    sha_raw = hashlib.sha256(est.read_bytes()).hexdigest()
    cids, lib = load_library()
    if mode == "check":
        rng = np.random.default_rng(2222)
        for s in ("SD1", "SD2", "SD3", "ST1", "ST2"):
            c = draw(rng, s, cids, lib)
            t, y = signal(c, lib)
            print(s, c["template"], y.shape, bool(np.all(np.isfinite(y))), round(c["alpha"], 4), round(c["f"], 4))
        print("estimator sha256 (raw, LF-normalized):", sha_raw, sha)
        return 0
    assert FROZEN_SHA in (sha, sha_raw), "frozen V3 estimator changed"
    rng = np.random.default_rng(20260913)
    rows, crashed, t0 = [], [], time.time()
    for s, n in COUNTS.items():
        for k in range(n):
            c = draw(rng, s, cids, lib)
            try:
                rows.append(row_of(c, evaluate(c, lib)))
            except Exception as exc:
                crashed.append({"stratum": s, "k": k, "error": f"{type(exc).__name__}: {exc}"})
            print(s, k, flush=True)
    df = pd.DataFrame(rows)
    g = gates(df)
    g["no_crash"] = not crashed
    g["pass"] = bool(g["pass"] and not crashed)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "sq1_synthetic.csv", index=False)
    (OUT / "sq1_synthetic.json").write_text(json.dumps({"seed": 20260913, "estimator_sha256": sha_raw, "gates": g, "crashed": crashed, "wall_s": round(time.time() - t0, 1)},
                                                       indent=1, default=float), encoding="utf-8")
    print(json.dumps(g, indent=1, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
