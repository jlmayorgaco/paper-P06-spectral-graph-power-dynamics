# ruff: noqa: E501  -- diagnostic tables kept on one line
"""Diagnostics of the V2 G3a failure (post hoc; NOT gates; G3a stays FAIL, later gates stay NOT RUN).

DV2-1  The frozen matrix_pencil on omega - 1 (the preregistered G3a ringdown sub-check) for EMT and
       reference, with its internal resolved flag and residual; and a fixed-order (2) least-squares
       damped-sinusoid fit of omega - 1 over [1.2, 10] s as an independent read of the physical mode.
DV2-2  Sanity of the frozen estimator (emt_estimator.classify) on synthetic multichannel signals with
       known modes (the IEEE-39 use case of EMT05-EMT18, which v1 and v2 never reached).
Writes results/EMTV2/G3a/G3a_diagnostics.json. Uses only G3a runs and synthetic data (no G3b/G4 run).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))

import EMTV2_gates as G  # noqa: E402
import numpy as np  # noqa: E402
from emt_estimator import classify, hilbert_envelope, matrix_pencil  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402


def mode_fit(t, y, t0=1.2, t1=10.0):
    m = (t >= t0) & (t <= t1)
    tt, yy = t[m] - t0, y[m]
    a0 = np.abs(yy).max()

    def res(p):
        a, al, f, ph, c = p
        return a * np.exp(al * tt) * np.cos(2 * np.pi * f * tt + ph) + c - yy

    best = None
    for f0 in (0.8, 1.0, 1.2):
        for ph0 in (0.0, 1.5, 3.0):
            r = least_squares(res, [a0, -0.2, f0, ph0, 0.0], x_scale="jac", max_nfev=5000)
            if best is None or r.cost < best.cost:
                best = r
    a, al, f, ph, c = best.x
    rel = float(np.linalg.norm(best.fun) / np.linalg.norm(yy))
    return {"alpha": float(al), "freq": float(abs(f)), "rel_residual": rel}


def main() -> int:
    out = {"DV2-1": {}, "DV2-2": []}
    for cid in G.SG_IDS:
        ref = np.load(G.REFS / f"{cid}_qs.npz")
        t, x, s, fin, s0 = G.sg_run(cid, True, 50e-6)
        for lab, (tt, w) in {"emt": (t, x[:, 1] - 1.0), "ref": (ref["t"], ref["x"][:, 1] - 1.0)}.items():
            mp = matrix_pencil(tt, w[None, :], 1e-3)
            out["DV2-1"][f"{cid}_{lab}"] = {"frozen_matrix_pencil": {k: mp[k] for k in ("resolved", "alpha", "freq", "order", "resid")},
                                             "fixed_order_fit": mode_fit(tt, w)}
    rng = np.random.default_rng(7)
    t = np.arange(0, 30.001, 1e-3)
    for alpha_true, f_true in ((0.127, 0.622), (-0.144, 0.637), (-0.0174, 0.71), (0.0036, 0.705)):
        chans = []
        for _ in range(20):
            a1, a2 = rng.uniform(0.2, 1.0), rng.uniform(0.2, 1.0)
            env = np.where(t >= 1.0, 1.0, 0.0)
            y = env * (a1 * np.exp(alpha_true * (t - 1.0)) * np.cos(2 * np.pi * f_true * (t - 1.0) + rng.uniform(0, 6.28))
                       + a2 * np.exp(-0.6 * (t - 1.0)) * np.cos(2 * np.pi * 0.95 * (t - 1.0) + rng.uniform(0, 6.28)))
            chans.append(y)
        row = {"alpha_true": alpha_true, "f_true": f_true}
        for lab, fn in (("A_matrix_pencil", matrix_pencil), ("B_hilbert", hilbert_envelope)):
            try:
                r = fn(t, np.array(chans), 1e-3)
                row[lab] = {k: r[k] for k in ("resolved", "alpha", "freq") if k in r}
            except Exception as exc:  # the frozen estimator may fail numerically; record, do not repair
                row[lab] = {"error": f"{type(exc).__name__}: {exc}"}
        try:
            row["classify"] = classify(t, np.array(chans), 1e-3)["verdict"]
        except Exception as exc:
            row["classify"] = f"ERROR {type(exc).__name__}"
        out["DV2-2"].append(row)
    (G.OUT / "G3a" / "G3a_diagnostics.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print(json.dumps(out, indent=1, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
