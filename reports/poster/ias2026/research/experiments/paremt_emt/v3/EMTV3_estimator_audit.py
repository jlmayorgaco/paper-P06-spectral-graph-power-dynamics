# ruff: noqa: E501  -- audit tables kept on one line
"""V3-0 audit of the frozen v1 estimator A (experiments/paremt_emt/emt_estimator.py; NOT modified).

Re-executes the arithmetic of emt_estimator.matrix_pencil step by step on the inputs where it was
seen to fail (the four V2 G3a evaluations and the four DV2-2 synthetic cases, seed 7 as in
EMTV2_G3a_diagnostics.py), and records the intermediate quantities that explain the failure:
model order, largest |z|, the dynamic range of the full-record Vandermonde z^k, the numerical rank
that lstsq keeps, the reconstruction residual, and the energy weights exp(2 Re(s) T) of the
selected modes. Writes results/EMTV3/V3-0/estimator_audit.json.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "v2"))

import emt_estimator as E  # noqa: E402
import numpy as np  # noqa: E402

OUT = HERE.parents[2] / "results" / "EMTV3" / "V3-0"


def instrumented(t, chans, dt_in):
    tt, x, _ = E._prep(t, chans)
    fac = int(round(0.1 / dt_in))
    y = E._decimate(x, fac)
    dt = dt_in * fac
    n = y.shape[1]
    L = n // 3
    H = np.vstack([np.array([yc[i : i + L + 1] for i in range(n - L)]) for yc in y])
    _, sv, vh = np.linalg.svd(H, full_matrices=False)
    order = int(np.clip(np.sum(sv / sv[0] >= E.SV_REL), 2, E.MAX_ORDER))
    v = vh[:order].T
    z = np.linalg.eigvals(np.linalg.pinv(v[:-1]) @ v[1:])
    s = np.log(z.astype(complex)) / dt
    out = {"n_samples": int(n), "order": order, "max_abs_z": float(np.abs(z).max()),
           "n_modes_abs_z_gt_1": int(np.sum(np.abs(z) > 1.0)),
           "log10_max_abs_z_power_n": float((n - 1) * np.log10(np.abs(z).max()))}
    with np.errstate(all="ignore"):
        vand = z[None, :] ** np.arange(n)[:, None]
        finite = bool(np.all(np.isfinite(vand)))
        out["vandermonde_finite"] = finite
        if finite:
            colnorm = np.linalg.norm(vand, axis=0)
            out["log10_column_norm_range"] = float(np.log10(colnorm.max() / colnorm.min()))
            try:
                res, _, rank, svals = np.linalg.lstsq(vand, y.T, rcond=None)
                out["lstsq_rank_kept"] = int(rank)
                out["lstsq_rank_full"] = int(order)
                recon = (vand @ res).real.T
                out["resid"] = float(np.linalg.norm(y - recon) / np.linalg.norm(y))
                T = n * dt
                integ = np.where(np.abs(2 * s.real) > 1e-12, (np.exp(2 * s.real * T) - 1) / (2 * s.real), T)
                energy = np.sum(np.abs(res) ** 2, axis=1) * integ
                f = np.abs(s.imag) / (2 * np.pi)
                band = (f >= E.BAND[0]) & (f <= E.BAND[1]) & (s.imag >= 0)
                rows = []
                for i in np.flatnonzero(band):
                    rows.append({"alpha": float(s[i].real), "f": float(f[i]), "abs_residue": float(np.sqrt(np.sum(np.abs(res[i]) ** 2))),
                                 "log10_energy_weight_exp2alphaT": float(2 * s[i].real * T / np.log(10)), "energy": float(energy[i])})
                out["in_band_modes"] = sorted(rows, key=lambda r: -r["energy"])[:6]
            except np.linalg.LinAlgError as exc:
                out["lstsq_error"] = f"LinAlgError: {exc}"
    try:
        mp = E.matrix_pencil(t, chans, dt_in)
        out["frozen_output"] = {k: mp[k] for k in ("resolved", "alpha", "freq", "order", "resid")}
    except Exception as exc:
        out["frozen_output"] = {"error": f"{type(exc).__name__}: {exc}"}
    return out


def main() -> int:
    import EMTV2_gates as G

    OUT.mkdir(parents=True, exist_ok=True)
    src = HERE.parent / "emt_estimator.py"
    res = {"estimator_file": "experiments/paremt_emt/emt_estimator.py", "estimator_sha256": hashlib.sha256(src.read_bytes()).hexdigest(),
           "G3a": {}, "DV2-2": []}
    with np.errstate(all="ignore"):
        for cid in G.SG_IDS:
            ref = np.load(G.REFS / f"{cid}_qs.npz")
            t, x, s, fin, s0 = G.sg_run(cid, True, 50e-6)
            res["G3a"][f"{cid}_emt"] = instrumented(t, (x[:, 1] - 1.0)[None, :], 1e-3)
            res["G3a"][f"{cid}_ref"] = instrumented(ref["t"], (ref["x"][:, 1] - 1.0)[None, :], 1e-3)
        rng = np.random.default_rng(7)
        t = np.arange(0, 30.001, 1e-3)
        for alpha_true, f_true in ((0.127, 0.622), (-0.144, 0.637), (-0.0174, 0.71), (0.0036, 0.705)):
            chans = []
            for _ in range(20):
                a1, a2 = rng.uniform(0.2, 1.0), rng.uniform(0.2, 1.0)
                env = np.where(t >= 1.0, 1.0, 0.0)
                chans.append(env * (a1 * np.exp(alpha_true * (t - 1.0)) * np.cos(2 * np.pi * f_true * (t - 1.0) + rng.uniform(0, 6.28))
                                    + a2 * np.exp(-0.6 * (t - 1.0)) * np.cos(2 * np.pi * 0.95 * (t - 1.0) + rng.uniform(0, 6.28))))
            res["DV2-2"].append({"alpha_true": alpha_true, "f_true": f_true, **instrumented(t, np.array(chans), 1e-3)})
    (OUT / "estimator_audit.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    print(json.dumps(res, indent=1, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
