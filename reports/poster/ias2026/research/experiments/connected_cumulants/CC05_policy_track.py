# ruff: noqa: E501  -- table labels kept on one line
"""CC05 - Phase 5: connected cumulants along the frozen F7B witness-contraction line.

Preregistered (connected_cumulants_prereg_v1.yaml, de09db3b). All 418 frozen points of
the line t = 0.8515625, k = h = 1 (labels and kappa read from F7B_points.csv.gz).
Evaluation at the label-independent critical eigenvalue lambda_V(g) of the full core.

H5: witness contraction corresponds to a change of the dominant connected support
    D(g) = argmax_{B subseteq V, |B| >= 2} delta_B(g),  delta_B = |chi_B F(V minus B)| / |dF(V)/ds|.
    SUPPORTED iff |D| = kappa at >= 80 % and D in H(g) at >= 80 % of non-empty points.
"""

from __future__ import annotations

import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from _cc import (
    CORE,
    F7B_POINTS,
    PATH_T,
    Realization,
    band_critical,
    label,
    out_dir,
    write_json,
)

from _overnight import pin_blas_threads
from ibr_cycles.cycles.connected import (
    characteristic_values,
    connected_share,
    cumulants,
    subsets,
)

OUT = out_dir("CC05")
H_S = 1e-5
TRACKED = ("30+33+35+37", "30+33+35", "30+33")


def f_core(r, s):
    q = r.q_matrix(s)[0]
    vals = characteristic_values(q, [2] * 4)
    return {frozenset(CORE[i] for i in k): v for k, v in vals.items()}


def point_task(g):
    pin_blas_threads()
    r = Realization(CORE, (g, 1.0, PATH_T, 1.0))
    lam = band_critical(r.vertex_spectrum(CORE))
    f = f_core(r, lam)
    fp, fm = f_core(r, lam + H_S), f_core(r, lam - H_S)
    full = frozenset(CORE)
    dfs = (fp[full] - fm[full]) / (2 * H_S)
    chi = cumulants(f, CORE)
    delta = {
        label(b): float(abs(chi[frozenset(b)] * f[full - frozenset(b)]) / abs(dfs))
        for b in subsets(CORE)
        if len(b) >= 2
    }
    axis = f_core(r, complex(0.0, lam.imag))
    chi_axis = cumulants(axis, CORE)
    raw = {label(b): abs(chi_axis[frozenset(b)]) for b in subsets(CORE) if len(b) >= 2}
    out = {
        "g": g,
        "lambda_V_re": lam.real,
        "lambda_V_hz": lam.imag / (2 * np.pi),
        "abs_F_V_at_lambda": abs(f[full]),
        "nu_V": connected_share(chi, CORE),
        "D_delta": max(delta, key=delta.get),
        "D_raw_axis": max(raw, key=raw.get),
    }
    for b in delta:
        out[f"delta_{b}"] = delta[b]
        out[f"rawaxis_{b}"] = raw[b]
    return out


def main(argv) -> int:
    started = time.time()
    pts = pd.read_csv(
        F7B_POINTS, usecols=["g", "k", "t", "h", "status", "label", "kappa"]
    )
    line = pts[
        (pts.t == PATH_T) & (pts.k == 1.0) & (pts.h == 1.0) & (pts.status == "OK")
    ].sort_values("g")
    with Pool(16) as pool:
        rows = pool.map(point_task, line.g.tolist(), chunksize=2)
    df = pd.DataFrame(rows)
    df["label"] = line.label.to_numpy()
    df["kappa"] = line.kappa.to_numpy()
    nonempty = df.label != "EMPTY"
    hyperedges = df.label.str.split("|")
    df["size_match"] = df.D_delta.str.count(r"\+") + 1 == df.kappa
    df["D_in_H"] = [d in h for d, h in zip(df.D_delta, hyperedges, strict=True)]
    df["size_match_raw"] = df.D_raw_axis.str.count(r"\+") + 1 == df.kappa
    df["raw_in_H"] = [d in h for d, h in zip(df.D_raw_axis, hyperedges, strict=True)]
    df.to_csv(OUT / "CC05_track.csv", index=False)
    a = float(df.loc[nonempty, "size_match"].mean())
    b = float(df.loc[nonempty, "D_in_H"].mean())
    verdict = "SUPPORTED" if (a >= 0.8 and b >= 0.8) else "REJECTED"
    segments = []
    prev = None
    for row in df.itertuples():
        key = (row.label, row.D_delta)
        if key != prev:
            segments.append(
                {
                    "g_start": row.g,
                    "label": row.label,
                    "kappa": row.kappa,
                    "D_delta": row.D_delta,
                    "D_raw_axis": row.D_raw_axis,
                }
            )
            prev = key
    pd.DataFrame(segments).to_csv(OUT / "CC05_segments.csv", index=False)
    summary = {
        "n_points": int(len(df)),
        "n_nonempty": int(nonempty.sum()),
        "fraction_size_match": a,
        "fraction_D_in_H": b,
        "fraction_size_match_raw_axis": float(
            df.loc[nonempty, "size_match_raw"].mean()
        ),
        "fraction_raw_in_H": float(df.loc[nonempty, "raw_in_H"].mean()),
        "H5_verdict": verdict,
        "D_counts_nonempty": df.loc[nonempty, "D_delta"].value_counts().to_dict(),
        "max_abs_F_V_at_lambda": float(df.abs_F_V_at_lambda.max()),
        "elapsed_s": round(time.time() - started, 1),
    }
    write_json(OUT / "CC05_summary.json", summary)
    pd.set_option("display.width", 250)
    print(pd.DataFrame(segments).to_string(index=False))
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
