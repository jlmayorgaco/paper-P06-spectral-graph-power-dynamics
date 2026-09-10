"""F8B - excitation fleet MEAN versus machine-to-machine HETEROGENEITY.

F7C varied heterogeneity at fixed GEOMETRIC mean time constant and F7B varied
the mean at native heterogeneity; neither separates the two, and "the mean" has
more than one physical meaning. Here, for a heterogeneity amplitude h and a
target mean m,

    TE_i(m, h) = c(m, h) * (TE_i / G)^h ,   c chosen so that MEAN(TE(m, h)) = m

with MEAN one of

    geometric   exp(mean(log TE))            (what F7C held fixed)
    arithmetic  mean(TE)                     (mean exciter time constant)
    harmonic    1 / mean(1 / TE)             (mean exciter BANDWIDTH, 1/TE)

so every (m, h) pair is an exactly matched case: design A varies h at fixed m,
design B varies m at fixed h, everything else (KA, network, dispatch, PSS, core)
fixed. The same is done for the gain KA at native time constants.

Evaluation uses the exact parametric assembly of A_red (F7, validated against
the direct path) generalised to per-machine KA and TE; a direct-path check is
run on a random subset of the grid.
"""

from __future__ import annotations

import json
import math
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd

from _bootstrap import RESULTS
from _f7_common import LABELS, SUBSETS, fast_model, in_gamma, native_ka, native_te
from _overnight import choose_workers, pin_blas_threads
from _postfreeze import PostFreezeExperiment
from _v2c_common import BAND_HZ, CORE
from ibr_cycles.diagnosis.composability import (
    hypergraph_label,
    hypergraph_order,
    incompatibility_hypergraph,
)
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case
from ibr_cycles.models.ieee39_devices import ConverterParameters

NAME = "F8B_mean_vs_heterogeneity"
OUT = RESULTS / "F8B"
MEANS = ("geometric", "arithmetic", "harmonic")
H_GRID = np.round(np.linspace(0.0, 2.0, 21), 3)
M_SCALE = np.round(
    np.linspace(0.6, 2.2, 33), 4
)  # target mean / native mean of that type
G_VALUES = (0.0, 0.05, 0.2)
DIRECT_CHECKS = 120
SEED = 20260912


def fleet_mean(values: np.ndarray, kind: str) -> float:
    if kind == "geometric":
        return float(np.exp(np.mean(np.log(values))))
    if kind == "arithmetic":
        return float(np.mean(values))
    if kind == "harmonic":
        return float(1.0 / np.mean(1.0 / values))
    raise ValueError(kind)


def matched(
    native: dict[int, float], kind: str, scale: float, h: float
) -> dict[int, float]:
    buses = sorted(native)
    base = np.array([native[b] for b in buses])
    g = fleet_mean(base, "geometric")
    shape = (base / g) ** h
    target = scale * fleet_mean(base, kind)
    c = target / fleet_mean(shape, kind)
    return {b: float(c * s) for b, s in zip(buses, shape, strict=True)}


def matrix(members, g: float, ka: dict, te: dict) -> np.ndarray:
    a0, delta_g, efd, _, _ = fast_model().parts[members]
    a = a0 + g * delta_g
    for i, bus, u, w in efd:
        a[i] = (ka[bus] / te[bus]) * u + (1.0 / te[bus]) * w
    return a


def classify(g, ka, te, *, direct=False) -> dict:
    counts, row = {}, {}
    for members, _label in zip(SUBSETS, LABELS, strict=True):
        if direct:
            scaling = {
                b: {"ka": ka[b] / native_ka()[b], "ta": te[b] / native_te()[b]}
                for b in native_te()
            }
            case = solve_case(
                ReplacementPlan.of({b: 1.0 for b in members}),
                converter=ConverterParameters(
                    voltage_control=True, voltage_gain=g, voltage_leak=0.05
                ),
                machine_bus_scaling=scaling,
            )
            values = np.linalg.eigvals(case.system.A)
        else:
            values = np.linalg.eigvals(matrix(members, g, ka, te))
        dyn = values[np.abs(values) > 1e-3]
        if not members and dyn.real.max() >= 0.0:
            return {"label": "BASE_UNSTABLE"}
        counts[frozenset(members)] = int(np.count_nonzero(in_gamma(values)))
        if members == CORE:
            f = np.abs(values.imag) / (2 * math.pi)
            band = (
                (values.imag >= 0)
                & (np.abs(values) > 1e-3)
                & (f >= BAND_HZ[0])
                & (f <= BAND_HZ[1])
            )
            row["flagship_band_alpha"] = float(values.real[band].max())
    edges = incompatibility_hypergraph(counts)
    row.update(label=hypergraph_label(edges), kappa=hypergraph_order(edges))
    return row


def _task(task):
    variable, kind, g, scale, h = task
    ka, te = native_ka(), native_te()
    if variable == "TE":
        te = matched(te, kind, scale, h)
    else:
        ka = matched(ka, kind, scale, h)
    row = classify(g, ka, te)
    return {
        "variable": variable,
        "mean": kind,
        "g": g,
        "m_scale": scale,
        "h": h,
        **row,
        "te_spread": float(np.ptp(list(te.values()))),
        "ka_spread": float(np.ptp(list(ka.values()))),
    }


def _direct_task(task):
    variable, kind, g, scale, h = task
    ka, te = native_ka(), native_te()
    if variable == "TE":
        te = matched(te, kind, scale, h)
    else:
        ka = matched(ka, kind, scale, h)
    return classify(g, ka, te, direct=True)["label"]


def _init():
    pin_blas_threads()
    fast_model()


def boundary_vs_h(frame: pd.DataFrame) -> pd.DataFrame:
    """Flagship stability boundary m*(h) at each fixed h, by interpolation in m."""

    rows = []
    for (var, kind, g, h), grp in frame.groupby(["variable", "mean", "g", "h"]):
        grp = grp.dropna(subset=["flagship_band_alpha"]).sort_values("m_scale")
        a, m = grp.flagship_band_alpha.to_numpy(), grp.m_scale.to_numpy()
        roots = [
            float(m[i] - a[i] * (m[i + 1] - m[i]) / (a[i + 1] - a[i]))
            for i in range(len(a) - 1)
            if np.sign(a[i]) != np.sign(a[i + 1])
        ]
        rows.append(
            {
                "variable": var,
                "mean": kind,
                "g": g,
                "h": h,
                "flagship_boundaries_in_m": roots,
            }
        )
    return pd.DataFrame(rows)


def anova(frame: pd.DataFrame) -> pd.DataFrame:
    """Flagship-margin variance on the (m, h) grid: mean, heterogeneity, interaction."""

    rows = []
    for (var, kind, g), grp in frame.groupby(["variable", "mean", "g"]):
        pivot = grp.pivot_table(
            index="m_scale", columns="h", values="flagship_band_alpha"
        )
        y = pivot.to_numpy()
        if np.isnan(y).any():
            y = np.where(np.isnan(y), np.nanmean(y), y)
        grand = y.mean()
        row_eff = y.mean(axis=1, keepdims=True) - grand
        col_eff = y.mean(axis=0, keepdims=True) - grand
        inter = y - grand - row_eff - col_eff
        total = ((y - grand) ** 2).sum()
        rows.append(
            {
                "variable": var,
                "mean": kind,
                "g": g,
                "share_mean": float((row_eff**2).sum() * y.shape[1] / total),
                "share_heterogeneity": float((col_eff**2).sum() * y.shape[0] / total),
                "share_interaction": float((inter**2).sum() / total),
            }
        )
    return pd.DataFrame(rows)


def label_changes(frame: pd.DataFrame) -> pd.DataFrame:
    """Design A: H changes along h at fixed m. Design B: along m at fixed h."""

    rows = []
    for (var, kind, g), grp in frame.groupby(["variable", "mean", "g"]):
        ok = grp[grp.label != "BASE_UNSTABLE"]
        a_lines = ok.groupby("m_scale").label.nunique()
        b_lines = ok.groupby("h").label.nunique()
        rows.append(
            {
                "variable": var,
                "mean": kind,
                "g": g,
                "A_fixed_mean_lines": int(len(a_lines)),
                "A_lines_where_h_changes_H": int((a_lines > 1).sum()),
                "B_fixed_h_lines": int(len(b_lines)),
                "B_lines_where_mean_changes_H": int((b_lines > 1).sum()),
                "distinct_H": int(ok.label.nunique()),
                "kappa_values": sorted(int(k) for k in ok.kappa.dropna().unique()),
                "base_unstable_points": int((grp.label == "BASE_UNSTABLE").sum()),
            }
        )
    return pd.DataFrame(rows)


def main(argv) -> int:
    pin_blas_threads()
    OUT.mkdir(parents=True, exist_ok=True)
    tasks = [
        ("TE", kind, g, s, h)
        for kind in MEANS
        for g in G_VALUES
        for s in M_SCALE
        for h in H_GRID
    ]
    tasks += [
        ("KA", kind, 0.0, s, h) for kind in MEANS for s in M_SCALE for h in H_GRID
    ]
    workers = choose_workers()
    exp = PostFreezeExperiment(
        name=NAME,
        question="Is heterogeneity itself causal?",
        config={
            "means": MEANS,
            "h_grid": H_GRID.tolist(),
            "m_scale": M_SCALE.tolist(),
            "g": G_VALUES,
            "direct_checks": DIRECT_CHECKS,
            "seed": SEED,
        },
        workers=workers,
    )
    started = time.time()
    with Pool(workers, initializer=_init) as pool:
        rows = pool.map(_task, tasks, chunksize=16)
        frame = pd.DataFrame(rows)
        rng = np.random.default_rng(SEED)
        pick = rng.choice(len(tasks), DIRECT_CHECKS, replace=False)
        direct = pool.map(_direct_task, [tasks[i] for i in pick], chunksize=1)
    check = frame.iloc[pick].assign(direct_label=direct)
    check["agree"] = check.label == check.direct_label
    frame.to_csv(OUT / "F8B_grid.csv", index=False)
    check.to_csv(OUT / "F8B_direct_check.csv", index=False)
    changes, shares, bounds = label_changes(frame), anova(frame), boundary_vs_h(frame)
    changes.to_csv(OUT / "F8B_label_changes.csv", index=False)
    shares.to_csv(OUT / "F8B_variance_shares.csv", index=False)
    bounds.to_csv(OUT / "F8B_boundary_vs_h.csv", index=False)
    summary = {
        "points": len(frame),
        "direct_agreement": float(check.agree.mean()),
        "direct_checks": int(len(check)),
    }
    (OUT / "F8B_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    exp.finish("COMPUTED", **summary, elapsed_s=time.time() - started)
    pd.set_option("display.width", 250)
    print(changes.to_string(index=False))
    print(shares.to_string(index=False))
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
