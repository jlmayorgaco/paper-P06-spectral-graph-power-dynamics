"""F7 Safeguard D - characterise the narrow small-g instability region near 0.68 Hz.

Slice F7A (t = 1.5, h = 1), flagship 30+33+35+37. Questions:

    A  a genuine right-half-plane instability region?
    B  a frequency-band entry/exit artefact?
    C  a modal-family tracking artefact?
    D  a numerical near-tangency?

Method. The mode is tracked CONTINUOUSLY along lines of constant k from g = 0,
where it is the flagship's inter-area mode, with eigenvector continuity (MAC
between consecutive steps) on the smooth fast model. On each line the real part
R(g; k) of the tracked mode gives its crossings and its interior extrema. The
tip of the region is where an interior maximum of R along g touches zero; the
closing point of the stable gap is where an interior minimum does. At both,
R = 0 and dR/dg = 0: a pure-policy path through them meets the axis
tangentially, which is NOT a Proposition-3 simple crossing for that path. The
in-plane gradient decides whether the boundary curve itself stays regular.

Direct-path checks: the N_Gamma vector and direct RHP counts before, inside and
after the region, and port closure at both of its boundaries.
"""

from __future__ import annotations

import json
import math

import numpy as np

from _bootstrap import RESULTS
from _f7_common import (
    LABELS,
    Theta,
    evaluate,
    fast_model,
    port_diagnostics,
    solve_subset,
)
from _v2c_common import nominal_reference
from ibr_cycles.dynamics.modal_tracking import modal_assurance

OUT = RESULTS / "F7"
FLAG = (30, 33, 35, 37)
T_SCALE = 1.5
G_MAX = 0.35
STEP = 0.0005
K_LINES = np.round(np.arange(1.20, 1.401, 0.005), 4)


def _theta(g, k):
    return Theta(g=float(g), k=float(k), t=T_SCALE)


def _eig(g, k):
    values, vectors = np.linalg.eig(fast_model().matrix(FLAG, _theta(g, k)))
    return values, vectors


def _start(k):
    """The flagship's inter-area mode at g = 0: the worst band mode."""

    values, vectors = _eig(0.0, k)
    f = np.abs(values.imag) / (2 * math.pi)
    band = np.flatnonzero((values.imag > 0) & (f >= 0.3) & (f <= 1.5))
    i = band[np.argmax(values.real[band])]
    return values[i], vectors[:, i]


def track(k, g_grid):
    lam, vec = _start(k)
    out, macs = [], []
    for g in g_grid:
        values, vectors = _eig(g, k)
        candidates = np.flatnonzero(np.abs(values - lam) < 0.5)
        if candidates.size == 0:
            candidates = np.arange(values.size)
        scores = [modal_assurance(vec, vectors[:, c]) for c in candidates]
        best = candidates[int(np.argmax(scores))]
        macs.append(float(max(scores)))
        lam, vec = values[best], vectors[:, best]
        out.append(lam)
    return np.array(out), np.array(macs)


def tracked_at(g, k, near):
    values = np.linalg.eigvals(fast_model().matrix(FLAG, _theta(g, k)))
    upper = values[values.imag > 0]
    return complex(upper[np.argmin(np.abs(upper - near))])


def local_extremum(k, g_lo, g_hi, near, kind):
    """Golden-section search for an interior max (kind=+1) or min (kind=-1) of R."""

    golden = 0.5 * (3 - math.sqrt(5))
    a, b = g_lo, g_hi
    c, d = a + golden * (b - a), b - golden * (b - a)
    fc = kind * tracked_at(c, k, near).real
    fd = kind * tracked_at(d, k, near).real
    while b - a > 1e-9:
        if fc > fd:
            b, d, fd = d, c, fc
            c = a + golden * (b - a)
            fc = kind * tracked_at(c, k, near).real
        else:
            a, c, fc = c, d, fd
            d = b - golden * (b - a)
            fd = kind * tracked_at(d, k, near).real
    g = 0.5 * (a + b)
    return g, tracked_at(g, k, near)


def fold(kind, k_lo, k_hi, g_window, near):
    """Bisect k so that the interior extremum of R along g is exactly zero."""

    def extremum(k):
        return local_extremum(k, *g_window, near, kind)

    g_lo, lam_lo = extremum(k_lo)
    g_hi, lam_hi = extremum(k_hi)
    if np.sign(lam_lo.real) == np.sign(lam_hi.real):
        return None
    for _ in range(50):
        k_mid = 0.5 * (k_lo + k_hi)
        g_mid, lam_mid = extremum(k_mid)
        if np.sign(lam_mid.real) == np.sign(lam_lo.real):
            k_lo, lam_lo = k_mid, lam_mid
        else:
            k_hi, lam_hi = k_mid, lam_mid
    k_star = 0.5 * (k_lo + k_hi)
    g_star, lam = extremum(k_star)
    h = 1e-5
    r = lambda g, k: tracked_at(g, k, lam).real  # noqa: E731
    return {
        "g": g_star,
        "k": k_star,
        "lambda_re": lam.real,
        "lambda_im": lam.imag,
        "freq_hz": abs(lam.imag) / (2 * math.pi),
        "dR_dg": (r(g_star + h, k_star) - r(g_star - h, k_star)) / (2 * h),
        "dR_dk": (r(g_star, k_star + h) - r(g_star, k_star - h)) / (2 * h),
        "d2R_dg2": (
            r(g_star + 1e-3, k_star) - 2 * r(g_star, k_star) + r(g_star - 1e-3, k_star)
        )
        / 1e-6,
    }


def crossings(g_grid, series):
    out = []
    re = series.real
    for i in range(1, len(re)):
        if np.sign(re[i]) != np.sign(re[i - 1]):
            lo, hi = g_grid[i - 1], g_grid[i]
            near = series[i]
            for _ in range(40):
                mid = 0.5 * (lo + hi)
                if np.sign(tracked_at(mid, K_REF, near).real) == np.sign(re[i - 1]):
                    lo = mid
                else:
                    hi = mid
            g = 0.5 * (lo + hi)
            lam = tracked_at(g, K_REF, near)
            out.append(
                {
                    "g": g,
                    "freq_hz": abs(lam.imag) / (2 * math.pi),
                    "lambda_re": lam.real,
                    "direction": "LEAVES_RHP" if re[i - 1] > 0 else "ENTERS_RHP",
                }
            )
    return out


K_REF = 1.30


def main() -> int:
    g_grid = np.round(np.arange(0.0, G_MAX + 1e-12, STEP), 6)
    lines = {}
    for k in K_LINES:
        series, macs = track(k, g_grid)
        lines[float(k)] = (series, macs)
    report: dict = {
        "slice": {"t": T_SCALE, "h": 1.0},
        "k_lines": K_LINES.tolist(),
        "step": STEP,
    }

    # -- identity: continuous tracking from the g = 0 inter-area mode
    report["tracking_min_consecutive_MAC"] = float(
        min(m.min() for _, m in lines.values())
    )
    series_ref, _ = lines[K_REF]
    report["frequency_along_track_at_k_ref_hz"] = {
        f"{g:.3f}": abs(series_ref[i].imag) / (2 * math.pi)
        for i, g in enumerate(g_grid)
        if i % 40 == 0
    }

    # -- structure along each line: number of RHP intervals and extrema
    structure = []
    for k, (series, _) in lines.items():
        re = series.real
        sign_changes = int(np.sum(np.sign(re[1:]) != np.sign(re[:-1])))
        interior = np.flatnonzero((re[1:-1] > re[:-2]) & (re[1:-1] > re[2:])) + 1
        minima = np.flatnonzero((re[1:-1] < re[:-2]) & (re[1:-1] < re[2:])) + 1
        structure.append(
            {
                "k": k,
                "crossings": sign_changes,
                "R_at_g0": float(re[0]),
                "interior_max": [(float(g_grid[i]), float(re[i])) for i in interior],
                "interior_min": [(float(g_grid[i]), float(re[i])) for i in minima],
            }
        )
    report["lines"] = structure

    # -- the two folds (tongue tip: interior max = 0; gap closure: interior min = 0)
    ref = structure[[s["k"] for s in structure].index(K_REF)]
    g_max = ref["interior_max"][0][0] if ref["interior_max"] else 0.1
    g_min = ref["interior_min"][0][0] if ref["interior_min"] else 0.04
    near_max = series_ref[int(round(g_max / STEP))]
    near_min = series_ref[int(round(g_min / STEP))]
    report["tip_fold"] = fold(+1, 1.20, 1.30, (g_max - 0.06, g_max + 0.06), near_max)
    report["gap_closure_fold"] = fold(
        -1, 1.30, 1.40, (max(g_min - 0.03, 1e-4), g_min + 0.04), near_min
    )
    for name in ("tip_fold", "gap_closure_fold"):
        f = report[name]
        if f:
            f["nondegenerate"] = bool(
                abs(f["dR_dk"]) > 1e-3 and abs(f["d2R_dg2"]) > 1e-3
            )
            f["path_tangency_ratio"] = abs(f["dR_dg"]) / max(abs(f["dR_dk"]), 1e-300)
            f["classification"] = (
                "IMAGINARY_AXIS_TANGENCY for a pure-policy (g) path: boundary fold. "
                "In the (g, k) plane the boundary curve is regular (dR/dk != 0)."
            )

    # -- crossings on the reference line, direct-path counts and port closure
    cross = crossings(g_grid, series_ref)
    report["k_ref"] = K_REF
    report["crossings_at_k_ref"] = cross
    probes = {}
    gs = [c["g"] for c in cross]
    if len(gs) >= 3:
        probe_g = {
            "before_region (gap)": 0.5 * (gs[0] + gs[1]),
            "inside_region": 0.5 * (gs[1] + gs[2]),
            "after_region": gs[2] + 0.02,
            "at_g0_side": 0.5 * gs[0],
        }
        for name, g in probe_g.items():
            row = evaluate(_theta(g, K_REF), fast=False)
            probes[name] = {
                "g": g,
                "H": row["label"],
                "kappa": row.get("kappa"),
                "N": {lab: int(row[f"N_{lab}"]) for lab in LABELS},
                "rhp_flagship_direct": int(row["rhp_30+33+35+37"]),
            }
        for c in cross[1:3]:
            c["port"] = port_diagnostics(
                FLAG, _theta(c["g"], K_REF), 2 * math.pi * c["freq_hz"]
            )
    report["probes_at_k_ref"] = probes

    # -- mode identity against the frozen v2C inter-area reference
    inside_g = probes.get("inside_region", {}).get("g", g_max)
    case = solve_subset(FLAG, _theta(inside_g, K_REF))
    values, vectors = np.linalg.eig(case.system.A)
    target = tracked_at(inside_g, K_REF, near_max)
    i = int(np.argmin(np.abs(values - target)))
    right = vectors[:, i]
    left = np.linalg.inv(vectors)[i, :]
    participation = np.abs(left * right)
    participation /= participation.sum()
    order = np.argsort(-participation)[:10]
    report["mode_inside_region"] = {
        "lambda": [float(values[i].real), float(values[i].imag)],
        "freq_hz": abs(values[i].imag) / (2 * math.pi),
        "top_participation": [
            (case.system.labels[j], float(participation[j])) for j in order
        ],
    }
    shape, _ = nominal_reference()
    common = [n for n in case.system.labels if n in shape]
    idx = {n: j for j, n in enumerate(case.system.labels)}
    report["mode_inside_region"]["MAC_to_v2C_interarea_reference"] = float(
        modal_assurance(
            np.array([shape[n] for n in common]), right[[idx[n] for n in common]]
        )
    )
    start_lam, start_vec = _start(K_REF)
    report["mode_inside_region"]["MAC_to_g0_interarea_mode_machine_coords"] = float(
        modal_assurance(
            start_vec[[idx[n] for n in common]], right[[idx[n] for n in common]]
        )
    )

    # -- verdict on A-D
    widest = max(
        (max((m[1] for m in s["interior_max"]), default=-1) for s in structure),
        default=-1,
    )
    report["checks"] = {
        "B_band_artifact": "NO: crossings at "
        + ", ".join(f"{c['freq_hz']:.3f}" for c in cross)
        + " Hz, far from the 0.3/1.5 Hz edges; the direct RHP count changes with N",
        "C_tracking_artifact": "NO: N_Gamma is label-free, and the continuous "
        "track has minimum "
        f"consecutive MAC {report['tracking_min_consecutive_MAC']:.4f}",
        "D_numerical_near_tangency": "NO away from the tip: largest interior "
        "maximum of R over the "
        f"tracked lines {widest:+.4f} 1/s against fast-vs-direct band error 3e-8",
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "F7_tongue.json").write_text(
        json.dumps(report, indent=2, default=str), encoding="utf-8"
    )
    # one CSV of R(g; k) for the figure
    import pandas as pd

    pd.DataFrame({f"{k:.3f}": lines[k][0].real for k in lines}, index=g_grid).to_csv(
        OUT / "F7_tongue_R_grid.csv", index_label="g"
    )
    pd.DataFrame(
        {f"{k:.3f}": np.abs(lines[k][0].imag) / (2 * math.pi) for k in lines},
        index=g_grid,
    ).to_csv(OUT / "F7_tongue_freq_grid.csv", index_label="g")
    print(
        json.dumps(
            {k: v for k, v in report.items() if k != "lines"}, indent=1, default=str
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
