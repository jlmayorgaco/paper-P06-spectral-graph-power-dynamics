"""F10 - the strongest feasible conventional baselines on the same F7/F8 cases.

Truth: N_Gamma of every core subset (F7 localized assembly, validated against the
direct path) and hence H_Gamma, kappa_Gamma.

B1  generalized Nyquist on the FULL portfolio: winding of det(I + M) on dGamma
    with the open-loop pole correction. One loop, one answer.
B2  generalized Nyquist on EVERY sub-loop: the same winding for every principal
    sub-matrix of I + M. This is our port-closure computation; it is included to
    make the equivalence explicit, not as a competitor.
B3  single-port impedance screening: at each replaced bus alone, the winding of
    its individual factor det(I + M_aa) (minor-loop / impedance-ratio criterion).
B4  modal participation: participation of each core machine's rotor states in
    the base inter-area mode, and a participation-ranked guess of the witness.
B5  first-order modal sensitivity: d lambda / d rho_a of the base inter-area
    mode by a 1 % partial replacement at each bus, extrapolated additively to
    every subset (lambda(S) = lambda_0 + sum_a dlambda_a).
B5b additive and pairwise extrapolation from exhaustive single and pair studies
    of the inter-area family envelope.
B6  SCR / gSCR of each subset, from the network alone.
B7  structured robustness: D-scaled upper bound of the structured singular
    value of M(i omega) for block structure diag(delta_a I_2) over the band;
    mu < 1 certifies that no fractional portfolio crosses the imaginary axis.

Questions (per point, and along one pure-policy line):
    A first destabilizing subset(s)  B all minimal subsets  C same kappa,
    different H  D policy-driven kappa changes  E boundary frequency
    F full-model evaluations required.
"""

from __future__ import annotations

import json
import math
import sys
import time
from itertools import combinations
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.optimize import minimize

from _bootstrap import RESULTS
from _f7_common import LABELS, LEAK, SUBSETS, Theta, evaluate, fast_model
from _overnight import choose_workers, pin_blas_threads
from _postfreeze import PostFreezeExperiment
from _v2c_common import BAND_HZ, CORE
from F7_report import kappa_of, points_file, raster, witnesses
from ibr_cycles.diagnosis.baselines import generalized_scr, nodal_metrics
from ibr_cycles.diagnosis.composability import (
    hypergraph_label,
    hypergraph_order,
    incompatibility_hypergraph,
    parse_hypergraph_label,
)
from ibr_cycles.models.ieee39_case import ReplacementPlan, build_dae, solve_case
from ibr_cycles.models.ieee39_devices import ConverterParameters
from ibr_cycles.models.ieee39_network import load_network, solve_power_flow
from ibr_cycles.models.port_core import PortCore

NAME = "F10_baselines"
OUT = RESULTS / "F10"
LO, HI = 2 * math.pi * BAND_HZ[0], 2 * math.pi * BAND_HZ[1]
SEED = 20260914


def _conv(theta):
    return ConverterParameters(
        voltage_control=True, voltage_gain=theta.g, voltage_leak=LEAK
    )


def _scaling(theta):
    return {"ka": theta.k, "ta": theta.t}


def label(counts):
    return hypergraph_label(
        incompatibility_hypergraph(
            {frozenset(s): int(n > 0) for s, n in counts.items()}
        )
    )


def band_modes(values):
    f = np.abs(values.imag) / (2 * math.pi)
    keep = (
        (values.imag >= 0)
        & (np.abs(values) > 1e-3)
        & (f >= BAND_HZ[0])
        & (f <= BAND_HZ[1])
    )
    return values[keep]


# ---------------------------------------------------------------- baselines ----


def truth(theta) -> dict:
    row = evaluate(theta, fast=True)
    counts = {s: int(row[f"N_{lab}"]) for s, lab in zip(SUBSETS, LABELS, strict=True)}
    alpha = {
        s: float(row[f"bandmax_{lab}"]) for s, lab in zip(SUBSETS, LABELS, strict=True)
    }
    return {
        "counts": counts,
        "alpha": alpha,
        "H": row["label"],
        "kappa": row.get("kappa"),
    }


def port_core(theta):
    base = build_dae(
        ReplacementPlan.of({}), converter=_conv(theta), machine_scaling=_scaling(theta)
    )
    repl = build_dae(
        ReplacementPlan.of({b: 1.0 for b in CORE}),
        converter=_conv(theta),
        machine_scaling=_scaling(theta),
    )
    return PortCore.from_daes(*base, *repl, CORE)


def b1_b2_b3(core: PortCore) -> dict:
    counts, info = core.counts(SUBSETS, omega_lo=LO, omega_hi=HI)
    n = {s: counts[s]["delta_N"] for s in SUBSETS}
    full = n[CORE]
    singles = {s: n[s] for s in SUBSETS if len(s) == 1}
    # frequency where the full-portfolio return difference is closest to -1
    omega = np.linspace(LO, HI, 241)
    m = core.m_matrix(omega * 1j)
    closure = np.min(np.abs(np.linalg.eigvals(m) + 1.0), axis=1)
    return {
        "B1_flagship_unstable": bool(full > 0),
        "B1_freq_hz": float(omega[np.argmin(closure)] / (2 * math.pi)),
        "B1_min_closure": float(closure.min()),
        "B2_counts": n,
        "B3_unsafe_singles": sorted(s[0] for s, v in singles.items() if v > 0),
        "contour_points": info["contour_points"],
    }


def b4_b5(theta) -> dict:
    """Base inter-area mode participation and first-order replacement sensitivity."""

    base = solve_case(
        ReplacementPlan.of({}), converter=_conv(theta), machine_scaling=_scaling(theta)
    )
    values, vectors = np.linalg.eig(base.system.A)
    band = [
        i
        for i in range(values.size)
        if values[i].imag > 0
        and BAND_HZ[0] <= abs(values[i].imag) / (2 * math.pi) <= BAND_HZ[1]
    ]
    i0 = max(band, key=lambda i: values[i].real)  # least-damped band mode
    lam0 = values[i0]
    right = vectors[:, i0]
    left = np.linalg.inv(vectors)[i0, :]
    part = np.abs(left * right)
    part /= part.sum()
    labels = base.system.labels
    pf = {
        b: float(
            sum(
                part[j]
                for j, n in enumerate(labels)
                if n in (f"delta_sg{b}", f"omega_sg{b}")
            )
        )
        for b in CORE
    }
    rho = 0.01
    dlam = {}
    for b in CORE:
        case = solve_case(
            ReplacementPlan.of({b: rho}),
            converter=_conv(theta),
            machine_scaling=_scaling(theta),
        )
        ev = np.linalg.eigvals(case.system.A)
        dlam[b] = (ev[np.argmin(np.abs(ev - lam0))] - lam0) / rho
    pred = {s: lam0 + sum(dlam[b] for b in s) for s in SUBSETS}
    counts = {
        s: int(
            pred[s].real > 0
            and BAND_HZ[0] <= abs(pred[s].imag) / (2 * math.pi) <= BAND_HZ[1]
        )
        for s in SUBSETS
    }
    ranking = sorted(CORE, key=lambda b: -pf[b])
    return {
        "B4_participation": pf,
        "B4_ranking": ranking,
        "B5_counts": counts,
        "B5_freq_flagship_hz": float(abs(pred[CORE].imag) / (2 * math.pi)),
        "base_mode": [float(lam0.real), float(lam0.imag)],
        "evaluations": 1 + len(CORE),
    }


def b5b(alpha: dict) -> dict:
    """Additive (singles) and pairwise extrapolation of the band envelope."""

    a0 = alpha[()]
    single = {b: alpha[(b,)] - a0 for b in CORE}
    pair = {
        (a, b): alpha[(a, b)] - a0 - single[a] - single[b]
        for a, b in combinations(CORE, 2)
    }
    add = {s: a0 + sum(single[b] for b in s) for s in SUBSETS}
    pw = {s: add[s] + sum(pair[p] for p in combinations(s, 2)) for s in SUBSETS}
    return {
        "B5b_additive": {s: int(v > 0) for s, v in add.items()},
        "B5b_pairwise": {
            s: int(v > 0) if len(s) > 2 else int(alpha[s] > 0) for s, v in pw.items()
        },
    }


def b6_scores() -> dict:
    net = load_network()
    flow = solve_power_flow(net)
    metrics = nodal_metrics(net, flow)
    return {
        s: {
            "gscr": float(generalized_scr(net, s, flow)),
            "min_scr": float(min(metrics[b].scr for b in s)),
        }
        for s in SUBSETS
        if s
    }


def b7_mu(core: PortCore) -> dict:
    """D-scaled upper bound of complex structured mu over the band."""

    omega = np.linspace(LO, HI, 61)
    m_all = core.m_matrix(omega * 1j)
    bound = []
    for m in m_all:

        def sigma(logd, m=m):
            d = np.repeat(np.exp(np.concatenate([[0.0], logd])), 2)
            return float(
                np.linalg.svd((d[:, None] * m) / d[None, :], compute_uv=False)[0]
            )

        res = minimize(
            sigma,
            np.zeros(len(CORE) - 1),
            method="Nelder-Mead",
            options={"xatol": 1e-4, "fatol": 1e-6, "maxiter": 400},
        )
        bound.append(res.fun)
    return {
        "B7_mu_upper": float(max(bound)),
        "B7_certifies_safe": bool(max(bound) < 1.0),
    }


# ------------------------------------------------------------------- scoring ----


def score(true_counts, pred_counts) -> dict:
    h_true, h_pred = label(true_counts), label(pred_counts)
    k_true = hypergraph_order(parse_hypergraph_label(h_true))
    k_pred = hypergraph_order(parse_hypergraph_label(h_pred))
    return {
        "H_exact": h_true == h_pred,
        "kappa_exact": k_true == k_pred,
        "witness_exact": witnesses(h_true) == witnesses(h_pred)
        if h_true != "EMPTY" or h_pred != "EMPTY"
        else True,
        "H_pred": h_pred,
    }


def point_task(task):
    name, theta_d = task
    theta = Theta(**theta_d)
    if evaluate(theta, fast=True)["label"] == "BASE_UNSTABLE":
        return {"point": name, **theta_d, "H_true": "BASE_UNSTABLE"}
    t = truth(theta)
    core = port_core(theta)
    g = b1_b2_b3(core)
    s = b4_b5(theta)
    bb = b5b(t["alpha"])
    mu = b7_mu(core)
    row = {
        "point": name,
        **theta_d,
        "H_true": t["H"],
        "kappa_true": t["kappa"],
        "B1_flagship_unstable": g["B1_flagship_unstable"],
        "B1_matches_flagship": g["B1_flagship_unstable"] == (t["counts"][CORE] > 0),
        "B1_freq_hz": g["B1_freq_hz"],
        "B1_min_closure": g["B1_min_closure"],
        "B3_unsafe_singles": "+".join(map(str, g["B3_unsafe_singles"])),
        "B4_ranking": "+".join(map(str, s["B4_ranking"])),
        "B7_mu_upper": mu["B7_mu_upper"],
        "B7_certifies_safe": mu["B7_certifies_safe"],
        "B5_freq_flagship_hz": s["B5_freq_flagship_hz"],
    }
    for key, counts in (
        ("B2", g["B2_counts"]),
        ("B5", s["B5_counts"]),
        ("B5b_add", bb["B5b_additive"]),
        ("B5b_pair", bb["B5b_pairwise"]),
    ):
        for k, v in score(t["counts"], counts).items():
            row[f"{key}_{k}"] = v
    # B3 predicts only singleton hyperedges
    single_true = sorted(s[0] for s in SUBSETS if len(s) == 1 and t["counts"][s] > 0)
    row["B3_singletons_exact"] = single_true == g["B3_unsafe_singles"]
    row["B3_recovers_H"] = bool(
        t["H"] == "|".join(map(str, single_true)) if single_true else t["H"] == "EMPTY"
    )
    # the flagship crossing frequency, when the truth flagship is near a boundary
    fv = band_modes(np.linalg.eigvals(fast_model().matrix(CORE, theta)))
    row["true_flagship_worst_freq_hz"] = float(
        abs(fv[np.argmax(fv.real)].imag) / (2 * math.pi)
    )
    return row


def same_kappa_pairs() -> dict:
    """Interior points of the two largest H regions sharing each kappa (frozen F7A)."""

    points = pd.read_csv(points_file("F7A"), low_memory=False)
    leaves = pd.read_csv(RESULTS / "F7" / "F7A_leaves.csv")
    grid = raster(points, leaves)
    xs = (
        points.groupby("I")
        .x.first()
        .reindex(range(grid.shape[0]))
        .interpolate()
        .to_numpy()
    )
    ys = (
        points.groupby("J")
        .y.first()
        .reindex(range(grid.shape[1]))
        .interpolate()
        .to_numpy()
    )
    out = {}
    for kappa in (1, 2, 3):
        labels = sorted(
            {s for s in set(grid.ravel()) if s and kappa_of(s) == kappa},
            key=lambda s: -(grid == s).sum(),
        )[:2]
        for n, lab in enumerate(labels):
            mask = grid == lab
            depth = ndimage.distance_transform_edt(mask)
            i, j = np.unravel_index(np.argmax(depth), depth.shape)
            out[f"K{kappa}_{n}"] = {
                "g": float(xs[i]),
                "k": float(ys[j]),
                "t": 1.5,
                "h": 1.0,
            }
    return out


def _init():
    pin_blas_threads()
    fast_model()


def main(argv) -> int:
    pin_blas_threads()
    OUT.mkdir(parents=True, exist_ok=True)
    pts = json.loads((RESULTS / "F8" / "F8_points.json").read_text(encoding="utf-8"))
    points = {n: {k: p[k] for k in ("g", "k", "t", "h")} for n, p in pts.items()}
    points.update(same_kappa_pairs())
    spots = pd.read_csv(RESULTS / "F7" / "F7A_spotcheck.csv")
    spots = spots[spots.direct_label != "BASE_UNSTABLE"].reset_index(drop=True)
    rng = np.random.default_rng(SEED)
    for i, r in enumerate(
        spots.iloc[rng.choice(len(spots), 40, replace=False)].itertuples()
    ):
        points[f"S{i:02d}"] = {"g": float(r.x), "k": float(r.y), "t": 1.5, "h": 1.0}
    line = {
        f"L{i:02d}": {"g": float(g), "k": 1.0, "t": 0.852, "h": 1.0}
        for i, g in enumerate(np.round(np.linspace(0.0, 0.35, 36), 4))
    }
    workers = choose_workers()
    exp = PostFreezeExperiment(
        name=NAME,
        question="Can conventional methods recover H and its policy dependence?",
        config={"points": len(points), "line": "F7B t=0.852, g 0-0.35"},
        workers=workers,
    )
    started = time.time()
    with Pool(workers, initializer=_init) as pool:
        rows = pool.map(
            point_task, list(points.items()) + list(line.items()), chunksize=1
        )
    frame = pd.DataFrame(rows)
    frame["set"] = np.where(
        frame.point.str.startswith("L"),
        "policy_line",
        np.where(
            frame.point.str.startswith("K"),
            "same_kappa_pairs",
            np.where(frame.point.str.startswith("S"), "random_F7A", "F8_points"),
        ),
    )
    scr = b6_scores()
    pd.DataFrame(
        [{"subset": "+".join(map(str, s)), **v} for s, v in scr.items()]
    ).to_csv(OUT / "F10_scr.csv", index=False)
    frame.to_csv(OUT / "F10_baselines.csv", index=False)
    frame = frame[frame.H_true != "BASE_UNSTABLE"]
    pts_only = frame[frame.set != "policy_line"]
    summary = {
        "points": int(len(pts_only)),
        "B1_flagship_status_correct": float(pts_only.B1_matches_flagship.mean()),
        "B2_H_exact": float(pts_only.B2_H_exact.mean()),
        "B3_singletons_exact": float(pts_only.B3_singletons_exact.mean()),
        "B3_recovers_H": float(pts_only.B3_recovers_H.mean()),
        "B5_sensitivity_H_exact": float(pts_only.B5_H_exact.mean()),
        "B5_sensitivity_kappa_exact": float(pts_only.B5_kappa_exact.mean()),
        "B5b_additive_H_exact": float(pts_only.B5b_add_H_exact.mean()),
        "B5b_additive_kappa_exact": float(pts_only.B5b_add_kappa_exact.mean()),
        "B5b_pairwise_H_exact": float(pts_only.B5b_pair_H_exact.mean()),
        "B5b_pairwise_kappa_exact": float(pts_only.B5b_pair_kappa_exact.mean()),
        "B7_certifies_safe_when_safe": float(
            pts_only[pts_only.H_true == "EMPTY"].B7_certifies_safe.mean()
        )
        if (pts_only.H_true == "EMPTY").any()
        else None,
        "B7_false_certification": int(
            (pts_only.B7_certifies_safe & (pts_only.H_true != "EMPTY")).sum()
        ),
        "B6_scr": "static: identical at every policy point by construction",
        "line_true_kappa": frame[frame.set == "policy_line"][
            ["g", "kappa_true"]
        ].to_dict("records"),
    }
    for method in ("B2", "B5", "B5b_add", "B5b_pair"):
        ln = frame[frame.set == "policy_line"]
        summary[f"line_{method}_kappa_exact"] = float(
            ln[f"{method}_kappa_exact"].mean()
        )
    (OUT / "F10_summary.json").write_text(
        json.dumps(summary, indent=2, default=str), encoding="utf-8"
    )
    exp.finish(
        "COMPUTED",
        **{k: v for k, v in summary.items() if k != "line_true_kappa"},
        elapsed_s=time.time() - started,
    )
    print(
        json.dumps(
            {k: v for k, v in summary.items() if k != "line_true_kappa"},
            indent=1,
            default=str,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
