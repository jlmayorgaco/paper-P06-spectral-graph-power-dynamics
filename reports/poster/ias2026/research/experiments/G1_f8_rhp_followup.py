"""Journal gate 1, follow-up - the F8 planning claims under whole-RHP safety.

G1_whole_rhp.py found that 14.6 % of the frozen F8 (point, intervention)
lattices change H when Gamma = RHP replaces the 0.3-1.5 Hz band. This script
explains every such difference and restates the planning claim under Gamma_RHP:

1. Mode anatomy. For every difference, the RHP modes of each hyperedge of H_RHP:
   class (APERIODIC / IN_BAND / BAND_EDGE_ARTIFACT / OSCILLATORY, thresholds of
   G1_whole_rhp.py), frequency, growth rate, and the device group that carries
   most of the participation (condenser at a retired bus, converter, surviving
   machine).
2. Which condenser configurations are RHP-safe, per point.
3. The electromagnetic-presence thresholds under Gamma_RHP: condenser rating
   swept 1e-3 -> 1 (per-unit services fixed, as in F8C) for the band-minimal
   condenser R_0000000 and for the two single-service additions R_0010000
   (damping 2 pu) and R_0001000 (dynamic EMFs); the rating at which H_RHP
   becomes and stays empty is bisected.

Frozen F8 results are read, never modified.
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
from _f7_common import SUBSETS, Theta
from _gates import GateExperiment
from _overnight import choose_workers, pin_blas_threads
from F8_service_attribution import R_FACTORS, r_configs, solve_config
from G1_whole_rhp import classify_mode, label_of, rhp_modes
from ibr_cycles.diagnosis.composability import parse_hypergraph_label
from ibr_cycles.models.ieee39_case import InfeasibleReplacement

OUT = RESULTS / "G1"
BAD = {"BASE_UNSTABLE", "INFEASIBLE"}
RATINGS = np.geomspace(1e-3, 1.0, 31)
PATH_CONFIGS = ("R_0000000", "R_0010000", "R_0001000")
ITER = 22


def _configs() -> dict:
    return {c["name"]: c for c in r_configs()}


def _group(label: str, members) -> str:
    name, _, dev = label.rpartition("_")
    if dev.startswith("gfl"):
        return "converter"
    bus = int(dev[2:])
    return "condenser" if bus in members else "surviving_machine"


def anatomy(values, vectors_r, vectors_l, labels, members) -> list[dict]:
    """RHP modes with their dominant participation group and state."""

    out = []
    for i in np.flatnonzero(
        (values.real > 0) & (np.abs(values) > 1e-3) & (values.imag >= 0)
    ):
        p = np.abs(vectors_r[:, i] * vectors_l[:, i])
        p = p / p.sum()
        groups: dict[str, float] = {}
        for lab, w in zip(labels, p, strict=True):
            groups[_group(lab, members)] = groups.get(_group(lab, members), 0.0) + w
        top = labels[int(np.argmax(p))]
        v = values[i]
        out.append(
            {
                "re": float(v.real),
                "f_hz": float(abs(v.imag) / (2 * math.pi)),
                "class": classify_mode(v),
                "dominant_group": max(groups, key=groups.get),
                "share_condenser": groups.get("condenser", 0.0),
                "share_converter": groups.get("converter", 0.0),
                "share_surviving": groups.get("surviving_machine", 0.0),
                "top_state": top,
            }
        )
    return out


def anatomy_task(task):
    point, theta_d, config, h_rhp, h_ia = task
    theta = Theta(**theta_d)
    rows = []
    for edge in parse_hypergraph_label(h_rhp):
        members = tuple(sorted(edge))
        case = solve_config(members, theta, config)
        a = case.system.A
        values, right = np.linalg.eig(a)
        left = np.linalg.inv(right).T
        for m in anatomy(values, right, left, case.system.labels, members):
            rows.append(
                {
                    "point": point,
                    "config": config["name"],
                    "H_IA": h_ia,
                    "H_RHP": h_rhp,
                    "hyperedge": "+".join(map(str, members)),
                    **m,
                }
            )
    return rows


def _em_config(config: dict, rating: float) -> dict:
    """Condenser at `rating` with per-unit services fixed at their 0.25 values."""

    cfg = json.loads(json.dumps(config))
    spec = cfg["condenser"]
    spec["S1_inertia"] = spec["S1_inertia"] * rating / 0.25
    spec["EM_rating"] = rating
    return cfg


def _h_rhp(theta, cfg) -> str:
    counts = {}
    for members in SUBSETS:
        try:
            v = np.linalg.eigvals(solve_config(members, theta, cfg).system.A)
        except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError):
            return "INFEASIBLE"
        if not members and v[np.abs(v) > 1e-3].real.max() >= 0:
            return "BASE_UNSTABLE"
        counts[frozenset(members)] = int(rhp_modes(v).size)
    return label_of(counts)


def path_task(task):
    point, theta_d, name, config, rating = task
    return {
        "point": point,
        "config": name,
        "rating": float(rating),
        "H_RHP": _h_rhp(Theta(**theta_d), _em_config(config, rating)),
    }


def threshold_task(task):
    """Bisect the last change to H_RHP = EMPTY between two grid ratings."""

    point, theta_d, name, config, lo, hi = task
    theta = Theta(**theta_d)
    for _ in range(ITER):
        mid = math.sqrt(lo * hi)
        if _h_rhp(theta, _em_config(config, mid)) == "EMPTY":
            hi = mid
        else:
            lo = mid
    return {
        "point": point,
        "config": name,
        "rating_empty": math.sqrt(lo * hi),
        "H_just_below": _h_rhp(theta, _em_config(config, lo)),
    }


def main(argv) -> int:
    pin_blas_threads()
    OUT.mkdir(parents=True, exist_ok=True)
    workers = choose_workers()
    exp = GateExperiment(
        name="G1_f8_rhp_followup",
        question="Do the F8 planning claims survive whole-RHP safety?",
        config={
            "ratings": RATINGS.tolist(),
            "path_configs": PATH_CONFIGS,
            "bisection_iterations": ITER,
        },
        workers=workers,
    )
    started = time.time()
    pts = json.loads((RESULTS / "F8" / "F8_points.json").read_text(encoding="utf-8"))
    thetas = {n: {k: p[k] for k in ("g", "k", "t", "h")} for n, p in pts.items()}
    configs = _configs()
    g1 = pd.read_csv(OUT / "G1_f8_rhp.csv")
    frozen = pd.read_csv(RESULTS / "F8" / "F8_interventions.csv", low_memory=False)
    g1 = g1.merge(
        frozen[["point", "config", *R_FACTORS]], on=["point", "config"], how="left"
    )
    diff = g1[
        (g1.H_RHP != g1.label)
        & ~g1.H_RHP.isin(BAD)
        & ~g1.label.isin(BAD)
        & g1.config.isin(configs)
    ]

    with Pool(workers, initializer=pin_blas_threads) as pool:
        modes = pd.DataFrame(
            [
                r
                for rows in pool.map(
                    anatomy_task,
                    [
                        (r.point, thetas[r.point], configs[r.config], r.H_RHP, r.label)
                        for r in diff.itertuples()
                    ],
                    chunksize=2,
                )
                for r in rows
            ]
        )
        modes.to_csv(OUT / "G1_f8_rhp_modes.csv", index=False)
        path = pd.DataFrame(
            pool.map(
                path_task,
                [
                    (p, thetas[p], n, configs[n], r)
                    for p in pts
                    for n in PATH_CONFIGS
                    for r in RATINGS
                ],
                chunksize=1,
            )
        )
        path.to_csv(OUT / "G1_em_presence_path_rhp.csv", index=False)
        jobs = []
        for (p, n), grp in path.groupby(["point", "config"]):
            grp = grp.sort_values("rating")
            empty = (grp.H_RHP == "EMPTY").to_numpy()
            if not empty[-1] or empty.all():
                continue
            j = int(np.flatnonzero(~empty)[-1])  # last non-empty grid rating
            jobs.append(
                (
                    p,
                    thetas[p],
                    n,
                    configs[n],
                    grp.rating.iloc[j],
                    grp.rating.iloc[j + 1],
                )
            )
        thresholds = (
            pd.DataFrame(pool.map(threshold_task, jobs, chunksize=1))
            if jobs
            else pd.DataFrame()
        )
    thresholds.to_csv(OUT / "G1_em_presence_thresholds_rhp.csv", index=False)

    rc = g1[g1.config.str.match(r"R_\d{7}$")]
    safe_by_point = (
        rc.groupby("point")
        .apply(lambda d: int((d.H_RHP == "EMPTY").sum()), include_groups=False)
        .to_dict()
    )
    ia_by_point = (
        rc.groupby("point")
        .apply(lambda d: int((d.label == "EMPTY").sum()), include_groups=False)
        .to_dict()
    )
    per_factor = {
        f: rc.groupby(f)
        .apply(lambda d: float((d.H_RHP == "EMPTY").mean()), include_groups=False)
        .to_dict()
        for f in R_FACTORS
    }
    unsafe_points = ["P1", "P2", "P3", "P4", "P_fold"]
    safe_everywhere = sorted(
        c
        for c, d in rc[rc.point.isin(unsafe_points[1:])].groupby("config")
        if (d.H_RHP == "EMPTY").all()
    )
    summary = {
        "differences": int(len(diff)),
        "difference_mode_classes": modes["class"].value_counts().to_dict()
        if len(modes)
        else {},
        "difference_dominant_group": modes.dominant_group.value_counts().to_dict()
        if len(modes)
        else {},
        "out_of_band_modes": modes[modes["class"] != "IN_BAND"]
        .groupby(["class", "dominant_group"])
        .agg(
            n=("re", "size"),
            f_min=("f_hz", "min"),
            f_max=("f_hz", "max"),
            re_max=("re", "max"),
        )
        .reset_index()
        .to_dict("records")
        if len(modes)
        else [],
        "R_configs_EMPTY_per_point_RHP": safe_by_point,
        "R_configs_EMPTY_per_point_IA": ia_by_point,
        "EMPTY_fraction_by_factor_level_RHP": per_factor,
        "R_configs_EMPTY_at_P2_P3_P4_Pfold_RHP": safe_everywhere,
        "em_path": path.pivot_table(
            index=["config", "rating"], columns="point", values="H_RHP", aggfunc="first"
        )
        .reset_index()
        .to_dict("records"),
        "em_thresholds": thresholds.to_dict("records"),
    }
    (OUT / "G1_f8_rhp_followup.json").write_text(
        json.dumps(summary, indent=1, default=str), encoding="utf-8"
    )
    exp.finish("COMPUTED", elapsed_s=time.time() - started)
    print(
        json.dumps(
            {k: v for k, v in summary.items() if k != "em_path"}, indent=1, default=str
        )
    )
    pd.set_option("display.width", 250)
    print(
        path.pivot_table(
            index=["config", "rating"], columns="point", values="H_RHP", aggfunc="first"
        ).to_string()
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
