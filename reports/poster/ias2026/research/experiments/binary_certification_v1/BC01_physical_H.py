"""BC00-B / BC01: H and kappa regenerated with the physical criterion, no disc.

For every audited (point, subset) of configs/binary_certification_v1/BC00_config.yaml:
the rotation symmetry is removed by T1, the verified neutral-frequency mode is
deflated as a declared exact eigenvector (and kept in the ledger), and the rest
is classified STABLE / UNSTABLE / BOUNDARY_OR_UNRESOLVED. No eigenvalue is
removed by magnitude.

    H_phys = minimal UNSTABLE subsets. It is exact only when no subset is
    BOUNDARY_OR_UNRESOLVED. Otherwise the unresolved subsets are listed and every
    hyperedge whose minimality depends on them is marked.

Compared against the frozen G1 labels (H_RHP with the |s| > 1e-3 disc).
Label changes are corrections of the criterion, reported as such.

Usage: python BC01_physical_H.py {ieee39|kundur|ieee68|all}
Deviation from the preregistration: the F7 points with g <= 0.005 (13 158
points) use the exact F7 parametric assembly (_f7_common.FastModel, validated
0/1 153 against the direct path in F7) instead of the direct path, for cost.
The error matrix is assembled from the |dA| of its three anchor solves. The 200
random points per map and the F8 and tongue points stay on the direct path.
"""

from __future__ import annotations

import json
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from _bc import WORKERS, BCExperiment, out_dir, write_json

from _f7_common import SUBSETS as SUB39  # noqa: E402
from _f7_common import Theta, fast_model, solve_subset  # noqa: E402
from _overnight import pin_blas_threads  # noqa: E402
from F8_service_attribution import r_configs, solve_config  # noqa: E402
from F12_kundur import SUBSETS as SUBK  # noqa: E402
from F12_kundur import solve as k_solve  # noqa: E402
from G3_ieee68 import G_GRID, K_GRID  # noqa: E402
from G3_ieee68 import OUT as G3_OUT
from G3_ieee68 import SUBSETS as SUB68  # noqa: E402
from G3_ieee68 import solve as s68  # noqa: E402
from ibr_cycles.certification.classify import BOUNDARY, STABLE, UNSTABLE  # noqa: E402
from ibr_cycles.certification.physical import (  # noqa: E402
    classify_physical,
    physical_matrices,
    physical_report,
)
from ibr_cycles.certification.symmetry import (  # noqa: E402
    frequency_partner,
    rotation_generator,
)
from ibr_cycles.diagnosis.composability import hypergraph_label  # noqa: E402

OUT = out_dir("BC01")
SEED = 20260911
RES = out_dir("BC01").parent.parent  # research/results


# ------------------------------------------------------------------ helpers --


def h_physical(statuses: dict) -> dict:
    """H from per-subset statuses; exact only without unresolved subsets."""

    base = statuses[frozenset()]
    if base != STABLE:
        return {"H_phys": f"BASE_{base}", "exact": False, "unresolved": ""}
    unsafe = [s for s, v in statuses.items() if v == UNSTABLE]
    unresolved = [s for s, v in statuses.items() if v == BOUNDARY]
    minimal = [s for s in unsafe if not any(r < s for r in unsafe)]
    doubtful = [e for e in minimal if any(u < e for u in unresolved)]
    label = hypergraph_label(minimal)
    sizes = [len(e) for e in minimal]
    kappa = min(sizes) if sizes else -1
    return {
        "H_phys": label,
        "kappa_phys": kappa,
        "exact": not unresolved,
        "n_unresolved": len(unresolved),
        "unresolved": "|".join("+".join(map(str, sorted(u))) for u in unresolved),
        "edges_with_unresolved_subsets": len(doubtful),
    }


def _row(tag, point, members, rep) -> dict:
    return {
        **point,
        "tag": tag,
        "subset": "+".join(map(str, members)) or "BASE",
        "status_physical": rep.status_physical,
        "status_rest": rep.status_rest,
        "n_positive_rest": rep.n_positive_rest,
        "n_positive_real_rest": rep.n_positive_real_rest,
        "abscissa_rest": rep.abscissa_rest,
        "d_axis_rest": rep.d_axis_rest,
        "threshold": rep.threshold,
        "neutral_verified": rep.neutral_verified,
        "n_dead": rep.n_dead,
        "a_rot": rep.identities.get("a_rot"),
        "a_jordan": rep.identities.get("a_jordan"),
        "nz_full_1e-3": rep.near_zero_full["|lambda|<=0.001"],
        "nz_rest_1e-2": rep.near_zero_rest["|lambda|<=0.01"],
        "nz_rest_1e-3": rep.near_zero_rest["|lambda|<=0.001"],
        "positive_rest": ";".join(f"{v:.4e}" for v in rep.positive_rest),
        "near_axis_rest": ";".join(f"{v:.4e}" for v in rep.near_axis_rest),
    }


# ------------------------------------------------------------------ IEEE-39 --

_R_NONE = None
_ANCH: dict = {}


def _init39():
    global _R_NONE
    pin_blas_threads()
    _R_NONE = {c["name"]: c for c in r_configs()}["R_none"]


def task39_direct(task):
    tag, point = task
    theta = Theta(point["g"], point["k"], point["t"], point["h"])
    rows = []
    for m in SUB39:
        case = solve_config(m, theta, _R_NONE)
        rows.append(_row(tag, point, m, physical_report(case)))
    return rows


def _anchors39():
    """Per subset: labels, r_x, partner, |dA| of the three FastModel anchors."""

    if _ANCH:
        return _ANCH
    for m in SUB39:
        errs, first = [], None
        for th in (Theta(0.0), Theta(1.0), Theta(0.0, k=2.0)):
            case = solve_subset(m, th)
            a, d, _ = physical_matrices(case)
            errs.append(np.abs(d))
            first = first or case
        r_x, _ = rotation_generator(first.dae, first.equilibrium.z)
        _ANCH[m] = (
            first.dae.labels,
            r_x,
            frequency_partner(first.dae),
            errs[0] + errs[1] + errs[2],
        )
    return _ANCH


def _init39_fast():
    _init39()
    fast_model()
    _anchors39()


def task39_fast(task):
    tag, point = task
    theta = Theta(point["g"], point["k"], point["t"], point["h"])
    fm, anch = fast_model(), _anchors39()
    rows = []
    for m in SUB39:
        labels, r_x, partner, err = anch[m]
        a = fm.matrix(m, theta)
        rows.append(
            _row(tag, point, m, classify_physical(a, err, labels, r_x, partner))
        )
    return rows


def ieee39(pool_workers: int) -> pd.DataFrame:
    frozen = pd.read_csv(RES / "G1" / "G1_ieee39_points.csv.gz")
    f8 = json.loads((RES / "F8" / "F8_points.json").read_text(encoding="utf-8"))
    direct = [
        (
            "F8",
            {
                "point": n,
                "map": "F8",
                "g": p["g"],
                "k": p["k"],
                "t": p["t"],
                "h": p["h"],
            },
        )
        for n, p in f8.items()
    ]
    direct += [
        (
            "TONGUE",
            {
                "point": f"k1.30_g{g}",
                "map": "TONGUE",
                "g": g,
                "k": 1.30,
                "t": 1.5,
                "h": 1.0,
            },
        )
        for g in (0.020, 0.042, 0.102, 0.173)
    ]
    rng = np.random.default_rng(SEED)
    fast = []
    for name, grp in frozen.groupby("map"):
        pick = grp.iloc[rng.choice(len(grp), size=200, replace=False)]
        for idx, r in pick.iterrows():
            direct.append(
                (
                    "RANDOM",
                    {
                        "point": f"{name}:{idx}",
                        "map": name,
                        "g": r.g,
                        "k": r.k,
                        "t": r.t,
                        "h": r.h,
                        "H_RHP_G1": r.H_RHP,
                        "H_IA": r.H_IA,
                    },
                )
            )
        for idx, r in grp[grp.g <= 0.005].iterrows():
            fast.append(
                (
                    "G_LE_0.005",
                    {
                        "point": f"{name}:{idx}",
                        "map": name,
                        "g": r.g,
                        "k": r.k,
                        "t": r.t,
                        "h": r.h,
                        "H_RHP_G1": r.H_RHP,
                        "H_IA": r.H_IA,
                    },
                )
            )
    rows = []
    with Pool(pool_workers, initializer=_init39) as pool:
        for out in pool.imap_unordered(task39_direct, direct, chunksize=2):
            rows += out
    with Pool(pool_workers, initializer=_init39_fast) as pool:
        for out in pool.imap_unordered(task39_fast, fast, chunksize=16):
            rows += out
    return pd.DataFrame(rows)


# ------------------------------------------------------------------- Kundur --


def taskK(task):
    point = task
    th = {"g": point["g"], "k": point["k"], "t": point["t"]}
    rows = []
    for m in SUBK:
        rows.append(_row("NODE", point, m, physical_report(k_solve(m, th))))
    return rows


def kundur(pool_workers: int) -> pd.DataFrame:
    nodes = pd.read_csv(RES / "G1" / "G1_kundur_nodes.csv")
    tasks = [
        {
            "point": f"{r.map}:{r.i}:{r.j}",
            "map": r.map,
            "g": r.g,
            "k": r.k,
            "t": r.t,
            "H_RHP_G1": r.H_RHP,
            "H_IA": r.H_IA,
        }
        for r in nodes.itertuples()
    ]
    rows = []
    with Pool(pool_workers, initializer=pin_blas_threads) as pool:
        for out in pool.imap_unordered(taskK, tasks, chunksize=8):
            rows += out
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ IEEE-68 --

_P68: dict = {}


def _init68():
    pin_blas_threads()
    for m in SUB68:
        mats = []
        first = None
        for g, k in ((0.0, 1.0), (1.0, 1.0), (0.0, 2.0)):
            case = s68(m, g, k)
            a, d, _ = physical_matrices(case)
            mats.append((a, np.abs(d)))
            first = first or case
        r_x, _ = rotation_generator(first.dae, first.equilibrium.z)
        _P68[m] = (mats, first.dae.labels, r_x, frequency_partner(first.dae))


def task68(task):
    point = task
    g, k = point["g"], point["k"]
    rows = []
    for m in SUB68:
        (a00, e00), (a10, e10), (a02, e02) = _P68[m][0]
        labels, r_x, partner = _P68[m][1:]
        a = a00 + g * (a10 - a00) + (k - 1.0) * (a02 - a00)
        err = e00 + e10 + e02
        rows.append(
            _row("NODE", point, m, classify_physical(a, err, labels, r_x, partner))
        )
    return rows


def ieee68(pool_workers: int) -> pd.DataFrame:
    nodes = pd.read_csv(G3_OUT / "G3_nodes.csv")
    j1 = int(np.argmin(np.abs(K_GRID - 1.0)))
    pick = nodes[(nodes.j == j1) | (nodes.i == 0)]
    rng = np.random.default_rng(SEED)
    rest = nodes.drop(pick.index)
    pick = pd.concat([pick, rest.iloc[rng.choice(len(rest), size=60, replace=False)]])
    tasks = [
        {
            "point": f"G68A:{r.i}:{r.j}",
            "map": "G68A",
            "g": float(G_GRID[r.i]),
            "k": float(K_GRID[r.j]),
            "H_RHP_G1": r.H_RHP,
            "H_IA": r.H_IA,
        }
        for r in pick.itertuples()
    ]
    rows = []
    with Pool(min(pool_workers, 8), initializer=_init68) as pool:
        for out in pool.imap_unordered(task68, tasks, chunksize=2):
            rows += out
    return pd.DataFrame(rows)


# --------------------------------------------------------------------- main --


def summarize(frame: pd.DataFrame, name: str) -> tuple[pd.DataFrame, dict]:
    points = []
    for point, grp in frame.groupby("point"):
        statuses = {
            frozenset(int(b) for b in s.split("+")) if s != "BASE" else frozenset(): v
            for s, v in zip(grp.subset, grp.status_rest, strict=True)
        }
        info = h_physical(statuses)
        first = grp.iloc[0]
        points.append(
            {
                "point": point,
                "map": first.get("map"),
                "tag": first.tag,
                "g": first.g,
                "k": first.k,
                **info,
                "H_RHP_G1": first.get("H_RHP_G1", np.nan),
                "H_IA": first.get("H_IA", np.nan),
                "neutral_all": bool(grp.neutral_verified.all()),
                "max_a_rot": float(grp.a_rot.max()),
                "max_a_jordan": float(grp.a_jordan.max()),
                "positive_real_small": bool(
                    (
                        (grp.status_rest == UNSTABLE)
                        & (grp.n_positive_real_rest > 0)
                        & (grp.abscissa_rest < 1e-3)
                    ).any()
                ),
            }
        )
    pts = pd.DataFrame(points)
    has_g1 = pts.H_RHP_G1.notna()
    cmp_ = pts[has_g1]
    summary = {
        "benchmark": name,
        "points": int(len(pts)),
        "subset_cases": int(len(frame)),
        "status_rest_counts": frame.status_rest.value_counts().to_dict(),
        "status_physical_counts": frame.status_physical.value_counts().to_dict(),
        "neutral_verified_fraction": float(frame.neutral_verified.mean()),
        "max_rotation_residual": float(frame.a_rot.max()),
        "max_jordan_residual": float(frame.a_jordan.max()),
        "points_exact": int(pts.exact.sum()),
        "points_with_unresolved": int((~pts.exact).sum()),
        "compared_with_G1": int(len(cmp_)),
        "H_phys_equals_G1": int((cmp_.H_phys == cmp_.H_RHP_G1).sum()),
        "label_changes": cmp_[cmp_.H_phys != cmp_.H_RHP_G1][
            ["point", "tag", "g", "k", "H_RHP_G1", "H_phys", "exact", "unresolved"]
        ]
        .head(60)
        .to_dict("records"),
        "points_with_small_positive_real_mode": int(pts.positive_real_small.sum()),
    }
    return pts, summary


def main(argv) -> int:
    which = argv[0] if argv else "all"
    exp = BCExperiment(
        name=f"BC01_physical_H_{which}",
        question="Does removing the disc around zero change H or kappa?",
        config={"benchmarks": which, "seed": SEED, "workers": WORKERS},
    )
    started = time.time()
    runs = {"ieee39": ieee39, "kundur": kundur, "ieee68": ieee68}
    report = {}
    for name, fn in runs.items():
        if which not in (name, "all"):
            continue
        frame = fn(WORKERS)
        frame.to_csv(OUT / f"BC01_{name}_subsets.csv.gz", index=False)
        pts, summary = summarize(frame, name)
        pts.to_csv(OUT / f"BC01_{name}_points.csv", index=False)
        write_json(OUT / f"BC01_{name}_summary.json", summary)
        report[name] = {k: v for k, v in summary.items() if k != "label_changes"}
        print(json.dumps(summary, indent=1, default=str)[:4000], flush=True)
    exp.finish("COMPUTED", elapsed_s=time.time() - started, **report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
