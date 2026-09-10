"""F7 - policy-dependent minimal spectral incompatibility hypergraph.

Three declared two-dimensional slices of theta = (g, k, t, h), see _f7_common:

    F7A   reactive-policy gain g  x  excitation-gain scale k           (t = 1.5, h = 1)
    F7B   reactive-policy gain g  x  excitation time-constant scale t  (k = 1,   h = 1)
    F7C   reactive-policy gain g  x  heterogeneity amplitude h         (k = 1,   t = 1)

Every point solves every subset of the core and is labelled by H_Gamma.

Refinement (Safeguard C). Every cell is sampled at its four corners, four edge
midpoints and centre. A cell is homogeneous only if all nine points share the
H_Gamma label, kappa and the full vector of protected counts N_Gamma(S), no
subset has a band mode within AXIS_TOL of the imaginary axis, and no unstable
mode is within EDGE_TOL_HZ of a band edge. Otherwise it is split, down to the
declared depth. A cell whose corners agree while an interior point differs is
recorded as ISLAND_OR_TONGUE. Hanging nodes evaluated by a neighbour reopen a
leaf they contradict.

Boundaries. On every finest-level lattice edge whose unsafe-family indicators
differ, each flipping subset is bisected (fast path) to its crossing and
classified (Safeguard B); imaginary-axis crossings are split into TRANSVERSAL
and TANGENCY and get direct-path port-closure diagnostics.

Validation. An independent random spot-check through the plotted domain is
evaluated on the DIRECT path and compared with the label the map predicts.

A 4-D Cartesian grid is never formed.
"""

from __future__ import annotations

import json
import math
import sys
import time
from dataclasses import dataclass, field
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.stats import beta

from _bootstrap import RESULTS
from _f7_common import (
    LABELS,
    LEAK,
    SUBSETS,
    Theta,
    eigenvalues,
    evaluate,
    fast_model,
    locate_event,
    port_diagnostics,
    te_geometric_mean,
)
from _overnight import Experiment, choose_workers, pin_blas_threads
from ibr_cycles.diagnosis.composability import (
    hypergraph_label,
    hypergraph_move,
    hypergraph_order,
    incompatibility_hypergraph,
)

NAME = "F7_policy_hypergraph"
OUT = RESULTS / "F7"
AXIS_TOL = 0.01
EDGE_TOL_HZ = 0.02
SPOT_CHECKS = 400
SPOT_SEED = 20260911
PASS_UPPER = 0.02
"""Spot-check pass rule, frozen before the maps were run: the 95 percent
Clopper-Pearson upper bound of the label-error rate inside homogeneous leaves
must be below 2 percent, and every mismatch is listed."""

G_COARSE = (
    0.0,
    0.01,
    0.02,
    0.03,
    0.04,
    0.05,
    0.06,
    0.08,
    0.10,
    0.12,
    0.14,
    0.17,
    0.20,
    0.25,
    0.30,
    0.40,
    0.50,
    0.65,
    0.80,
    1.00,
)
BAD = {"BASE_UNSTABLE", "INFEASIBLE"}


@dataclass(frozen=True)
class MapSpec:
    name: str
    y_name: str
    y_coarse: tuple[float, ...]
    fixed: dict = field(default_factory=dict)
    depth: int = 4
    x_name: str = "g"
    x_coarse: tuple[float, ...] = G_COARSE

    @property
    def unit(self) -> int:
        """Lattice units per coarse interval; the finest cell is two units."""
        return 2 ** (self.depth + 1)

    def _coordinate(self, grid, index):
        i0, r = divmod(index, self.unit)
        if i0 >= len(grid) - 1:
            return float(grid[-1])
        return float(grid[i0] + (r / self.unit) * (grid[i0 + 1] - grid[i0]))

    def _index(self, grid, value) -> float:
        i0 = int(
            np.clip(np.searchsorted(grid, value, side="right") - 1, 0, len(grid) - 2)
        )
        return (i0 + (value - grid[i0]) / (grid[i0 + 1] - grid[i0])) * self.unit

    def xy(self, key):
        return self._coordinate(self.x_coarse, key[0]), self._coordinate(
            self.y_coarse, key[1]
        )

    def lattice(self, x, y):
        return self._index(self.x_coarse, x), self._index(self.y_coarse, y)

    def theta_xy(self, x, y) -> Theta:
        return Theta(
            **{
                **{"g": 0.0, "k": 1.0, "t": 1.0, "h": 1.0},
                **self.fixed,
                self.x_name: x,
                self.y_name: y,
            }
        )

    def theta(self, key) -> Theta:
        return self.theta_xy(*self.xy(key))

    @property
    def scales(self) -> dict[str, float]:
        return {
            self.x_name: self.x_coarse[-1] - self.x_coarse[0],
            self.y_name: self.y_coarse[-1] - self.y_coarse[0],
        }


MAPS = {
    "F7A": MapSpec(
        "F7A", "k", tuple(np.round(np.arange(0.5, 2.31, 0.1), 3)), {"t": 1.5, "h": 1.0}
    ),
    "F7B": MapSpec(
        "F7B",
        "t",
        tuple(np.round(np.arange(0.5, 3.01, 0.125), 4)),
        {"k": 1.0, "h": 1.0},
    ),
    "F7C": MapSpec(
        "F7C", "h", tuple(np.round(np.arange(0.0, 2.01, 0.1), 3)), {"k": 1.0, "t": 1.0}
    ),
}


def _init():
    pin_blas_threads()
    fast_model()


def _fast(theta):
    return evaluate(theta, fast=True)


def _direct(theta):
    return evaluate(theta, fast=False)


# ------------------------------------------------------------- refinement ----


def _nine(cell):
    i, j, s = cell
    h = s // 2
    return [(i + a * h, j + b * h) for a in range(3) for b in range(3)]


def _corners(cell):
    i, j, s = cell
    return [(i, j), (i + s, j), (i, j + s), (i + s, j + s)]


def _closed(cell, points):
    i, j, s = cell
    return [
        (a, b)
        for a in range(i, i + s + 1)
        for b in range(j, j + s + 1)
        if (a, b) in points
    ]


def _near(row) -> bool:
    return row["axis_gap"] < AXIS_TOL or row["edge_gap_hz"] < EDGE_TOL_HZ


def refine_map(spec: MapSpec, pool, log) -> tuple[dict, list, list]:
    u = spec.unit
    nx, ny = len(spec.x_coarse), len(spec.y_coarse)
    points: dict = {}

    def run(keys):
        keys = [k for k in dict.fromkeys(keys) if k not in points]
        if not keys:
            return
        rows = pool.map(_fast, [spec.theta(k) for k in keys], chunksize=16)
        for key, row in zip(keys, rows, strict=True):
            row["I"], row["J"] = key
            points[key] = row

    pending = [(i * u, j * u, u) for i in range(nx - 1) for j in range(ny - 1)]
    leaves: dict = {}
    islands: list = []
    passes = 0
    while pending:
        run([k for cell in pending for k in _nine(cell)])
        split = []
        for cell in pending:
            keys = _closed(cell, points)
            signatures = {points[k]["signature"] for k in keys}
            near = any(_near(points[k]) for k in keys)
            homogeneous = len(signatures) == 1 and not near
            if (
                len(signatures) > 1
                and len({points[k]["signature"] for k in _corners(cell)}) == 1
            ):
                islands.append(
                    {
                        "I": cell[0],
                        "J": cell[1],
                        "size": cell[2],
                        "corner_label": points[cell[:2]]["label"],
                        "labels_inside": "/".join(
                            sorted({points[k]["label"] for k in keys})
                        ),
                    }
                )
            if homogeneous or cell[2] == 2:
                # at the finest size the near-axis monitor has done its job; the
                # leaf is mixed only if its labels actually differ
                leaves[cell] = len(signatures) == 1
                continue
            h = cell[2] // 2
            i, j = cell[:2]
            split += [(i, j, h), (i + h, j, h), (i, j + h, h), (i + h, j + h, h)]
        passes += 1
        # hanging nodes evaluated by a neighbour can contradict an accepted leaf
        reopened = [
            c
            for c, hom in leaves.items()
            if hom
            and c[2] > 2
            and len({points[k]["signature"] for k in _closed(c, points)}) > 1
        ]
        for c in reopened:
            del leaves[c]
        pending = split + reopened
        log(
            f"{spec.name}: pass {passes}: {len(split)} split, "
            f"{len(reopened)} reopened, "
            f"{len(points)} points, {len(leaves)} leaves"
        )
    return points, [(c, hom) for c, hom in leaves.items()], islands


# ------------------------------------------------------------ boundary edges --


def _unsafe(row) -> tuple:
    return tuple(int(row[f"N_{lab}"] > 0) for lab in LABELS)


def boundary_edges(points, leaves) -> list:
    edges = set()
    for (i, j, s), hom in leaves:
        if hom:
            continue
        for a in range(i, i + s + 1):
            for b in range(j, j + s + 1):
                for nb in ((a + 1, b), (a, b + 1)):
                    if (
                        nb[0] > i + s
                        or nb[1] > j + s
                        or nb not in points
                        or (a, b) not in points
                    ):
                        continue
                    p, q = points[(a, b)], points[nb]
                    bad = p["label"] in BAD or q["label"] in BAD
                    if (bad and p["label"] != q["label"]) or (
                        not bad and _unsafe(p) != _unsafe(q)
                    ):
                        edges.add(((a, b), nb))
    return sorted(edges)


def _edge_task(task) -> list[dict]:
    spec, key_a, key_b, row_a, row_b = task
    theta_a, theta_b = spec.theta(key_a), spec.theta(key_b)
    axis = spec.x_name if key_a[1] == key_b[1] else spec.y_name
    common = {
        "map": spec.name,
        "I_a": key_a[0],
        "J_a": key_a[1],
        "I_b": key_b[0],
        "J_b": key_b[1],
        "edge_axis": axis,
        "label_a": row_a["label"],
        "label_b": row_b["label"],
    }
    if row_a["label"] in BAD or row_b["label"] in BAD:
        return [{**common, "kind": "ADMISSIBILITY", **_base_event(theta_a, theta_b)}]

    ind_a, ind_b = _unsafe(row_a), _unsafe(row_b)
    flips = [s for s, x, y in zip(SUBSETS, ind_a, ind_b, strict=True) if s and x != y]
    events = [
        {
            **locate_event(
                m, theta_a, theta_b, scales=spec.scales, axes=(spec.x_name, spec.y_name)
            ),
            "subset": "+".join(map(str, m)),
            "subset_size": len(m),
        }
        for m in flips
    ]
    state = {frozenset(s): x for s, x in zip(SUBSETS, ind_a, strict=True)}
    current = incompatibility_hypergraph(state)
    out = []
    ordered = sorted(events, key=lambda e: e.get("s_star", 0.0))
    for n, event in enumerate(ordered):
        members = frozenset(int(m) for m in event["subset"].split("+"))
        state[members] = 1 - state[members]
        after = incompatibility_hypergraph(state)
        simultaneous = bool(
            n > 0
            and abs(event.get("s_star", 0) - ordered[n - 1].get("s_star", 9)) < 1e-6
        )
        event.update(
            {
                **common,
                "kind": "SUBSET_CROSSING",
                "H_before": hypergraph_label(current),
                "H_after": hypergraph_label(after),
                "kappa_before": hypergraph_order(current),
                "kappa_after": hypergraph_order(after),
                "changes_H": hypergraph_label(current) != hypergraph_label(after),
                "move": hypergraph_move(current, after)["kind"],
                "simultaneous": simultaneous,
            }
        )
        if simultaneous:
            event["boundary_type"] = "CORNER_OR_MULTIPLE"
        if str(event.get("boundary_type", "")).startswith("IMAGINARY_AXIS"):
            theta_star = Theta(
                event["g_star"], event["k_star"], event["t_star"], event["h_star"]
            )
            try:
                event.update(
                    port_diagnostics(
                        tuple(sorted(members)), theta_star, abs(event["lambda_im"])
                    )
                )
            except (np.linalg.LinAlgError, ValueError) as error:
                event["port_error"] = str(error)[:80]
        current = after
        out.append(event)
    if hypergraph_label(current) != row_b["label"]:
        for event in out:
            event["replay_mismatch"] = True
    return out


def _base_event(theta_a, theta_b) -> dict:
    """Bisect the base-case stability boundary and report its crossing mode."""

    def unstable(theta):
        values = eigenvalues((), theta)
        dynamic = values[np.abs(values) > 1e-3]
        return bool(dynamic.real.max() >= 0.0), dynamic

    ua, va = unstable(theta_a)
    ub, vb = unstable(theta_b)
    if ua == ub:
        return {
            "status": "INFEASIBILITY_EDGE",
            "boundary_type": "DAE_OR_EQUILIBRIUM_FAILURE",
        }
    lo, hi = 0.0, 1.0
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        um, vm = unstable(theta_a.lerp(theta_b, mid))
        if um == ua:
            lo, va = mid, vm
        else:
            hi, vb = mid, vm
    worst = va if ua else vb
    lam = worst[np.argmax(worst.real)]
    theta_star = theta_a.lerp(theta_b, 0.5 * (lo + hi))
    return {
        "status": "LOCATED",
        "subset": "BASE",
        "boundary_type": "BASE_STABILITY_BOUNDARY",
        **{f"{k}_star": v for k, v in theta_star.as_dict().items()},
        "freq_star_hz": float(abs(lam.imag) / (2 * np.pi)),
        "lambda_re": float(lam.real),
        "lambda_im": float(abs(lam.imag)),
    }


# ------------------------------------------------------------------ areas ----

WEIGHTS = {
    (0, 0): 1,
    (0, 2): 1,
    (2, 0): 1,
    (2, 2): 1,
    (1, 0): 2,
    (0, 1): 2,
    (2, 1): 2,
    (1, 2): 2,
    (1, 1): 4,
}


def areas(spec: MapSpec, points, leaves) -> pd.DataFrame:
    rows = []
    for (i, j, s), hom in leaves:
        x0, y0 = spec.xy((i, j))
        x1, y1 = spec.xy((i + s, j + s))
        area = (x1 - x0) * (y1 - y0)
        h = s // 2
        for (a, b), w in WEIGHTS.items():
            key = (i + a * h, j + b * h)
            rows.append(
                {
                    "label": points[key]["label"],
                    "area": area * w / 16.0,
                    "boundary_cell": not hom,
                }
            )
    frame = pd.DataFrame(rows)
    return frame.groupby(["label", "boundary_cell"], as_index=False).area.sum()


# ------------------------------------------------------------- spot-check ----


def predicted_label(spec, leaves, points, x, y):
    """(leaf label or BOUNDARY_CELL, homogeneous?, nearest-node label)."""

    fi, fj = spec.lattice(x, y)
    for (i, j, s), hom in leaves:
        if i <= fi <= i + s and j <= fj <= j + s:
            h = s // 2
            nodes = [(i + a * h, j + b * h) for a in range(3) for b in range(3)]
            nearest = min(nodes, key=lambda k: (k[0] - fi) ** 2 + (k[1] - fj) ** 2)
            if hom:
                return points[(i, j)]["label"], True, points[nearest]["label"]
            return "BOUNDARY_CELL", False, points[nearest]["label"]
    return "OUTSIDE", False, "OUTSIDE"


def spot_check(spec: MapSpec, points, leaves, pool, rng) -> pd.DataFrame:
    x0, x1 = spec.x_coarse[0], spec.x_coarse[-1]
    y0, y1 = spec.y_coarse[0], spec.y_coarse[-1]
    half = SPOT_CHECKS // 2
    xs = np.concatenate(
        [
            rng.uniform(x0, x1, half),
            rng.uniform(math.sqrt(x0), math.sqrt(x1), SPOT_CHECKS - half) ** 2,
        ]
    )
    ys = rng.uniform(y0, y1, SPOT_CHECKS)
    thetas = [spec.theta_xy(float(x), float(y)) for x, y in zip(xs, ys, strict=True)]
    direct = pool.map(_direct, thetas, chunksize=2)
    fast = pool.map(_fast, thetas, chunksize=8)
    rows = []
    for n, (x, y, d, f) in enumerate(zip(xs, ys, direct, fast, strict=True)):
        label, homogeneous, nearest = predicted_label(spec, leaves, points, x, y)
        rows.append(
            {
                "x": x,
                "y": y,
                "sampling": "uniform" if n < half else "sqrt_g",
                "predicted": label,
                "in_homogeneous_leaf": homogeneous,
                "nearest_node_label": nearest,
                "nearest_match": nearest == d["label"],
                "direct_label": d["label"],
                "fast_label": f["label"],
                "match": label == d["label"],
                "fast_equals_direct": f["label"] == d["label"],
                **{
                    c: d.get(c) for c in d if c.startswith(("nx_", "zeros_", "gzcond_"))
                },
                "direct_axis_gap": d.get("axis_gap"),
            }
        )
    return pd.DataFrame(rows)


def clopper_pearson_upper(errors: int, n: int, level: float = 0.95) -> float:
    if n == 0:
        return float("nan")
    return 1.0 if errors >= n else float(beta.ppf(level, errors + 1, n - errors))


# ------------------------------------------------------------------- main ----


def main(argv) -> int:
    pin_blas_threads()
    names = [a for a in argv if a in MAPS] or list(MAPS)
    OUT.mkdir(parents=True, exist_ok=True)
    workers = choose_workers()
    experiment = Experiment(
        name=NAME,
        question="Does physical policy change the minimal spectral "
        "incompatibility hypergraph?",
        config={
            "maps": {
                n: {
                    "y": MAPS[n].y_name,
                    "y_range": [MAPS[n].y_coarse[0], MAPS[n].y_coarse[-1]],
                    "y_step": MAPS[n].y_coarse[1] - MAPS[n].y_coarse[0],
                    "fixed": MAPS[n].fixed,
                    "depth": MAPS[n].depth,
                }
                for n in names
            },
            "g_coarse": list(G_COARSE),
            "regulator": {"kp_v": 2.0, "ki_v": 20.0, "leak_rad_s": LEAK},
            "te_geometric_mean_s": te_geometric_mean(),
            "gamma": "Re s > 0 and 0.3 <= |Im s|/2pi <= 1.5 Hz, "
            "one representative per pair",
            "tested_family": "power set of the core 30, 33, 35, 37",
            "refinement": {
                "rule": "nine-point homogeneity on label, kappa and N vector",
                "axis_tol": AXIS_TOL,
                "edge_tol_hz": EDGE_TOL_HZ,
            },
            "spot_check": {
                "n_per_map": SPOT_CHECKS,
                "seed": SPOT_SEED,
                "path": "direct",
                "pass_rule": f"95% Clopper-Pearson upper bound < {PASS_UPPER}",
            },
            "operating_point": "matched dispatch; invariant along every coordinate",
        },
        workers=workers,
    )
    started = time.time()
    summary = {}
    rng = np.random.default_rng(SPOT_SEED)
    with Pool(processes=workers, initializer=_init) as pool:
        for name in names:
            spec = MAPS[name]
            t0 = time.time()
            points, leaves, islands = refine_map(spec, pool, experiment.note)
            frame = pd.DataFrame(points.values())
            frame["x"] = [spec.xy((r["I"], r["J"]))[0] for r in points.values()]
            frame["y"] = [spec.xy((r["I"], r["J"]))[1] for r in points.values()]
            frame.to_csv(OUT / f"{name}_points.csv", index=False)
            pd.DataFrame(
                [
                    {"I": c[0], "J": c[1], "size": c[2], "homogeneous": h}
                    for c, h in leaves
                ]
            ).to_csv(OUT / f"{name}_leaves.csv", index=False)
            pd.DataFrame(islands).to_csv(OUT / f"{name}_islands.csv", index=False)
            areas(spec, points, leaves).to_csv(OUT / f"{name}_areas.csv", index=False)
            edges = boundary_edges(points, leaves)
            experiment.note(
                f"{name}: {len(points)} points, {len(islands)} island/tongue cells, "
                f"{len(edges)} boundary edges"
            )
            print(flush=True)
            tasks = [(spec, a, b, points[a], points[b]) for a, b in edges]
            results = pool.map(_edge_task, tasks, chunksize=4)
            events = pd.DataFrame([e for group in results for e in group])
            events.to_csv(OUT / f"{name}_events.csv", index=False)
            spots = spot_check(spec, points, leaves, pool, rng)
            spots.to_csv(OUT / f"{name}_spotcheck.csv", index=False)
            scored = spots[spots.in_homogeneous_leaf]
            errors = int((~scored.match).sum())
            summary[name] = {
                "points": len(points),
                "leaves": len(leaves),
                "islands": len(islands),
                "boundary_edges": len(edges),
                "events": len(events),
                "spot_checks": len(spots),
                "spot_scored": len(scored),
                "spot_errors": errors,
                "spot_boundary_cells": int((~spots.in_homogeneous_leaf).sum()),
                "spot_nearest_node_errors_all": int((~spots.nearest_match).sum()),
                "spot_error_upper95": clopper_pearson_upper(errors, len(scored)),
                "fast_direct_disagreements": int((~spots.fast_equals_direct).sum()),
                "elapsed_s": round(time.time() - t0, 1),
            }
            summary[name]["spot_check_pass"] = bool(
                summary[name]["spot_error_upper95"] < PASS_UPPER
            )
            experiment.note(f"{name}: {json.dumps(summary[name])}")
            print(flush=True)
    (OUT / "F7_run_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    experiment.finish(
        "COMPUTED"
        if all(v["spot_check_pass"] for v in summary.values())
        else "SPOT_CHECK_FAILED",
        maps=summary,
        elapsed_s=time.time() - started,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
