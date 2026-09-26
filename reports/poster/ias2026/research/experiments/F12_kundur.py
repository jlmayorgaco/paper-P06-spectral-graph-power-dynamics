"""F12 - policy-dependent composability on a second documented benchmark (Kundur).

Implements configs/kundur/F12_preregistration.yaml exactly (committed before
this script produced any composability number): candidates {2, 3, 4}, the F7
converter and leaky Q/V policy, band 0.3-1.5 Hz, maps K12A (g x k) and K12B
(g x t) on 61 x 61 nodes, bisection on every grid edge whose unsafe-family
indicators differ, Safeguard B classification, port closure at imaginary-axis
crossings, and 200 random off-grid points per map on the direct path.
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

from _bootstrap import RESULTS, ROOT
from _overnight import choose_workers, pin_blas_threads
from _postfreeze import PostFreezeExperiment
from ibr_cycles.diagnosis.composability import (
    hypergraph_label,
    hypergraph_move,
    hypergraph_order,
    incompatibility_hypergraph,
)
from ibr_cycles.models.ieee39_case import (
    InfeasibleReplacement,
    ReplacementPlan,
    solve_case,
)
from ibr_cycles.models.ieee39_devices import ConverterParameters
from ibr_cycles.models.ieee39_network import load_network
from ibr_cycles.models.port_admittance import build_action_space

NAME = "F12_kundur"
OUT = RESULTS / "F12"
NETWORK = (ROOT / "configs" / "kundur" / "kundur_network.json").resolve()
PREREG = ROOT / "configs" / "kundur" / "F12_preregistration.yaml"
CANDIDATES = (2, 3, 4)
SUBSETS = [
    tuple(c) for k in range(len(CANDIDATES) + 1) for c in combinations(CANDIDATES, k)
]
BAND = (0.3, 1.5)
LEAK = 0.05
N = 61
G_GRID = np.linspace(0.0, 1.0, N) ** 2
Y_GRID = np.linspace(0.5, 2.0, N)
MAPS = {"K12A": "k", "K12B": "t"}
SPOTS = 200
SEED = 20260915
ITER = 30


def net():
    return load_network(NETWORK)


def theta_of(map_name, g, y):
    return {
        "g": float(g),
        "k": float(y) if MAPS[map_name] == "k" else 1.0,
        "t": float(y) if MAPS[map_name] == "t" else 1.0,
    }


def solve(members, th):
    return solve_case(
        ReplacementPlan.of({b: 1.0 for b in members}),
        network=net(),
        converter=ConverterParameters(
            voltage_control=True, voltage_gain=th["g"], voltage_leak=LEAK
        ),
        machine_scaling={"ka": th["k"], "ta": th["t"]},
    )


def in_gamma(v):
    f = np.abs(v.imag) / (2 * math.pi)
    return (
        (v.real > 0)
        & (v.imag >= 0)
        & (np.abs(v) > 1e-3)
        & (f >= BAND[0])
        & (f <= BAND[1])
    )


def evaluate(th) -> dict:
    row, counts = dict(th), {}
    for members in SUBSETS:
        try:
            values = np.linalg.eigvals(solve(members, th).system.A)
        except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError) as error:
            return {**row, "label": "INFEASIBLE", "reason": str(error)[:60]}
        name = "+".join(map(str, members)) or "BASE"
        row[f"N_{name}"] = int(np.count_nonzero(in_gamma(values)))
        f = np.abs(values.imag) / (2 * math.pi)
        band = (
            (values.imag >= 0)
            & (np.abs(values) > 1e-3)
            & (f >= BAND[0])
            & (f <= BAND[1])
        )
        row[f"bandmax_{name}"] = (
            float(values.real[band].max()) if band.any() else float("-inf")
        )
        counts[frozenset(members)] = row[f"N_{name}"]
        if not members:
            dyn = values[np.abs(values) > 1e-3]
            row["base_abscissa"] = float(dyn.real.max())
            if row["base_abscissa"] >= 0:
                return {**row, "label": "BASE_UNSTABLE"}
    edges = incompatibility_hypergraph(counts)
    return {**row, "label": hypergraph_label(edges), "kappa": hypergraph_order(edges)}


def _node(task):
    map_name, i, j = task
    return {
        "map": map_name,
        "i": i,
        "j": j,
        **evaluate(theta_of(map_name, G_GRID[i], Y_GRID[j])),
    }


def _indicator(row):
    return tuple(
        int(row.get(f"N_{'+'.join(map(str, s)) or 'BASE'}", 0) > 0) for s in SUBSETS
    )


def _lerp(a, b, s):
    return {k: a[k] + s * (b[k] - a[k]) for k in a}


def _violation(v):
    f = abs(v.imag) / (2 * math.pi)
    parts = [
        p
        for p, c in (
            ("IMAGINARY_AXIS", v.real <= 0),
            ("LOWER_BAND_EDGE", f < BAND[0]),
            ("UPPER_BAND_EDGE", f > BAND[1]),
        )
        if c
    ]
    return "+".join(parts) or "NONE"


def _edge(task):
    map_name, ta, tb, ra, rb = task
    out = []
    ia, ib = _indicator(ra), _indicator(rb)
    state = {frozenset(s): x for s, x in zip(SUBSETS, ia, strict=True)}
    for members, xa, xb in zip(SUBSETS, ia, ib, strict=True):
        if not members or xa == xb:
            continue
        lo, hi = 0.0, 1.0
        va = np.linalg.eigvals(solve(members, ta).system.A)
        vb = np.linalg.eigvals(solve(members, tb).system.A)
        for _ in range(ITER):
            mid = 0.5 * (lo + hi)
            vm = np.linalg.eigvals(solve(members, _lerp(ta, tb, mid)).system.A)
            if (np.count_nonzero(in_gamma(vm)) > 0) == bool(xa):
                lo, va = mid, vm
            else:
                hi, vb = mid, vm
        unsafe, safe = (va, vb) if xa else (vb, va)
        upper = safe[safe.imag >= 0]
        kind, star = "UNRESOLVED", None
        for lam in unsafe[in_gamma(unsafe)]:
            partner = upper[np.argmin(np.abs(upper - lam))]
            if _violation(partner) != "NONE":
                kind, star = _violation(partner), 0.5 * (lam + partner)
                break
        before = incompatibility_hypergraph(state)
        state[frozenset(members)] = xb
        after = incompatibility_hypergraph(state)
        row = {
            "map": map_name,
            "witness": "+".join(map(str, members)),
            "boundary_type": kind,
            **{f"{k}_star": v for k, v in _lerp(ta, tb, 0.5 * (lo + hi)).items()},
            "H_before": hypergraph_label(before),
            "H_after": hypergraph_label(after),
            "changes_H": hypergraph_label(before) != hypergraph_label(after),
            "move": hypergraph_move(before, after)["kind"],
            "edge_axis": "g" if ta["g"] != tb["g"] else MAPS[map_name],
        }
        if star is not None:
            row.update(
                freq_hz=float(abs(star.imag) / (2 * math.pi)),
                lambda_re=float(star.real),
            )
        if kind == "IMAGINARY_AXIS":
            th = _lerp(ta, tb, 0.5 * (lo + hi))
            space = build_action_space(solve((), th), solve(members, th), members)
            split = space.split(complex(0.0, abs(star.imag)))
            ref = np.median(
                [
                    abs(complex(space.split(complex(0, w))["full"]))
                    for w in np.linspace(
                        2 * math.pi * BAND[0], 2 * math.pi * BAND[1], 25
                    )
                ]
            )
            row.update(
                port_visible=bool(abs(complex(split["full"])) / ref < 1e-3),
                closure_distance=float(
                    np.min(np.abs(np.asarray(split["q_eigenvalues"]) + 1.0))
                )
                if len(members) > 1
                else float("nan"),
                individual_min=float(min(abs(d) for d in split["diagonals"])),
            )
        out.append(row)
    return out


def _spot(task):
    map_name, g, y = task
    return {"map": map_name, "g": g, "y": y, **evaluate(theta_of(map_name, g, y))}


def main(argv) -> int:
    pin_blas_threads()
    OUT.mkdir(parents=True, exist_ok=True)
    import hashlib

    prereg_sha = hashlib.sha256(PREREG.read_bytes()).hexdigest()
    workers = choose_workers()
    exp = PostFreezeExperiment(
        name=NAME,
        question="Is policy-dependent composability reproduced on Kundur?",
        config={
            "preregistration_sha256": prereg_sha,
            "candidates": CANDIDATES,
            "grid": N,
            "maps": MAPS,
        },
        workers=workers,
    )
    started = time.time()
    rng = np.random.default_rng(SEED)
    with Pool(workers, initializer=pin_blas_threads) as pool:
        nodes = pd.DataFrame(
            pool.map(
                _node,
                [(m, i, j) for m in MAPS for i in range(N) for j in range(N)],
                chunksize=8,
            )
        )
        tasks = []
        for m in MAPS:
            grid = {
                (r["i"], r["j"]): r for r in nodes[nodes["map"] == m].to_dict("records")
            }
            for (i, j), ra in grid.items():
                for nb in ((i + 1, j), (i, j + 1)):
                    if nb not in grid:
                        continue
                    rb = grid[nb]
                    bad = {"BASE_UNSTABLE", "INFEASIBLE"}
                    if (
                        ra["label"] in bad
                        or rb["label"] in bad
                        or _indicator(ra) == _indicator(rb)
                    ):
                        continue
                    tasks.append(
                        (
                            m,
                            theta_of(m, G_GRID[i], Y_GRID[j]),
                            theta_of(m, G_GRID[nb[0]], Y_GRID[nb[1]]),
                            ra,
                            rb,
                        )
                    )
        events = pd.DataFrame(
            [e for group in pool.map(_edge, tasks, chunksize=2) for e in group]
        )
        spot_tasks = [
            (m, float(rng.uniform(0, 1) ** 2), float(rng.uniform(0.5, 2.0)))
            for m in MAPS
            for _ in range(SPOTS)
        ]
        spots = pd.DataFrame(pool.map(_spot, spot_tasks, chunksize=4))
    nodes.to_csv(OUT / "F12_nodes.csv", index=False)
    events.to_csv(OUT / "F12_events.csv", index=False)
    # nearest-node label for every spot
    near = []
    for r in spots.itertuples():
        i = int(np.argmin(np.abs(G_GRID - r.g)))
        j = int(np.argmin(np.abs(Y_GRID - r.y)))
        near.append(
            nodes[(nodes["map"] == r.map) & (nodes.i == i) & (nodes.j == j)].label.iloc[
                0
            ]
        )
    spots["nearest_node_label"] = near
    spots["match"] = spots.label == spots.nearest_node_label
    spots.to_csv(OUT / "F12_spotcheck.csv", index=False)

    summary = {"preregistration_sha256": prereg_sha}
    for m in MAPS:
        ok = nodes[
            (nodes["map"] == m) & ~nodes.label.isin(["BASE_UNSTABLE", "INFEASIBLE"])
        ]
        ev = events[events["map"] == m] if len(events) else events
        lines = ok.groupby("j").label.nunique()
        pure_policy = (
            ev[
                (ev.edge_axis == "g")
                & ev.changes_H
                & (ev.boundary_type == "IMAGINARY_AXIS")
            ]
            if len(ev)
            else ev
        )
        summary[m] = {
            "nodes": int((nodes["map"] == m).sum()),
            "base_unstable": int(
                (nodes[(nodes["map"] == m)].label == "BASE_UNSTABLE").sum()
            ),
            "kappa_values": sorted(int(k) for k in ok.kappa.dropna().unique()),
            "distinct_H": int(ok.label.nunique()),
            "H_counts": ok.label.value_counts().to_dict(),
            "policy_lines_with_2plus_H": int((lines > 1).sum()),
            "policy_lines": int(len(lines)),
            "events": int(len(ev)),
            "boundary_types": ev.boundary_type.value_counts().to_dict()
            if len(ev)
            else {},
            "pure_policy_imaginary_axis_H_changes": int(len(pure_policy)),
            "port_visible": int(
                ev.get("port_visible", pd.Series(dtype=bool))
                .fillna(False)
                .astype(bool)
                .sum()
            )
            if len(ev)
            else 0,
            "closure_max": float(ev.closure_distance.max())
            if len(ev) and "closure_distance" in ev
            else None,
            "spot_mismatch": int((~spots[spots["map"] == m].match).sum()),
        }
        summary[m]["decision"] = (
            "REPRODUCED"
            if summary[m]["pure_policy_imaginary_axis_H_changes"] > 0
            else "NOT_REPRODUCED"
        )
    (OUT / "F12_summary.json").write_text(
        json.dumps(summary, indent=2, default=str), encoding="utf-8"
    )
    exp.finish(
        "COMPUTED",
        **{k: v for k, v in summary.items() if k != "preregistration_sha256"},
        elapsed_s=time.time() - started,
    )
    print(json.dumps(summary, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
