"""F7 post-processing: the questions of the post-F7 report.

Reads results/F7/{map}_points.csv, _leaves.csv, _areas.csv, _events.csv,
_islands.csv, _spotcheck.csv and writes results/F7/F7_summary.json plus the
tables quoted in docs/F7_POLICY_HYPERGRAPH.md. Nothing here re-solves a case.
"""

from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.stats import beta

from _bootstrap import RESULTS
from _f7_common import LABELS
from ibr_cycles.diagnosis.composability import (
    hypergraph_order,
    maximal_free_sets,
    parse_hypergraph_label,
)

F7 = RESULTS / "F7"
CORE = (30, 33, 35, 37)
BAD = {"BASE_UNSTABLE", "INFEASIBLE"}
Y = {"F7A": "k", "F7B": "t", "F7C": "h"}
CHAIN = ["30+33+35+37", "30+33+35", "30+33"]
DETECT = 1e-4


def points_file(name: str):
    """The map point table: the working CSV if present, else the frozen .csv.gz.

    The freeze stores the three large point tables gzip-compressed (lossless;
    SHA-256 of the uncompressed bytes in results/F7/SHA256SUMS_points.txt).
    """

    plain = F7 / f"{name}_points.csv"
    return plain if plain.exists() else F7 / f"{name}_points.csv.gz"


def kappa_of(label: str) -> int:
    if label in BAD:
        return -2
    return hypergraph_order(parse_hypergraph_label(label))


def witnesses(label: str) -> str:
    if label in BAD:
        return "-"
    edges = parse_hypergraph_label(label)
    k = hypergraph_order(edges)
    return (
        "|".join("+".join(map(str, sorted(e))) for e in edges if len(e) == k) or "none"
    )


def free_sets(label: str) -> str:
    return "|".join(
        "+".join(map(str, sorted(s)))
        for s in maximal_free_sets(parse_hypergraph_label(label), CORE)
    )


def upper95(errors: int, n: int) -> float:
    if n == 0:
        return float("nan")
    return 1.0 if errors >= n else float(beta.ppf(0.95, errors + 1, n - errors))


# ------------------------------------------------------------------ raster ----


def raster(points: pd.DataFrame, leaves: pd.DataFrame) -> np.ndarray:
    """Label of every lattice node: homogeneous leaves fill their rectangle."""

    ni, nj = int(points.I.max()) + 1, int(points.J.max()) + 1
    grid = np.full((ni, nj), "", dtype=object)
    lookup = {
        (int(i), int(j)): lab
        for i, j, lab in zip(points.I, points.J, points.label, strict=True)
    }
    for i, j, s, hom in leaves[["I", "J", "size", "homogeneous"]].itertuples(
        index=False
    ):
        if hom:
            grid[i : i + s + 1, j : j + s + 1] = lookup[(i, j)]
    for (i, j), lab in lookup.items():
        grid[i, j] = lab
    return grid


EIGHT = np.ones((3, 3), dtype=int)
SUBSTANTIAL = 25
"""A component counts as substantial when it holds at least this many lattice
nodes; smaller ones are slivers along nearly coincident boundaries, reported
but not counted as fragmentation."""


def _component_row(name, kappa, mask, total):
    lab, n = ndimage.label(mask, structure=EIGHT)
    sizes = sorted(np.bincount(lab.ravel())[1:].tolist(), reverse=True) if n else []
    return {
        "label": name,
        "kappa": kappa,
        "components": int(n),
        "substantial_components": int(sum(s >= SUBSTANTIAL for s in sizes)),
        "component_sizes": sizes[:12],
        "node_fraction": float(mask.sum() / total),
    }


def components(grid: np.ndarray) -> pd.DataFrame:
    """Connected components of every H region and every kappa level (8-connectivity)."""

    rows = []
    total = np.isin(grid, list(BAD) + [""], invert=True).sum()
    for label in sorted(set(grid.ravel()) - {""} - BAD):
        rows.append(_component_row(label, kappa_of(label), grid == label, total))
    kappa_grid = np.vectorize(lambda s: kappa_of(s) if s else -3)(grid)
    for k in sorted(set(kappa_grid.ravel()) - {-3, -2}):
        rows.append(_component_row(f"KAPPA={k}", k, kappa_grid == k, total))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ tables ----


def label_table(
    points: pd.DataFrame, areas: pd.DataFrame, comps: pd.DataFrame
) -> pd.DataFrame:
    stable = areas[~areas.label.isin(BAD)]
    total = stable.area.sum()
    by = stable.groupby("label").area.sum()
    boundary = stable[stable.boundary_cell].groupby("label").area.sum()
    ok = points[~points.label.isin(BAD)]
    defect = ok.groupby("label").closure_defect.apply(
        lambda s: int(s.fillna("").ne("").sum())
    )
    ncomp = comps.set_index("label").substantial_components
    rows = []
    for label in by.index:
        rows.append(
            {
                "label": label,
                "kappa": kappa_of(label),
                "n_edges": len(parse_hypergraph_label(label)),
                "kappa_witnesses": witnesses(label),
                "maximal_safe_portfolios": free_sets(label),
                "area_fraction": float(by[label] / total),
                "area_fraction_in_boundary_cells": float(
                    boundary.get(label, 0.0) / total
                ),
                "components": int(ncomp.get(label, 0)),
                "points": int((ok.label == label).sum()),
                "points_with_closure_defect": int(defect.get(label, 0)),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["kappa", "area_fraction"], ascending=[True, False]
    )


def event_summary(events: pd.DataFrame) -> dict:
    sub = events[events.kind == "SUBSET_CROSSING"].copy()
    sub["changes_H"] = sub.changes_H.astype(bool)
    moving = sub[sub.changes_H]
    imag = sub[sub.boundary_type.astype(str).str.startswith("IMAGINARY_AXIS")]
    multi = imag[imag.subset_size > 1]
    single = imag[imag.subset_size == 1]
    rhp_change = (sub.rhp_to - sub.rhp_from).abs()
    port = {}
    if len(imag):
        port = {
            "imaginary_axis_events": int(len(imag)),
            "port_visible": int(imag.port_visible.astype(bool).sum()),
            "port_invisible": int((~imag.port_visible.astype(bool)).sum()),
            "multi_port_events": int(len(multi)),
            "multi_port_detected_closure_lt_1e-4": int(
                (multi.closure_distance < DETECT).sum()
            ),
            "closure_distance_max": float(multi.closure_distance.max())
            if len(multi)
            else None,
            "closure_distance_median": float(multi.closure_distance.median())
            if len(multi)
            else None,
            "closure_min_refined_max": float(multi.closure_min_refined.max())
            if len(multi)
            else None,
            "freq_port_error_max_hz": float(
                (multi.freq_port_hz - multi.freq_star_hz).abs().max()
            )
            if len(multi)
            else None,
            "single_port_events": int(len(single)),
            "single_port_individual_factor_max": float(single.port_min_individual.max())
            if len(single)
            else None,
            "factorization_residual_max": float(imag.factorization_residual.max()),
            "solve_residual_max": float(imag.closure_solve_residual.max()),
            "direct_re_at_star_max_abs": float(imag.direct_re_at_star.abs().max()),
            "port_errors": int(
                imag.get("port_error", pd.Series(dtype=str)).notna().sum()
            ),
        }
    return {
        "subset_crossings": int(len(sub)),
        "H_changing_crossings": int(len(moving)),
        "boundary_types_all": sub.boundary_type.value_counts().to_dict(),
        "boundary_types_H_changing": moving.boundary_type.value_counts().to_dict(),
        "moves": moving.move.value_counts().to_dict(),
        "H_witness_subsets": moving.subset.value_counts().to_dict(),
        "unresolved": int((sub.boundary_type == "UNRESOLVED").sum()),
        "replay_mismatch": int(
            sub.get("replay_mismatch", pd.Series(dtype=object)).notna().sum()
        ),
        "rhp_change_by_type": {
            str(t): sorted(
                rhp_change[sub.boundary_type == t]
                .dropna()
                .astype(int)
                .unique()
                .tolist()
            )
            for t in sub.boundary_type.dropna().unique()
        },
        "freq_range_hz_by_type": {
            str(t): [float(g.freq_star_hz.min()), float(g.freq_star_hz.max())]
            for t, g in moving.groupby("boundary_type")
        },
        "tangencies": moving[moving.boundary_type == "IMAGINARY_AXIS_TANGENCY"][
            [
                "subset",
                "g_star",
                "k_star",
                "t_star",
                "h_star",
                "freq_star_hz",
                "fold_distance_edges",
                "H_before",
                "H_after",
            ]
        ].to_dict("records")[:30],
        "admissibility_events": int((events.kind == "ADMISSIBILITY").sum()),
        "admissibility_types": events[events.kind == "ADMISSIBILITY"]
        .boundary_type.value_counts()
        .to_dict(),
        "port": port,
    }


def transitions_table(events: pd.DataFrame) -> pd.DataFrame:
    sub = events[(events.kind == "SUBSET_CROSSING") & events.changes_H.astype(bool)]
    grouped = sub.groupby(
        ["H_before", "H_after", "subset", "boundary_type", "move", "edge_axis"]
    )
    return (
        grouped.agg(
            n=("subset", "size"),
            kappa_before=("kappa_before", "first"),
            kappa_after=("kappa_after", "first"),
            freq_min=("freq_star_hz", "min"),
            freq_max=("freq_star_hz", "max"),
            closure_max=("closure_distance", "max"),
        )
        .reset_index()
        .sort_values("n", ascending=False)
    )


def axis_moves(events: pd.DataFrame, axis: str) -> dict:
    """kappa transitions along one coordinate, oriented in the increasing direction."""

    sub = events[
        (events.kind == "SUBSET_CROSSING")
        & events.changes_H.astype(bool)
        & (events.edge_axis == axis)
    ]
    inf = lambda k: "inf" if k == -1 else str(int(k))  # noqa: E731
    pairs = (
        (sub.kappa_before.map(inf) + "->" + sub.kappa_after.map(inf))
        .value_counts()
        .to_dict()
    )
    unsafe_to_safe = int(((sub.kappa_before != -1) & (sub.kappa_after == -1)).sum())
    safe_to_unsafe = int(((sub.kappa_before == -1) & (sub.kappa_after != -1)).sum())
    return {
        "events": int(len(sub)),
        "kappa_pairs_increasing_coordinate": pairs,
        "unsafe_to_safe_increasing": unsafe_to_safe,
        "safe_to_unsafe_increasing": safe_to_unsafe,
        "H_changes_same_kappa": int((sub.kappa_before == sub.kappa_after).sum()),
    }


def lines(points: pd.DataFrame, along: str, across: str) -> pd.DataFrame:
    rows = []
    for value, group in points.groupby(across):
        group = group.sort_values(along)
        labels = group.label.tolist()
        seq = [lab for i, lab in enumerate(labels) if i == 0 or lab != labels[i - 1]]
        kap = [kappa_of(lab) for lab in seq]
        kseq = [k for i, k in enumerate(kap) if i == 0 or k != kap[i - 1]]
        rows.append(
            {
                across: value,
                "n_points": len(group),
                "H_sequence": " > ".join(seq),
                "kappa_sequence": " > ".join(
                    "inf" if k == -1 else ("U" if k == -2 else str(k)) for k in kseq
                ),
                "n_labels": len(set(labels) - BAD),
            }
        )
    return pd.DataFrame(rows)


def chain_check(points: pd.DataFrame, along: str) -> pd.DataFrame:
    """Along ``along`` at each coarse g: does the witness chain 4 > 3 > 2 appear?"""

    rows = []
    unit = 32
    for g, group in points[points.I % unit == 0].groupby("g"):
        group = points[points.g == g].sort_values(along)
        group = group[~group.label.isin(BAD)]
        wit = [witnesses(lab) for lab in group.label]
        seq = [w for i, w in enumerate(wit) if i == 0 or w != wit[i - 1]]
        position = []
        for target in CHAIN:
            hit = [i for i, w in enumerate(seq) if target in w.split("|")]
            position.append(hit[0] if hit else None)
        found = all(p is not None for p in position)
        forward = found and position == sorted(position)
        backward = found and position == sorted(position, reverse=True)
        rows.append(
            {
                "g": g,
                "kappa_witness_sequence": " > ".join(seq),
                "contains_4_3_2_chain_in_order": forward or backward,
                "chain_direction": "increasing"
                if forward
                else ("decreasing" if backward else ""),
            }
        )
    return pd.DataFrame(rows)


def flagship_paths(points: pd.DataFrame, y: str) -> pd.DataFrame:
    rows = []
    for value, group in points[~points.label.isin(BAD)].groupby(y):
        group = group.sort_values("g")
        if group.g.min() > 0 or len(group) < 3:
            continue
        flag = group["N_30+33+35+37"].to_numpy() > 0
        switches = int(np.sum(flag[1:] != flag[:-1]))
        rows.append(
            {
                y: value,
                "flagship_unsafe_at_g0": bool(flag[0]),
                "flagship_switches_along_g": switches,
                "first_g_flagship_safe": float(group.g[~flag].min())
                if (~flag).any()
                else None,
                "last_g_flagship_unsafe": float(group.g[flag].max())
                if flag.any()
                else None,
                "first_g_H_empty": float(group.g[group.label == "EMPTY"].min())
                if (group.label == "EMPTY").any()
                else None,
                "H_at_g0": group.iloc[0]["label"],
                "H_at_gmax": group.iloc[-1]["label"],
            }
        )
    return pd.DataFrame(rows)


def safeguard_a(spots: pd.DataFrame) -> dict:
    ok = spots[~spots.direct_label.isin(BAD)]
    return {
        "direct_points": int(len(ok)),
        "state_dimension_constant_per_subset": all(
            ok[f"nx_{lab}"].nunique() == 1 for lab in LABELS if f"nx_{lab}" in ok
        ),
        "artificial_zero_eigenvalues": int(
            sum(ok[f"zeros_{lab}"].sum() for lab in LABELS if f"zeros_{lab}" in ok)
        ),
        "gz_condition_max_relative_spread": float(
            max(
                np.ptp(ok[f"gzcond_{lab}"]) / ok[f"gzcond_{lab}"].median()
                for lab in LABELS
                if f"gzcond_{lab}" in ok
            )
        ),
    }


def main(argv) -> int:
    names = [a for a in argv if a in Y] or [
        n for n in Y if (F7 / f"{n}_events.csv").exists()
    ]
    summary = {}
    for name in names:
        y = Y[name]
        points = pd.read_csv(points_file(name), low_memory=False)
        leaves = pd.read_csv(F7 / f"{name}_leaves.csv")
        areas = pd.read_csv(F7 / f"{name}_areas.csv")
        events = pd.read_csv(F7 / f"{name}_events.csv", low_memory=False)
        islands = (
            pd.read_csv(F7 / f"{name}_islands.csv")
            if (F7 / f"{name}_islands.csv").stat().st_size > 2
            else pd.DataFrame()
        )
        spots = pd.read_csv(F7 / f"{name}_spotcheck.csv")
        ok = points[~points.label.isin(BAD)]
        grid = raster(points, leaves)
        comps = components(grid)
        comps.to_csv(F7 / f"{name}_components.csv", index=False)
        table = label_table(points, areas, comps)
        table.to_csv(F7 / f"{name}_labels.csv", index=False)
        transitions_table(events).to_csv(F7 / f"{name}_transitions.csv", index=False)
        lines(points, "g", y).to_csv(F7 / f"{name}_lines_along_g.csv", index=False)
        lines(points, y, "g").to_csv(F7 / f"{name}_lines_along_{y}.csv", index=False)
        paths = flagship_paths(points, y)
        paths.to_csv(F7 / f"{name}_flagship_paths.csv", index=False)
        chain = chain_check(points, y)
        chain.to_csv(F7 / f"{name}_chain_check.csv", index=False)
        scored = spots[spots.in_homogeneous_leaf]
        err = int((~scored.match).sum())
        near_err = int((~spots.nearest_match).sum())
        entry = {
            "points": int(len(points)),
            "stable_base_points": int(len(ok)),
            "kappa_values": sorted(int(v) for v in ok.kappa.unique()),
            "distinct_hypergraphs": int(ok.label.nunique()),
            "hypergraphs_per_kappa": {
                int(k): int(v) for k, v in ok.groupby("kappa").label.nunique().items()
            },
            "fragmented_labels": comps[
                (comps.substantial_components > 1)
                & ~comps.label.str.startswith("KAPPA")
            ][
                ["label", "components", "substantial_components", "component_sizes"]
            ].to_dict("records"),
            "kappa_components": comps[comps.label.str.startswith("KAPPA")][
                [
                    "kappa",
                    "components",
                    "substantial_components",
                    "component_sizes",
                    "node_fraction",
                ]
            ].to_dict("records"),
            "points_with_closure_defect": int(
                ok.closure_defect.fillna("").ne("").sum()
            ),
            "closure_defect_sets": sorted(
                set(";".join(ok.closure_defect.dropna()).split(";")) - {""}
            ),
            "island_or_tongue_cells": int(len(islands)),
            "island_examples": islands.head(10).to_dict("records")
            if len(islands)
            else [],
            "events": event_summary(events),
            "policy_moves_along_g": axis_moves(events, "g"),
            f"moves_along_{y}": axis_moves(events, y),
            "flagship_lines": int(len(paths)),
            "flagship_lines_unsafe_at_g0": int(paths.flagship_unsafe_at_g0.sum()),
            "flagship_lines_made_safe_by_g": int(
                paths.first_g_flagship_safe.notna().sum()
            ),
            "flagship_lines_with_2plus_switches": int(
                (paths.flagship_switches_along_g >= 2).sum()
            ),
            "chain_4_3_2_in_order_g_columns": int(
                chain.contains_4_3_2_chain_in_order.sum()
            ),
            "chain_g_columns": int(len(chain)),
            "spot_check": {
                "n": int(len(spots)),
                "in_homogeneous_leaves": int(len(scored)),
                "errors": err,
                "error_upper95": upper95(err, len(scored)),
                "boundary_cell_points": int((~spots.in_homogeneous_leaf).sum()),
                "nearest_node_errors_all_points": near_err,
                "nearest_node_error_upper95": upper95(near_err, len(spots)),
                "fast_vs_direct_label_disagreements": int(
                    (~spots.fast_equals_direct).sum()
                ),
                "mismatches": spots[~spots.match & spots.in_homogeneous_leaf][
                    ["x", "y", "predicted", "direct_label"]
                ].to_dict("records"),
            },
            "safeguard_A_direct_spot_points": safeguard_a(spots),
        }
        summary[name] = entry
        print(f"===== {name} =====")
        print(json.dumps(entry, indent=1, default=str)[:6000])
        print(table.to_string(index=False))
    (F7 / "F7_summary.json").write_text(
        json.dumps(summary, indent=2, default=str), encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
