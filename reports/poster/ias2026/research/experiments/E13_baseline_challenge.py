"""E13 - PHASE D. Conventional baseline challenge (R3).

Question
    Do conventional screening metrics predict which replacement portfolios fail,
    and if so, does the interaction framework add anything?

Independent variables
    every portfolio already evaluated by E12.

Predictors under test
    minimum and mean SCR over the replaced buses, generalized SCR of the
    replaced set, maximum multi-infeed interaction factor, replaced MW, replaced
    inertia fraction, and two order-truncated predictions built from the census
    itself: the additive and the pairwise reconstruction of the spectral
    abscissa.

Target
    instability of the portfolio.

Why the truncated reconstructions are included
    They are the strongest possible "lower-order screening" a planner could
    build from exhaustive single and pair studies. If they predicted the outcome,
    no higher-order theory would be needed, and that is the honest bar to clear.

Status of the result
    NUMERICAL OBSERVATION at one operating point, one converter tuning, rho = 1.

Usage
    python experiments/E13_baseline_challenge.py
"""

from __future__ import annotations

import argparse
from itertools import combinations

import numpy as np
import pandas as pd

from _bootstrap import RESULTS
from ibr_cycles.cycles.mobius import decompose
from ibr_cycles.diagnosis.baselines import (
    generalized_scr,
    interaction_factors,
    nodal_metrics,
)
from ibr_cycles.io.manifest import Manifest
from ibr_cycles.models.ieee39_network import load_network, solve_power_flow

EXPERIMENT = "E13_baseline_challenge"
CENSUS = "E12_compatibility_census_portfolios.csv"

PREDICTORS = (
    ("min_scr", -1.0),
    ("mean_scr", -1.0),
    ("gscr", -1.0),
    ("max_miif", +1.0),
    ("replaced_mw", +1.0),
    ("replaced_inertia_fraction", +1.0),
    ("additive_prediction", +1.0),
    ("pairwise_prediction", +1.0),
)


def roc_auc(scores: np.ndarray, labels: np.ndarray) -> float:
    """Mann-Whitney AUC; 0.5 means the predictor carries no information."""

    positive = scores[labels]
    negative = scores[~labels]
    if positive.size == 0 or negative.size == 0:
        return float("nan")
    order = np.argsort(np.concatenate([positive, negative]), kind="mergesort")
    ranks = np.empty(order.size, dtype=float)
    ranks[order] = np.arange(1, order.size + 1)
    values = np.concatenate([positive, negative])
    for value in np.unique(values):
        mask = values == value
        ranks[mask] = ranks[mask].mean()
    rank_sum = ranks[: positive.size].sum()
    return float(
        (rank_sum - positive.size * (positive.size + 1) / 2)
        / (positive.size * negative.size)
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="PHASE D baseline challenge")
    parser.add_argument("--seed", type=int, default=20260909)
    args = parser.parse_args()

    network = load_network()
    power_flow = solve_power_flow(network)
    metrics = nodal_metrics(network, power_flow)
    candidates = network.replacement_candidates
    miif = interaction_factors(network, candidates)
    position = {bus: i for i, bus in enumerate(candidates)}

    census = pd.read_csv(RESULTS / "tables" / CENSUS)
    census["members_tuple"] = census["members"].apply(
        lambda s: () if s == "BASE" else tuple(int(b) for b in s.split("+"))
    )
    abscissa = {
        frozenset(map(str, row.members_tuple)): float(row.spectral_abscissa)
        for row in census.itertuples()
    }
    total_inertia = sum(
        network.machine_on_system_base(b)["M"] for b in network.generator_buses
    )

    rows = []
    for row in census.itertuples():
        members = row.members_tuple
        if len(members) < 2:
            continue
        names = tuple(map(str, members))
        decomposition = decompose(lambda s: abscissa[frozenset(s)], names)
        scr = [metrics[b].scr for b in members]
        pairs = [
            miif[position[a], position[b]] for a, b in combinations(members, 2)
        ] or [0.0]
        rows.append(
            {
                "members": row.members,
                "size": len(members),
                "unstable": bool(row.unstable),
                "spectral_abscissa": row.spectral_abscissa,
                "min_scr": min(scr),
                "mean_scr": float(np.mean(scr)),
                "gscr": generalized_scr(network, members, power_flow),
                "max_miif": max(pairs),
                "replaced_mw": row.replaced_mw,
                "replaced_inertia_fraction": sum(
                    network.machine_on_system_base(b)["M"] for b in members
                )
                / total_inertia,
                "additive_prediction": decomposition.truncated(names, 1),
                "pairwise_prediction": decomposition.truncated(names, 2),
            }
        )
    table = pd.DataFrame(rows)

    results = []
    for size in sorted(table["size"].unique()):
        subset = table[table["size"] == size]
        labels = subset.unstable.to_numpy()
        if labels.all() or not labels.any():
            continue
        for name, sign in PREDICTORS:
            results.append(
                {
                    "size": int(size),
                    "predictor": name,
                    "auc": roc_auc(sign * subset[name].to_numpy(float), labels),
                    "cases": len(subset),
                    "unstable": int(labels.sum()),
                }
            )
    auc_table = pd.DataFrame(results)

    # False negatives of the best conventional screen at the minimum failing order.
    minimum_order = int(table[table.unstable]["size"].min())
    critical = table[table["size"] == minimum_order].copy()
    conventional = [
        n
        for n, _ in PREDICTORS
        if n not in ("additive_prediction", "pairwise_prediction")
    ]
    screening = []
    for name in conventional:
        sign = dict(PREDICTORS)[name]
        scores = sign * critical[name].to_numpy(float)
        k = int(critical.unstable.sum())
        flagged = set(np.argsort(-scores)[:k])
        actual = set(np.flatnonzero(critical.unstable.to_numpy()))
        screening.append(
            {
                "predictor": name,
                "precision_at_k": len(flagged & actual) / max(k, 1),
                "recall_at_k": len(flagged & actual) / max(len(actual), 1),
                "false_negatives": len(actual - flagged),
                "k": k,
            }
        )
    screen_table = pd.DataFrame(screening)

    tables = RESULTS / "tables"
    table.to_csv(tables / f"{EXPERIMENT}_predictors.csv", index=False)
    auc_table.to_csv(tables / f"{EXPERIMENT}_auc.csv", index=False)
    screen_table.to_csv(tables / f"{EXPERIMENT}_screening.csv", index=False)

    manifest = Manifest(
        experiment=EXPERIMENT,
        seed=args.seed,
        config={"census": CENSUS, "predictors": [n for n, _ in PREDICTORS]},
    )
    manifest.finish(
        "SUCCESS",
        minimum_failing_order=minimum_order,
        auc_by_size=auc_table.to_dict("records"),
        screening_at_minimum_order=screen_table.to_dict("records"),
        best_conventional_auc_at_minimum_order=float(
            auc_table[
                (auc_table["size"] == minimum_order)
                & (auc_table.predictor.isin(conventional))
            ].auc.max()
        ),
        pairwise_auc_at_minimum_order=float(
            auc_table[
                (auc_table["size"] == minimum_order)
                & (auc_table.predictor == "pairwise_prediction")
            ].auc.iloc[0]
        ),
    )
    path = manifest.write(RESULTS / "manifests")

    print("minimum failing portfolio order: %d" % minimum_order)
    print()
    print("ROC AUC for predicting instability (0.5 = no information)")
    pivot = auc_table.pivot(index="predictor", columns="size", values="auc")
    print(pivot.to_string(float_format=lambda v: f"{v:7.3f}", na_rep="   -"))
    print()
    print(
        "screening at the minimum failing order (size %d, %d unstable of %d)"
        % (minimum_order, int(critical.unstable.sum()), len(critical))
    )
    print(screen_table.to_string(index=False, float_format=lambda v: f"{v:8.3f}"))
    print()
    print("the ten failing portfolios at the minimum order, against the screens:")
    view = critical[critical.unstable].sort_values("spectral_abscissa", ascending=False)
    print(
        view[
            [
                "members",
                "spectral_abscissa",
                "min_scr",
                "gscr",
                "max_miif",
                "replaced_mw",
                "additive_prediction",
                "pairwise_prediction",
            ]
        ].to_string(index=False, float_format=lambda v: f"{v:10.4f}")
    )
    print()
    print("for comparison, the ten SAFEST portfolios at the same order:")
    safe = critical[~critical.unstable].nsmallest(10, "spectral_abscissa")
    print(
        safe[
            [
                "members",
                "spectral_abscissa",
                "min_scr",
                "gscr",
                "max_miif",
                "replaced_mw",
            ]
        ].to_string(index=False, float_format=lambda v: f"{v:10.4f}")
    )
    print()
    print(f"manifest -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
