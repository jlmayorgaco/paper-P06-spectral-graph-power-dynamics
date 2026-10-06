"""Make requested but unperformed comparisons visibly incomplete, not invented."""
from __future__ import annotations

import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent

TABLES = {
    "TABLE_F5_MULTISTART.csv": (
        ["experiment", "start", "optimized_tau_crit_ms", "same_basin_as_best", "status", "reason"],
        [["multistart", "Z / N / previous best / 3 seeded random starts", "", "",
          "NOT_EXECUTED", "Full event-constrained optimization per start was not completed"]],
    ),
    "TABLE_F6_GENERIC_OPTIMIZER_COMPARISON.csv": (
        ["experiment", "solver", "optimized_tau_crit_ms", "full_model_evaluations", "wall_time_s", "status", "reason"],
        [["controlled_generic_comparison", "", "", "", "", "NOT_EXECUTED",
          "No full-model event-constrained generic comparison was completed"]],
    ),
    "TABLE_F7_PARETO_FRONTIER.csv": (
        ["experiment", "sigma0_s_inv", "optimized_tau_crit_ms", "active_family", "status", "reason"],
        [["nominal_vs_latency_pareto", "", "", "", "NOT_EXECUTED",
          "Only observed design scatter is available; no sigma0 sweep optimization was completed"]],
    ),
}

for name, (header, data) in TABLES.items():
    path = HERE / name
    if path.exists():
        raise FileExistsError(f"Refusing to replace existing result: {path}")
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(data)
    print(name, "NOT_EXECUTED")
