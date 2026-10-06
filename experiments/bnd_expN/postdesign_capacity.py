"""Post-design small-signal capacity extrapolation using frozen ExpG limits."""
from __future__ import annotations

import csv
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports" / "experiment_N"


def main() -> None:
    config = tomllib.loads((ROOT / "experiments" / "bnd_expG" / "configs" /
                            "DESIGN_SCENARIO_FROZEN.toml").read_text())
    with (OUT / "TABLE_N14_time_domain_validation.csv").open(newline="") as fh:
        source = list(csv.DictReader(fh))
    rows = []
    for r in source:
        rg = float(r["peak_RoCoF_per_MW"])
        fg = float(r["peak_frequency_per_MW"])
        rcap = float(config["rocof_limit_Hz_s"])/rg
        fcap = float(config["frequency_limit_Hz"])/fg
        rows.append({"event_bus": r["event_bus"],
                     "pulse_fraction": r["pulse_fraction"],
                     "unit_RoCoF_gain_Hz_s_per_MW": rg,
                     "unit_frequency_gain_Hz_per_MW": fg,
                     "ExpG_frozen_RoCoF_limit_Hz_s": config["rocof_limit_Hz_s"],
                     "ExpG_frozen_frequency_limit_Hz": config["frequency_limit_Hz"],
                     "delta_P_max_RoCoF_MW_linear_extrapolation": rcap,
                     "delta_P_max_frequency_MW_linear_extrapolation": fcap,
                     "delta_P_max_both_MW_linear_extrapolation": min(rcap, fcap),
                     "nonlinear_validation_at_capacity": "NOT_RUN",
                     "design_constraint_status": "POST_DESIGN_DIAGNOSTIC"})
    path = OUT / "TABLE_N14_disturbance_capacity.csv"
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    worst = min(rows, key=lambda r: r["delta_P_max_both_MW_linear_extrapolation"])
    print("WORST_TESTED_BUS", worst["event_bus"])
    print("DELTA_P_MAX_ROCOF_LINEAR", worst["delta_P_max_RoCoF_MW_linear_extrapolation"])
    print("DELTA_P_MAX_FREQUENCY_LINEAR", worst["delta_P_max_frequency_MW_linear_extrapolation"])


if __name__ == "__main__":
    main()
