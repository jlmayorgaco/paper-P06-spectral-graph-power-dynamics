"""Use an observed one-dimensional event secant to choose a guarded line step."""

from __future__ import annotations

import csv
import math
import tomllib
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
PARENT = "N_step1_multimode_lp"
FULL = f"{PARENT}_next_full"


def model(name: str) -> dict:
    return tomllib.loads((HERE / "designs" / f"{name}.toml").read_text())


def slack(name: str, folder: str) -> float:
    result = tomllib.loads((HERE / folder / name / "Q0_RESULT.toml").read_text())
    event = next(e for e in result["events"] if e["event"] == "bus16_plus100")
    assert event["complete"]
    return event["min_SG_actuator_fraction_slack"]


def main() -> None:
    parent_slack = slack(PARENT, "event_validations")
    full_slack = slack(FULL, "event_screens")
    assert parent_slack > 0.002 and full_slack < 0.002
    secant_slope = full_slack - parent_slack
    linear_boundary = (0.002 - parent_slack) / secant_slope
    # 5% directional buffer and a deterministic 1/16 line-search grid.
    scale = math.floor((0.95 * linear_boundary) * 16) / 16
    assert 0.5 < scale < 1.0
    base = model(PARENT)
    far = model(FULL)
    kp = np.asarray(base["Kp"]) * np.exp(scale * np.log(np.asarray(far["Kp"]) / np.asarray(base["Kp"])))
    ki = np.asarray(base["Ki"]) * np.exp(scale * np.log(np.asarray(far["Ki"]) / np.asarray(base["Ki"])))
    name = f"{PARENT}_next_{str(scale).replace('.', 'p')}"
    target = HERE / "designs" / f"{name}.toml"
    target.write_text(
        "rho = [" + ", ".join(repr(float(x)) for x in base["rho"]) + "]\n"
        + "Kp = [" + ", ".join(repr(float(x)) for x in kp) + "]\n"
        + "Ki = [" + ", ".join(repr(float(x)) for x in ki) + "]\n"
        + f'parent = "{PARENT}"\nline_search_scale = "{scale}"\n',
        encoding="utf-8",
    )
    row = {
        "parent": PARENT,
        "full_trial": FULL,
        "active_event": "bus16_plus100",
        "actuator_slack_parent": parent_slack,
        "actuator_slack_full": full_slack,
        "minimum_slack": 0.002,
        "observed_directional_secant_per_full_step": secant_slope,
        "linearized_max_scale": linear_boundary,
        "guard_fraction": 0.95,
        "grid_increment": 0.0625,
        "chosen_scale": scale,
        "chosen_candidate": name,
        "status": "LOCAL_EVENT_SECANT_PREDICTION_REQUIRES_EXACT_EVENT_CORRECTOR",
    }
    with (HERE / "L3_ACTIVE_EVENT_SECANT.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)
    print(name, "linear event scale", linear_boundary, "chosen", scale)


if __name__ == "__main__":
    main()
