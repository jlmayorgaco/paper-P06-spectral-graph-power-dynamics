"""Preselect gain-box probes before action-space reconstruction is evaluated."""

from __future__ import annotations

import json
import math
import random
import tomllib
from pathlib import Path


HERE = Path(__file__).resolve().parent


def main() -> None:
    out = HERE / "L1_FROZEN_GAIN_PROBES.json"
    if out.exists():
        raise SystemExit("L1 probes already frozen")
    m = json.loads((HERE / "EXPERIMENT_MANIFEST.json").read_text())
    z = tomllib.loads((HERE / "designs/Z_zero_delay_tuned.toml").read_text())
    n = tomllib.loads((HERE / "designs/N_nominal.toml").read_text())
    rng = random.Random(20261003)
    kp_min, kp_max = m["gain_bounds"]["Kp_min"], m["gain_bounds"]["Kp_max"]
    ki_min, ki_max = m["gain_bounds"]["Ki_min"], m["gain_bounds"]["Ki_max"]

    def logdraw(low: float, high: float) -> float:
        return math.exp(math.log(low) + rng.random() * math.log(high / low))

    probes = [
        {"id": "Z", "Kp": z["Kp"], "Ki": z["Ki"], "meaning": "stored zero-delay tuned"},
        {"id": "N", "Kp": n["Kp"], "Ki": n["Ki"], "meaning": "stored nominal"},
    ]
    probes += [
        {
            "id": f"random_{j:02d}",
            "Kp": [logdraw(kp_min, kp_max) for _ in range(10)],
            "Ki": [logdraw(ki_min, ki_max) for _ in range(10)],
            "meaning": "gain-box admissible only; zero-delay events not screened",
        }
        for j in range(20)
    ]
    state = {
        "random_seed": 20261003,
        "distribution": "independent log-uniform in unchanged frozen gain box",
        "rho": z["rho"],
        "tau_ms": [0.0, 20.0, 37.38, 39.38, 40.0],
        "s_points": [
            {"real": -0.05, "frequency_hz": 0.017},
            {"real": -0.05, "frequency_hz": 5.26},
            {"real": -0.05, "frequency_hz": 5.64},
        ],
        "probes": probes,
    }
    out.write_text(json.dumps(state, indent=2) + "\n")
    print(f"frozen {len(probes)} gain probes")


if __name__ == "__main__":
    main()
