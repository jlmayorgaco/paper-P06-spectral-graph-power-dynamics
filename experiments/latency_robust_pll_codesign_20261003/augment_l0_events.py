"""Join newly executed five-event metrics into the L0 spectral reproduction."""

from __future__ import annotations

import csv
import tomllib
from pathlib import Path


HERE = Path(__file__).resolve().parent


def main() -> None:
    path = HERE / "L0_REPRODUCTION.csv"
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert {r["design_id"] for r in rows} == {"Z", "N"}
    for row in rows:
        d = tomllib.loads(
            (HERE / "event_validations" / ("Z_zero_delay_tuned" if row["design_id"] == "Z" else "N_nominal") / "Q0_RESULT.toml").read_text()
        )
        events = d["events"]
        assert len(events) == 5 and all(e["complete"] for e in events)
        row["zero_delay_event_status"] = "PASS_ALL_FIVE" if all(e["pass"] for e in events) else "FAIL"
        row["equilibrium_residual_inf"] = d["equilibrium_residual_inf"]
        row["zero_delay_event_max_frequency_deviation_hz"] = max(e["F_peak_Hz"] for e in events)
        row["zero_delay_event_max_RoCoF_hz_s"] = max(e["RoCoF_peak_Hz_s"] for e in events)
        row["zero_delay_event_min_SG_actuator_slack"] = min(e["min_SG_actuator_fraction_slack"] for e in events)
        row["zero_delay_worst_frequency_event"] = max(events, key=lambda e: e["F_peak_Hz"])["event"]
        assert abs(float(row["alpha_zero_delay_s_inv"]) - d["critical_real_part_s_inv"]) < 1e-7
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print([(r["design_id"], r["zero_delay_event_status"], r["zero_delay_event_min_SG_actuator_slack"]) for r in rows])


if __name__ == "__main__":
    main()
