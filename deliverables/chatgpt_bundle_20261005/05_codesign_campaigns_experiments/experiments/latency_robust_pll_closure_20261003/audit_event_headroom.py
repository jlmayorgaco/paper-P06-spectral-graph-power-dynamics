"""Audit the realized step-13 bus-16 actuator event-size bracket."""
from __future__ import annotations

import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent


def first(path: Path):
    with path.open(newline="") as stream:
        return next(csv.DictReader(stream))


cases = [
    ("frozen", 100.0, HERE / "event_validations/full20_step13_tiny/Q0_EVENT_METRICS.csv"),
    ("plus_0p01_MW", 100.01, HERE / "event_screens/uncertainty_step13_event_plus0p01MW/Q0_EVENT_METRICS.csv"),
    ("plus_0p02_MW", 100.02, HERE / "event_screens/uncertainty_step13_event_plus0p02MW/Q0_EVENT_METRICS.csv"),
    ("plus_1_MW", 101.0, HERE / "event_screens/uncertainty_step13_event_plus1MW/Q0_EVENT_METRICS.csv"),
]
result = []
for name, load_mw, path in cases:
    with path.open(newline="") as stream:
        rec = next(r for r in csv.DictReader(stream) if r["event"] == "bus16_plus100")
    result.append(dict(case=name, load_MW=load_mw,
                       actuator_slack=float(rec["min_SG_actuator_fraction_slack"]),
                       frozen_limit=0.002,
                       excess_above_limit=float(rec["min_SG_actuator_fraction_slack"]) - 0.002,
                       event_guard_pass=rec["pass"],
                       scope="ONE_EVENT_ONLY_FIXED_STEP13_DESIGN"))
with (HERE / "TABLE_F13_EVENT_SIZE_HEADROOM.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(result[0]))
    writer.writeheader()
    writer.writerows(result)
print("step13 event headroom: +0.01 MW", result[1]["event_guard_pass"],
      "+0.02 MW", result[2]["event_guard_pass"])

fig, ax = plt.subplots(figsize=(7, 4.2), dpi=150)
near = result[:3]
ax.plot([r["load_MW"] for r in near], [r["actuator_slack"] for r in near],
        color="#087e8b", linewidth=2)
for r in near:
    ax.scatter(r["load_MW"], r["actuator_slack"], s=80,
               color="#087e8b" if r["event_guard_pass"] == "true" else "#dc5b40", zorder=3)
ax.axhline(0.002, color="#14213d", linestyle="--", label="Frozen actuator limit")
ax.set_xticks([100, 100.01, 100.02])
ax.set_xlabel("Bus-16 load event magnitude (MW)")
ax.set_ylabel("Minimum normalized SG actuator slack")
ax.set_title("A 10–20 kW event-size bracket at the step-13 design")
ax.legend(loc="upper right")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(HERE / "FIG_F7_EVENT_SIZE_HEADROOM.png", dpi=220)
