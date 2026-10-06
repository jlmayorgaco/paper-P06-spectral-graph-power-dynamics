"""Freeze sparse active-event directional probes near accepted N1 design."""

from __future__ import annotations

import csv
import math
import tomllib
from pathlib import Path


HERE = Path(__file__).resolve().parent
BASE = "N_step1_multimode_lp"
COORDS = [("Kp", 30), ("Kp", 35), ("Kp", 36), ("Kp", 37), ("Kp", 39), ("Ki", 30), ("Ki", 35), ("Ki", 36)]
STEP_LOG = 0.002


def write(path: Path, rho: list[float], kp: list[float], ki: list[float]) -> None:
    path.write_text(
        "rho = [" + ", ".join(repr(x) for x in rho) + "]\n"
        + "Kp = [" + ", ".join(repr(x) for x in kp) + "]\n"
        + "Ki = [" + ", ".join(repr(x) for x in ki) + "]\n",
        encoding="utf-8",
    )


def main() -> None:
    out = HERE / "L3_FROZEN_EVENT_GRADIENT_PROBES.csv"
    if out.exists():
        raise SystemExit("event probes already frozen")
    base = tomllib.loads((HERE / "designs" / f"{BASE}.toml").read_text())
    rows = [{"candidate_id": BASE, "gain_kind": "BASE", "bus": "", "delta_log_gain": 0.0, "path": f"designs/{BASE}.toml"}]
    for kind, bus in COORDS:
        kp, ki = list(base["Kp"]), list(base["Ki"])
        (kp if kind == "Kp" else ki)[bus - 30] *= math.exp(STEP_LOG)
        name = f"{BASE}_probe_{kind}{bus}_plus"
        path = HERE / "designs" / f"{name}.toml"
        write(path, base["rho"], kp, ki)
        rows.append({"candidate_id": name, "gain_kind": kind, "bus": bus, "delta_log_gain": STEP_LOG, "path": f"designs/{name}.toml"})
    with out.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print("frozen", len(rows), "bus16-plus100 actuator probes")


if __name__ == "__main__":
    main()
