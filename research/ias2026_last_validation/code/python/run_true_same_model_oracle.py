"""Run the frozen Python oracle for base and full V4 reconciliation."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

from campaign_root import campaign_root_from_argv


ROOT = campaign_root_from_argv()
sys.path.insert(0, str(Path(__file__).resolve().parent))
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case  # noqa: E402


OUT = ROOT / "raw" / "reconciliation"
OUT.mkdir(parents=True, exist_ok=True)
GAUGE_TOL = 1e-4


def run_case(label: str, mapping: dict[int, float]) -> None:
    case = solve_case(ReplacementPlan.of(mapping), tol=1e-9)
    values, vectors = np.linalg.eig(case.system.A)
    np.savez(
        OUT / f"python_{label}.npz",
        state=case.equilibrium.x,
        algebraic=case.equilibrium.z,
        matrix=case.system.A,
        eigenvalues=values,
        eigenvectors=vectors,
    )
    np.savetxt(OUT / f"python_{label}_state.csv", case.equilibrium.x, delimiter=",")
    np.savetxt(OUT / f"python_{label}_algebraic.csv", case.equilibrium.z, delimiter=",")
    with (OUT / f"python_{label}_spectrum.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["index", "real", "imag", "frequency_hz", "magnitude", "transverse"])
        for index, value in enumerate(values):
            writer.writerow([index, value.real, value.imag, abs(value.imag) / (2.0 * np.pi), abs(value), abs(value) >= GAUGE_TOL])
    with (OUT / "python_equilibrium_summary.csv").open("a", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        if handle.tell() == 0:
            writer.writerow(["case", "states", "norm_f", "norm_g", "gz_condition", "alpha_full", "alpha_transverse", "critical_frequency_hz"])
        transverse = values[np.abs(values) >= GAUGE_TOL]
        critical = transverse[np.argmax(transverse.real)] if transverse.size else values[np.argmax(values.real)]
        writer.writerow([label, case.n_states, case.equilibrium.norm_f, case.equilibrium.norm_g, case.gz_condition, values.real.max(), transverse.real.max() if transverse.size else np.nan, abs(critical.imag) / (2.0 * np.pi)])
    print(f"PYTHON_TRUE_SAME_MODEL_PASS case={label} states={case.n_states} alpha={values.real.max():.6e}")


run_case("base", {})
run_case("v4", {30: 1.0, 33: 1.0, 35: 1.0, 37: 1.0})

target_buses = (30, 33, 35, 37)
with (OUT / "python_true_same_model_v4_census.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle)
    writer.writerow(["portfolio", "replaced_count", "states", "alpha_full", "alpha_transverse", "critical_frequency_hz", "stable"])
    for mask in range(16):
        buses = tuple(target_buses[index] for index in range(4) if (mask >> index) & 1)
        label = "census_" + ("none" if not buses else "+".join(str(bus) for bus in buses))
        run_case(label, {bus: 1.0 for bus in buses})
        values = np.load(OUT / f"python_{label}.npz")["eigenvalues"]
        transverse = values[np.abs(values) >= GAUGE_TOL]
        critical = transverse[np.argmax(transverse.real)]
        writer.writerow(["none" if not buses else "+".join(str(bus) for bus in buses), len(buses), len(np.load(OUT / f"python_{label}.npz")["state"]), values.real.max(), transverse.real.max(), abs(critical.imag) / (2.0 * np.pi), transverse.real.max() < 0.0])
