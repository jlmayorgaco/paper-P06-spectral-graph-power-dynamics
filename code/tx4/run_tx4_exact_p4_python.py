"""Run the exact frozen TX4 P4/GFL11 Python reference for the 16-case census."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

from campaign_root import campaign_root_from_argv


ROOT = campaign_root_from_argv()
SOURCE = ROOT / "research" / "ias2026_last_validation" / "code" / "python"
sys.path.insert(0, str(SOURCE))
from ibr_cycles.dynamics.linearize import central_difference_jacobians  # noqa: E402
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case  # noqa: E402
from ibr_cycles.models.ieee39_devices import ConverterParameters  # noqa: E402


OUT = ROOT / "raw" / "reconciliation"
OUT.mkdir(parents=True, exist_ok=True)
GAUGE_TOL = 1e-4
P4_G = 0.03625
P4_K = 1.425
P4_T = 1.5
P4_H = 1.0
TARGET_BUSES = (30, 33, 35, 37)

P4_CONVERTER = ConverterParameters(
    voltage_control=True,
    voltage_gain=P4_G,
    voltage_leak=0.05,
)
P4_MACHINE_SCALING = {"ka": P4_K, "ta": P4_T}


def run_case(label: str, mapping: dict[int, float]) -> dict[str, float | bool]:
    plan = ReplacementPlan.of(mapping, q_policy="matched", device="gfl")
    case = solve_case(
        plan,
        converter=P4_CONVERTER,
        machine_scaling=P4_MACHINE_SCALING,
        tol=1e-9,
    )
    values, vectors = np.linalg.eig(case.system.A)
    jac = central_difference_jacobians(
        case.dae, case.equilibrium.x, case.equilibrium.z, case.equilibrium.theta
    )
    voltage_modes = -np.linalg.solve(jac.gz, jac.gx @ vectors)
    voltages = case.dae.voltages(case.equilibrium.z)
    np.savez(
        OUT / f"python_{label}.npz",
        state=case.equilibrium.x,
        algebraic=case.equilibrium.z,
        matrix=case.system.A,
        eigenvalues=values,
        eigenvectors=vectors,
        voltage_modes=voltage_modes,
        voltages=voltages,
        state_labels=np.asarray(case.system.labels, dtype=str),
    )
    np.savetxt(OUT / f"python_{label}_state.csv", case.equilibrium.x, delimiter=",")
    np.savetxt(OUT / f"python_{label}_algebraic.csv", case.equilibrium.z, delimiter=",")
    with (OUT / f"python_{label}_spectrum.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["index", "real", "imag", "frequency_hz", "magnitude", "transverse"])
        for index, value in enumerate(values):
            writer.writerow([index, value.real, value.imag, abs(value.imag) / (2.0 * np.pi), abs(value), abs(value) >= GAUGE_TOL])
    transverse = values[np.abs(values) >= GAUGE_TOL]
    critical = transverse[np.argmax(transverse.real)] if transverse.size else values[np.argmax(values.real)]
    result: dict[str, float | bool] = {
        "states": float(case.n_states),
        "norm_f": float(case.equilibrium.norm_f),
        "norm_g": float(case.equilibrium.norm_g),
        "gz_condition": float(case.gz_condition),
        "alpha_full": float(values.real.max()),
        "alpha_transverse": float(transverse.real.max()) if transverse.size else float("nan"),
        "critical_frequency_hz": float(abs(critical.imag) / (2.0 * np.pi)),
        "stable": bool(transverse.real.max() < 0.0) if transverse.size else False,
    }
    return result


summary_path = OUT / "python_equilibrium_summary.csv"
if summary_path.exists():
    summary_path.unlink()

with summary_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle)
    writer.writerow(["case", "states", "norm_f", "norm_g", "gz_condition", "alpha_full", "alpha_transverse", "critical_frequency_hz", "stable"])
    for label, mapping in (
        ("base", {}),
        ("v4", {bus: 1.0 for bus in TARGET_BUSES}),
    ):
        result = run_case(label, mapping)
        writer.writerow([label, result["states"], result["norm_f"], result["norm_g"], result["gz_condition"], result["alpha_full"], result["alpha_transverse"], result["critical_frequency_hz"], result["stable"]])
        print(f"PYTHON_TX4_EXACT_P4_PASS case={label} states={int(result['states'])} alpha={result['alpha_full']:.6e} transverse={result['alpha_transverse']:.6e}")

with (OUT / "python_true_same_model_v4_census.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle)
    writer.writerow(["portfolio", "replaced_count", "states", "norm_f", "norm_g", "gz_condition", "alpha_full", "alpha_transverse", "critical_frequency_hz", "stable"])
    for mask in range(16):
        buses = tuple(TARGET_BUSES[index] for index in range(4) if (mask >> index) & 1)
        label = "census_" + ("none" if not buses else "+".join(str(bus) for bus in buses))
        result = run_case(label, {bus: 1.0 for bus in buses})
        writer.writerow([
            "none" if not buses else "+".join(str(bus) for bus in buses),
            len(buses), result["states"], result["norm_f"], result["norm_g"], result["gz_condition"],
            result["alpha_full"], result["alpha_transverse"], result["critical_frequency_hz"], result["stable"],
        ])
        print(f"PYTHON_TX4_EXACT_P4_PASS case={label} states={int(result['states'])} alpha={result['alpha_full']:.6e} transverse={result['alpha_transverse']:.6e}")
