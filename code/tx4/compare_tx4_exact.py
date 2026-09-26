"""Reconcile the exact Python and Julia P4/GFL11 results."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

from campaign_root import campaign_root_from_argv

ROOT = campaign_root_from_argv()
RAW = ROOT / "raw" / "reconciliation"
RESULTS = ROOT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)
GAUGE_TOL = 1e-4
TARGET_BUSES = (30, 33, 35, 37)
PORTFOLIOS = []
for mask in range(16):
    buses = tuple(TARGET_BUSES[index] for index in range(4) if (mask >> index) & 1)
    PORTFOLIOS.append(("none" if not buses else "+".join(map(str, buses)), buses))


def spectrum(path: Path) -> np.ndarray:
    data = np.genfromtxt(path, delimiter=",", names=True)
    return np.asarray(data["real"], float) + 1j * np.asarray(data["imag"], float)


def modes(path_real: Path, path_imag: Path) -> np.ndarray:
    real = np.loadtxt(path_real, delimiter=",")
    imag = np.loadtxt(path_imag, delimiter=",")
    return real[:, 1:].T + 1j * imag[:, 1:].T


def voltages_from_z(z: np.ndarray) -> np.ndarray:
    return z[0::2] + 1j * z[1::2]


def angle_delta(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return np.angle(np.exp(1j * (np.angle(left) - np.angle(right))))


def mac(left: np.ndarray, right: np.ndarray) -> float:
    numerator = abs(np.vdot(left, right)) ** 2
    denominator = float(np.vdot(left, left).real * np.vdot(right, right).real)
    return float(numerator / denominator) if denominator else 0.0


def load_case(portfolio: str) -> dict[str, np.ndarray]:
    label = "census_" + portfolio
    py = np.load(RAW / f"python_{label}.npz")
    return {
        "py_values": py["eigenvalues"],
        "py_voltage_modes": py["voltage_modes"],
        "py_state": py["state"],
        "py_z": py["algebraic"],
        "jl_values": spectrum(RAW / f"julia_{label}_spectrum.csv"),
        "jl_voltage_modes": modes(RAW / f"julia_{label}_voltage_modes_real.csv", RAW / f"julia_{label}_voltage_modes_imag.csv"),
        "jl_state": np.loadtxt(RAW / f"julia_{label}_state.csv"),
        "jl_z": np.loadtxt(RAW / f"julia_{label}_algebraic.csv"),
    }


def metrics(portfolio: str, buses: tuple[int, ...]) -> dict[str, float | str | bool]:
    data = load_case(portfolio)
    py_values = data["py_values"]
    jl_values = data["jl_values"]
    assignment_cost = abs(py_values[:, None] - jl_values[None, :])
    rows, cols = linear_sum_assignment(assignment_cost)
    mapping = {int(row): int(col) for row, col in zip(rows, cols, strict=True)}
    py_transverse = np.flatnonzero(abs(py_values) >= GAUGE_TOL)
    jl_transverse = np.flatnonzero(abs(jl_values) >= GAUGE_TOL)
    py_critical = int(py_transverse[np.argmax(py_values[py_transverse].real)])
    jl_critical = mapping[py_critical]
    py_alpha = float(py_values[py_transverse].real.max())
    jl_alpha = float(jl_values[jl_transverse].real.max())
    py_critical_freq = float(abs(py_values[py_critical].imag) / (2.0 * np.pi))
    jl_critical_freq = float(abs(jl_values[jl_critical].imag) / (2.0 * np.pi))
    py_v = voltages_from_z(data["py_z"])
    jl_v = voltages_from_z(data["jl_z"])
    critical_mac = mac(data["py_voltage_modes"][:, py_critical], data["jl_voltage_modes"][:, jl_critical])
    transverse_set = set(map(int, jl_transverse))
    transverse_macs = [
        mac(data["py_voltage_modes"][:, int(row)], data["jl_voltage_modes"][:, mapping[int(row)]])
        for row in py_transverse
        if int(row) in mapping and mapping[int(row)] in transverse_set
    ]
    py_stable = bool(py_alpha < 0.0)
    jl_stable = bool(jl_alpha < 0.0)
    return {
        "portfolio": portfolio,
        "replaced_count": len(buses),
        "python_states": len(data["py_state"]),
        "julia_states": len(data["jl_state"]),
        "state_count_match": len(data["py_state"]) == len(data["jl_state"]),
        "python_alpha": py_alpha,
        "julia_alpha": jl_alpha,
        "alpha_diff": abs(py_alpha - jl_alpha),
        "python_frequency_hz": py_critical_freq,
        "julia_frequency_hz": jl_critical_freq,
        "frequency_diff_hz": abs(py_critical_freq - jl_critical_freq),
        "python_stable": py_stable,
        "julia_stable": jl_stable,
        "verdict_identical": py_stable == jl_stable,
        "critical_voltage_mac": critical_mac,
        "mean_transverse_voltage_mac": float(np.mean(transverse_macs)) if transverse_macs else float("nan"),
        "state_max_abs_diff": float(np.max(np.abs(data["py_state"] - data["jl_state"]))),
        "algebraic_max_abs_diff": float(np.max(np.abs(data["py_z"] - data["jl_z"]))),
        "voltage_magnitude_max_diff": float(np.max(np.abs(np.abs(py_v) - np.abs(jl_v)))),
        "voltage_angle_max_diff": float(np.max(np.abs(angle_delta(py_v, jl_v)))),
        "max_eigenvalue_pair_diff": float(np.max(assignment_cost[rows, cols])),
    }


crosscode = [metrics(portfolio, buses) for portfolio, buses in PORTFOLIOS]
with (RESULTS / "TX4_EXACT_P4_V4_CROSSCODE.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(crosscode[0]))
    writer.writeheader()
    writer.writerows(crosscode)

h4 = next(row for row in crosscode if row["portfolio"] == "30+33+35+37")
with (RESULTS / "TX4_EXACT_P4_MODE_MATCH.csv").open("w", newline="", encoding="utf-8") as handle:
    fields = ["portfolio", "critical_voltage_mac", "mean_transverse_voltage_mac", "frequency_diff_hz", "max_eigenvalue_pair_diff"]
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    writer.writerow({field: h4[field] for field in fields})

old_gfl10 = ROOT / "research" / "ias2026_last_validation" / "raw" / "reconciliation" / "julia_v4_state.csv"
state_rows = [
    {"case": "old_gfl10_H4", "model": "prior independent Julia GFL10", "states": int(np.loadtxt(old_gfl10).size)},
    {"case": "exact_gfl11_H4", "model": "this branch exact Julia GFL11", "states": int(h4["julia_states"])},
]
with (RESULTS / "TX4_EXACT_P4_STATE_COUNTS.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["case", "model", "states"])
    writer.writeheader()
    writer.writerows(state_rows)

base = next(row for row in crosscode if row["portfolio"] == "none")
with (RESULTS / "TX4_EXACT_P4_ALLSG_PARITY.csv").open("w", newline="", encoding="utf-8") as handle:
    fields = ["case", "python_states", "julia_states", "alpha_diff", "frequency_diff_hz", "state_max_abs_diff", "algebraic_max_abs_diff", "voltage_magnitude_max_diff", "voltage_angle_max_diff", "verdict_identical", "pass"]
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    base["case"] = "all_sg_base"
    base["pass"] = bool(base["state_count_match"] and base["verdict_identical"] and base["alpha_diff"] <= 1e-4 and base["frequency_diff_hz"] <= 1e-4 and base["voltage_magnitude_max_diff"] <= 1e-6 and base["voltage_angle_max_diff"] <= 1e-6)
    writer.writerow({field: base.get(field) for field in fields})

device_root = ROOT / "raw" / "tx4_exact_p4_julia"
with (device_root / "device_parity_python.csv").open(encoding="utf-8") as handle:
    device_py = {row["quantity"]: float(row["value"]) for row in csv.DictReader(handle)}
with (device_root / "device_parity_julia.csv").open(encoding="utf-8") as handle:
    device_jl = {row["quantity"]: float(row["value"]) for row in csv.DictReader(handle)}
device_rows = [{"quantity": key, "python": device_py[key], "julia": device_jl[key], "abs_diff": abs(device_py[key] - device_jl[key])} for key in device_py]
with (RESULTS / "TX4_GFL11_DEVICE_PARITY.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["quantity", "python", "julia", "abs_diff"])
    writer.writeheader()
    writer.writerows(device_rows)

proper = [row for row in crosscode if row["portfolio"] != "30+33+35+37"]
device_pass = max(row["abs_diff"] for row in device_rows) <= 1e-12
gates = [
    {"claim_id": "C1", "status": "PASS" if device_pass else "FAIL", "evidence": "TX4_GFL11_DEVICE_PARITY.csv"},
    {"claim_id": "C2", "status": "PASS", "evidence": "TX4_P4_SG_POLICY_AUDIT.md"},
    {"claim_id": "C3", "status": "PASS" if all(row["verdict_identical"] for row in crosscode) else "FAIL", "evidence": "TX4_EXACT_P4_V4_CROSSCODE.csv"},
    {"claim_id": "C4", "status": "PASS" if h4["julia_states"] == 86 and h4["python_states"] == 86 else "FAIL", "evidence": "TX4_EXACT_P4_STATE_COUNTS.csv"},
    {"claim_id": "C5", "status": "PASS" if h4["alpha_diff"] <= 1e-4 and h4["frequency_diff_hz"] <= 1e-4 else "FAIL", "evidence": "TX4_EXACT_P4_V4_CROSSCODE.csv"},
    {"claim_id": "C6", "status": "PASS" if h4["critical_voltage_mac"] >= 0.95 else "FAIL", "evidence": "TX4_EXACT_P4_MODE_MATCH.csv"},
    {"claim_id": "C7", "status": "PASS" if len(proper) == 15 and all(row["python_stable"] and row["julia_stable"] for row in proper) else "FAIL", "evidence": "TX4_EXACT_P4_V4_CROSSCODE.csv"},
    {"claim_id": "C8", "status": "PASS", "evidence": "TX4_EXACT_P4_JULIA_DEVIATIONS.md"},
]
with (RESULTS / "TX4_EXACT_P4_CLAIMS.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["claim_id", "status", "evidence"])
    writer.writeheader()
    writer.writerows(gates)

headline = {
    "branch": "research/tx4-exact-p4-julia-reproduction",
    "h4": h4,
    "proper_subsets_stable_python": sum(bool(row["python_stable"]) for row in proper),
    "proper_subsets_stable_julia": sum(bool(row["julia_stable"]) for row in proper),
    "all_16_verdicts_identical": all(bool(row["verdict_identical"]) for row in crosscode),
    "device_parity_max_abs_diff": max(row["abs_diff"] for row in device_rows),
    "all_sg_parity_pass": bool(base["pass"]),
    "gates": gates,
}
(RESULTS / "TX4_EXACT_P4_HEADLINE.json").write_text(json.dumps(headline, indent=2, default=str) + "\n", encoding="utf-8")

print(json.dumps(headline, indent=2, default=str))
