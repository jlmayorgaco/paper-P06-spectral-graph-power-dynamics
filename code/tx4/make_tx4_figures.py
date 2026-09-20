"""Create the four preregistered TX4 validation figures."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import linear_sum_assignment

from campaign_root import campaign_root_from_argv

ROOT = campaign_root_from_argv()
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

with (RESULTS / "TX4_EXACT_P4_V4_CROSSCODE.csv").open(encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle))
labels = [row["portfolio"] for row in rows]
py_alpha = np.array([float(row["python_alpha"]) for row in rows])
jl_alpha = np.array([float(row["julia_alpha"]) for row in rows])
py_freq = np.array([float(row["python_frequency_hz"]) for row in rows])
jl_freq = np.array([float(row["julia_frequency_hz"]) for row in rows])

plt.figure(figsize=(7.2, 5.1))
plt.scatter(py_alpha, jl_alpha, s=32, color="#0f766e")
lo = min(py_alpha.min(), jl_alpha.min())
hi = max(py_alpha.max(), jl_alpha.max())
plt.plot([lo, hi], [lo, hi], "k--", linewidth=1, label="y=x")
for index, label in enumerate(labels):
    if label in {"30+33+35+37", "none"}:
        plt.annotate(label, (py_alpha[index], jl_alpha[index]), xytext=(5, 5), textcoords="offset points", fontsize=8)
plt.xlabel("Python alpha_perp (s$^{-1}$)")
plt.ylabel("Julia alpha_perp (s$^{-1}$)")
plt.title("F1. Exact P4 alpha census")
plt.grid(alpha=0.25)
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig(FIGURES / "TX4_F1_ALPHA_CROSSCODE.png", dpi=220)
plt.close()

positions = np.arange(len(labels))
width = 0.38
plt.figure(figsize=(10.0, 5.1))
plt.bar(positions - width / 2, py_freq, width, label="Python", color="#2563eb")
plt.bar(positions + width / 2, jl_freq, width, label="Julia", color="#f97316")
plt.xticks(positions, labels, rotation=55, ha="right", fontsize=8)
plt.ylabel("Critical frequency (Hz)")
plt.title("F2. Exact P4 critical-frequency census")
plt.grid(axis="y", alpha=0.25)
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig(FIGURES / "TX4_F2_FREQUENCY_CROSSCODE.png", dpi=220)
plt.close()

raw = ROOT / "raw" / "reconciliation"
py = np.load(raw / "python_census_30+33+35+37.npz")
jl_spec = np.genfromtxt(raw / "julia_census_30+33+35+37_spectrum.csv", delimiter=",", names=True)
jl_values = jl_spec["real"] + 1j * jl_spec["imag"]
py_modes = py["voltage_modes"]
jl_real = np.loadtxt(raw / "julia_census_30+33+35+37_voltage_modes_real.csv", delimiter=",")
jl_imag = np.loadtxt(raw / "julia_census_30+33+35+37_voltage_modes_imag.csv", delimiter=",")
jl_modes = jl_real[:, 1:].T + 1j * jl_imag[:, 1:].T
py_transverse = np.flatnonzero(abs(py["eigenvalues"]) >= 1e-4)
jl_transverse = np.flatnonzero(abs(jl_values) >= 1e-4)
py_critical = int(py_transverse[np.argmax(py["eigenvalues"][py_transverse].real)])
cost = abs(py["eigenvalues"][:, None] - jl_values[None, :])
assignment = linear_sum_assignment(cost)
mapping = {int(row): int(col) for row, col in zip(*assignment, strict=True)}
jl_critical = mapping[py_critical]
py_voltage_mode = py_modes[:, py_critical][0::2] + 1j * py_modes[:, py_critical][1::2]
jl_voltage_mode = jl_modes[:, jl_critical][0::2] + 1j * jl_modes[:, jl_critical][1::2]
py_shape = abs(py_voltage_mode) / abs(py_voltage_mode).max()
jl_shape = abs(jl_voltage_mode) / abs(jl_voltage_mode).max()
plt.figure(figsize=(9.0, 4.8))
bus_positions = np.arange(1, 40)
plt.plot(bus_positions, py_shape, "o-", ms=3, label="Python", color="#2563eb")
plt.plot(bus_positions, jl_shape, "s--", ms=3, label="Julia", color="#f97316")
plt.xlabel("IEEE-39 bus")
plt.ylabel("Normalized voltage-mode magnitude")
plt.title("F3. H4 critical voltage-mode shape (MAC >= 0.95 gate)")
plt.xticks(np.arange(1, 40, 3))
plt.grid(alpha=0.25)
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig(FIGURES / "TX4_F3_H4_VOLTAGE_MODE.png", dpi=220)
plt.close()

old_states = int(np.loadtxt(ROOT / "research" / "ias2026_last_validation" / "raw" / "reconciliation" / "julia_v4_state.csv").size)
exact_states = int(np.loadtxt(raw / "julia_census_30+33+35+37_state.csv").size)
plt.figure(figsize=(6.6, 4.7))
plt.bar(["Prior GFL10", "Exact GFL11"], [old_states, exact_states], color=["#94a3b8", "#0f766e"])
plt.ylabel("Dynamic state count, H4")
plt.title("F4. State-count correction")
for x, value in enumerate([old_states, exact_states]):
    plt.text(x, value + 0.8, str(value), ha="center", va="bottom", fontweight="bold")
plt.ylim(0, max(exact_states + 8, 20))
plt.grid(axis="y", alpha=0.25)
plt.tight_layout()
plt.savefig(FIGURES / "TX4_F4_STATE_COUNTS.png", dpi=220)
plt.close()

print("TX4_FIGURES_PASS")
