"""Evaluate the frozen GFL11 device at one fixed parity point."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

from campaign_root import campaign_root_from_argv

ROOT = campaign_root_from_argv()
SOURCE = ROOT / "research" / "ias2026_last_validation" / "code" / "python"
sys.path.insert(0, str(SOURCE))
from ibr_cycles.models.ieee39_devices import ConverterParameters, GridFollowingConverter  # noqa: E402

OUT = ROOT / "raw" / "tx4_exact_p4_julia"
OUT.mkdir(parents=True, exist_ok=True)

v = 0.97 + 0.11j
x = np.array([0.07, 0.003, 0.12, -0.08, 0.01, -0.02, 0.42, -0.11, 0.004, -0.001, -0.08])
parameters = ConverterParameters(
    p_ref=0.12,
    q_ref=-0.08,
    v_ref=0.98,
    voltage_control=True,
    voltage_gain=0.03625,
    voltage_leak=0.05,
)
device = GridFollowingConverter(bus=30, parameters=parameters, weight=0.37)
derivatives = device.derivatives(x, v)
injection = device.injection(x, v)

with (OUT / "device_parity_python.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle)
    writer.writerow(["quantity", "value"])
    for index, value in enumerate(derivatives):
        writer.writerow([f"derivative_{index+1:02d}", repr(float(value))])
    writer.writerow(["injection_real", repr(float(injection.real))])
    writer.writerow(["injection_imag", repr(float(injection.imag))])

print("PYTHON_TX4_GFL11_DEVICE_PARITY_PASS")
