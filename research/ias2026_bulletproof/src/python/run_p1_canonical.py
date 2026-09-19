"""Run the P1 oracle cases through the frozen canonical Python class.

The source imported here is the byte-for-byte snapshot in raw/gfl11/canonical_source.
This runner is deliberately separate from the independent equation transcription in
run_p1_gfl11.py so the evidence can distinguish specification parity from canonical
implementation parity.
"""

from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np


CAMPAIGN = Path(__file__).resolve().parents[2]
RAW = CAMPAIGN / "raw" / "gfl11"
SNAPSHOT = RAW / "canonical_source"
SOURCE = SNAPSHOT / "ieee39_devices.py"
ORACLE = RAW / "p1_python_oracle.csv"


def load_canonical():
    spec = importlib.util.spec_from_file_location("canonical_ieee39_devices", SOURCE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load canonical source: {SOURCE}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def err(a: float, b: float) -> float:
    return abs(a - b) / max(1.0, abs(a), abs(b))


def main() -> int:
    if not ORACLE.is_file():
        raise FileNotFoundError(ORACLE)
    canonical = load_canonical()
    rows = list(csv.DictReader(ORACLE.open(encoding="utf-8", newline="")))
    output: list[dict[str, float | int | str]] = []
    max_derivative = 0.0
    max_current = 0.0
    max_power = 0.0
    max_initialization = 0.0
    for row in rows:
        v = float(row["v"])
        angle = float(row["a"])
        p0 = float(row["p0"])
        q0 = float(row["q0"])
        weight = float(row["w"])
        gain = float(row["g"])
        x = np.array([float(row[f"x{i}"]) for i in range(11)], dtype=float)
        expected_xeq = np.array([float(row[f"x_eq{i}"]) for i in range(11)], dtype=float)
        expected_f = np.array([float(row[f"f{i}"]) for i in range(11)], dtype=float)

        params = canonical.ConverterParameters(
            p_ref=float(row["p_ref"]),
            q_ref=float(row["q_ref"]),
            v_ref=float(row["v_ref"]),
            voltage_control=True,
            voltage_gain=gain,
            voltage_leak=0.05,
        )
        device = canonical.GridFollowingConverter(bus=39, parameters=params, weight=weight)
        terminal = v * np.exp(1j * angle)
        actual_f = np.asarray(device.derivatives(x, terminal), dtype=float)
        current = complex(device.injection(x, terminal))
        power = terminal * np.conj(current)
        tuned, actual_xeq = device.initialize(terminal, p0 + 1j * q0)
        actual_xeq = np.asarray(actual_xeq, dtype=float)

        derivative_error = max(err(float(a), float(b)) for a, b in zip(actual_f, expected_f))
        current_error = max(err(current.real, float(row["i_r"])), err(current.imag, float(row["i_i"])))
        power_error = max(err(power.real, float(row["P_system"])), err(power.imag, float(row["Q_system"])))
        initialization_error = max(err(float(a), float(b)) for a, b in zip(actual_xeq, expected_xeq))
        max_derivative = max(max_derivative, derivative_error)
        max_current = max(max_current, current_error)
        max_power = max(max_power, power_error)
        max_initialization = max(max_initialization, initialization_error)
        output.append(
            {
                "case": int(row["case"]),
                "derivative_max_relative_error": derivative_error,
                "terminal_current_max_relative_error": current_error,
                "system_base_power_max_relative_error": power_error,
                "initialization_max_relative_error": initialization_error,
                "canonical_n_states": tuned.n_states,
                "canonical_i_r": current.real,
                "canonical_i_i": current.imag,
                "canonical_P_system": power.real,
                "canonical_Q_system": power.imag,
            }
        )

    out_path = RAW / "canonical_gfl_response.csv"
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)

    manifest_path = RAW / "canonical_source_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["canonical_evaluation"] = {
        "runner": "src/python/run_p1_canonical.py",
        "oracle_input": "raw/gfl11/p1_python_oracle.csv",
        "cases": len(rows),
        "class": "GridFollowingConverter",
        "function_set": ["derivatives", "injection", "initialize"],
        "n_states": 11,
        "max_derivative_relative_error": max_derivative,
        "max_terminal_current_relative_error": max_current,
        "max_system_base_power_relative_error": max_power,
        "max_initialization_relative_error": max_initialization,
        "status": "PASS" if max(max_derivative, max_current, max_power, max_initialization) < 1e-12 else "FAIL",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    passed = manifest["canonical_evaluation"]["status"] == "PASS"
    print(
        "P1_CANONICAL_PASS" if passed else "P1_CANONICAL_FAIL",
        f"cases={len(rows)} max_derivative={max_derivative:.3e}",
        f"max_current={max_current:.3e} max_power={max_power:.3e} max_init={max_initialization:.3e}",
    )
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
