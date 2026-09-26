# ruff: noqa: E501  -- export records kept on one line
"""EMT-OP - canonical TX4 operating points and reference data for the EMT adapter (tx3-analysis).

Committed with the EMT preregistration, before any TX4 EMT run. Writes
results/EMT_PRED/tx4_operating_points.json with:
  - the canonical power-flow bus voltages (complex, per bus) of the frozen network and of the
    12 frozen holdout-line networks scaled by 1.5 (EMT12); matched dispatch makes the power
    flow identical for every portfolio, condenser, governor and parameter draw;
  - generator complex power per generator bus and the loads (system base);
  - the canonical 60-Hz Ybus (for the EMT01 comparison).
The EMT side reads device parameters from the configs directly (configs/ias2026/*.json).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

import numpy as np  # noqa: E402
from EMT_phasor_predictions import HOLDOUT, scaled_network  # noqa: E402

import _bootstrap  # noqa: E402,F401
from ibr_cycles.models.ieee39_network import (  # noqa: E402
    load_network,
    solve_power_flow,
)

RESEARCH = HERE.parents[1]
OUT = RESEARCH / "results" / "EMT_PRED" / "tx4_operating_points.json"


def c(z):
    return [float(np.real(z)), float(np.imag(z))]


def point(net):
    flow = solve_power_flow(net)
    assert flow.converged, "power flow"
    v = {str(b): c(flow.voltages[net.position(b)]) for b in net.bus_idx}
    gen = {
        str(b): c(flow.injection(b, net.ybus) + net.loads.get(b, 0j))
        for b in net.generator_buses
    }
    return {"voltages": v, "generation": gen, "max_mismatch": float(flow.max_mismatch)}


def main() -> int:
    net = load_network()
    data = {
        "system_base_mva": 100.0,
        "frequency_hz": 60.0,
        "bus_order": [int(b) for b in net.bus_idx],
        "loads": {str(b): c(s) for b, s in net.loads.items()},
        "ybus_real": np.real(net.ybus).tolist(),
        "ybus_imag": np.imag(net.ybus).tolist(),
        "base": point(net),
        "lines_x1p5": {},
    }
    for li in json.loads(HOLDOUT.read_text())["lines"]:
        data["lines_x1p5"][str(int(li))] = point(scaled_network(int(li), 1.5))
    OUT.write_text(json.dumps(data, indent=1), encoding="utf-8")
    print("wrote", OUT.name, "base mismatch", data["base"]["max_mismatch"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
