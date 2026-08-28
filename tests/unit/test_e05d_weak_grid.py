from __future__ import annotations

import sys
from pathlib import Path

import andes
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
E05D = ROOT / "experiments" / "tx3" / "E05D_weak_grid"
sys.path[:0] = [str(ROOT), str(E05D)]

from weak_grid import (  # noqa: E402
    CASE,
    EXPECTED_STRESSED_LINES,
    line_branch_records,
    short_circuit_ratios_from_ybus,
)


def test_external_grid_branch_rule_is_exact_and_topological() -> None:
    system = andes.load(str(CASE), setup=True, no_output=True)
    records = line_branch_records(system)
    included = tuple(record["line_id"] for record in records if record["included"])
    excluded = {record["line_id"]: record["reason"] for record in records if not record["included"]}
    assert included == EXPECTED_STRESSED_LINES
    assert len(records) == 46 and len(included) == 34 and len(excluded) == 12
    assert excluded["Line_45"] == "excluded_A8_branch_and_local_GFL37_step_up_interface"
    assert excluded["Line_44"] == "excluded_local_GFL_terminal_step_up_interface"
    assert excluded["Line_46"] == "excluded_local_GFL_terminal_step_up_interface"


def test_scr_matches_one_bus_analytical_thevenin_network() -> None:
    # Source bus 0 is grounded; POI bus 1 sees Z_th = 0.01+j0.10 pu.
    impedance = 0.01 + 0.10j
    admittance = 1.0 / impedance
    ybus = np.asarray([[admittance, -admittance], [-admittance, admittance]], dtype=complex)
    result = short_circuit_ratios_from_ybus(
        ybus, poi_indices=[1], source_indices=[0], converter_ratings_system_pu=[2.0], poi_voltages_pu=[1.0]
    )
    assert np.isclose(result[1], 1.0 / (abs(impedance) * 2.0), rtol=1e-13, atol=0.0)


def test_scr_scales_inversely_with_series_impedance_without_shunts() -> None:
    impedance = 0.02 + 0.20j
    admittance = 1.0 / impedance
    ybus = np.asarray([[admittance, -admittance], [-admittance, admittance]], dtype=complex)
    base = short_circuit_ratios_from_ybus(ybus, [1], [0], [1.0])[1]
    weakened = short_circuit_ratios_from_ybus(ybus / 2.0, [1], [0], [1.0])[1]
    assert np.isclose(weakened, 0.5 * base, rtol=1e-13, atol=0.0)
