"""The IEEE-39 network must match an independent ANDES solution."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from ibr_cycles.models.ieee39_network import load_network, solve_power_flow

CONFIG = Path(__file__).resolve().parents[1] / "configs" / "ias2026"


@pytest.fixture(scope="module")
def reference():
    return json.loads(
        (CONFIG / "ieee39_andes_powerflow_reference.json").read_text(encoding="utf-8")
    )


def test_ybus_matches_andes(reference):
    network = load_network()
    lines_only = np.load(CONFIG / "ieee39_ybus_andes.npy")
    shunts = np.zeros_like(lines_only)
    for bus, value in ((4, 1.0), (5, 2.0)):
        position = network.position(bus)
        shunts[position, position] += 1j * value
    assert np.abs(network.ybus - (lines_only + shunts)).max() < 1e-9


def test_power_flow_matches_andes(reference):
    solution = solve_power_flow(load_network())
    assert solution.converged
    assert solution.max_mismatch < 1e-10
    assert np.abs(np.abs(solution.voltages) - np.array(reference["v"])).max() < 1e-5
    assert np.abs(np.angle(solution.voltages) - np.array(reference["a"])).max() < 1e-5


def test_replacement_candidates_exclude_the_system_equivalent():
    network = load_network()
    assert network.slack_bus == 39
    assert network.replacement_candidates == (30, 31, 32, 33, 34, 35, 36, 37, 38)
    assert network.machine_on_system_base(39)["M"] > 1000.0
