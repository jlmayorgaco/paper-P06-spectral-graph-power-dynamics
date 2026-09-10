"""The v2 dispatcher must conserve power, respect ratings and be reproducible.

These run before any v2 campaign. The v1 sampler had none of them and shipped a
dispatch that violated active-power limits on four buses.
"""

from __future__ import annotations

import numpy as np
import pytest

from ibr_cycles.models.ieee39_network import load_network, solve_power_flow
from ibr_cycles.uncertainty.sampling import (
    LOSS_ALLOWANCE,
    audit_limits,
    sample_operating_point,
    stratified_grid,
)

CORE = (30, 33, 35, 37)
GRID = [
    (0.90, 1.00, 0.75),
    (1.00, 1.00, 1.00),
    (1.05, 0.95, 0.90),
    (1.10, 1.05, 1.00),
    (0.85, 1.10, 0.70),
]


def draw(load, reactive, availability, seed=0):
    return sample_operating_point(
        load_network(),
        active_load=load,
        reactive_load=reactive,
        availability=availability,
        availability_buses=CORE,
        rng=np.random.default_rng(seed),
    )


@pytest.mark.parametrize("load,reactive,availability", GRID)
def test_active_limits_are_respected(load, reactive, availability):
    point = draw(load, reactive, availability)
    network = load_network()
    for bus, spec in point.network.pv.items():
        assert spec["p"] <= network.pv[bus]["pmax"] + 1e-9, bus
        assert spec["p"] >= network.pv[bus]["pmin"] - 1e-9, bus


@pytest.mark.parametrize("load,reactive,availability", GRID)
def test_power_is_conserved(load, reactive, availability):
    """Scheduled generation plus the slack estimate must cover load and losses."""

    point = draw(load, reactive, availability)
    if not point.dispatchable:
        pytest.skip("sample is not dispatchable, which is a reported outcome")
    assert abs(point.power_balance_residual_pu) < 1e-8
    assert abs(point.unserved_pu) < 1e-6


@pytest.mark.parametrize("load,reactive,availability", GRID)
def test_availability_caps_the_resource_not_the_schedule(load, reactive, availability):
    network = load_network()
    point = draw(load, reactive, availability)
    for bus in CORE:
        cap = network.pv[bus]["pmax"] * availability
        assert point.network.pv[bus]["p"] <= min(network.pv[bus]["p"], cap) + 1e-9


def test_the_slack_is_not_the_default_balancer():
    """At nominal conditions the slack must stay near its scheduled output."""

    network = load_network()
    point = draw(1.00, 1.00, 1.00)
    nominal_slack = (
        solve_power_flow(network).injection(network.slack_bus, network.ybus)
        + network.loads.get(network.slack_bus, 0j)
    ).real
    assert abs(point.slack_estimate_pu - nominal_slack) < 0.35 * abs(nominal_slack)


def test_sampling_is_reproducible():
    first = draw(1.03, 0.98, 0.88, seed=7)
    second = draw(1.03, 0.98, 0.88, seed=7)
    assert first.dispatched == second.dispatched
    assert first.slack_estimate_pu == second.slack_estimate_pu


def test_stratified_grid_fills_every_cell():
    loads = np.array([0.9, 1.0, 1.1])
    availability = np.array([0.8, 0.9, 1.0])
    rng = np.random.default_rng(3)
    points = stratified_grid(loads, availability, 5, rng)
    assert len(points) == 4 * 5
    for load, reactive, avail in points:
        assert 0.9 <= load <= 1.1
        assert 0.8 <= avail <= 1.0
        assert 0.90 <= reactive <= 1.10


def test_limit_audit_reports_reactive_binding():
    network = load_network()
    audit = audit_limits(network, solve_power_flow(network))
    assert audit.slack_within_limits
    assert isinstance(audit.reactive_binding, tuple)


def test_loss_allowance_is_declared():
    assert 0.0 < LOSS_ALLOWANCE < 0.1
