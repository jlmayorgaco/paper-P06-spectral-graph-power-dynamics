"""Gate 3: the documented 68-bus model reproduces the source report.

Checked against Singh & Pal (2013), Table 1 (power flow) and Table 4 (modes with
PSS on G1-G12). These tests fail on a transcription error in
configs/ieee68/ieee68_network.json or an error in the documented equations.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case
from ibr_cycles.models.ieee39_devices import ConverterParameters
from ibr_cycles.models.ieee39_network import load_network, solve_power_flow

NETWORK = (
    Path(__file__).resolve().parents[1] / "configs" / "ieee68" / "ieee68_network.json"
)
TABLE1 = {
    1: (1.045, -8.9563),
    13: (1.011, -28.6539),
    17: (0.9499, -36.0269),
    39: (0.9915, -39.2902),
    41: (0.9996, 9.4272),
    64: (0.8367, -8.377),
}
INTERAREA = [(33.537, 0.314), (3.621, 0.52), (9.625, 0.591), (3.381, 0.779)]


@pytest.fixture(scope="module")
def base():
    return solve_case(ReplacementPlan.of({}), network=load_network(NETWORK))


def test_power_flow_reproduces_table1():
    net = load_network(NETWORK)
    pf = solve_power_flow(net)
    assert pf.converged
    for bus, (v, angle) in TABLE1.items():
        assert abs(abs(pf.at(bus)) - v) < 1e-3
        assert abs(math.degrees(np.angle(pf.at(bus))) - angle) < 0.05


def test_every_device_is_at_equilibrium(base):
    f = base.dae.f(base.equilibrium.x, base.equilibrium.z, {})
    g = base.dae.g(base.equilibrium.x, base.equilibrium.z, {})
    assert np.abs(f).max() < 1e-8
    assert np.abs(g).max() < 1e-8


def test_interarea_modes_reproduce_table4(base):
    values = np.linalg.eigvals(base.system.A)
    osc = values[(np.abs(values) > 1e-3) & (values.imag > 0)]
    freq = osc.imag / (2 * math.pi)
    zeta = -osc.real / np.abs(osc)
    for z, f in INTERAREA:
        i = int(np.argmin(np.hypot((freq - f) / 0.05, (zeta - z / 100) / 0.01)))
        assert abs(freq[i] - f) <= 0.005
        assert abs(100 * zeta[i] - z) <= 0.05


def test_base_is_stable_with_two_reference_zeros(base):
    values = np.linalg.eigvals(base.system.A)
    assert np.count_nonzero(np.abs(values) <= 1e-3) == 2
    assert values[np.abs(values) > 1e-3].real.max() < 0


def test_equilibrium_inside_documented_limits(base):
    margins = [s.device.limit_margin for s in base.dae.slots]
    assert min(margins) > 0


def test_policy_and_gain_enter_linearly():
    """The exact-assembly premise of G3: A(g, k) is affine in g and in k."""

    net = load_network(NETWORK)

    def a(g, k):
        return solve_case(
            ReplacementPlan.of({9: 1.0, 6: 1.0}),
            network=net,
            converter=ConverterParameters(
                voltage_control=True, voltage_gain=g, voltage_leak=0.05
            ),
            machine_scaling={"ka": k},
        ).system.A

    a00, a10, a02 = a(0.0, 1.0), a(1.0, 1.0), a(0.0, 2.0)
    direct = a(0.37, 1.6)
    assembled = a00 + 0.37 * (a10 - a00) + 0.6 * (a02 - a00)
    assert np.abs(direct - assembled).max() <= 1e-9 * np.abs(direct).max()


def test_gain_scale_keeps_the_operating_point():
    net = load_network(NETWORK)
    one = solve_case(ReplacementPlan.of({}), network=net, machine_scaling={"ka": 1.0})
    two = solve_case(ReplacementPlan.of({}), network=net, machine_scaling={"ka": 2.0})
    assert np.abs(one.equilibrium.z - two.equilibrium.z).max() < 1e-10
