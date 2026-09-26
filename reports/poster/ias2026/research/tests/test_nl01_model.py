"""NL01: the L1 matrix model reproduces L0, and the model respects its physics.

Every test here is on the frozen L0 code or on synthetic branches; none adjusts
a benchmark.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case
from ibr_cycles.models.ieee39_devices import ConverterParameters
from ibr_cycles.models.ieee39_network import load_network, solve_power_flow
from ibr_cycles.nonlinear.model import PhasorModel
from ibr_cycles.nonlinear.network import (
    interleave_permutation,
    load_branch_network,
    realify,
)

RESEARCH = Path(__file__).resolve().parents[1]
NETS = {
    "ieee39": RESEARCH / "configs" / "ias2026" / "ieee39_network.json",
    "kundur": RESEARCH / "configs" / "kundur" / "kundur_network.json",
    "ieee68": RESEARCH / "configs" / "ieee68" / "ieee68_network.json",
}
CONV = ConverterParameters(voltage_control=True, voltage_gain=0.05, voltage_leak=0.05)


@pytest.mark.parametrize("name", list(NETS))
def test_l1_branch_ybus_equals_l0(name):
    l1 = load_branch_network(NETS[name])
    l0 = load_network(NETS[name])
    assert np.abs(l1.ybus - l0.ybus).max() < 1e-12 * np.abs(l0.ybus).max()


@pytest.mark.parametrize("name", list(NETS))
def test_branch_losses_with_taps(name):
    l1 = load_branch_network(NETS[name])
    pf = solve_power_flow(load_network(NETS[name]))
    flows = l1.branch_flows(pf.voltages)
    assert np.abs(flows["loss"] - flows["loss_series"]).max() < 1e-10
    # KCL: sum of branch injections at each bus equals Ybus V minus shunts
    injected = l1.cf.T @ flows["I_f"] + l1.ct.T @ flows["I_t"]
    assert (
        np.abs(injected + l1.y_shunt * pf.voltages - l1.ybus @ pf.voltages).max()
        < 1e-10
    )


def test_synthetic_branches_with_phase_shifters():
    rng = np.random.default_rng(3)
    for _ in range(40):
        r, x = rng.uniform(0.001, 0.05), rng.uniform(0.01, 0.3)
        b = rng.uniform(0, 0.5)
        t = rng.uniform(0.9, 1.1) * np.exp(1j * rng.uniform(-0.3, 0.3))
        y = 1 / complex(r, x)
        yff, yft, ytf, ytt = (
            (y + 0.5j * b) / abs(t) ** 2,
            -y / np.conj(t),
            -y / t,
            y + 0.5j * b,
        )
        vf, vt = (
            rng.uniform(0.9, 1.1) * np.exp(1j * rng.uniform(-1, 1)),
            rng.uniform(0.9, 1.1) * np.exp(1j * rng.uniform(-1, 1)),
        )
        i_f, i_t = yff * vf + yft * vt, ytf * vf + ytt * vt
        loss = (vf * np.conj(i_f) + vt * np.conj(i_t)).real
        assert abs(loss - r * abs(y * (vf / t - vt)) ** 2) < 1e-13
        if abs(np.angle(t)) > 1e-3:
            assert abs(yft - ytf) > 1e-6  # no symmetry to impose


def test_realification_orders_are_one_permutation():
    rng = np.random.default_rng(4)
    m = rng.standard_normal((5, 5)) + 1j * rng.standard_normal((5, 5))
    v = rng.standard_normal(5) + 1j * rng.standard_normal(5)
    stacked = np.concatenate([v.real, v.imag])
    p = interleave_permutation(5)
    out_s = realify(m, "stacked") @ stacked
    out_i = realify(m, "interleaved") @ stacked[p]
    assert np.allclose(out_s[p], out_i)
    assert np.allclose(out_s[:5] + 1j * out_s[5:], m @ v)


@pytest.fixture(scope="module")
def flagship():
    return solve_case(
        ReplacementPlan.of({b: 1.0 for b in (30, 33, 35, 37)}), converter=CONV
    )


def test_equilibrium_and_device_power(flagship):
    model = PhasorModel(flagship)
    x, z = model.initialize_from_pf()
    assert np.abs(model.residual_f(x, z).values).max() < 1e-9
    assert np.abs(model.residual_g(x, z).values).max() < 1e-9
    out = model.outputs(x, z)
    names = dict(zip(out.names, out.values, strict=True))
    pf = flagship.dae.power_flow
    net = flagship.dae.network
    for bus in (30, 33):
        gen = pf.injection(bus, net.ybus) + net.loads.get(bus, 0j)
        assert abs(names[f"P_gfl{bus}"] - gen.real) < 1e-9
        assert abs(names[f"Q_gfl{bus}"] - gen.imag) < 1e-9


def test_l1_kcl_equals_l0_residual(flagship):
    """g of L0 equals the L1 matrix KCL: Ybus V - Cdev I_dev + I_load."""

    dae = flagship.dae
    l1 = load_branch_network(NETS["ieee39"])
    rng = np.random.default_rng(5)
    x = flagship.equilibrium.x + 1e-3 * rng.standard_normal(dae.n_x)
    z = flagship.equilibrium.z + 1e-3 * rng.standard_normal(dae.n_z)
    v = dae.voltages(z)
    inj = np.zeros(dae.network.n_bus, complex)
    for slot in dae.slots:
        pos = dae.network.position(slot.bus)
        inj[pos] += slot.device.injection(x[slot.start : slot.stop], complex(v[pos]))
    for bus, load in dae.network.loads.items():
        pos = dae.network.position(bus)
        inj[pos] -= np.conj(load) / np.conj(v[pos])
    res = l1.ybus @ v - inj
    g0 = dae.g(x, z, {})
    assert np.abs(g0[0::2] + 1j * g0[1::2] - res).max() < 1e-12


def _rest_spectrum(case):
    from ibr_cycles.certification.symmetry import (
        deflate_eigenvector,
        frequency_partner,
        quotient,
        rotation_generator,
    )

    r_x, _ = rotation_generator(case.dae, case.equilibrium.z)
    q = quotient(case.system.A, r_x)
    v = q.z.T @ frequency_partner(case.dae).w
    return np.linalg.eigvals(deflate_eigenvector(q.a_q, v / np.linalg.norm(v)).a_q)


def test_bus_reordering_does_not_change_the_spectrum(tmp_path):
    payload = json.loads(NETS["kundur"].read_text(encoding="utf-8"))
    reordered = dict(payload)
    reordered["buses"] = list(reversed(payload["buses"]))
    path = tmp_path / "kundur_reversed.json"
    path.write_text(json.dumps(reordered), encoding="utf-8")
    kw = {"converter": CONV, "machine_scaling": {"ka": 1.25, "ta": 1.0}}
    a = solve_case(
        ReplacementPlan.of({2: 1.0}), network=load_network(NETS["kundur"]), **kw
    )
    b = solve_case(ReplacementPlan.of({2: 1.0}), network=load_network(path), **kw)
    # compare the physical spectrum: the rotation and its Jordan partner (an
    # exact double zero split by finite-difference noise) are removed exactly
    ea = np.sort_complex(_rest_spectrum(a))
    eb = np.sort_complex(_rest_spectrum(b))
    assert np.max(np.abs(ea - eb)) < 1e-6


def test_input_perturbation_changes_residuals_only_where_declared(flagship):
    model = PhasorModel(flagship)
    x, z = model.initialize_from_pf()
    u = np.zeros(model.n_u)
    k = model.input_names.index("dP_load_16")
    u[k] = 0.01
    dg = model.residual_g(x, z, u).values - model.residual_g(x, z).values
    pos = flagship.dae.network.position(16)
    nonzero = np.flatnonzero(np.abs(dg) > 1e-12)
    assert set(nonzero) <= {2 * pos, 2 * pos + 1}
    assert (
        np.abs(model.residual_f(x, z, u).values - model.residual_f(x, z).values).max()
        == 0
    )
