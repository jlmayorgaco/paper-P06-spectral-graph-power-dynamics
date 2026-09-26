"""BC02 adapter contracts on the Kundur model (small, fast) and IEEE-39 identity."""

from __future__ import annotations

import numpy as np
import pytest

from ibr_cycles.certification.adapter import (
    NOT_APPLICABLE,
    UNKNOWN,
    Ieee39Adapter,
    KundurAdapter,
)
from ibr_cycles.certification.binary import all_vertices, reduced, stacked
from ibr_cycles.models.ieee39_devices import OMEGA_B

POLICY = {"g": 0.08, "k": 1.25, "t": 1.0}


@pytest.fixture(scope="module")
def kundur():
    return KundurAdapter()


def test_explicit_statuses(kundur):
    assert kundur.equilibrate((2,), POLICY, uncertainty={"load": 1.1}) == NOT_APPLICABLE
    assert kundur.operator_error_bound(cell=None) == UNKNOWN


def test_symmetry_and_jordan_partner_from_the_equations(kundur):
    case = kundur.equilibrate((2, 3), POLICY)
    lin = kundur.linearize_full(case)
    r_x, r_z, partner = kundur.symmetry_generators(case)
    scale = max(np.abs(lin["fx"]).max(), np.abs(lin["gz"]).max())
    assert np.abs(lin["fx"] @ r_x + lin["fz"] @ r_z).max() < 1e-8 * scale
    assert np.abs(lin["gx"] @ r_x + lin["gz"] @ r_z).max() < 1e-8 * scale
    a = case.system.A
    assert partner.reason == "DERIVED"
    assert np.abs(a @ partner.w - OMEGA_B * r_x).max() < 1e-8 * np.abs(a).max()


def test_port_determinant_identity(kundur):
    case = kundur.equilibrate((2,), POLICY)
    t_of, h_of = kundur.port_model(case)
    lin = kundur.linearize_full(case)
    n = lin["fx"].shape[0]
    for s in (0.3 + 2j, 1.0 + 0.1j):
        p = np.block([[s * np.eye(n) - lin["fx"], -lin["fz"]], [lin["gx"], lin["gz"]]])
        assert np.isclose(np.linalg.det(p), h_of(s) * np.linalg.det(t_of(s)), rtol=1e-8)


def test_common_realization_is_affine_and_vertex_equivalent(kundur):
    cr = kundur.affine_binary_model(POLICY)
    jac = {}
    for delta, s in all_vertices(cr.candidates):
        j, residual, _ = cr.jacobian(delta, "C2")
        assert residual < 1e-9  # one equilibrium for every vertex
        jac[s] = (stacked(j), reduced(j))
    j0 = jac[()][0]
    for s, (jm, _) in jac.items():
        pred = j0 + sum(jac[(b,)][0] - j0 for b in s)
        assert np.linalg.norm(jm - pred) <= 1e-12 * np.linalg.norm(jm)
    # vertex equivalence: filler adds -1's, the rest is the real portfolio
    real = np.sort_complex(
        np.linalg.eigvals(kundur.equilibrate((2, 4), POLICY).system.A)
    )
    common = np.linalg.eigvals(jac[(2, 4)][1])
    extra = common[np.isclose(common, -1.0, atol=1e-7)]
    rest = np.sort_complex(common[~np.isclose(common, -1.0, atol=1e-7)])
    assert rest.size + extra.size == common.size
    assert (
        np.max(np.abs(np.sort_complex(real[~np.isclose(real, -1.0, atol=1e-7)]) - rest))
        < 1e-4
    )


def test_ieee39_adapter_matches_direct_path():
    from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case

    adapter = Ieee39Adapter()
    case = adapter.equilibrate((30,), {"g": 0.0, "k": 1.0, "t": 1.0})
    direct = solve_case(ReplacementPlan.of({30: 1.0}), **adapter._kwargs({"g": 0.0}))
    assert np.allclose(case.system.A, direct.system.A)
