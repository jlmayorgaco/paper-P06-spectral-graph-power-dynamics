from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
from scipy.sparse import csc_matrix, diags


ROOT = Path(__file__).resolve().parents[2]
ENGINE = ROOT / "experiments" / "tx3" / "E05_discovery_holdout" / "mechanism_consequence.py"
SPEC = importlib.util.spec_from_file_location("mechanism_consequence", ENGINE)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_reduced_state_matrix_matches_schur_complement() -> None:
    fx = np.array([[-1.0, 0.3], [0.2, -2.0]])
    fy = np.array([[0.4], [-0.2]])
    gx = np.array([[0.5, 0.1]])
    gy = np.array([[-3.0]])
    jacobian = np.block([[fx, fy], [gx, gy]])
    mass = np.array([2.0, 4.0, 0.0])
    point = MODULE.OperatorPoint(
        csc_matrix(jacobian),
        diags(mass, format="csc"),
        "SUCCESS",
        {"dynamic_state_count": 2, "algebraic_variable_count": 1},
    )
    expected = (fx - fy @ np.linalg.solve(gy, gx)) / mass[:2, None]
    np.testing.assert_allclose(MODULE.reduced_state_matrix(point), expected, rtol=1e-13, atol=1e-13)


def test_biorthogonal_assignment_recovers_permuted_modes() -> None:
    reference_values = np.array([-0.1 + 1.0j, -0.2 + 2.0j])
    identity = np.eye(2, dtype=complex)
    reference = MODULE.ModeSet(reference_values, identity, identity)
    permutation = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
    target = MODULE.ModeSet(reference_values[::-1], permutation, permutation)
    mapping = MODULE.track_modes(
        reference,
        target,
        reference_band_hz=(0.1, 1.0),
        target_band_hz=(0.1, 1.0),
    )
    assert mapping[0][0] == 1
    assert mapping[1][0] == 0
    assert mapping[0][1] == 1.0
    assert mapping[1][1] == 1.0


def test_subset_enumeration_is_complete_and_ordered() -> None:
    assert MODULE.all_subsets((2, 7, 8)) == [
        (),
        (2,),
        (7,),
        (8,),
        (2, 7),
        (2, 8),
        (7, 8),
        (2, 7, 8),
    ]
