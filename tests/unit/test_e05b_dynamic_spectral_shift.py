from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "experiments" / "tx3" / "E05B_dynamic_spectral_shift" / "dynamic_spectral_shift.py"
SPEC = importlib.util.spec_from_file_location("dynamic_spectral_shift", MODULE)
assert SPEC is not None and SPEC.loader is not None
SHIFT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SHIFT)


def test_connected_resolvent_and_contour_moments_match_scalar_residues() -> None:
    spectra = {
        (): np.asarray([-1.0 + 2.0j]),
        (1,): np.asarray([-1.1 + 2.1j]),
        (2,): np.asarray([-0.9 + 1.8j]),
        (1, 2): np.asarray([-0.8 + 2.4j]),
    }
    exact = SHIFT.exact_contour_moments(spectra, (1, 2), -1.0 + 2.0j, 1.0)
    expected = spectra[(1, 2)][0] - spectra[(1,)][0] - spectra[(2,)][0] + spectra[()][0]
    assert exact["single_pole_each_vertex"]
    assert exact["mu0"] == 0
    assert abs(exact["mu1"] - expected) < 1e-14
    mu0, mu1 = SHIFT.numerical_contour_moments(spectra, (1, 2), -1.0 + 2.0j, 1.0, 256)
    assert abs(mu0) < 1e-12
    assert abs(mu1 - expected) < 1e-12


def test_pole_logdet_is_characteristic_determinant_logabs() -> None:
    poles = np.asarray([-1.0 + 2.0j, -1.0 - 2.0j, -3.0 + 0.0j])
    frequencies = np.asarray([0.1, 1.0, 10.0])
    actual = SHIFT.pole_logdet(poles, frequencies)
    expected = []
    for frequency in frequencies:
        matrix = 2j * np.pi * frequency * np.eye(3) - np.diag(poles)
        expected.append(np.linalg.slogdet(matrix)[1])
    assert np.allclose(actual, expected, atol=1e-13, rtol=0.0)


def test_vertex_and_coalition_enumeration_is_complete() -> None:
    assert len(SHIFT.all_vertices()) == 93
    assert len(SHIFT.all_coalitions()) == 84
    assert sum(len(SHIFT.mobius_vertices(value)) for value in SHIFT.all_coalitions()) == 560
