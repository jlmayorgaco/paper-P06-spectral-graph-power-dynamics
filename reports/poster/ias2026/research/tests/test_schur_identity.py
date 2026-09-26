"""EXACT IDENTITY: the Schur determinant factorization of the reduced system."""

from __future__ import annotations

import numpy as np
import pytest

from ibr_cycles.models import toy5
from ibr_cycles.models.toy5_case import ACTION_NAMES, PROBE_POINTS
from ibr_cycles.reduction.residues import expand_self_energy
from ibr_cycles.reduction.schur import schur_determinant_residual
from ibr_cycles.reduction.self_energy import SelfEnergy


def test_schur_determinant_identity(toy):
    residual = schur_determinant_residual(toy.a0, toy.partition, PROBE_POINTS)
    assert residual.max_relative_error < 1e-10


def test_schur_identity_holds_for_every_portfolio(toy):
    for members in [(), ("A",), ("B", "C"), ACTION_NAMES]:
        residual = schur_determinant_residual(
            toy.state_matrix(members), toy.partition, PROBE_POINTS
        )
        assert residual.max_relative_error < 1e-10


def test_angle_level_self_energy_is_exact(toy):
    """det(sI - A) = (s + a) det(T_delta(s)) with Sigma_delta = k s/(s+a) b b^T."""

    for s in PROBE_POINTS:
        full = np.linalg.det(s * np.eye(5) - toy.a0)
        factored = (s + toy.parameters.a) * np.linalg.det(
            toy5.angle_operator(toy.parameters, s)
        )
        assert abs(full - factored) / max(abs(full), 1e-300) < 1e-12


def test_self_energy_pole_expansion(toy):
    sigma = SelfEnergy.from_matrix(toy.a0, toy.partition)
    expansion = expand_self_energy(sigma)
    assert expansion.diagonalizable
    assert expansion.max_reconstruction_error < 1e-12
    assert len(expansion.poles) == 1
    assert expansion.poles[0].pole == pytest.approx(-toy.parameters.a)
