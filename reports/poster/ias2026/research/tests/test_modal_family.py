"""Subspace metrics and the family envelope, on cases with known answers."""

from __future__ import annotations

import numpy as np
import pytest

from ibr_cycles.dynamics.modal_family import (
    ModalFamily,
    descendant_family,
    principal_cosines,
    subspace_similarity,
)
from ibr_cycles.dynamics.modes import eigen_analysis

BAND = (0.3, 1.5)


def test_identical_spans_have_similarity_one():
    rng = np.random.default_rng(1)
    a = rng.normal(size=(9, 3)) + 1j * rng.normal(size=(9, 3))
    assert subspace_similarity(a, a) == pytest.approx(1.0, abs=1e-12)
    assert subspace_similarity(a, a @ rng.normal(size=(3, 3))) == pytest.approx(
        1.0, abs=1e-10
    )


def test_orthogonal_spans_have_similarity_zero():
    basis = np.eye(6, dtype=complex)
    assert subspace_similarity(basis[:, :2], basis[:, 2:4]) == pytest.approx(0.0, abs=1e-12)


def test_similarity_survives_a_swap_of_the_members():
    """The point of the whole module: order and mixing must not matter."""

    rng = np.random.default_rng(2)
    a = rng.normal(size=(8, 2)) + 1j * rng.normal(size=(8, 2))
    swapped = a[:, ::-1]
    mixed = a @ np.array([[0.6, 0.8], [-0.8, 0.6]], dtype=complex)
    assert subspace_similarity(a, swapped) == pytest.approx(1.0, abs=1e-10)
    assert subspace_similarity(a, mixed) == pytest.approx(1.0, abs=1e-10)


def test_a_subspace_inside_a_larger_one_scores_one():
    basis = np.eye(5, dtype=complex)
    assert subspace_similarity(basis[:, :1], basis[:, :3]) == pytest.approx(1.0, abs=1e-12)


def test_principal_cosines_match_a_known_angle():
    angle = 0.4
    a = np.array([[1.0], [0.0]], dtype=complex)
    b = np.array([[np.cos(angle)], [np.sin(angle)]], dtype=complex)
    assert principal_cosines(a, b)[0] == pytest.approx(np.cos(angle), abs=1e-12)


def _rotation(frequency_hz, real):
    omega = 2.0 * np.pi * frequency_hz
    return np.array([[real, omega], [-omega, real]])


def test_family_collects_both_descendants_and_takes_the_envelope():
    """Two oscillators in the band that share a shape, plus one that does not."""

    block = np.zeros((6, 6))
    block[0:2, 0:2] = _rotation(0.57, +0.05)
    block[2:4, 2:4] = _rotation(0.66, -0.60)
    block[4:6, 4:6] = _rotation(1.20, -0.10)
    mix = np.eye(6)
    mix[0, 2] = mix[2, 0] = 0.9  # make the first two share a direction
    matrix = mix @ block @ np.linalg.inv(mix)

    spectrum = eigen_analysis(matrix)
    reference = np.array([1.0, 0.0, 1.0, 0.0, 0.0, 0.0], dtype=complex)
    reference /= np.linalg.norm(reference)
    family = descendant_family(
        reference, list(range(6)), spectrum, band_hz=BAND, threshold=0.30
    )

    assert family is not None
    assert family.size == 2
    assert family.alpha == pytest.approx(0.05, abs=1e-9)
    assert family.frequency_worst_hz == pytest.approx(0.57, abs=1e-6)
    assert all(f < 0.7 for f in family.frequencies_hz)  # 1.20 Hz is excluded


def test_family_reduces_to_the_single_branch_when_there_is_only_one():
    block = np.zeros((4, 4))
    block[0:2, 0:2] = _rotation(0.60, -0.20)
    block[2:4, 2:4] = _rotation(0.90, -0.30)
    spectrum = eigen_analysis(block)
    # the exact shape of the 0.60 Hz mode: [[a, w], [-w, a]] has eigenvector [1, i]
    reference = np.array([1.0, 1.0j, 0.0, 0.0], dtype=complex) / np.sqrt(2.0)
    family = descendant_family(
        reference, list(range(4)), spectrum, band_hz=BAND, threshold=0.80
    )
    assert family is not None and family.size == 1
    assert family.alpha == pytest.approx(-0.20, abs=1e-9)


def test_empty_family_returns_none_rather_than_a_fallback():
    block = np.zeros((2, 2))
    block[0:2, 0:2] = _rotation(0.60, -0.20)
    spectrum = eigen_analysis(block)
    orthogonal = np.array([0.0, 0.0], dtype=complex)
    assert (
        descendant_family(
            orthogonal, [0, 1], spectrum, band_hz=BAND, threshold=0.80
        )
        is None
    )


def test_envelope_is_never_below_any_member():
    block = np.zeros((4, 4))
    block[0:2, 0:2] = _rotation(0.55, +0.11)
    block[2:4, 2:4] = _rotation(0.65, -0.44)
    spectrum = eigen_analysis(block)
    reference = np.array([1.0, 0.0, 1.0, 0.0], dtype=complex)
    reference /= np.linalg.norm(reference)
    family = descendant_family(
        reference, list(range(4)), spectrum, band_hz=BAND, threshold=0.10
    )
    assert isinstance(family, ModalFamily)
    assert family.alpha == max(family.reals)
    assert family.alpha >= family.best.real
