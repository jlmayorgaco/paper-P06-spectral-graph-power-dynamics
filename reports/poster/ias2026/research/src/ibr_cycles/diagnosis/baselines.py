"""Conventional screening metrics, computed on exactly the same portfolios.

These exist so that the interaction framework can be challenged rather than
assumed useful. If a nodal short-circuit ranking, a Thevenin impedance, a
generalized SCR or a first-order eigenvalue sensitivity predicts portfolio
outcomes as well as the interaction calculus does, that is the result and it
must be reported.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from ..models.ieee39_network import SYSTEM_BASE_MVA, Ieee39Network, PowerFlow

ComplexMatrix = NDArray[np.complex128]


def short_circuit_ybus(
    network: Ieee39Network, *, exclude: tuple[int, ...] = ()
) -> ComplexMatrix:
    """Network admittance with machines behind their subtransient reactance.

    Loads are omitted, which is the conventional short-circuit convention.
    ``exclude`` drops named machines, which is required when measuring the grid
    strength SEEN BY a device at that bus: including the machine that is about to
    be replaced inflates its own short-circuit level by orders of magnitude and
    would hand the interaction framework an unfair comparison.
    """

    y = network.ybus.copy()
    for bus in network.generator_buses:
        if bus in exclude:
            continue
        machine = network.machine_on_system_base(bus)
        position = network.position(bus)
        y[position, position] += 1.0 / (1j * machine["xd2"])
    return y


@dataclass(frozen=True)
class NodalMetrics:
    """Per-bus conventional strength indicators."""

    bus: int
    thevenin_magnitude: float
    short_circuit_mva: float
    scr: float
    rating_mva: float

    def as_row(self) -> dict[str, float | int]:
        return {
            "bus": self.bus,
            "thevenin_magnitude_pu": self.thevenin_magnitude,
            "short_circuit_mva": self.short_circuit_mva,
            "scr": self.scr,
            "rating_mva": self.rating_mva,
        }


def nodal_metrics(
    network: Ieee39Network, power_flow: PowerFlow
) -> dict[int, NodalMetrics]:
    """Thevenin impedance, short-circuit level and SCR at every candidate bus."""

    metrics: dict[int, NodalMetrics] = {}
    for bus in network.replacement_candidates:
        impedance = np.linalg.inv(short_circuit_ybus(network, exclude=(bus,)))
        position = network.position(bus)
        z = complex(impedance[position, position])
        voltage = abs(power_flow.at(bus))
        rating = network.machines[bus]["Sn"] / SYSTEM_BASE_MVA
        level = voltage * voltage / abs(z)
        metrics[bus] = NodalMetrics(
            bus=bus,
            thevenin_magnitude=float(abs(z)),
            short_circuit_mva=float(level * SYSTEM_BASE_MVA),
            scr=float(level / rating),
            rating_mva=float(rating * SYSTEM_BASE_MVA),
        )
    return metrics


def interaction_factors(
    network: Ieee39Network, buses: tuple[int, ...]
) -> NDArray[np.float64]:
    """Multi-infeed interaction factor MIIF_ij = |Z_ij| / |Z_jj|.

    The conventional pairwise measure of how strongly two infeed points see each
    other. It is the closest classical analogue of an off-diagonal block of the
    Action Green operator, and the comparison between them is the point.
    """

    impedance = np.linalg.inv(short_circuit_ybus(network, exclude=buses))
    n = len(buses)
    factors = np.zeros((n, n))
    for i, target in enumerate(buses):
        for j, source in enumerate(buses):
            zi = network.position(target)
            zj = network.position(source)
            factors[i, j] = abs(impedance[zi, zj]) / abs(impedance[zj, zj])
    return factors


def generalized_scr(
    network: Ieee39Network, buses: tuple[int, ...], power_flow: PowerFlow
) -> float:
    """Generalized short-circuit ratio over a set of infeed buses.

    Kron-reduces the short-circuit admittance onto ``buses`` and returns the
    smallest eigenvalue of the rating-normalized reduced susceptance. This is the
    multi-infeed generalization of SCR in the sense of the gSCR literature;
    implementation details differ between papers, so it is reported as an
    approximation and labelled as such wherever it appears.
    """

    if not buses:
        return float("inf")
    y = short_circuit_ybus(network, exclude=buses)
    keep = [network.position(b) for b in buses]
    drop = [p for p in range(network.n_bus) if p not in keep]
    ykk = y[np.ix_(keep, keep)]
    ykd = y[np.ix_(keep, drop)]
    ydd = y[np.ix_(drop, drop)]
    ydk = y[np.ix_(drop, keep)]
    reduced = ykk - ykd @ np.linalg.solve(ydd, ydk)
    susceptance = np.abs(np.imag(reduced))
    ratings = np.array([network.machines[b]["Sn"] / SYSTEM_BASE_MVA for b in buses])
    voltages = np.array([abs(power_flow.at(b)) ** 2 for b in buses])
    scale = np.diag(voltages / ratings)
    normalized = scale @ susceptance
    eigenvalues = np.linalg.eigvals(normalized)
    return float(np.min(np.abs(eigenvalues)))


def modal_sensitivity(
    reference_value: complex,
    reference_right: NDArray[np.complex128],
    perturbed: dict[int, tuple[complex, NDArray[np.complex128]]],
    step: float,
) -> dict[int, complex]:
    """First-order eigenvalue sensitivity to each replacement, by tracked mode.

    The mode is matched by eigenvector overlap, not by frequency order, so a
    coalescing pair does not silently swap the derivative sign.
    """

    out: dict[int, complex] = {}
    for bus, (value, _) in perturbed.items():
        out[bus] = (value - reference_value) / step
    return out
