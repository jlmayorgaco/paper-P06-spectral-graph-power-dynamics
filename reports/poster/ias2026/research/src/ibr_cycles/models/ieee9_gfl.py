"""Modified WSCC/IEEE 9-bus case with one grid-following converter.

Topology and branch data are the standard nine-bus case (MATPOWER ``case9``
convention, loads at buses 5, 7 and 9). Bus 1 is a stiff reference, bus 2 keeps
a classical synchronous machine, and the machine at bus 3 is replaced by a
grid-following converter with an explicit SRF PLL, power measurement filters,
outer P and Q loops, inner current loops and an L filter.

    12 differential states + 16 algebraic voltage variables = 28 variables

Every controller equation is written out rather than lumped, because the whole
mechanism under study lives in the hidden controller states.

.. warning::

   The eigenvalue envelope quoted in the original specification could not be
   certified: the four reference scripts it cites are not in this repository and
   the textual description fixes neither the dq sign convention nor the machine
   and power-flow data to the six decimals quoted. This module is an independent
   implementation from the stated equations. Its numbers are reproducible from
   this source, and they are NOT a reproduction of that envelope.
   The one convention the specification does pin down is the PLL parameterization:
   action C maps ``(wn, zeta) = (12, 0.15)`` to ``kp = 3.6``, ``ki = 144``, which
   forces ``ki = wn^2`` and ``kp = 2 zeta wn`` with the loop written directly in
   rad/s on ``v_q``. That is the convention used here.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import root

OMEGA_B = 2.0 * np.pi * 60.0
N_BUS = 9
SLACK_BUS = 1
SG_BUS = 2
GFL_BUS = 3

#: (from, to, r, x, total line charging b) on a 100 MVA base.
BRANCHES: tuple[tuple[int, int, float, float, float], ...] = (
    (1, 4, 0.0000, 0.0576, 0.000),
    (4, 5, 0.0170, 0.0920, 0.158),
    (5, 6, 0.0390, 0.1700, 0.358),
    (3, 6, 0.0000, 0.0586, 0.000),
    (6, 7, 0.0119, 0.1008, 0.209),
    (7, 8, 0.0085, 0.0720, 0.149),
    (8, 2, 0.0000, 0.0625, 0.000),
    (8, 9, 0.0320, 0.1610, 0.306),
    (9, 4, 0.0100, 0.0850, 0.176),
)

#: Constant PQ loads, per unit on a 100 MVA base.
LOADS: dict[int, complex] = {5: 0.90 + 0.30j, 7: 1.00 + 0.35j, 9: 1.25 + 0.50j}

#: Power-flow scheduling: slack voltage, and (P, V) for each PV bus.
SLACK_VOLTAGE = 1.04
PV_SCHEDULE: dict[int, tuple[float, float]] = {2: (1.63, 1.025), 3: (0.85, 1.025)}

STATE_LABELS: tuple[str, ...] = (
    "delta_sg",
    "dw_sg",
    "theta_pll",
    "x_pll",
    "p_filt",
    "q_filt",
    "x_p",
    "x_q",
    "i_d",
    "i_q",
    "x_id",
    "x_iq",
)

#: Hidden coordinates for the self-energy: everything inside the converter
#: controller. The retained coordinates are the machine and the converter
#: currents, i.e. what the network actually sees.
RETAINED_STATES: tuple[str, ...] = ("delta_sg", "dw_sg", "i_d", "i_q")

ALGEBRAIC_BUSES: tuple[int, ...] = (2, 3, 4, 5, 6, 7, 8, 9)


@dataclass(frozen=True)
class Ieee9Parameters:
    """Machine, converter and network parameters, all per unit on 100 MVA."""

    # synchronous machine at bus 2, classical model
    h_sg: float = 6.4
    d_sg: float = 2.0
    xdp_sg: float = 0.1198
    # PLL, written directly in rad/s on v_q: ki = wn^2, kp = 2 zeta wn
    kp_pll: float = 53.0
    ki_pll: float = 1400.0
    # power measurement filter
    tau_p: float = 0.03
    # outer loops
    kp_p: float = 0.20
    ki_p: float = 8.0
    kp_q: float = 0.20
    ki_q: float = 8.0
    # inner current loop
    kp_i: float = 0.25
    ki_i: float = 6.0
    # converter filter
    xf: float = 0.15
    rf: float = 0.01
    # references, filled from the power flow
    p_ref: float = 0.85
    q_ref: float = 0.0
    pm_sg: float = 1.63
    e_sg: float = 1.0

    @staticmethod
    def pll_from_design(wn: float, zeta: float) -> dict[str, float]:
        return {"ki_pll": wn * wn, "kp_pll": 2.0 * zeta * wn}


#: Stress actions. Each maps parameter name -> absolute replacement value.
ACTIONS: dict[str, dict[str, float]] = {
    "A": {"d_sg": 0.0},
    "B": {"kp_p": 2.0},
    "C": Ieee9Parameters.pll_from_design(12.0, 0.15),
}


def apply_actions(base: Ieee9Parameters, members: tuple[str, ...]) -> Ieee9Parameters:
    """Actions replace parameter values; overlapping targets are rejected."""

    updates: dict[str, float] = {}
    for name in members:
        if name not in ACTIONS:
            raise KeyError(f"unknown action {name!r}")
        for key, value in ACTIONS[name].items():
            if key in updates and updates[key] != value:
                raise ValueError(f"actions disagree on parameter {key!r}")
            updates[key] = value
    return replace(base, **updates)


def build_ybus() -> NDArray[np.complex128]:
    """Bus admittance matrix including line charging."""

    y = np.zeros((N_BUS, N_BUS), dtype=np.complex128)
    for src, dst, r, x, b in BRANCHES:
        i, j = src - 1, dst - 1
        series = 1.0 / complex(r, x)
        y[i, i] += series + 0.5j * b
        y[j, j] += series + 0.5j * b
        y[i, j] -= series
        y[j, i] -= series
    return y


@dataclass(frozen=True)
class PowerFlow:
    """A converged AC operating point with its mismatch evidence."""

    voltages: NDArray[np.complex128]
    max_mismatch: float
    converged: bool

    @property
    def min_voltage(self) -> float:
        return float(np.abs(self.voltages).min())

    @property
    def max_voltage(self) -> float:
        return float(np.abs(self.voltages).max())

    def injection(self, bus: int, ybus: NDArray[np.complex128]) -> complex:
        current = ybus[bus - 1, :] @ self.voltages
        return complex(self.voltages[bus - 1] * np.conj(current))


def solve_power_flow(
    ybus: NDArray[np.complex128] | None = None, *, tol: float = 1e-12
) -> PowerFlow:
    """Newton AC power flow: bus 1 slack, buses 2 and 3 PV, the rest PQ."""

    y = build_ybus() if ybus is None else ybus
    pq_buses = tuple(b for b in range(2, N_BUS + 1) if b not in PV_SCHEDULE)

    def unpack(w: NDArray[np.float64]) -> NDArray[np.complex128]:
        angles = np.zeros(N_BUS)
        magnitudes = np.full(N_BUS, 1.0)
        magnitudes[0] = SLACK_VOLTAGE
        for bus, (_, v) in PV_SCHEDULE.items():
            magnitudes[bus - 1] = v
        angles[1:] = w[: N_BUS - 1]
        for offset, bus in enumerate(pq_buses):
            magnitudes[bus - 1] = w[N_BUS - 1 + offset]
        return magnitudes * np.exp(1j * angles)

    def mismatch(w: NDArray[np.float64]) -> NDArray[np.float64]:
        v = unpack(w)
        s = v * np.conj(y @ v)
        scheduled = np.zeros(N_BUS, dtype=np.complex128)
        for bus, (p, _) in PV_SCHEDULE.items():
            scheduled[bus - 1] += p
        for bus, load in LOADS.items():
            scheduled[bus - 1] -= load
        residual = list(np.real(s - scheduled)[1:])
        residual += [float(np.imag(s - scheduled)[bus - 1]) for bus in pq_buses]
        return np.array(residual)

    guess = np.concatenate([np.zeros(N_BUS - 1), np.ones(len(pq_buses))])
    solution = root(mismatch, guess, method="hybr", tol=tol)
    voltages = unpack(solution.x)
    residual = float(np.abs(mismatch(solution.x)).max())
    return PowerFlow(
        voltages=voltages,
        max_mismatch=residual,
        converged=bool(solution.success and residual < 1e-8),
    )


@dataclass
class Ieee9GflModel:
    """The semi-explicit DAE of the modified nine-bus case."""

    parameters: Ieee9Parameters
    ybus: NDArray[np.complex128] = field(default_factory=build_ybus)
    n_x: int = 12
    n_z: int = 16

    @property
    def labels(self) -> tuple[str, ...]:
        return STATE_LABELS

    def _voltages(self, z: NDArray[np.float64]) -> NDArray[np.complex128]:
        v = np.zeros(N_BUS, dtype=np.complex128)
        v[SLACK_BUS - 1] = complex(SLACK_VOLTAGE, 0.0)
        for index, bus in enumerate(ALGEBRAIC_BUSES):
            v[bus - 1] = complex(z[2 * index], z[2 * index + 1])
        return v

    def _dq(self, v: complex, theta: float) -> tuple[float, float]:
        rotated = v * np.exp(-1j * theta)
        return float(rotated.real), float(rotated.imag)

    def _sg_current(self, x: NDArray[np.float64], v2: complex) -> complex:
        p = self.parameters
        internal = p.e_sg * np.exp(1j * x[0])
        return (internal - v2) / (1j * p.xdp_sg)

    def f(
        self, x: NDArray[np.float64], z: NDArray[np.float64], theta: dict[str, float]
    ) -> NDArray[np.float64]:
        p = replace(self.parameters, **theta) if theta else self.parameters
        v = self._voltages(z)
        v2, v3 = v[SG_BUS - 1], v[GFL_BUS - 1]

        internal = p.e_sg * np.exp(1j * x[0])
        i_sg = (internal - v2) / (1j * p.xdp_sg)
        pe = float(np.real(internal * np.conj(i_sg)))

        theta_pll = float(x[2])
        v_d, v_q = self._dq(v3, theta_pll)
        i_d, i_q = float(x[8]), float(x[9])

        power = v_d * i_d + v_q * i_q
        reactive = v_q * i_d - v_d * i_q

        id_ref = p.kp_p * (p.p_ref - x[4]) + x[6]
        # With the d axis aligned to the terminal voltage, injected reactive
        # power is Q = v_q i_d - v_d i_q, so dQ/di_q = -v_d < 0. The q-axis
        # loop therefore carries an explicit inversion; without it the PI acts
        # as positive feedback and the base case shows a spurious real pole
        # near +7.6 rad/s that has nothing to do with any grid interaction.
        iq_ref = -(p.kp_q * (p.q_ref - x[5]) + x[7])
        e_d = v_d + p.kp_i * (id_ref - i_d) + x[10] - p.xf * i_q
        e_q = v_q + p.kp_i * (iq_ref - i_q) + x[11] + p.xf * i_d

        return np.array(
            [
                OMEGA_B * x[1],
                (p.pm_sg - pe - p.d_sg * x[1]) / (2.0 * p.h_sg),
                p.kp_pll * v_q + x[3],
                p.ki_pll * v_q,
                (power - x[4]) / p.tau_p,
                (reactive - x[5]) / p.tau_p,
                p.ki_p * (p.p_ref - x[4]),
                p.ki_q * (p.q_ref - x[5]),
                (OMEGA_B / p.xf) * (e_d - v_d - p.rf * i_d + p.xf * i_q),
                (OMEGA_B / p.xf) * (e_q - v_q - p.rf * i_q - p.xf * i_d),
                p.ki_i * (id_ref - i_d),
                p.ki_i * (iq_ref - i_q),
            ]
        )

    def g(
        self, x: NDArray[np.float64], z: NDArray[np.float64], theta: dict[str, float]
    ) -> NDArray[np.float64]:
        v = self._voltages(z)
        current = self.ybus @ v
        injection = np.zeros(N_BUS, dtype=np.complex128)
        injection[SG_BUS - 1] = self._sg_current(x, v[SG_BUS - 1])
        injection[GFL_BUS - 1] = complex(x[8], x[9]) * np.exp(1j * float(x[2]))
        for bus, load in LOADS.items():
            injection[bus - 1] -= np.conj(load) / np.conj(v[bus - 1])
        residual = current - injection
        out = np.zeros(self.n_z)
        for index, bus in enumerate(ALGEBRAIC_BUSES):
            out[2 * index] = float(residual[bus - 1].real)
            out[2 * index + 1] = float(residual[bus - 1].imag)
        return out


def initialize(
    power_flow: PowerFlow, base: Ieee9Parameters
) -> tuple[Ieee9Parameters, NDArray[np.float64], NDArray[np.float64]]:
    """Derive machine and converter setpoints from the AC operating point."""

    ybus = build_ybus()
    v = power_flow.voltages
    v2, v3 = v[SG_BUS - 1], v[GFL_BUS - 1]

    s2 = power_flow.injection(SG_BUS, ybus)
    i2 = np.conj(s2 / v2)
    internal = v2 + 1j * base.xdp_sg * i2
    delta_sg = float(np.angle(internal))
    e_sg = float(np.abs(internal))
    pm = float(np.real(internal * np.conj(i2)))

    s3 = power_flow.injection(GFL_BUS, ybus)
    i3 = np.conj(s3 / v3)
    theta_pll = float(np.angle(v3))
    current_dq = i3 * np.exp(-1j * theta_pll)
    i_d, i_q = float(current_dq.real), float(current_dq.imag)
    v_d, v_q = float(np.abs(v3)), 0.0

    p_ref = v_d * i_d + v_q * i_q
    q_ref = v_q * i_d - v_d * i_q
    x_id = base.rf * i_d
    x_iq = base.rf * i_q
    x_q0 = -i_q

    parameters = replace(base, p_ref=p_ref, q_ref=q_ref, pm_sg=pm, e_sg=e_sg)
    x0 = np.array(
        [delta_sg, 0.0, theta_pll, 0.0, p_ref, q_ref, i_d, x_q0, i_d, i_q, x_id, x_iq]
    )
    z0 = np.concatenate([[v[bus - 1].real, v[bus - 1].imag] for bus in ALGEBRAIC_BUSES])
    return parameters, x0, z0
