"""Dynamic devices for the IEEE-39 replacement study.

Two device families share one interface, so that a bus can carry either, or a
weighted parallel combination of both:

    derivatives(x, v_xy)      internal dynamics
    injection(x, v_xy)        current injected into the network, system pu
    initialize(v, s)          state and setpoint initialization from the AC point

Both are written per unit on their OWN rating and then scaled to the system base
by a single rating factor. That is what makes the replacement operator exactly
linear in the replacement fraction: a device at rating fraction ``w`` has system
base admittance ``w * Y_unit(s)``, so a bus carrying ``(1-rho)`` of a machine and
``rho`` of a converter has

    Y_bus(s, rho) = (1-rho) Y_SG(s) + rho Y_GFL(s)

and therefore ``Delta T(s, rho) = rho [Y_GFL(s) - Y_SG(s)]`` exactly, not as an
approximation.

MACHINE FIDELITY, declared once and not varied between experiments:

    retained    4th-order two-axis (delta, omega, Eq', Ed'), first-order AVR,
                two-state IEEEST power-input PSS (washout plus lag)
    NOT retained
                the IEEEX1 self-excited exciter, its saturation, and the TGOV1N
                governor.

The IEEEX1 exciter is dropped deliberately. Six of the ten machines in the source
case carry ``KE < 0``, and ANDES' own eigenanalysis of that case reports six
unstable real modes near ``+1.0`` together with an ill-conditioned algebraic
Jacobian (rcond about 1e-20). Reproducing that pathology would make base-case
eligibility impossible to assess, so the regulator is retained in reduced,
well-posed form with the case gain and exciter time constant. The governor is
dropped because its reheat time constant of 2.1 s lies outside the 0.2-10 Hz band
under study and its damping contribution ``Dt`` is zero in this case.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Protocol

import numpy as np
from numpy.typing import NDArray

OMEGA_B = 2.0 * np.pi * 60.0


class Device(Protocol):
    """A bus-connected dynamic device in system per unit."""

    bus: int
    weight: float

    @property
    def labels(self) -> tuple[str, ...]: ...

    @property
    def n_states(self) -> int: ...

    def derivatives(
        self, x: NDArray[np.float64], v: complex
    ) -> NDArray[np.float64]: ...

    def injection(self, x: NDArray[np.float64], v: complex) -> complex: ...


def _rotate_to_dq(v: complex, delta: float) -> tuple[float, float]:
    """Machine convention: d leads q by 90 degrees behind the rotor angle."""

    return (
        float(v.real * np.sin(delta) - v.imag * np.cos(delta)),
        float(v.real * np.cos(delta) + v.imag * np.sin(delta)),
    )


def _rotate_to_network(d: float, q: float, delta: float) -> complex:
    return complex(
        d * np.sin(delta) + q * np.cos(delta), -d * np.cos(delta) + q * np.sin(delta)
    )


@dataclass(frozen=True)
class MachineParameters:
    """Two-axis machine with a first-order AVR and a power-input PSS."""

    ra: float
    xd: float
    xq: float
    xd1: float
    xq1: float
    td10: float
    tq10: float
    m: float
    d: float
    ka: float
    ta: float
    pss_gain: float
    pss_washout: float
    pss_wash_lag: float
    pss_lag: float
    # setpoints, filled at initialization
    pm: float = 0.0
    vref: float = 1.0
    #: Service-isolation switches. ``avr_manual`` freezes the field voltage,
    #: which is manual excitation rather than a gain of zero: driving ``ka`` to
    #: zero would let the field decay instead of holding it.
    avr_manual: bool = False

    def with_services(
        self,
        *,
        inertia_scale: float = 1.0,
        pss_scale: float = 1.0,
        avr_scale: float = 1.0,
        avr_manual: bool = False,
    ) -> MachineParameters:
        """Scale individual synchronous services for attribution experiments."""

        from dataclasses import replace as _replace

        return _replace(
            self,
            m=self.m * inertia_scale,
            pss_gain=self.pss_gain * pss_scale,
            ka=self.ka * avr_scale,
            avr_manual=avr_manual,
        )


MACHINE_LABELS = ("delta", "omega", "eq1", "ed1", "efd", "pss_w", "pss_l")


@dataclass
class SynchronousMachine:
    """Rating-scaled synchronous machine.

    ``weight`` is the surviving rating fraction. Every extensive quantity
    (inertia, damping, admittance) scales with it and every intensive quantity
    (time constants, gains on per-unit signals) does not.
    """

    bus: int
    parameters: MachineParameters
    weight: float = 1.0

    @property
    def labels(self) -> tuple[str, ...]:
        return tuple(f"{name}_sg{self.bus}" for name in MACHINE_LABELS)

    @property
    def n_states(self) -> int:
        return 7

    def _currents(
        self, x: NDArray[np.float64], v: complex
    ) -> tuple[float, float, float, float]:
        p = self.parameters
        delta = float(x[0])
        vd, vq = _rotate_to_dq(v, delta)
        determinant = p.ra * p.ra + p.xd1 * p.xq1
        rd = float(x[3]) - vd
        rq = float(x[2]) - vq
        id_ = (p.ra * rd + p.xq1 * rq) / determinant
        iq = (-p.xd1 * rd + p.ra * rq) / determinant
        return id_, iq, vd, vq

    def derivatives(self, x: NDArray[np.float64], v: complex) -> NDArray[np.float64]:
        p = self.parameters
        id_, iq, vd, vq = self._currents(x, v)
        power = vd * id_ + vq * iq + p.ra * (id_ * id_ + iq * iq)
        terminal = float(abs(v))
        stabilizer = p.pss_gain * float(x[6])
        washout_out = p.pss_washout * (power - float(x[5])) / p.pss_wash_lag
        return np.array(
            [
                OMEGA_B * (float(x[1]) - 1.0),
                (p.pm - power - p.d * (float(x[1]) - 1.0)) / p.m,
                (float(x[4]) - float(x[2]) - (p.xd - p.xd1) * id_) / p.td10,
                (-float(x[3]) + (p.xq - p.xq1) * iq) / p.tq10,
                0.0
                if p.avr_manual
                else (p.ka * (p.vref + stabilizer - terminal) - float(x[4])) / p.ta,
                (power - float(x[5])) / p.pss_wash_lag,
                (washout_out - float(x[6])) / p.pss_lag,
            ]
        )

    def injection(self, x: NDArray[np.float64], v: complex) -> complex:
        id_, iq, _, _ = self._currents(x, v)
        return self.weight * _rotate_to_network(id_, iq, float(x[0]))

    def scaled(self, weight: float) -> SynchronousMachine:
        """Rating change: inertia and damping scale, per-unit reactances do not.

        The device stays on its own base, so its parameters are untouched; only
        the current it injects into the system base is scaled.
        """

        return SynchronousMachine(
            bus=self.bus, parameters=self.parameters, weight=weight
        )

    def initialize(
        self, v: complex, s: complex
    ) -> tuple[SynchronousMachine, NDArray[np.float64]]:
        """Derive rotor angle, internal emfs and setpoints from the AC point.

        ``s`` is the injection on the SYSTEM base; it is converted to the device
        base by the rating weight before the machine equations are applied.
        """

        p = self.parameters
        current = np.conj(s / self.weight / v)
        internal = v + (p.ra + 1j * p.xq) * current
        delta = float(np.angle(internal))
        vd, vq = _rotate_to_dq(v, delta)
        id_, iq = _rotate_to_dq(current, delta)
        eq1 = vq + p.ra * iq + p.xd1 * id_
        ed1 = vd + p.ra * id_ - p.xq1 * iq
        efd = eq1 + (p.xd - p.xd1) * id_
        power = vd * id_ + vq * iq + p.ra * (id_ * id_ + iq * iq)
        terminal = float(abs(v))
        vref = terminal if p.avr_manual else terminal + efd / p.ka
        tuned = replace(p, pm=power, vref=vref)
        state = np.array([delta, 1.0, eq1, ed1, efd, power, 0.0])
        return replace(self, parameters=tuned), state


@dataclass(frozen=True)
class ConverterParameters:
    """Grid-following converter, per unit on its own rating."""

    kp_pll: float = 53.0
    ki_pll: float = 1400.0
    tau_p: float = 0.03
    kp_p: float = 0.20
    ki_p: float = 8.0
    kp_q: float = 0.20
    ki_q: float = 8.0
    kp_i: float = 0.25
    ki_i: float = 6.0
    xf: float = 0.15
    rf: float = 0.01
    p_ref: float = 0.0
    q_ref: float = 0.0
    #: Outer voltage regulator. When enabled the reactive reference is produced
    #: by a PI on the terminal voltage instead of being a fixed setpoint, which
    #: is what a real plant controller does. It adds one state.
    voltage_control: bool = False
    kp_v: float = 2.0
    ki_v: float = 20.0
    v_ref: float = 1.0
    #: Reactive-policy coordinate for continuation studies (F7). ``voltage_gain``
    #: scales both regulator gains; ``voltage_leak`` bleeds the integrator back to
    #: the dispatched reactive setpoint at that rate (rad/s). With a leak, gain 0
    #: is exact fixed-Q plus one decoupled stable pole at ``-voltage_leak``; with
    #: no leak, gain 0 leaves a decoupled marginal integrator, a zero eigenvalue
    #: and a continuum of equilibria. The defaults reproduce the plain PI exactly.
    voltage_gain: float = 1.0
    voltage_leak: float = 0.0

    @staticmethod
    def pll_from_design(wn: float, zeta: float) -> dict[str, float]:
        return {"ki_pll": wn * wn, "kp_pll": 2.0 * zeta * wn}


CONVERTER_LABELS = (
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
VOLTAGE_LABEL = "x_v"


@dataclass
class GridFollowingConverter:
    """Battery-free PV grid-following converter, rating-scaled.

    Internal per-unit quantities are on the device rating; ``weight`` converts
    the injected current to the system base. The reactive loop carries the
    explicit sign inversion established in the IEEE-9 gate: with the d axis
    aligned to the terminal voltage, injected reactive power is
    ``Q = v_q i_d - v_d i_q``, so ``dQ/di_q < 0``.
    """

    bus: int
    parameters: ConverterParameters
    weight: float = 1.0

    @property
    def labels(self) -> tuple[str, ...]:
        names = CONVERTER_LABELS
        if self.parameters.voltage_control:
            names = (*names, VOLTAGE_LABEL)
        return tuple(f"{name}_gfl{self.bus}" for name in names)

    @property
    def n_states(self) -> int:
        return 11 if self.parameters.voltage_control else 10

    def derivatives(self, x: NDArray[np.float64], v: complex) -> NDArray[np.float64]:
        p = self.parameters
        theta = float(x[0])
        rotated = v * np.exp(-1j * theta)
        v_d, v_q = float(rotated.real), float(rotated.imag)
        i_d, i_q = float(x[6]), float(x[7])
        power = v_d * i_d + v_q * i_q
        reactive = v_q * i_d - v_d * i_q
        if p.voltage_control:
            error = p.v_ref - float(abs(v))
            if p.voltage_gain != 1.0:
                error = p.voltage_gain * error
            q_command = p.kp_v * error + float(x[10])
        else:
            error = 0.0
            q_command = p.q_ref
        id_ref = p.kp_p * (p.p_ref - float(x[2])) + float(x[4])
        iq_ref = -(p.kp_q * (q_command - float(x[3])) + float(x[5]))
        e_d = v_d + p.kp_i * (id_ref - i_d) + float(x[8]) - p.xf * i_q
        e_q = v_q + p.kp_i * (iq_ref - i_q) + float(x[9]) + p.xf * i_d
        derivatives = [
            p.kp_pll * v_q + float(x[1]),
            p.ki_pll * v_q,
            (power - float(x[2])) / p.tau_p,
            (reactive - float(x[3])) / p.tau_p,
            p.ki_p * (p.p_ref - float(x[2])),
            p.ki_q * (q_command - float(x[3])),
            (OMEGA_B / p.xf) * (e_d - v_d - p.rf * i_d + p.xf * i_q),
            (OMEGA_B / p.xf) * (e_q - v_q - p.rf * i_q - p.xf * i_d),
            p.ki_i * (id_ref - i_d),
            p.ki_i * (iq_ref - i_q),
        ]
        if p.voltage_control:
            integral = p.ki_v * error
            if p.voltage_leak:
                integral -= p.voltage_leak * (float(x[10]) - p.q_ref)
            derivatives.append(integral)
        return np.array(derivatives)

    def injection(self, x: NDArray[np.float64], v: complex) -> complex:
        return (
            self.weight * complex(float(x[6]), float(x[7])) * np.exp(1j * float(x[0]))
        )

    def initialize(
        self, v: complex, s: complex
    ) -> tuple[GridFollowingConverter, NDArray[np.float64]]:
        p = self.parameters
        current = np.conj(s / self.weight / v)
        theta = float(np.angle(v))
        dq = current * np.exp(-1j * theta)
        i_d, i_q = float(dq.real), float(dq.imag)
        v_d = float(abs(v))
        p_ref = v_d * i_d
        q_ref = -v_d * i_q
        tuned = replace(p, p_ref=p_ref, q_ref=q_ref, v_ref=v_d)
        values = [theta, 0.0, p_ref, q_ref, i_d, -i_q, i_d, i_q, p.rf * i_d, p.rf * i_q]
        if p.voltage_control:
            # The regulator integrator holds the steady reactive command, so the
            # voltage error is zero at the operating point.
            values.append(q_ref)
        return replace(self, parameters=tuned), np.array(values)


@dataclass
class StaticInjection:
    """Ideal constant P/Q injection with no internal dynamics.

    The negative control for controller mediation. It replaces a machine by the
    same complex power injection the machine was producing, contributing nothing
    to the state vector and no frequency-dependent port admittance beyond the
    algebraic constant-power characteristic. If a portfolio is unstable with
    grid-following converters and stable with this device at the same operating
    point, the mechanism is carried by the converter dynamics and not by the loss
    of synchronous machines alone.
    """

    bus: int
    injection_pu: complex
    weight: float = 1.0
    characteristic: str = "power"
    voltage_pu: complex = 1.0 + 0j

    @property
    def labels(self) -> tuple[str, ...]:
        return ()

    @property
    def n_states(self) -> int:
        return 0

    def derivatives(self, x: NDArray[np.float64], v: complex) -> NDArray[np.float64]:
        return np.zeros(0)

    def injection(self, x: NDArray[np.float64], v: complex) -> complex:
        """Injected current under the chosen static characteristic.

        Two bracketing references are provided. ``power`` is an ideal constant
        P/Q source and is the harshest: constant-power injection has a negative
        incremental impedance and is destabilizing on its own, so a portfolio
        that fails under it has not been shown to fail because of converter
        control. ``impedance`` is the most benign static reference. Reporting
        both keeps the control from being either a strawman or a free pass.

        A constant-current variant is deliberately absent: with the injection
        angle following the bus voltage it leaves the voltage magnitude
        undetermined and the equilibrium degenerate.
        """

        if self.characteristic == "impedance":
            admittance = np.conj(self.injection_pu) / (abs(self.voltage_pu) ** 2)
            return admittance * v
        return np.conj(self.injection_pu) / np.conj(v)

    def initialize(
        self, v: complex, s: complex
    ) -> tuple[StaticInjection, NDArray[np.float64]]:
        tuned = replace(self, injection_pu=complex(s), voltage_pu=complex(v))
        return tuned, np.zeros(0)
