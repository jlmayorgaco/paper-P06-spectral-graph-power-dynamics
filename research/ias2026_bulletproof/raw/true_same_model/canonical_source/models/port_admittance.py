"""Port-space nonlinear eigenvalue operator for the replacement study.

The full DAE gives the ground-truth spectrum. This module gives the same system
as an operator on bus voltages only,

    T(s) = Y_net - sum_i E_i [ D_i + C_i (sI - A_i)^-1 B_i ] E_i^T

where each device contributes a local 2x2 port admittance in rectangular
coordinates. Loads contribute a constant 2x2 block. Zeros of ``det T(s)`` are the
system modes that are observable at the buses; modes confined to a device and
invisible at its terminals are absent by construction, and this module reports
which ones those are rather than hiding them.

A replacement at one bus is then a rank-at-most-2 update of ``T(s)`` at that bus,
so a four-replacement portfolio lives in an 8x8 action space regardless of how
many states the devices carry.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .ieee39_case import Ieee39Dae

ComplexMatrix = NDArray[np.complex128]
_EPS = float(np.finfo(np.float64).eps)


def _step(value: float, scale: float) -> float:
    return float(np.cbrt(_EPS) * max(abs(value), scale))


@dataclass(frozen=True)
class DevicePort:
    """State-space port model of one device: input dV_xy, output dI_xy."""

    bus: int
    kind: str
    a: ComplexMatrix
    b: ComplexMatrix
    c: ComplexMatrix
    d: ComplexMatrix

    @property
    def n_states(self) -> int:
        return int(self.a.shape[0])

    def admittance(self, s: complex) -> ComplexMatrix:
        """Port contribution ``D + C (sI - A)^-1 B`` as seen from the bus."""

        if self.n_states == 0:
            return self.d.astype(np.complex128)
        identity = np.eye(self.n_states, dtype=np.complex128)
        return self.d + self.c @ np.linalg.solve(s * identity - self.a, self.b)


def linearize_device(slot, x: NDArray[np.float64], v: complex) -> DevicePort:
    """Central-difference port linearization of a single device."""

    device = slot.device
    n = device.n_states
    state = np.asarray(x[slot.start : slot.stop], dtype=np.float64)

    a = np.zeros((n, n))
    c = np.zeros((2, n))
    for j in range(n):
        h = _step(float(state[j]), 1.0)
        up, down = state.copy(), state.copy()
        up[j] += h
        down[j] -= h
        a[:, j] = (device.derivatives(up, v) - device.derivatives(down, v)) / (2.0 * h)
        delta = device.injection(up, v) - device.injection(down, v)
        c[0, j] = delta.real / (2.0 * h)
        c[1, j] = delta.imag / (2.0 * h)

    b = np.zeros((n, 2))
    d = np.zeros((2, 2))
    for j, direction in enumerate((1.0 + 0j, 1j)):
        h = _step(abs(v), 1.0)
        up, down = v + h * direction, v - h * direction
        if n:
            b[:, j] = (
                device.derivatives(state, up) - device.derivatives(state, down)
            ) / (2.0 * h)
        delta = device.injection(state, up) - device.injection(state, down)
        d[0, j] = delta.real / (2.0 * h)
        d[1, j] = delta.imag / (2.0 * h)

    return DevicePort(
        bus=slot.bus,
        kind=slot.kind,
        a=a.astype(np.complex128),
        b=b.astype(np.complex128),
        c=c.astype(np.complex128),
        d=d.astype(np.complex128),
    )


def load_admittance(load: complex, v: complex) -> ComplexMatrix:
    """Constant-power load, linearized as a 2x2 injected-current block."""

    d = np.zeros((2, 2))
    for j, direction in enumerate((1.0 + 0j, 1j)):
        h = _step(abs(v), 1.0)
        up, down = v + h * direction, v - h * direction
        delta = (-np.conj(load) / np.conj(up)) - (-np.conj(load) / np.conj(down))
        d[0, j] = delta.real / (2.0 * h)
        d[1, j] = delta.imag / (2.0 * h)
    return d.astype(np.complex128)


@dataclass(frozen=True)
class PortOperator:
    """The bus-space operator T(s) of one solved replacement case."""

    ybus_real: ComplexMatrix
    ports: tuple[DevicePort, ...]
    load_blocks: dict[int, ComplexMatrix]
    bus_index: dict[int, int]
    n_bus: int

    @property
    def dimension(self) -> int:
        return 2 * self.n_bus

    def bus_admittance(self, s: complex, bus: int) -> ComplexMatrix:
        """Total device port admittance at one bus."""

        total = np.zeros((2, 2), dtype=np.complex128)
        for port in self.ports:
            if port.bus == bus:
                total = total + port.admittance(s)
        return total

    def evaluate(self, s: complex) -> ComplexMatrix:
        t = self.ybus_real.copy()
        for bus, position in self.bus_index.items():
            block = np.zeros((2, 2), dtype=np.complex128)
            for port in self.ports:
                if port.bus == bus:
                    block = block + port.admittance(s)
            if bus in self.load_blocks:
                block = block + self.load_blocks[bus]
            if np.any(block):
                rows = slice(2 * position, 2 * position + 2)
                t[rows, rows] -= block
        return t

    def log_determinant(self, s: complex) -> complex:
        sign, magnitude = np.linalg.slogdet(self.evaluate(s))
        return complex(np.log(sign) + magnitude)

    def determinant(self, s: complex) -> complex:
        return complex(np.linalg.det(self.evaluate(s)))

    def smallest_singular_value(self, s: complex) -> float:
        return float(np.linalg.svd(self.evaluate(s), compute_uv=False)[-1])


def _real_ybus(ybus: ComplexMatrix) -> ComplexMatrix:
    """Expand a complex nodal admittance into rectangular 2x2 blocks."""

    n = ybus.shape[0]
    out = np.zeros((2 * n, 2 * n), dtype=np.complex128)
    for i in range(n):
        for j in range(n):
            y = ybus[i, j]
            out[2 * i, 2 * j] = y.real
            out[2 * i, 2 * j + 1] = -y.imag
            out[2 * i + 1, 2 * j] = y.imag
            out[2 * i + 1, 2 * j + 1] = y.real
    return out


def build_port_operator(
    dae: Ieee39Dae, x: NDArray[np.float64], z: NDArray[np.float64]
) -> PortOperator:
    """Assemble the port operator of a solved case at its equilibrium."""

    voltages = dae.voltages(z)
    network = dae.network
    ports = tuple(
        linearize_device(slot, x, complex(voltages[network.position(slot.bus)]))
        for slot in dae.slots
    )
    loads = {
        bus: load_admittance(load, complex(voltages[network.position(bus)]))
        for bus, load in network.loads.items()
    }
    return PortOperator(
        ybus_real=_real_ybus(network.ybus),
        ports=ports,
        load_blocks=loads,
        bus_index={bus: network.position(bus) for bus in network.bus_idx},
        n_bus=network.n_bus,
    )


@dataclass(frozen=True)
class PortActionSpace:
    """The replacement action operator of a portfolio, in port coordinates.

    With ``E_a`` the selector of the two port coordinates of bus ``a``,

        T_S(s) = T_0(s) + sum_a E_a dY_a(s) E_a^T
        K_ab(s) = E_a^T T_0(s)^-1 E_b
        M_ab(s) = dY_a(s) K_ab(s)

    and ``det(T_S)/det(T_0) = det(I + M)``.

    GAUGE. The port coordinates at each bus are a choice. Under an invertible
    per-port change ``E_a -> E_a S_a`` the update must transform as
    ``dY_a -> S_a^-1 dY_a S_a^-T`` to leave the operator unchanged, and then

        K_ab -> S_a^T K_ab S_b,      M -> S^-1 M S

    with ``S = blockdiag(S_a)``. M changes by SIMILARITY, and it is
    block-conformal, so each diagonal block transforms within itself. Hence
    ``det(I+M)``, every ``det(I+M_aa)``, their product, ``det(I+Q)`` and the
    eigenvalues of ``Q`` are gauge invariant, while ``sigma_min(I+Q)``, ``||Q||``
    and the entries of ``Q`` are not. Only the invariants may be reported.

    This is a stronger position than the abstract low-rank factorization allows:
    fixing the ports removes the freedom to redistribute between ``U`` and ``V``
    that would otherwise make the self/interaction split arbitrary.
    """

    t0: PortOperator
    ts: PortOperator
    buses: tuple[int, ...]
    bases: tuple[ComplexMatrix, ...] | None = None

    @property
    def order(self) -> int:
        return len(self.buses)

    @property
    def dimension(self) -> int:
        return 2 * self.order

    def _basis(self, index: int) -> ComplexMatrix:
        if self.bases is None:
            return np.eye(2, dtype=np.complex128)
        return self.bases[index]

    def selector(self) -> ComplexMatrix:
        n = self.t0.dimension
        u = np.zeros((n, self.dimension), dtype=np.complex128)
        for k, bus in enumerate(self.buses):
            p = self.t0.bus_index[bus]
            u[2 * p : 2 * p + 2, 2 * k : 2 * k + 2] = self._basis(k)
        return u

    def update(self, s: complex) -> ComplexMatrix:
        """Block-diagonal ``dY_a(s)`` in the current port basis."""

        block = np.zeros((self.dimension, self.dimension), dtype=np.complex128)
        for k, bus in enumerate(self.buses):
            raw = self.t0.bus_admittance(s, bus) - self.ts.bus_admittance(s, bus)
            basis = self._basis(k)
            inverse = np.linalg.inv(basis)
            block[2 * k : 2 * k + 2, 2 * k : 2 * k + 2] = inverse @ raw @ inverse.T
        return block

    def k(self, s: complex) -> ComplexMatrix:
        u = self.selector()
        return u.T @ np.linalg.solve(self.t0.evaluate(s), u)

    def m(self, s: complex) -> ComplexMatrix:
        return self.update(s) @ self.k(s)

    def block(self, matrix: ComplexMatrix, a: int, b: int) -> ComplexMatrix:
        return matrix[2 * a : 2 * a + 2, 2 * b : 2 * b + 2]

    def split(self, s: complex) -> dict[str, object]:
        m = self.m(s)
        identity = np.eye(self.dimension, dtype=np.complex128)
        total = identity + m
        self_block = np.zeros_like(total)
        individual = 1.0 + 0.0j
        diagonals = []
        for k in range(self.order):
            block = self.block(total, k, k)
            self_block[2 * k : 2 * k + 2, 2 * k : 2 * k + 2] = block
            value = complex(np.linalg.det(block))
            diagonals.append(value)
            individual *= value
        q = np.linalg.solve(self_block, total) - identity
        eigenvalues = np.linalg.eigvals(q)
        return {
            "full": complex(np.linalg.det(total)),
            "individual": individual,
            "diagonals": tuple(diagonals),
            "collective": complex(np.linalg.det(identity + q)),
            "q_eigenvalues": np.sort_complex(eigenvalues),
            "spectral_radius_q": float(np.max(np.abs(eigenvalues))),
            "closest_to_minus_one": complex(
                eigenvalues[int(np.argmin(np.abs(eigenvalues + 1.0)))]
            ),
            "sigma_min_i_plus_q": float(
                np.linalg.svd(identity + q, compute_uv=False)[-1]
            ),
            "norm_q": float(np.linalg.norm(q, 2)),
        }

    def polar_bases(self, voltages: dict[int, complex]) -> tuple[ComplexMatrix, ...]:
        """Magnitude and angle instead of rectangular coordinates.

        A physically motivated admissible port basis, not a random one: the
        Jacobian of (vx, vy) with respect to (|V|, theta) at the operating point.
        """

        out = []
        for bus in self.buses:
            v = voltages[bus]
            magnitude, angle = abs(v), float(np.angle(v))
            out.append(
                np.array(
                    [
                        [np.cos(angle), -magnitude * np.sin(angle)],
                        [np.sin(angle), magnitude * np.cos(angle)],
                    ],
                    dtype=np.complex128,
                )
            )
        return tuple(out)

    def with_bases(self, bases) -> PortActionSpace:
        return PortActionSpace(
            t0=self.t0,
            ts=self.ts,
            buses=self.buses,
            bases=tuple(np.asarray(b, dtype=np.complex128) for b in bases),
        )


def build_action_space(base_case, flagship_case, buses) -> PortActionSpace:
    """Port action space of one replacement portfolio against the base case."""

    return PortActionSpace(
        t0=build_port_operator(
            base_case.dae, base_case.equilibrium.x, base_case.equilibrium.z
        ),
        ts=build_port_operator(
            flagship_case.dae, flagship_case.equilibrium.x, flagship_case.equilibrium.z
        ),
        buses=tuple(buses),
    )


@dataclass(frozen=True)
class ClosureMargin:
    """Imaginary-axis port closure margin over a frozen frequency band.

        m = min over omega in the band of min over mu in spectrum(Q(j omega))
            of abs(mu + 1)

    Why the imaginary axis. Evaluating the closure at an eigenvalue of the
    REPLACED system returns zero identically, because that is what being an
    eigenvalue means. Evaluating it at an eigenvalue of the BASE system is
    inadmissible for the opposite reason: the operator is built on the resolvent
    of the base system and is singular there. The imaginary axis sits away from
    both spectra and is where an oscillatory stability boundary actually lives,
    so ``m -> 0`` is the natural signature of approaching such a boundary.

    ``omega_port`` is the frequency attaining the minimum. Near a boundary it
    should approach the imaginary part of the tracked mode, and that agreement is
    a stronger statement than the margin alone.
    """

    margin: float
    omega_port: float
    frequency_port_hz: float
    worst_condition: float
    samples: int

    @property
    def well_conditioned(self) -> bool:
        return bool(np.isfinite(self.worst_condition) and self.worst_condition < 1e10)


def imaginary_axis_closure(
    space: PortActionSpace,
    *,
    band_hz: tuple[float, float] = (0.3, 1.5),
    samples: int = 121,
) -> ClosureMargin:
    """Minimise the port closure distance over a frozen frequency band."""

    frequencies = np.linspace(band_hz[0], band_hz[1], samples)
    best = np.inf
    best_omega = float("nan")
    worst_condition = 0.0
    for frequency in frequencies:
        omega = 2.0 * np.pi * float(frequency)
        s = complex(0.0, omega)
        operator = space.t0.evaluate(s)
        worst_condition = max(worst_condition, float(np.linalg.cond(operator)))
        try:
            eigenvalues = space.split(s)["q_eigenvalues"]
        except (np.linalg.LinAlgError, ValueError):
            continue
        value = float(np.min(np.abs(np.asarray(eigenvalues) + 1.0)))
        if value < best:
            best = value
            best_omega = omega
    return ClosureMargin(
        margin=float(best),
        omega_port=float(best_omega),
        frequency_port_hz=float(best_omega / (2.0 * np.pi)),
        worst_condition=float(worst_condition),
        samples=int(samples),
    )


@dataclass(frozen=True)
class ClosureSample:
    """The closure margin at one frequency, with the numerics that produced it."""

    omega: float
    margin: float
    condition: float
    solve_residual: float
    """Backward residual ``||T0 X - U|| / (||T0|| ||X||)`` of the K(s) solve.

    ``K`` is formed by a linear solve against ``T_0``, never by an explicit
    inverse. This is the quantity that says whether a small margin is a property
    of the operator or of the solve.
    """

    @property
    def frequency_hz(self) -> float:
        return float(self.omega / (2.0 * np.pi))


def closure_at(space: PortActionSpace, omega: float) -> ClosureSample:
    """Closure margin on the imaginary axis at one frequency, with diagnostics."""

    s = complex(0.0, float(omega))
    operator = space.t0.evaluate(s)
    selector = space.selector()
    solution = np.linalg.solve(operator, selector)
    scale = float(np.linalg.norm(operator, 2)) * float(np.linalg.norm(solution, 2))
    residual = float(np.linalg.norm(operator @ solution - selector, 2))
    eigenvalues = np.asarray(space.split(s)["q_eigenvalues"])
    return ClosureSample(
        omega=float(omega),
        margin=float(np.min(np.abs(eigenvalues + 1.0))),
        condition=float(np.linalg.cond(operator)),
        solve_residual=residual / scale if scale > 0.0 else float("nan"),
    )


def closure_profile(
    space: PortActionSpace,
    *,
    band_hz: tuple[float, float],
    samples: int,
) -> tuple[ClosureSample, ...]:
    """The closure margin across a frozen frequency band."""

    frequencies = np.linspace(band_hz[0], band_hz[1], samples)
    out = []
    for frequency in frequencies:
        try:
            out.append(closure_at(space, 2.0 * np.pi * float(frequency)))
        except (np.linalg.LinAlgError, ValueError):
            continue
    return tuple(out)


def refine_closure(
    space: PortActionSpace,
    omega: float,
    *,
    half_width: float,
    tolerance: float = 1e-10,
) -> ClosureSample:
    """Golden-section refinement of the margin around a grid minimum.

    The grid minimum is an upper bound on the true band minimum. This says how
    much of the reported margin is grid resolution rather than the operator.
    """

    golden = 0.5 * (3.0 - np.sqrt(5.0))
    low, high = omega - half_width, omega + half_width
    left = low + golden * (high - low)
    right = high - golden * (high - low)
    f_left, f_right = closure_at(space, left), closure_at(space, right)
    while high - low > tolerance * max(abs(omega), 1.0):
        if f_left.margin <= f_right.margin:
            high, right, f_right = right, left, f_left
            left = low + golden * (high - low)
            f_left = closure_at(space, left)
        else:
            low, left, f_left = left, right, f_right
            right = high - golden * (high - low)
            f_right = closure_at(space, right)
    return f_left if f_left.margin <= f_right.margin else f_right
