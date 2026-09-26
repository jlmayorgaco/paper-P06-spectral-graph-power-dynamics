"""Reduced port-core evaluation of every replacement subset at once.

With the operating point shared by every subset (the matched dispatch), replacing
the machines of a subset S changes the bus operator only at S's ports:

    T_S(s) = T_0(s) + U_S dY_S(s) U_S^T,     dY_a = Y_machine,a - Y_converter,a

so, with K(s) = U^T T_0(s)^-1 U over ALL candidate ports and M(s) = dY(s) K(s),

    det T_S(s) / det T_0(s) = det( I + M(s)[S, S] ),

a principal minor of ONE matrix. This is the return difference of the portfolio
interconnection; its eigenvalues reaching -1 on the imaginary axis is the
classical generalized Nyquist condition, which is not claimed as new here.

Counting. ``det T`` has zeros at the port-visible system eigenvalues and poles at
the open-loop device eigenvalues (each device with its terminal voltage held).
The argument principle on the boundary of the truncated band region

    Gamma_R = { 0 < Re s < sigma_R,  omega_lo <= Im s <= omega_hi }

therefore gives, for every subset,

    N_Gamma(S) = N_Gamma(0) + wind(det(I + M_SS)) + P_conv(S) - P_mach(S),

P the open-loop device poles inside Gamma_R at S's buses. Modes hidden inside a
device (uncontrollable or unobservable at its port) are invisible to this count
by construction; the comparison with the full eigenproblem measures that.

Valid ONLY where the operating point is common to every subset. Where a policy
re-equilibrates per subset, T_S is not a low-rank update of T_0 and this module
must not be used.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import numpy as np
from numpy.typing import NDArray

from .port_admittance import DevicePort, PortOperator, build_port_operator

Complex = NDArray[np.complex128]


def device_admittance(port: DevicePort, s: Complex) -> Complex:
    """``D + C (sI - A)^-1 B`` for an array of s, shape (len(s), 2, 2)."""

    s = np.asarray(s, dtype=np.complex128)
    if port.n_states == 0:
        return np.broadcast_to(port.d, (s.size, 2, 2)).copy()
    eye = np.eye(port.n_states, dtype=np.complex128)
    pencil = s[:, None, None] * eye - port.a
    x = np.linalg.solve(pencil, np.broadcast_to(port.b, (s.size, *port.b.shape)))
    return port.d + port.c @ x


@dataclass
class PortCore:
    base: PortOperator
    replaced: PortOperator
    buses: tuple[int, ...]

    def __post_init__(self) -> None:
        n = self.base.dimension
        const = self.base.ybus_real.copy()
        for bus, block in self.base.load_blocks.items():
            p = self.base.bus_index[bus]
            const[2 * p : 2 * p + 2, 2 * p : 2 * p + 2] -= block
        self._const = const
        self._selector = np.zeros((n, 2 * len(self.buses)), dtype=np.complex128)
        for k, bus in enumerate(self.buses):
            p = self.base.bus_index[bus]
            self._selector[2 * p : 2 * p + 2, 2 * k : 2 * k + 2] = np.eye(2)

    @classmethod
    def from_daes(
        cls, base_dae, base_x, base_z, replaced_dae, replaced_x, replaced_z, buses
    ):
        return cls(
            build_port_operator(base_dae, base_x, base_z),
            build_port_operator(replaced_dae, replaced_x, replaced_z),
            tuple(buses),
        )

    # --------------------------------------------------------------- operators
    def t0(self, s: Complex) -> Complex:
        s = np.asarray(s, dtype=np.complex128)
        t = np.broadcast_to(self._const, (s.size, *self._const.shape)).copy()
        for port in self.base.ports:
            p = self.base.bus_index[port.bus]
            t[:, 2 * p : 2 * p + 2, 2 * p : 2 * p + 2] -= device_admittance(port, s)
        return t

    def _bus_y(self, operator: PortOperator, bus: int, s: Complex) -> Complex:
        out = np.zeros((s.size, 2, 2), dtype=np.complex128)
        for port in operator.ports:
            if port.bus == bus:
                out += device_admittance(port, s)
        return out

    def delta_y(self, s: Complex) -> Complex:
        s = np.asarray(s, dtype=np.complex128)
        m = len(self.buses)
        out = np.zeros((s.size, 2 * m, 2 * m), dtype=np.complex128)
        for k, bus in enumerate(self.buses):
            out[:, 2 * k : 2 * k + 2, 2 * k : 2 * k + 2] = self._bus_y(
                self.base, bus, s
            ) - self._bus_y(self.replaced, bus, s)
        return out

    def m_matrix(self, s: Complex) -> Complex:
        s = np.asarray(s, dtype=np.complex128)
        u = self._selector
        x = np.linalg.solve(self.t0(s), np.broadcast_to(u, (s.size, *u.shape)))
        k = np.swapaxes(u, 0, 1)[None] @ x
        return self.delta_y(s) @ k

    def subset_dets(self, s: Complex, subsets) -> Complex:
        """det(I + M[S,S]) for every subset, shape (len(s), len(subsets))."""

        m = self.m_matrix(s)
        index = {b: k for k, b in enumerate(self.buses)}
        out = np.empty((m.shape[0], len(subsets)), dtype=np.complex128)
        for j, subset in enumerate(subsets):
            rows = [r for b in subset for r in (2 * index[b], 2 * index[b] + 1)]
            if not rows:
                out[:, j] = 1.0
                continue
            block = m[:, rows][:, :, rows]
            out[:, j] = np.linalg.det(np.eye(len(rows)) + block)
        return out

    # ------------------------------------------------------------------ poles
    def device_poles(self, region) -> dict[int, tuple[int, int]]:
        """(machine poles, converter poles) inside Gamma_R at every candidate bus."""

        sigma, lo, hi = region
        out = {}
        for bus in self.buses:
            counts = []
            for operator in (self.base, self.replaced):
                n = 0
                for port in operator.ports:
                    if port.bus != bus or port.n_states == 0:
                        continue
                    ev = np.linalg.eigvals(port.a)
                    n += int(
                        np.count_nonzero(
                            (ev.real > 0)
                            & (ev.real < sigma)
                            & (ev.imag >= lo)
                            & (ev.imag <= hi)
                        )
                    )
                counts.append(n)
            out[bus] = (counts[0], counts[1])
        return out

    # --------------------------------------------------------------- counting
    def counts(
        self,
        subsets,
        *,
        omega_lo: float,
        omega_hi: float,
        sigma_r: float = 50.0,
        max_phase: float = 0.35,
        max_passes: int = 14,
    ):
        """N_Gamma(S) - N_Gamma(0) for every subset by the argument principle."""

        corners = [
            complex(0, omega_lo),
            complex(sigma_r, omega_lo),
            complex(sigma_r, omega_hi),
            complex(0, omega_hi),
            complex(0, omega_lo),
        ]
        path = [
            np.linspace(a, b, 48, endpoint=False)
            for a, b in zip(corners[:-1], corners[1:], strict=True)
        ]
        s = np.concatenate(path + [np.array([corners[0]])])
        d = self.subset_dets(s, subsets)
        for _ in range(max_passes):
            step = np.angle(d[1:] / d[:-1])
            bad = np.flatnonzero(np.abs(step).max(axis=1) > max_phase)
            if bad.size == 0:
                break
            mids = 0.5 * (s[bad] + s[bad + 1])
            dm = self.subset_dets(mids, subsets)
            s = np.insert(s, bad + 1, mids)
            d = np.insert(d, bad + 1, dm, axis=0)
        step = np.angle(d[1:] / d[:-1])
        winding = step.sum(axis=0) / (2 * np.pi)
        poles = self.device_poles((sigma_r, omega_lo, omega_hi))
        out = {}
        for j, subset in enumerate(subsets):
            w = float(winding[j])
            correction = sum(poles[b][1] - poles[b][0] for b in subset)
            out[tuple(subset)] = {
                "winding": w,
                "winding_residual": abs(w - round(w)),
                "delta_N": int(round(w)) + correction,
                "pole_correction": correction,
            }
        return out, {
            "contour_points": int(s.size),
            "max_final_phase_step": float(np.abs(step).max()),
        }


def all_subsets(buses) -> list[tuple[int, ...]]:
    return [tuple(c) for k in range(len(buses) + 1) for c in combinations(buses, k)]
