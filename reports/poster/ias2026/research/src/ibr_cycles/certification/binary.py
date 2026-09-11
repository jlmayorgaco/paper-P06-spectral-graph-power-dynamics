"""BC02: representations of the binary replacement family.

Three different statements, kept apart:

A. For a FIXED portfolio S, affinity of A_red in the continuous control
   parameters (F11/G3 exact assembly). Says nothing about the indicators delta.
B. Port locality at a common equilibrium: T_S(s) = T_0(s) + sum_{i in S} dT_i(s),
   where dT_i is supported on the port of bus i.
C. A common dynamic realization that is affine in the indicators:
   J(delta) = J_0 + sum_i delta_i J_i on one fixed state dimension.

C is constructed here, not assumed. Every candidate bus carries BOTH devices:
the machine with weight (1 - delta_i) Sn/100 and the converter with weight
delta_i * rating. The absent device is kept at its intensive (per-unit)
equilibrium, which is the same at every vertex under matched dispatch. Weights
multiply only the injections, so the DAE Jacobian [f_x f_z; g_x g_z] is affine
in delta. The absent device's own rows remain:

  C1 (ghost)   its true dynamics, driven by the bus voltage, injecting nothing.
               Its spectrum is added to every vertex where it is absent.
  C2 (filler)  its rows replaced by -(x - x0). The same algebra, delta-affine;
               it adds eigenvalue -1 with known multiplicity.

The REDUCED matrix A(delta) = f_x - f_z g_z^-1 g_x is in general NOT affine,
because g_z depends on delta. The checks here measure how far it is.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from itertools import combinations

import numpy as np

from ..dynamics.linearize import central_difference_jacobians
from ..models.ieee39_case import ReplacementPlan, solve_case


class FillerDevice:
    """Stable filler with the state dimension of a replaced device: x' = -(x - x0)."""

    def __init__(self, device, x0):
        self.bus = device.bus
        self.weight = 0.0
        self._labels = tuple(device.labels)
        self._x0 = np.array(x0, dtype=np.float64)

    @property
    def labels(self):
        return self._labels

    @property
    def n_states(self):
        return len(self._labels)

    def derivatives(self, x, v):
        return -(np.asarray(x) - self._x0)

    def injection(self, x, v):
        return 0j


@dataclass
class CommonRealization:
    case: object  # the solved rho = 0.5 case (both devices present at every candidate)
    candidates: tuple[int, ...]
    slot_index: dict  # bus -> (sg slot position, gfl slot position)

    @property
    def x0(self):
        return self.case.equilibrium.x

    @property
    def z0(self):
        return self.case.equilibrium.z

    def dae_at(self, delta: dict[int, int], variant: str = "C1"):
        dae = self.case.dae
        slots = list(dae.slots)
        for bus in self.candidates:
            d = int(delta.get(bus, 0))
            i_sg, i_gfl = self.slot_index[bus]
            for pos, factor in ((i_sg, 2.0 * (1 - d)), (i_gfl, 2.0 * d)):
                slot = slots[pos]
                device = slot.device
                if factor == 0.0 and variant == "C2":
                    new = FillerDevice(device, self.x0[slot.start : slot.stop])
                else:
                    new = replace(device, weight=device.weight * factor)
                slots[pos] = replace(slot, device=new, weight=slot.weight * factor)
        return replace(dae, slots=tuple(slots))

    def jacobian(self, delta: dict[int, int], variant: str = "C1"):
        dae = self.dae_at(delta, variant)
        jac = central_difference_jacobians(dae, self.x0, self.z0, {})
        residual = max(
            float(np.abs(dae.f(self.x0, self.z0, {})).max()),
            float(np.abs(dae.g(self.x0, self.z0, {})).max()),
        )
        return jac, residual, dae

    def absent_blocks(self, delta: dict[int, int]):
        """(bus, kind, eigenvalues) of every absent device's own dynamics at fixed v."""

        dae = self.case.dae
        out = []
        jac = central_difference_jacobians(dae, self.x0, self.z0, {})
        for bus in self.candidates:
            d = int(delta.get(bus, 0))
            i_sg, i_gfl = self.slot_index[bus]
            pos = i_gfl if d == 0 else i_sg
            slot = dae.slots[pos]
            block = jac.fx[slot.start : slot.stop, slot.start : slot.stop]
            out.append((bus, slot.kind, np.linalg.eigvals(block)))
        return out


def build_common_realization(candidates, **solve_kwargs) -> CommonRealization:
    plan = ReplacementPlan.of({b: 0.5 for b in candidates})
    case = solve_case(plan, **solve_kwargs)
    index = {}
    for bus in candidates:
        sg = [
            k for k, s in enumerate(case.dae.slots) if s.bus == bus and s.kind == "sg"
        ]
        gfl = [
            k for k, s in enumerate(case.dae.slots) if s.bus == bus and s.kind == "gfl"
        ]
        index[bus] = (sg[0], gfl[0])
    return CommonRealization(case=case, candidates=tuple(candidates), slot_index=index)


def stacked(jac) -> np.ndarray:
    return np.block([[jac.fx, jac.fz], [jac.gx, jac.gz]])


def reduced(jac) -> np.ndarray:
    return jac.fx - jac.fz @ np.linalg.solve(jac.gz, jac.gx)


def all_vertices(candidates):
    for r in range(len(candidates) + 1):
        for s in combinations(candidates, r):
            yield {b: (1 if b in s else 0) for b in candidates}, s
