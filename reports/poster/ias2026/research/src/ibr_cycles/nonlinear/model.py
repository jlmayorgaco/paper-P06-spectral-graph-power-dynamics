"""NL01: one model object that generates residuals, derivatives, ports and TDS.

It wraps the L0 DAE of a solved case (src/ibr_cycles/models/ieee39_case.py:
Ieee39Dae, the same code the frozen experiments ran). It does not re-implement
the equations, so there is one formulation for TDS, Jacobians and ports.

Contract (every method returns arrays with names, units and the model version):

    residual_f(x, z, u, p)   F(x, z, u) with M F = f; in L0 M = I (each device
                             already returns the normalized derivative)
    residual_g(x, z, u, p)   KCL residual, interleaved real coordinates
    mass(x, z, p)            identity (L0), declared, not inferred
    outputs(x, z, u)         |V| per bus, machine speed, device P and Q
    initialize_from_pf(...)  the frozen solve_case path (AC power flow + devices)
    symmetry_generators()    rotation R and neutral-frequency partner w
    active_set(x, z)         L0 has no limiters: the empty set, declared ABSENT
    guard_values(x, z)       none in L0 (ABSENT)
    reset_map(...)           NOT_APPLICABLE in L0 (no discrete events)

Inputs u (declared physical perturbations, all zero at the operating point):
    u = [dP_load(bus b) for each load bus] + [dQ_load(bus b)] + [dPm(machine)]
    in system per unit; they enter L0 by editing the load and the machine's
    mechanical power, exactly as the frozen TDS did (E32/G2).
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from numpy.typing import NDArray

from ..certification.symmetry import frequency_partner, rotation_generator
from . import MODEL_VERSION

ABSENT = "ABSENT"
NOT_APPLICABLE = "NOT_APPLICABLE"


@dataclass(frozen=True)
class Named:
    values: NDArray
    names: tuple[str, ...]
    units: tuple[str, ...]
    version: str = MODEL_VERSION


def _state_unit(label: str) -> str:
    base = label.rsplit("_", 1)[0]
    return {
        "delta": "rad",
        "theta_pll": "rad",
        "omega": "pu speed (68-bus: slip)",
        "x_pll": "rad/s",
    }.get(base, "pu (device base)")


class PhasorModel:
    """L0/L1 phasor DAE of one solved portfolio, with declared physical inputs."""

    def __init__(self, case):
        self.case = case
        self.dae = case.dae
        net = self.dae.network
        self.load_buses = tuple(sorted(net.loads))
        self.machines = tuple(
            (k, s.bus) for k, s in enumerate(self.dae.slots) if s.kind == "sg"
        )
        self.input_names = (
            tuple(f"dP_load_{b}" for b in self.load_buses)
            + tuple(f"dQ_load_{b}" for b in self.load_buses)
            + tuple(f"dPm_sg{b}" for _, b in self.machines)
        )
        self.n_u = len(self.input_names)

    # ------------------------------------------------------------ inputs --
    def _with_inputs(self, u):
        """A DAE copy with u applied (loads and mechanical powers)."""

        if u is None or not np.any(u):
            return self.dae
        nl = len(self.load_buses)
        net = self.dae.network
        loads = dict(net.loads)
        for k, b in enumerate(self.load_buses):
            loads[b] = loads[b] + complex(u[k], u[nl + k])
        new_net = replace(net, loads=loads)
        slots = list(self.dae.slots)
        for j, (pos, _bus) in enumerate(self.machines):
            dpm = float(u[2 * nl + j])
            if dpm:
                slot = slots[pos]
                dev = slot.device
                p = dev.parameters
                if hasattr(p, "pm"):
                    newp = replace(p, pm=p.pm + dpm / dev.weight)
                else:  # 68-bus machine: mechanical torque on the machine base
                    newp = replace(p, tm=p.tm + dpm / dev.weight)
                slots[pos] = replace(slot, device=replace(dev, parameters=newp))
        return replace(self.dae, network=new_net, slots=tuple(slots))

    # --------------------------------------------------------- contracts --
    def residual_f(self, x, z, u=None, p=None) -> Named:
        f = self._with_inputs(u).f(x, z, {})
        return Named(f, self.dae.labels, tuple(_state_unit(n) for n in self.dae.labels))

    def residual_g(self, x, z, u=None, p=None) -> Named:
        g = self._with_inputs(u).g(x, z, {})
        names = tuple(
            f"{c}_bus{b}"
            for b in self.dae.network.bus_idx
            for c in ("KCL_re", "KCL_im")
        )
        return Named(g, names, ("pu current (system base)",) * len(names))

    def mass(self, x=None, z=None, p=None) -> Named:
        n = self.dae.n_x
        return Named(np.eye(n), self.dae.labels, ("dimensionless",) * n)

    def outputs(self, x, z, u=None) -> Named:
        dae = self._with_inputs(u)
        v = dae.voltages(z)
        vals, names, units = [], [], []
        for k, b in enumerate(dae.network.bus_idx):
            vals.append(abs(v[k]))
            names.append(f"Vmag_bus{b}")
            units.append("pu")
        for slot in dae.slots:
            pos = dae.network.position(slot.bus)
            i = slot.device.injection(x[slot.start : slot.stop], complex(v[pos]))
            s = v[pos] * np.conj(i)
            vals += [s.real, s.imag]
            names += [f"P_{slot.kind}{slot.bus}", f"Q_{slot.kind}{slot.bus}"]
            units += ["pu (system base)"] * 2
        return Named(np.array(vals), tuple(names), tuple(units))

    def initialize_from_pf(self):
        return self.case.equilibrium.x, self.case.equilibrium.z

    def symmetry_generators(self):
        r_x, r_z = rotation_generator(self.dae, self.case.equilibrium.z)
        return {"R_x": r_x, "R_z": r_z, "partner": frequency_partner(self.dae)}

    def active_set(self, x=None, z=None):
        return {"status": ABSENT, "active": ()}

    def guard_values(self, x=None, z=None) -> Named:
        return Named(np.zeros(0), (), ())

    def reset_map(self, *args, **kwargs):
        return NOT_APPLICABLE
