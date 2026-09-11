"""ieee39_governed_documented_v1: the frozen IEEE-39 DAE plus its DOCUMENTED governors.

A separately versioned model (Decision 2). The frozen benchmark is not modified:
this module wraps a solved frozen case and adds, to every synchronous machine
still present, the TGOV1N turbine governor whose parameters are documented in
the SAME source workbook (data/raw/ieee39_full.xlsx, sheet TGOV1N;
configs/ias2026/ieee39_governed_documented_v1.json). Nothing is tuned; no
damping is added (GENROU D = 0 in the source).

TGOV1N (andes.models.governor.tgov1.TGOV1N), machine base (Tn = Sn):

    pd  = pref - (omega - 1) / R
    T1 dx1/dt = pd - x1                    (anti-windup limits [VMIN, VMAX])
    T3 dx2/dt = x1 - x2
    y   = (T2 / T3) (x1 - x2) + x2
    Pm  = y - Dt (omega - 1),              pref = Pm0

The valve limits are not active at any equilibrium used here (x1 = Pm0 in
[VMIN, VMAX], checked); small-signal analysis does not see them. A condenser
(no active power) carries no governor.

The rotation symmetry is unchanged (governors act on speed deviations, not on
angles): A R_x = 0 still holds. The uniform-frequency shift is no longer a
Jordan partner: the governors restore the common frequency, so the structural
center subspace is C = span{R_x}, of dimension one.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

MODEL_VERSION = "ieee39_governed_documented_v1"
CONFIG = (
    Path(__file__).resolve().parents[3]
    / "configs"
    / "ias2026"
    / "ieee39_governed_documented_v1.json"
)


@dataclass(frozen=True)
class Tgov1n:
    r: float
    t1: float
    t2: float
    t3: float
    dt: float
    vmax: float
    vmin: float


def documented_governors(path: Path = CONFIG) -> dict[int, Tgov1n]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return {
        int(g["bus"]): Tgov1n(g["R"], g["T1"], g["T2"], g["T3"], g["Dt"], g["VMAX"], g["VMIN"])
        for g in payload["governors"]
    }


@dataclass
class GovernedMachine:
    base: object  # SynchronousMachine of the frozen model
    gov: Tgov1n
    pref: float

    @property
    def bus(self) -> int:
        return self.base.bus

    @property
    def weight(self) -> float:
        return self.base.weight

    @property
    def labels(self) -> tuple[str, ...]:
        b = self.base.bus
        return (*self.base.labels, f"gov_lag_sg{b}", f"gov_ll_sg{b}")

    @property
    def n_states(self) -> int:
        return self.base.n_states + 2

    def mechanical(self, x) -> float:
        n0 = self.base.n_states
        omega = float(x[1])
        x1, x2 = float(x[n0]), float(x[n0 + 1])
        y = (self.gov.t2 / self.gov.t3) * (x1 - x2) + x2
        return y - self.gov.dt * (omega - 1.0)

    def derivatives(self, x, v):
        n0 = self.base.n_states
        omega = float(x[1])
        x1, x2 = float(x[n0]), float(x[n0 + 1])
        pd = self.pref - (omega - 1.0) / self.gov.r
        dev = replace(self.base, parameters=replace(self.base.parameters, pm=self.mechanical(x)))
        return np.concatenate(
            [
                dev.derivatives(np.asarray(x[:n0]), v),
                [(pd - x1) / self.gov.t1, (x1 - x2) / self.gov.t3],
            ]
        )

    def injection(self, x, v):
        return self.base.injection(np.asarray(x[: self.base.n_states]), v)


def govern(case, governors: dict[int, Tgov1n] | None = None):
    """(dae, x0, z0) of the governed version of a solved frozen case."""

    from ..models.ieee39_case import DeviceSlot

    governors = governors or documented_governors()
    dae = case.dae
    x_old = case.equilibrium.x
    slots, states, cursor = [], [], 0
    for slot in dae.slots:
        dev = slot.device
        xs = x_old[slot.start : slot.stop]
        params = getattr(dev, "parameters", None)
        is_generator = slot.kind == "sg" and params is not None and hasattr(params, "pm")
        if is_generator and abs(params.pm) > 1e-9 and slot.bus in governors:
            gov = governors[slot.bus]
            pref = float(params.pm)
            if not (gov.vmin <= pref <= gov.vmax):
                raise ValueError(f"valve limit active at bus {slot.bus}: Pm0 = {pref:.3f}")
            dev = GovernedMachine(dev, gov, pref)
            xs = np.concatenate([xs, [pref, pref]])
        slots.append(replace(slot, device=dev, start=cursor, stop=cursor + xs.size))
        states.append(xs)
        cursor += xs.size
    new = replace(dae, slots=tuple(slots))
    new.__post_init__()
    x0 = np.concatenate(states)
    z0 = case.equilibrium.z.copy()
    return new, x0, z0
