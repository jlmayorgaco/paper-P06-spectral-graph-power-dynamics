"""Operating-point sampling with a physical dispatch.

The v1 held-out campaign scaled every generator by the load multiplier. That is
not a dispatch: buses 31, 32, 35 and 36 sit within 0.04 to 0.25 pu of their
active-power limits, so a uniform scaling pushed them past those limits and 29 %
of samples were rejected for a reason that belongs to the sampler and not to the
power system.

The dispatch here does four things the v1 sampler did not.

1. Photovoltaic availability caps the RESOURCE, not the schedule: a plant of
   capacity ``pmax`` produces ``min(schedule, pmax * availability)``.
2. The load change and the photovoltaic shortfall are absorbed by the machines
   that have headroom, in proportion to it, with ITERATIVE SATURATION: a unit
   that hits its limit drops out of the pool and its share is redistributed.
3. The slack takes only the final residual and the losses. It is not a default
   balancing generator, and its own active limit is checked.
4. Reactive limits are audited. Correcting the active dispatch and then building
   high-load points whose real infeasibility is reactive would be an easy way to
   fool ourselves, so every generator that reaches ``qmax`` or ``qmin`` in the
   solved power flow is recorded.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from numpy.typing import NDArray

from ..models.ieee39_network import Ieee39Network, PowerFlow

#: Active-power losses are not known before the power flow is solved; this
#: allowance is used only to pre-screen a sample as dispatchable.
LOSS_ALLOWANCE = 0.03


@dataclass(frozen=True)
class OperatingPoint:
    """One sampled operating point and the evidence that it is dispatchable."""

    network: Ieee39Network
    active_load: float
    reactive_load: float
    availability: float
    dispatched: dict[int, float]
    saturated: tuple[int, ...]
    curtailed: tuple[int, ...]
    unserved_pu: float
    slack_estimate_pu: float
    slack_min_pu: float
    slack_max_pu: float

    @property
    def dispatchable(self) -> bool:
        """Whether the fleet can meet this load within its active limits."""

        return (
            abs(self.unserved_pu) < 1e-6
            and self.slack_min_pu - 1e-9 <= self.slack_estimate_pu
            and self.slack_estimate_pu <= self.slack_max_pu + 1e-9
        )

    @property
    def total_load_pu(self) -> float:
        return float(sum(load.real for load in self.network.loads.values()))

    @property
    def total_generation_pu(self) -> float:
        return float(
            sum(spec["p"] for spec in self.network.pv.values())
            + self.slack_estimate_pu
        )

    @property
    def power_balance_residual_pu(self) -> float:
        """Scheduled generation minus load, before losses are known."""

        return self.total_generation_pu - self.total_load_pu * (1.0 + LOSS_ALLOWANCE)


def sample_operating_point(
    base: Ieee39Network,
    *,
    active_load: float,
    reactive_load: float,
    availability: float,
    availability_buses: tuple[int, ...],
    rng: np.random.Generator,
    load_scatter: float = 0.05,
    dispatch_jitter: float = 0.10,
) -> OperatingPoint:
    """Build one operating point with a saturating headroom-proportional redispatch."""

    loads: dict[int, complex] = {}
    for bus, load in base.loads.items():
        jitter = 1.0 + rng.uniform(-load_scatter, load_scatter)
        loads[bus] = complex(
            load.real * active_load * jitter, load.imag * reactive_load * jitter
        )
    load_change = sum(load.real for load in loads.values()) - sum(
        load.real for load in base.loads.values()
    )

    dispatched = {bus: spec["p"] for bus, spec in base.pv.items()}
    curtailed: list[int] = []
    shortfall = 0.0
    for bus in availability_buses:
        if bus not in dispatched:
            continue
        capped = min(dispatched[bus], base.pv[bus]["pmax"] * availability)
        if capped < dispatched[bus] - 1e-12:
            curtailed.append(bus)
        shortfall += dispatched[bus] - capped
        dispatched[bus] = capped

    # Iterative saturation: a unit at its limit leaves the pool and its share is
    # redistributed among the units that still have room.
    pool = [b for b in base.pv if b not in availability_buses]
    remaining = load_change + shortfall
    saturated: list[int] = []
    weights = {b: 1.0 + rng.uniform(-dispatch_jitter, dispatch_jitter) for b in pool}
    for _ in range(len(pool) + 2):
        if abs(remaining) < 1e-10 or not pool:
            break
        if remaining > 0:
            room = {b: max(base.pv[b]["pmax"] - dispatched[b], 0.0) for b in pool}
        else:
            room = {b: max(dispatched[b] - base.pv[b]["pmin"], 0.0) for b in pool}
        total = sum(room[b] * weights[b] for b in pool)
        if total <= 1e-12:
            break
        step = remaining
        hit: list[int] = []
        for bus in pool:
            share = room[bus] * weights[bus] / total
            proposed = dispatched[bus] + share * step
            limited = float(
                np.clip(proposed, base.pv[bus]["pmin"], base.pv[bus]["pmax"])
            )
            remaining -= limited - dispatched[bus]
            dispatched[bus] = limited
            if abs(limited - proposed) > 1e-12:
                hit.append(bus)
        for bus in hit:
            pool.remove(bus)
            saturated.append(bus)

    pv = {bus: {**spec, "p": dispatched[bus]} for bus, spec in base.pv.items()}
    total_load = sum(load.real for load in loads.values())
    slack_estimate = total_load * (1.0 + LOSS_ALLOWANCE) - sum(
        spec["p"] for spec in pv.values()
    )
    return OperatingPoint(
        network=replace(base, loads=loads, pv=pv),
        active_load=active_load,
        reactive_load=reactive_load,
        availability=availability,
        dispatched=dispatched,
        saturated=tuple(sorted(saturated)),
        curtailed=tuple(sorted(curtailed)),
        unserved_pu=float(remaining) if abs(remaining) > 1e-10 else 0.0,
        slack_estimate_pu=float(slack_estimate),
        slack_min_pu=0.0,
        slack_max_pu=float(base.slack_pmax),
    )


@dataclass(frozen=True)
class LimitAudit:
    """Which generators reach a limit in the solved power flow."""

    slack_active_pu: float
    slack_within_limits: bool
    reactive_binding: tuple[tuple[int, str], ...]

    @property
    def clean(self) -> bool:
        return self.slack_within_limits and not self.reactive_binding


def audit_limits(
    network: Ieee39Network, power_flow: PowerFlow, *, tolerance: float = 1e-4
) -> LimitAudit:
    """Audit the solved point against the active and reactive limits.

    The power flow holds voltage at the regulating buses without enforcing
    reactive limits, so a point can be active-feasible and reactive-infeasible.
    Recording which units bind is the difference between knowing that and not.
    """

    binding: list[tuple[int, str]] = []
    for bus, spec in network.pv.items():
        reactive = (
            power_flow.injection(bus, network.ybus) + network.loads.get(bus, 0j)
        ).imag
        if reactive > spec["qmax"] - tolerance:
            binding.append((bus, "qmax"))
        elif reactive < spec["qmin"] + tolerance:
            binding.append((bus, "qmin"))
    slack = (
        power_flow.injection(network.slack_bus, network.ybus)
        + network.loads.get(network.slack_bus, 0j)
    ).real
    return LimitAudit(
        slack_active_pu=float(slack),
        slack_within_limits=bool(0.0 <= slack <= network.slack_pmax + tolerance),
        reactive_binding=tuple(binding),
    )


def stratified_grid(
    load_strata: NDArray[np.float64],
    availability_strata: NDArray[np.float64],
    per_cell: int,
    rng: np.random.Generator,
) -> list[tuple[float, float, float]]:
    """Latin-hypercube fill inside every (load, availability) cell.

    Stratifying rather than sampling the whole envelope is the point of v2A: the
    v1 campaign averaged over a range in which the effect changes sign, so its
    marginal statistics could not see the interaction.
    """

    points: list[tuple[float, float, float]] = []
    for i in range(len(load_strata) - 1):
        for j in range(len(availability_strata) - 1):
            unit = (rng.permutation(per_cell) + rng.random(per_cell)) / per_cell
            other = (rng.permutation(per_cell) + rng.random(per_cell)) / per_cell
            reactive = rng.uniform(0.90, 1.10, size=per_cell)
            for k in range(per_cell):
                load = load_strata[i] + unit[k] * (load_strata[i + 1] - load_strata[i])
                availability = availability_strata[j] + other[k] * (
                    availability_strata[j + 1] - availability_strata[j]
                )
                points.append((float(load), float(reactive[k]), float(availability)))
    return points
