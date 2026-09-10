"""Reactive-power policies for the engineering re-equilibration study.

The mechanism-isolation campaigns hold the operating point fixed by giving the
converter the reactive output the machine was producing. That is the right
control for asking whether the effect is caused by the replacement rather than by
a moved power flow, and it is the wrong thing to ship in a paper on its own: a
reviewer will ask what happens when the plant is dispatched the way a real plant
is dispatched.

Four policies, all at full replacement:

``matched``
    the frozen mechanism-isolation reference. The bus injection is unchanged, so
    the AC solution is identical to the base case.

``unity_pf``
    the converter carries no reactive power at all. The bus becomes PQ and the
    operating point genuinely moves.

``vreg_pf095``
    the converter regulates the bus to its pre-replacement magnitude, with a
    reactive capability of +-0.3287 of its rating, the range of a plant rated to
    0.95 power factor. A plant that would need more is clamped to the limit and
    runs on a fixed reactive command from there.

``vreg_scap``
    the same regulation, with capability set instead by apparent-power
    saturation, ``Q_max = sqrt(S^2 - P^2)`` on the plant rating. A plant carrying
    its full megawatts has no reactive headroom left, which is the honest version
    of the constraint.

Both regulating policies enforce a limit. Neither allows unlimited reactive
power. The capability constants are frozen in
``configs/ias2026/overnight_policies.yaml`` before any policy case was solved.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from .ieee39_case import (
    SYSTEM_BASE_MVA,
    InfeasibleReplacement,
    ReplacementCase,
    ReplacementPlan,
    solve_case,
)
from .ieee39_devices import ConverterParameters
from .ieee39_network import Ieee39Network, load_network, solve_power_flow

#: Reactive capability of a plant rated to 0.95 power factor, per unit of rating.
POWER_FACTOR_CAPABILITY = 0.3287
#: Iterations allowed while switching regulating plants onto their limit.
MAX_CLAMP_PASSES = 12

POLICIES = ("matched", "unity_pf", "vreg_pf095", "vreg_scap")


@dataclass(frozen=True)
class PolicyOutcome:
    """A solved policy case plus the reactive evidence a reviewer will ask for."""

    case: ReplacementCase
    policy: str
    members: tuple[int, ...]
    q_required_pu: dict[int, float]
    q_capability_pu: dict[int, float]
    utilisation: dict[int, float]
    """Apparent power carried over the plant rating, per replaced bus."""
    saturated: tuple[int, ...]
    voltage_regulating: tuple[int, ...]
    clamp_passes: int

    @property
    def worst_utilisation(self) -> float:
        return max(self.utilisation.values(), default=0.0)

    @property
    def worst_q_ratio(self) -> float:
        ratios = [
            abs(self.q_required_pu[b]) / self.q_capability_pu[b]
            for b in self.q_required_pu
            if self.q_capability_pu.get(b, 0.0) > 0.0
        ]
        return max(ratios, default=0.0)


def capability(policy: str, rating_pu: float, active_pu: float) -> float:
    """Reactive capability of one plant, per unit on the system base."""

    if policy == "vreg_pf095":
        return POWER_FACTOR_CAPABILITY * rating_pu
    if policy == "vreg_scap":
        return float(np.sqrt(max(rating_pu**2 - active_pu**2, 0.0)))
    raise ValueError(f"{policy} has no reactive capability rule")


def _with_pv(network: Ieee39Network, pv: dict[int, dict[str, float]]) -> Ieee39Network:
    return replace(network, pv=pv)


def solve_policy(
    members, policy: str, *, network: Ieee39Network | None = None
) -> PolicyOutcome:
    """Solve one subset under one reactive policy, at full replacement."""

    if policy not in POLICIES:
        raise ValueError(f"unknown policy {policy}")
    net = network or load_network()
    members = tuple(sorted(int(b) for b in members))
    plan_mapping = {bus: 1.0 for bus in members}

    if policy in ("matched", "unity_pf"):
        case = solve_case(ReplacementPlan.of(plan_mapping, q_policy=policy))
        flow = case.dae.power_flow
        required, caps, use = {}, {}, {}
        for bus in members:
            generation = flow.injection(bus, net.ybus) + net.loads.get(bus, 0j)
            rating = net.machines[bus]["Sn"] / SYSTEM_BASE_MVA
            required[bus] = float(generation.imag)
            caps[bus] = float(rating)
            use[bus] = float(abs(generation) / rating)
        return PolicyOutcome(
            case=case,
            policy=policy,
            members=members,
            q_required_pu=required,
            q_capability_pu=caps,
            utilisation=use,
            saturated=(),
            voltage_regulating=() if policy != "matched" else (),
            clamp_passes=0,
        )

    # Regulating policies. The replaced bus stays voltage-controlled at its own
    # pre-replacement setpoint until it asks for more reactive power than the
    # plant has, at which point it is switched onto the limit and re-solved.
    clamped: dict[int, float] = {}
    passes = 0
    for passes in range(1, MAX_CLAMP_PASSES + 1):
        pv = {bus: dict(spec) for bus, spec in net.pv.items()}
        for bus, value in clamped.items():
            pv[bus] = {**pv[bus], "q": value}
            pv[bus].pop("v", None)
        candidate = _with_pv(net, pv)
        flow = solve_power_flow(candidate)
        if not flow.converged:
            raise InfeasibleReplacement(
                f"power flow did not converge for {policy} on {members}"
            )
        changed = False
        for bus in members:
            generation = flow.injection(bus, net.ybus) + net.loads.get(bus, 0j)
            rating = net.machines[bus]["Sn"] / SYSTEM_BASE_MVA
            limit = capability(policy, rating, float(generation.real))
            if bus in clamped:
                continue
            if abs(generation.imag) > limit:
                clamped[bus] = float(np.sign(generation.imag) * limit)
                changed = True
        if not changed:
            break
    else:
        raise InfeasibleReplacement(
            f"reactive limit switching did not settle for {policy} on {members}"
        )

    pv = {bus: dict(spec) for bus, spec in net.pv.items()}
    for bus, value in clamped.items():
        pv[bus] = {**pv[bus], "q": value}
        pv[bus].pop("v", None)
    prepared = _with_pv(net, pv)

    regulating = tuple(bus for bus in members if bus not in clamped)
    converters = {
        bus: ConverterParameters(voltage_control=bus in regulating) for bus in members
    }
    case = solve_case(
        ReplacementPlan.of(plan_mapping, q_policy="matched"),
        network=prepared,
        converter=ConverterParameters(voltage_control=False),
        converters=converters,
    )

    flow = case.dae.power_flow
    required, caps, use = {}, {}, {}
    for bus in members:
        generation = flow.injection(bus, net.ybus) + net.loads.get(bus, 0j)
        rating = net.machines[bus]["Sn"] / SYSTEM_BASE_MVA
        required[bus] = float(generation.imag)
        caps[bus] = float(capability(policy, rating, float(generation.real)))
        use[bus] = float(abs(generation) / rating)
    return PolicyOutcome(
        case=case,
        policy=policy,
        members=members,
        q_required_pu=required,
        q_capability_pu=caps,
        utilisation=use,
        saturated=tuple(sorted(clamped)),
        voltage_regulating=regulating,
        clamp_passes=passes,
    )
