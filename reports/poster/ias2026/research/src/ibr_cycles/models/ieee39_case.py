"""Assembly of the IEEE-39 SG to PV-GFL replacement DAE.

A replacement plan assigns a fraction ``rho`` to each candidate generator bus.
The bus then carries ``(1-rho)`` of its synchronous machine and ``rho`` of a
battery-free PV grid-following converter, both on their own rating base and
scaled into the system base by their rating share.

Two reactive-power policies are supported and MUST be reported separately:

``matched``
    the converter supplies its rating share of the reactive power the machine
    was producing. The bus injection is unchanged, so the AC power flow solution
    is IDENTICAL to the base case at every ``rho``. This is the mechanism
    isolation campaign; the frozen operating point is exact, not approximate.

``unity_pf``
    the converter operates at unity power factor, as a battery-free PV plant
    normally would. For ``rho < 1`` the surviving machine must carry all of the
    reactive output on a reduced rating, which is often infeasible and is
    rejected explicitly. At ``rho = 1`` the bus becomes a PQ bus with zero
    reactive injection and the operating point genuinely moves. This is the
    engineering replacement campaign.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np
from numpy.typing import NDArray

from ..dynamics.dae import LinearSystem
from ..dynamics.equilibrium import Equilibrium, solve_equilibrium
from ..dynamics.linearize import central_difference_jacobians, reduce_index_one
from .ieee39_devices import (
    ConverterParameters,
    GridFollowingConverter,
    MachineParameters,
    StaticInjection,
    SynchronousMachine,
)
from .ieee39_network import Ieee39Network, PowerFlow, load_network, solve_power_flow

SYSTEM_BASE_MVA = 100.0
#: Largest apparent-power loading a device may carry, on its own rating.
LOADING_LIMIT = 1.0
#: Rating fraction below which a device is not instantiated at all.
RATING_FLOOR = 1e-6


@dataclass(frozen=True)
class ReplacementPlan:
    """Replacement fractions per candidate bus plus the reactive policy."""

    rho: tuple[tuple[int, float], ...]
    q_policy: str = "matched"
    device: str = "gfl"
    machine_damping: float | None = None
    condenser: tuple[tuple[int, float], ...] = ()
    condenser_services: tuple[tuple[str, float], ...] = ()

    @classmethod
    def of(
        cls,
        mapping: dict[int, float],
        q_policy: str = "matched",
        device: str = "gfl",
        machine_damping: float | None = None,
        condenser: dict[int, float] | None = None,
        condenser_services: dict[str, float] | None = None,
    ) -> ReplacementPlan:
        """Build a plan.

        ``device`` selects what replaces the machine: ``gfl`` for the dynamic
        converter, ``static`` for an ideal constant P/Q injection with no
        dynamics, which is the negative control for controller mediation.

        ``machine_damping`` overrides the damping coefficient of every surviving
        machine, on the machine base. The source case has ``D = 0`` throughout,
        so a multiplicative "increase damping" control would be vacuous; this is
        an absolute addition and is reported as such.
        """

        return cls(
            rho=tuple(
                sorted((int(b), float(r)) for b, r in mapping.items() if r > 0.0)
            ),
            q_policy=q_policy,
            device=device,
            machine_damping=machine_damping,
            condenser=tuple(
                sorted((int(b), float(v)) for b, v in (condenser or {}).items())
            ),
            condenser_services=tuple(
                sorted(
                    (str(k), float(v)) for k, v in (condenser_services or {}).items()
                )
            ),
        )

    @property
    def mapping(self) -> dict[int, float]:
        return {bus: value for bus, value in self.rho}

    @property
    def members(self) -> tuple[int, ...]:
        return tuple(bus for bus, _ in self.rho)

    @property
    def size(self) -> int:
        return len(self.rho)

    @property
    def label(self) -> str:
        if not self.rho:
            body = "BASE"
        else:
            body = "+".join(f"{bus}@{value:g}" for bus, value in self.rho)
        tags = [] if self.device == "gfl" else [self.device]
        if self.machine_damping is not None:
            tags.append(f"D={self.machine_damping:g}")
        if self.condenser:
            tags.append("SC=" + ",".join(f"{b}:{v:g}" for b, v in self.condenser))
        for name, value in self.condenser_services:
            tags.append(f"{name}={value:g}")
        return body + ("[" + ",".join(tags) + "]" if tags else "")


@dataclass(frozen=True)
class DeviceSlot:
    """One instantiated device with its state slice in the global vector."""

    device: object
    start: int
    stop: int
    kind: str
    bus: int
    weight: float
    loading: float


@dataclass
class Ieee39Dae:
    """Semi-explicit DAE of the replacement case, in rectangular bus voltages."""

    network: Ieee39Network
    power_flow: PowerFlow
    slots: tuple[DeviceSlot, ...]
    plan: ReplacementPlan
    n_x: int = 0
    n_z: int = 0
    labels: tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        self.n_x = sum(slot.stop - slot.start for slot in self.slots)
        self.n_z = 2 * self.network.n_bus
        self.labels = tuple(name for slot in self.slots for name in slot.device.labels)

    def voltages(self, z: NDArray[np.float64]) -> NDArray[np.complex128]:
        return z[0::2] + 1j * z[1::2]

    def f(
        self, x: NDArray[np.float64], z: NDArray[np.float64], theta: dict[str, float]
    ) -> NDArray[np.float64]:
        v = self.voltages(z)
        out = np.empty(self.n_x)
        for slot in self.slots:
            bus_position = self.network.position(slot.bus)
            out[slot.start : slot.stop] = slot.device.derivatives(
                x[slot.start : slot.stop], complex(v[bus_position])
            )
        return out

    def g(
        self, x: NDArray[np.float64], z: NDArray[np.float64], theta: dict[str, float]
    ) -> NDArray[np.float64]:
        v = self.voltages(z)
        injection = np.zeros(self.network.n_bus, dtype=np.complex128)
        for slot in self.slots:
            position = self.network.position(slot.bus)
            injection[position] += slot.device.injection(
                x[slot.start : slot.stop], complex(v[position])
            )
        if self.network.load_model == "impedance":
            # constant shunt impedance fixed by the power-flow voltage (68-bus)
            v0 = self.power_flow.voltages
            for bus, load in self.network.loads.items():
                position = self.network.position(bus)
                injection[position] -= (
                    np.conj(load) / abs(v0[position]) ** 2 * v[position]
                )
        else:
            for bus, load in self.network.loads.items():
                position = self.network.position(bus)
                injection[position] -= np.conj(load) / np.conj(v[position])
        residual = self.network.ybus @ v - injection
        out = np.empty(self.n_z)
        out[0::2] = residual.real
        out[1::2] = residual.imag
        return out


def _machine_parameters(
    network: Ieee39Network, bus: int, payload: dict, damping: float | None = None
) -> MachineParameters:
    raw = network.machines[bus]
    avr = payload["avr_by_bus"][bus]
    pss = payload["pss_by_bus"][bus]
    return MachineParameters(
        ra=raw["ra"],
        xd=raw["xd"],
        xq=raw["xq"],
        xd1=raw["xd1"],
        xq1=raw["xq1"],
        td10=raw["Td10"],
        tq10=raw["Tq10"],
        m=raw["M"],
        d=raw["D"] if damping is None else float(damping),
        ka=avr["KA"],
        ta=avr["TE"],
        avr_tatb=float(avr.get("TATB", 1.0)),
        avr_tb=float(avr.get("TB", 0.0)),
        pss_gain=pss["KS"],
        pss_washout=pss["T5"],
        pss_wash_lag=pss["T6"],
        pss_lag=pss["T4"],
    )


@lru_cache(maxsize=8)
def _controller_payload(path: str | None = None) -> dict:
    import json
    from pathlib import Path

    from .ieee39_network import CONFIG

    source = Path(path) if path else CONFIG
    data = json.loads(source.read_text(encoding="utf-8"))
    order = [int(row["bus"]) for row in data["machines"]]
    return {
        "avr_by_bus": {bus: row for bus, row in zip(order, data["avr"], strict=True)},
        "pss_by_bus": {bus: row for bus, row in zip(order, data["pss"], strict=True)},
    }


@lru_cache(maxsize=2)
def _ieee68_parameters(path: str) -> dict:
    import json
    from pathlib import Path

    from .ieee68_devices import Ieee68MachineParameters

    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return {
        int(bus): Ieee68MachineParameters.from_rows(m, a, p)
        for bus, m, a, p in zip(
            data["machine_order"],
            data["machines"],
            data["avr"],
            data["pss"],
            strict=True,
        )
    }


def _ieee68_machine(
    network: Ieee39Network, bus: int, weight: float, scaling: dict[str, float] | None
):
    """Documented 68-bus generator; ``scaling['ka']`` is the regulator-gain scale."""

    from dataclasses import replace as _replace

    from .ieee68_devices import Ieee68Machine

    parameters = _ieee68_parameters(network.config_path)[bus]
    unknown = set(scaling or {}) - {"ka"}
    if unknown:
        raise ValueError(
            f"68-bus machines accept only the 'ka' coordinate, not {sorted(unknown)}"
        )
    if scaling and "ka" in scaling:
        parameters = _replace(parameters, gain_scale=float(scaling["ka"]))
    return Ieee68Machine(bus=bus, parameters=parameters, weight=weight)


class InfeasibleReplacement(RuntimeError):
    """Raised when a plan asks a device to exceed its own rating."""


SCALABLE_MACHINE_FIELDS = (
    "m",
    "d",
    "xd1",
    "xq1",
    "ka",
    "ta",
    "pss_gain",
    "pss_washout",
    "pss_wash_lag",
    "pss_lag",
)


def _override_machine(
    parameters: MachineParameters, overrides: dict[str, float] | None
) -> MachineParameters:
    """Set named machine parameters to ABSOLUTE values on every machine.

    Scaling is the right tool for an uncertainty study, where each machine keeps
    its own value and is perturbed around it. A parameter STUDY needs the same
    value everywhere, which scaling cannot express when the machines differ: the
    excitation time constants in this case run from 0.25 to 0.50 s.
    """

    if not overrides:
        return parameters
    unknown = set(overrides) - set(SCALABLE_MACHINE_FIELDS)
    if unknown:
        raise ValueError(f"cannot override machine fields {sorted(unknown)}")
    from dataclasses import replace as _replace

    return _replace(parameters, **{k: float(v) for k, v in overrides.items()})


def _scale_machine(
    parameters: MachineParameters, scaling: dict[str, float] | None
) -> MachineParameters:
    """Multiply named machine parameters, for uncertainty studies.

    Omitting ``scaling`` reproduces the earlier behaviour exactly. Only the
    fields in ``SCALABLE_MACHINE_FIELDS`` may be scaled; anything else raises,
    so a typo in a sampling script cannot silently do nothing.
    """

    if not scaling:
        return parameters
    unknown = set(scaling) - set(SCALABLE_MACHINE_FIELDS)
    if unknown:
        raise ValueError(f"cannot scale machine fields {sorted(unknown)}")
    from dataclasses import replace as _replace

    return _replace(
        parameters,
        **{
            name: getattr(parameters, name) * float(factor)
            for name, factor in scaling.items()
        },
    )


def build_dae(
    plan: ReplacementPlan,
    *,
    network: Ieee39Network | None = None,
    converter: ConverterParameters | None = None,
    converters: dict[int, ConverterParameters] | None = None,
    machine_scaling: dict[str, float] | None = None,
    machine_services: dict[str, float] | None = None,
    machine_overrides: dict[str, float] | None = None,
    machine_bus_scaling: dict[int, dict[str, float]] | None = None,
) -> tuple[Ieee39Dae, NDArray[np.float64], NDArray[np.float64]]:
    """Instantiate the devices, solve the AC point and initialize every state.

    ``converters`` optionally overrides the converter parameters at individual
    buses. It is needed when one plant in a portfolio is at its reactive
    capability limit and must run on a fixed reactive command while the others
    still regulate voltage. Omitting it reproduces the earlier behaviour exactly.

    ``machine_bus_scaling`` multiplies named parameters of the machine at one bus,
    after the fleet-wide ``machine_scaling``. It is what a heterogeneity study
    needs: one scalar cannot move machines with different data by different
    factors. Omitting it reproduces the earlier behaviour exactly.
    """

    net = network or load_network()
    payload = _controller_payload(net.config_path or None)
    rho = plan.mapping
    base_converter = converter or ConverterParameters()

    pv_overrides: dict[int, dict[str, float]] = {}
    if plan.q_policy == "unity_pf":
        base = solve_power_flow(net)
        for bus, value in rho.items():
            if value >= 1.0 - RATING_FLOOR:
                gen = base.injection(bus, net.ybus) + net.loads.get(bus, 0j)
                pv_overrides[bus] = {"p": float(gen.real), "q": 0.0}
    power_flow = solve_power_flow(net, pv_overrides=pv_overrides or None)
    if not power_flow.converged:
        raise InfeasibleReplacement(f"power flow did not converge for {plan.label}")

    slots: list[DeviceSlot] = []
    states: list[NDArray[np.float64]] = []
    cursor = 0
    condensers = dict(plan.condenser)
    services = dict(plan.condenser_services)
    for bus in net.generator_buses:
        voltage = power_flow.at(bus)
        generation = power_flow.injection(bus, net.ybus) + net.loads.get(bus, 0j)
        rating = net.machines[bus]["Sn"] / SYSTEM_BASE_MVA
        fraction = float(rho.get(bus, 0.0))

        condenser_fraction = float(condensers.get(bus, 0.0))
        if condenser_fraction > RATING_FLOOR and fraction > RATING_FLOOR:
            # Synchronous condenser: the machine keeps its full dynamic apparatus
            # at the declared rating but produces no active power, and the
            # converter carries every displaced megawatt. This is the physically
            # interpretable version of "retain the synchronous service without
            # retaining the synchronous generation".
            # ``q_share`` (F8, default 1) splits the bus reactive output between
            # the condenser and the converter; the bus injection is unchanged.
            q_share = float(services.get("q_share", 1.0))
            machine_share = complex(0.0, q_share * generation.imag)
            converter_share = complex(
                generation.real, (1.0 - q_share) * generation.imag
            )
            entries = (
                ("sg", machine_share, condenser_fraction),
                ("gfl", converter_share, fraction),
            )
        else:
            if plan.q_policy == "matched":
                machine_share = complex((1.0 - fraction) * generation)
                converter_share = complex(fraction * generation)
            else:
                machine_share = complex(
                    (1.0 - fraction) * generation.real, generation.imag
                )
                converter_share = complex(fraction * generation.real, 0.0)
                if fraction >= 1.0 - RATING_FLOOR:
                    machine_share = 0j
                    converter_share = complex(generation)
            entries = (
                ("sg", machine_share, 1.0 - fraction),
                ("gfl", converter_share, fraction),
            )

        for kind, share, weight_fraction in entries:
            if weight_fraction <= RATING_FLOOR:
                continue
            weight = weight_fraction * rating
            if kind == "gfl" and net.converter_loading > 0.0:
                # 68-bus rule (G3 preregistration): the machine MVA column is a
                # per-unit base, not a rating, so the converter is rated at
                # |S_gen| / converter_loading.
                weight = weight_fraction * abs(generation) / net.converter_loading
            loading = float(abs(share) / weight)
            if loading > LOADING_LIMIT and net.enforce_ratings:
                raise InfeasibleReplacement(
                    f"{kind} at bus {bus} would run at {loading:.3f} of its rating "
                    f"under plan {plan.label}"
                )
            if kind == "sg" and net.machine_model == "ieee68_subtransient":
                if (
                    condenser_fraction > RATING_FLOOR
                    or machine_services
                    or machine_overrides
                ):
                    raise ValueError(
                        "services and condensers are not defined for the 68-bus model"
                    )
                device = _ieee68_machine(net, bus, weight, machine_scaling)
            elif kind == "sg":
                parameters = _scale_machine(
                    _override_machine(
                        _machine_parameters(net, bus, payload, plan.machine_damping),
                        machine_overrides,
                    ),
                    machine_scaling,
                )
                if machine_bus_scaling and bus in machine_bus_scaling:
                    parameters = _scale_machine(parameters, machine_bus_scaling[bus])
                if machine_services:
                    # Physically meaningful service ablation on EVERY surviving
                    # machine, the same switch the condenser study already uses.
                    # Excitation is frozen with avr_manual rather than by setting
                    # the gain to zero, which would let the field decay instead of
                    # holding it.
                    parameters = parameters.with_services(
                        inertia_scale=machine_services.get("inertia", 1.0),
                        pss_scale=machine_services.get("pss", 1.0),
                        avr_scale=machine_services.get("avr", 1.0),
                        avr_manual=bool(machine_services.get("avr_manual", 0.0)),
                        avr_blend=machine_services.get("avr_blend"),
                        flux_blend=machine_services.get("flux_blend"),
                        damping=machine_services.get("damping"),
                    )
                if condenser_fraction > RATING_FLOOR and fraction > RATING_FLOOR:
                    parameters = parameters.with_services(
                        inertia_scale=services.get("inertia", 1.0),
                        pss_scale=services.get("pss", 1.0),
                        avr_scale=services.get("avr", 1.0),
                        avr_manual=bool(services.get("avr_manual", 0.0)),
                        avr_blend=services.get("avr_blend"),
                        flux_blend=services.get("flux_blend"),
                        damping=services.get("damping"),
                    )
                device = SynchronousMachine(
                    bus=bus, parameters=parameters, weight=weight
                )
            elif plan.device.startswith("static"):
                characteristic = (
                    plan.device.split("_", 1)[1] if "_" in plan.device else "power"
                )
                device = StaticInjection(
                    bus=bus,
                    injection_pu=0j,
                    weight=weight,
                    characteristic=characteristic,
                )
            else:
                device = GridFollowingConverter(
                    bus=bus,
                    parameters=(converters or {}).get(bus, base_converter),
                    weight=weight,
                )
            device, state = device.initialize(voltage, share)
            slots.append(
                DeviceSlot(
                    device=device,
                    start=cursor,
                    stop=cursor + device.n_states,
                    kind=kind,
                    bus=bus,
                    weight=weight,
                    loading=loading,
                )
            )
            states.append(state)
            cursor += device.n_states

    dae = Ieee39Dae(network=net, power_flow=power_flow, slots=tuple(slots), plan=plan)
    x0 = np.concatenate(states)
    z0 = np.empty(dae.n_z)
    z0[0::2] = power_flow.voltages.real
    z0[1::2] = power_flow.voltages.imag
    return dae, x0, z0


@dataclass(frozen=True)
class ReplacementCase:
    """A solved replacement case: equilibrium, reduced matrix and evidence."""

    plan: ReplacementPlan
    dae: Ieee39Dae
    equilibrium: Equilibrium
    system: LinearSystem
    gz_condition: float

    @property
    def n_states(self) -> int:
        return self.system.n

    @property
    def max_loading(self) -> float:
        return max((slot.loading for slot in self.dae.slots), default=0.0)

    @property
    def replaced_mw(self) -> float:
        return sum(
            slot.weight * SYSTEM_BASE_MVA
            for slot in self.dae.slots
            if slot.kind == "gfl"
        )


def solve_case(
    plan: ReplacementPlan,
    *,
    network: Ieee39Network | None = None,
    converter: ConverterParameters | None = None,
    converters: dict[int, ConverterParameters] | None = None,
    machine_scaling: dict[str, float] | None = None,
    machine_services: dict[str, float] | None = None,
    machine_overrides: dict[str, float] | None = None,
    machine_bus_scaling: dict[int, dict[str, float]] | None = None,
    tol: float = 1e-9,
) -> ReplacementCase:
    """Build, equilibrate and linearize one replacement plan."""

    dae, x0, z0 = build_dae(
        plan,
        network=network,
        converter=converter,
        converters=converters,
        machine_scaling=machine_scaling,
        machine_services=machine_services,
        machine_overrides=machine_overrides,
        machine_bus_scaling=machine_bus_scaling,
    )
    equilibrium = solve_equilibrium(dae, {}, x0, z0, tol=tol)
    if not equilibrium.ok:
        raise InfeasibleReplacement(
            f"equilibrium failed for {plan.label}: {equilibrium.status}"
        )
    jacobians = central_difference_jacobians(dae, equilibrium.x, equilibrium.z, {})
    system = reduce_index_one(jacobians, dae.labels)
    return ReplacementCase(
        plan=plan,
        dae=dae,
        equilibrium=equilibrium,
        system=system,
        gz_condition=jacobians.gz_condition,
    )
