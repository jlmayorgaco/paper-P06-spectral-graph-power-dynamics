"""IEEE/NE 39-bus network: frozen data, admittance matrix and AC power flow.

The network is imported once from the repository case ``data/raw/ieee39_full.xlsx``
and frozen into ``configs/ias2026/ieee39_network.json`` with the source hash, so
that no experiment depends on the ANDES environment being installed. ANDES is
used only as an independent cross-check of the admittance matrix and of the
power-flow solution; those references are stored beside the frozen data.

Base convention follows the source case: machine parameters are given on the
machine base ``Sn`` and converted to the 100 MVA system base as
``x_sys = x * 100/Sn`` and ``M_sys = M * Sn/100``.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import root

SYSTEM_BASE_MVA = 100.0
CONFIG = (
    Path(__file__).resolve().parents[3] / "configs" / "ias2026" / "ieee39_network.json"
)


@dataclass(frozen=True)
class Ieee39Network:
    """Frozen network topology, injections and the derived admittance matrix."""

    bus_idx: tuple[int, ...]
    ybus: NDArray[np.complex128]
    loads: dict[int, complex]
    pv: dict[int, dict[str, float]]
    slack_bus: int
    slack_voltage: float
    slack_angle: float
    slack_pmax: float
    machines: dict[int, dict[str, float]]
    source_sha256: str
    #: the frozen JSON this network was read from; selects the controller data
    config_path: str = ""
    #: Model switches of the frozen JSON's ``model`` block (IEEE 68-bus, Gate 3).
    #: The defaults are the IEEE-39 / Kundur conventions and change nothing there.
    machine_model: str = "two_axis"
    load_model: str = "power"
    enforce_ratings: bool = True
    #: converter rating = |S_gen| / converter_loading when set (68-bus rule)
    converter_loading: float = 0.0

    @property
    def n_bus(self) -> int:
        return len(self.bus_idx)

    def position(self, bus: int) -> int:
        return self.bus_idx.index(bus)

    @property
    def generator_buses(self) -> tuple[int, ...]:
        """Every bus carrying a synchronous machine, slack included."""

        return tuple(sorted(self.machines))

    @property
    def replacement_candidates(self) -> tuple[int, ...]:
        """Machines eligible for SG to PV-GFL replacement.

        The slack machine is excluded: in this case it carries ``M = 100`` on its
        own base, i.e. it represents the rest of the interconnection rather than
        a replaceable plant. Excluding it is a declared modelling decision, not a
        result-dependent choice.
        """

        return tuple(b for b in self.generator_buses if b != self.slack_bus)

    def machine_on_system_base(self, bus: int) -> dict[str, float]:
        """Convert one machine from its own base to the 100 MVA system base."""

        raw = self.machines[bus]
        ratio = SYSTEM_BASE_MVA / raw["Sn"]
        converted = {"Sn": raw["Sn"]}
        for name in ("ra", "xl", "xd", "xq", "xd1", "xq1", "xd2", "xq2"):
            converted[name] = raw[name] * ratio
        for name in ("M", "D"):
            converted[name] = raw[name] / ratio
        for name in ("Td10", "Td20", "Tq10", "Tq20"):
            converted[name] = raw[name]
        return converted


def _build_ybus(payload: dict) -> NDArray[np.complex128]:
    order = [int(b["idx"]) for b in payload["buses"]]
    index = {bus: position for position, bus in enumerate(order)}
    n = len(order)
    y = np.zeros((n, n), dtype=np.complex128)
    for line in payload["lines"]:
        if float(line["u"]) == 0.0:
            continue
        f, t = index[int(line["bus1"])], index[int(line["bus2"])]
        series = 1.0 / complex(line["r"], line["x"])
        charging = complex(line["g"], line["b"]) / 2.0
        m = float(line["tap"]) * np.exp(1j * float(line["phi"]))
        m2 = float(abs(m) ** 2)
        y[f, f] += (series + charging) / m2
        y[t, t] += series + charging
        y[f, t] += -series / np.conj(m)
        y[t, f] += -series / m
    for shunt in payload["shunts"]:
        position = index[int(shunt["bus"])]
        y[position, position] += complex(shunt["g"], shunt["b"])
    return y


@lru_cache(maxsize=1)
def load_network(path: Path = CONFIG) -> Ieee39Network:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    order = tuple(int(b["idx"]) for b in payload["buses"])
    loads: dict[int, complex] = {}
    for row in payload["loads"]:
        bus = int(row["bus"])
        loads[bus] = loads.get(bus, 0j) + complex(row["p0"], row["q0"])
    pv = {
        int(row["bus"]): {
            "p": float(row["p0"]),
            "v": float(row["v0"]),
            "sn": float(row["Sn"]),
            "qmax": float(row["qmax"]),
            "qmin": float(row["qmin"]),
            "pmax": float(row["pmax"]),
            "pmin": float(row["pmin"]),
        }
        for row in payload["pv"]
    }
    slack = payload["slack"][0]
    machines = {int(row["bus"]): dict(row) for row in payload["machines"]}
    model = payload.get("model", {})
    return Ieee39Network(
        bus_idx=order,
        ybus=_build_ybus(payload),
        loads=loads,
        pv=pv,
        slack_bus=int(slack["bus"]),
        slack_voltage=float(slack["v0"]),
        slack_angle=float(slack["a0"]),
        slack_pmax=float(slack["pmax"]),
        machines=machines,
        source_sha256=payload["source"]["sha256"],
        config_path=str(Path(path)),
        machine_model=str(model.get("machine_model", "two_axis")),
        load_model=str(model.get("load_model", "power")),
        enforce_ratings=bool(model.get("enforce_ratings", True)),
        converter_loading=float(model.get("converter_loading", 0.0)),
    )


@dataclass(frozen=True)
class PowerFlow:
    """A converged AC operating point with its acceptance evidence."""

    voltages: NDArray[np.complex128]
    bus_idx: tuple[int, ...]
    max_mismatch: float
    converged: bool
    iterations: int

    def at(self, bus: int) -> complex:
        return complex(self.voltages[self.bus_idx.index(bus)])

    @property
    def min_voltage(self) -> float:
        return float(np.abs(self.voltages).min())

    @property
    def max_voltage(self) -> float:
        return float(np.abs(self.voltages).max())

    def injection(self, bus: int, ybus: NDArray[np.complex128]) -> complex:
        position = self.bus_idx.index(bus)
        current = ybus[position, :] @ self.voltages
        return complex(self.voltages[position] * np.conj(current))


def solve_power_flow(
    network: Ieee39Network,
    *,
    pv_overrides: dict[int, dict[str, float]] | None = None,
    pq_overrides: dict[int, complex] | None = None,
    load_scale: complex | float = 1.0,
    tol: float = 1e-12,
) -> PowerFlow:
    """Newton AC power flow.

    ``pv_overrides`` retypes or rescales generator buses; a bus mapped to
    ``{"p": P, "q": Q}`` without ``"v"`` becomes a PQ injection, which is how a
    unity-power-factor PV plant replacing a voltage-regulating machine is
    represented. ``pq_overrides`` adds a fixed complex injection at a bus.
    """

    order = network.bus_idx
    index = {bus: position for position, bus in enumerate(order)}
    n = len(order)
    pv = {bus: dict(spec) for bus, spec in network.pv.items()}
    for bus, spec in (pv_overrides or {}).items():
        if spec is None:
            pv.pop(bus, None)
            continue
        merged = {**pv.get(bus, {}), **spec}
        if "q" in spec and "v" not in spec:
            # Retyping a voltage-controlled bus to PQ. Merging alone would leave
            # the old voltage setpoint in place and the bus would silently stay
            # PV, so the requested policy would never take effect.
            merged.pop("v", None)
        pv[bus] = merged

    scheduled = np.zeros(n, dtype=np.complex128)
    for bus, load in network.loads.items():
        scheduled[index[bus]] -= load * load_scale
    for bus, extra in (pq_overrides or {}).items():
        scheduled[index[bus]] += extra
    voltage_controlled: dict[int, float] = {}
    for bus, spec in pv.items():
        scheduled[index[bus]] += spec["p"]
        if "v" in spec:
            voltage_controlled[bus] = float(spec["v"])
        else:
            scheduled[index[bus]] += 1j * float(spec.get("q", 0.0))

    slack_position = index[network.slack_bus]
    pq_positions = [
        position
        for position, bus in enumerate(order)
        if bus != network.slack_bus and bus not in voltage_controlled
    ]
    angle_positions = [p for p in range(n) if p != slack_position]

    magnitudes0 = np.ones(n)
    magnitudes0[slack_position] = network.slack_voltage
    for bus, value in voltage_controlled.items():
        magnitudes0[index[bus]] = value

    def unpack(w: NDArray[np.float64]) -> NDArray[np.complex128]:
        angles = np.zeros(n)
        angles[slack_position] = network.slack_angle
        magnitudes = magnitudes0.copy()
        angles[angle_positions] = w[: len(angle_positions)]
        magnitudes[pq_positions] = w[len(angle_positions) :]
        return magnitudes * np.exp(1j * angles)

    def mismatch(w: NDArray[np.float64]) -> NDArray[np.float64]:
        v = unpack(w)
        s = v * np.conj(network.ybus @ v)
        delta = s - scheduled
        return np.concatenate(
            [np.real(delta)[angle_positions], np.imag(delta)[pq_positions]]
        )

    guess = np.concatenate([np.zeros(len(angle_positions)), magnitudes0[pq_positions]])
    solution = root(mismatch, guess, method="hybr", tol=tol)
    residual = float(np.abs(mismatch(solution.x)).max())
    return PowerFlow(
        voltages=unpack(solution.x),
        bus_idx=order,
        max_mismatch=residual,
        converged=bool(solution.success and residual < 1e-9),
        iterations=int(solution.get("nfev", -1)),
    )
