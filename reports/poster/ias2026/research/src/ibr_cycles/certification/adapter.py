"""BC02 research adapter: explicit contracts over the three benchmarks.

Each method either returns the object it names or an explicit status string
(NOT_APPLICABLE, UNKNOWN, INFEASIBLE). It never substitutes an approximation
for missing structure, and never pads dimensions with dead states or
integrators without the ledger that proves vertex equivalence.

    equilibrate(portfolio, policy, uncertainty)  -> ReplacementCase | "INFEASIBLE"
    linearize_full(case)       -> {"fx","fz","gx","gz","labels","bases"}
    symmetry_generators(case)  -> (R_x, R_z, FrequencyPartner), from the equations
    quotient_model(case)       -> PhysicalReport (spectrum preserved, zero ledger)
    port_model(case)           -> (T(s), h(s)) with det P = h det T
    affine_binary_model(policy)-> CommonRealization | "NOT_APPLICABLE"
    operator_error_bound(...)  -> "UNKNOWN" (no validated bound is implemented)

Policies: {"g", "k", "t", "h"}, the leaky Q/V gain and the excitation
coordinates of F7 (IEEE-39), F12 (Kundur) and G3 (IEEE-68). Only matched
dispatch is supported; the uncertainty argument must be None.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from ..dynamics.linearize import central_difference_jacobians
from ..models.ieee39_case import InfeasibleReplacement, ReplacementPlan, solve_case
from ..models.ieee39_devices import ConverterParameters
from ..models.ieee39_network import load_network
from .binary import build_common_realization
from .physical import physical_report
from .symmetry import frequency_partner, rotation_generator

NOT_APPLICABLE = "NOT_APPLICABLE"
UNKNOWN = "UNKNOWN"
INFEASIBLE = "INFEASIBLE"
LEAK = 0.05
_RESEARCH = Path(__file__).resolve().parents[3]


class BenchmarkAdapter:
    name = ""
    network_path: Path | None = None
    candidates: tuple[int, ...] = ()
    scaling_keys: tuple[str, ...] = ("ka", "ta")

    def _network(self):
        return load_network(self.network_path) if self.network_path else None

    def _kwargs(self, policy: dict) -> dict:
        conv = ConverterParameters(
            voltage_control=True,
            voltage_gain=float(policy.get("g", 0.0)),
            voltage_leak=LEAK,
        )
        scaling = {"ka": float(policy.get("k", 1.0))}
        if "ta" in self.scaling_keys:
            scaling["ta"] = float(policy.get("t", 1.0))
        kwargs = {"converter": conv, "machine_scaling": scaling}
        net = self._network()
        if net is not None:
            kwargs["network"] = net
        if float(policy.get("h", 1.0)) != 1.0:
            kwargs["machine_bus_scaling"] = self._heterogeneity(policy)
        return kwargs

    def _heterogeneity(self, policy):  # pragma: no cover - IEEE-39 only
        raise NotImplementedError(NOT_APPLICABLE)

    # -- contracts ---------------------------------------------------------
    def equilibrate(self, portfolio, policy: dict, uncertainty=None):
        if uncertainty is not None:
            return NOT_APPLICABLE
        try:
            return solve_case(
                ReplacementPlan.of({b: 1.0 for b in portfolio}), **self._kwargs(policy)
            )
        except (InfeasibleReplacement, ValueError):
            return INFEASIBLE

    def linearize_full(self, case) -> dict:
        jac = central_difference_jacobians(
            case.dae, case.equilibrium.x, case.equilibrium.z, {}
        )
        return {
            "fx": jac.fx,
            "fz": jac.fz,
            "gx": jac.gx,
            "gz": jac.gz,
            "labels": case.dae.labels,
            "bases": {"system_MVA": 100.0, "devices": "own rating, scaled by weight"},
        }

    def symmetry_generators(self, case):
        r_x, r_z = rotation_generator(case.dae, case.equilibrium.z)
        return r_x, r_z, frequency_partner(case.dae)

    def quotient_model(self, case):
        return physical_report(case)

    def port_model(self, case):
        lin = self.linearize_full(case)
        fx, fz, gx, gz = lin["fx"], lin["fz"], lin["gx"], lin["gz"]
        n = fx.shape[0]

        def t_of(s: complex):
            return gz + gx @ np.linalg.solve(s * np.eye(n) - fx, fz)

        def h_of(s: complex):
            return complex(np.linalg.det(s * np.eye(n) - fx))

        return t_of, h_of

    def affine_binary_model(self, policy: dict, candidates=None):
        try:
            return build_common_realization(
                tuple(candidates or self.candidates), **self._kwargs(policy)
            )
        except (InfeasibleReplacement, ValueError):
            return NOT_APPLICABLE

    def operator_error_bound(self, cell, *args, **kwargs):
        return UNKNOWN


class Ieee39Adapter(BenchmarkAdapter):
    name = "IEEE-39"
    candidates = (30, 33, 35, 37)

    def _heterogeneity(self, policy):
        from ..models.ieee39_network import load_network as _ln

        net = _ln()
        import json

        payload = json.loads(Path(net.config_path).read_text(encoding="utf-8"))
        te = {
            int(m["bus"]): a["TE"]
            for m, a in zip(payload["machines"], payload["avr"], strict=True)
        }
        mean = math.exp(sum(math.log(v) for v in te.values()) / len(te))
        h = float(policy["h"])
        return {b: {"ta": (mean / v) ** (1.0 - h)} for b, v in te.items()}


class KundurAdapter(BenchmarkAdapter):
    name = "Kundur"
    network_path = _RESEARCH / "configs" / "kundur" / "kundur_network.json"
    candidates = (2, 3, 4)


class Ieee68Adapter(BenchmarkAdapter):
    name = "IEEE-68"
    network_path = _RESEARCH / "configs" / "ieee68" / "ieee68_network.json"
    candidates = (3, 4, 6, 9)
    scaling_keys = ("ka",)
