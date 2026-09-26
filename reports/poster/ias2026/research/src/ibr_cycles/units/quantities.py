"""Unit data contract for replacement portfolios (Phase-I unit correction).

The legacy property ``ReplacementCase.replaced_mw`` returns

    sum over converter slots of weight * 100 MVA,

that is, the apparent-power RATING of the replacement converters in MVA. It is
not active power. This module defines the canonical quantities that every new
result must use instead, each with its own unit and provenance:

    replaced_sn_mva      converter apparent-power rating as built [MVA]
                         (numerically equal to the legacy replaced_mw)
    retired_sg_sn_mva    sum of rho_b * Sn_b of the retired machine share [MVA]
                         (NaN where Sn is a per-unit base, not a rating: IEEE-68)
    replaced_pg_mw       active power actually carried by the replacement
                         converters at the solved equilibrium [MW], MEASURED
                         from the device injections, never derived from MVA
    replaced_q_mvar      reactive power carried by the converters [Mvar]
    replaced_pmax_mw     sum of rho_b * Pmax_b of the retired machines [MW], only
                         where the source case documents Pmax (PMAX_DOCUMENTED);
                         NaN otherwise. It is the SYNCHRONOUS machine's active
                         limit, not a PV nameplate: PV nameplate is not modelled
    sg_sn_mva_left_at_replaced_buses
                         Sn of the machine share left at candidate buses plus
                         condenser ratings [MVA]
    condenser_sn_mva     condenser rating added at candidate buses [MVA]

No function in this module converts MVA to MW.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

SYSTEM_BASE_MVA = 100.0

#: Pmax is used only where the source case documents it as a machine limit.
#: IEEE-39: ANDES PV/Slack.pmax of data/raw/ieee39_full.xlsx (sha 9c2048dc...).
#: Kundur: the ANDES case carries placeholder limits (qmax = 99 pu, pmax 9/9/20/99
#: pu), which are not capabilities. IEEE-68: no Pmax in the source.
PMAX_DOCUMENTED = {
    "ieee39_network.json": True,
    "kundur_network.json": False,
    "ieee68_network.json": False,
}


@dataclass(frozen=True)
class ReplacementQuantities:
    replaced_sn_mva: float
    retired_sg_sn_mva: float
    replaced_pg_mw: float
    replaced_q_mvar: float
    replaced_pmax_mw: float
    sg_sn_mva_left_at_replaced_buses: float
    condenser_sn_mva: float
    pmax_documented: bool
    pg_mw_base_dispatch: float  # sum of rho_b * Pg_b of the base power flow (check)

    def as_dict(self) -> dict:
        return asdict(self)


def _documented(network) -> bool:
    return bool(PMAX_DOCUMENTED.get(Path(network.config_path).name, False))


def _pmax_pu(network, bus: int) -> float:
    if bus == network.slack_bus:
        return float(network.slack_pmax)
    return float(network.pv[bus]["pmax"]) if bus in network.pv else math.nan


def measure(case) -> ReplacementQuantities:
    """Canonical quantities of a SOLVED ReplacementCase, measured at equilibrium."""

    dae = case.dae
    net = dae.network
    x, z = case.equilibrium.x, case.equilibrium.z
    v = dae.voltages(z)
    rho = dict(case.plan.rho)
    condensers = dict(case.plan.condenser)
    sn_rating = pg = q = 0.0
    for slot in dae.slots:
        if slot.kind != "gfl":
            continue
        pos = net.position(slot.bus)
        s = complex(v[pos]) * np.conj(
            slot.device.injection(x[slot.start : slot.stop], complex(v[pos]))
        )
        pg += s.real * SYSTEM_BASE_MVA
        q += s.imag * SYSTEM_BASE_MVA
        sn_rating += slot.weight * SYSTEM_BASE_MVA
    is_rating = net.converter_loading <= 0.0  # IEEE-68: Sn is a per-unit base
    documented = _documented(net)
    retired_sn = pmax = base_pg = retained_sn = cond_sn = 0.0
    pf = dae.power_flow
    for bus, fraction in rho.items():
        machine_sn = float(net.machines[bus]["Sn"])
        retired_sn += fraction * machine_sn
        pmax += fraction * _pmax_pu(net, bus) * SYSTEM_BASE_MVA
        gen = pf.injection(bus, net.ybus) + net.loads.get(bus, 0j)
        base_pg += fraction * gen.real * SYSTEM_BASE_MVA
        c = float(condensers.get(bus, 0.0))
        if c > 0.0:
            cond_sn += c * machine_sn
            retained_sn += c * machine_sn
        else:
            retained_sn += (1.0 - fraction) * machine_sn
    return ReplacementQuantities(
        replaced_sn_mva=sn_rating,
        retired_sg_sn_mva=retired_sn if is_rating else math.nan,
        replaced_pg_mw=pg,
        replaced_q_mvar=q,
        replaced_pmax_mw=pmax if documented else math.nan,
        sg_sn_mva_left_at_replaced_buses=retained_sn if is_rating else math.nan,
        condenser_sn_mva=cond_sn,
        pmax_documented=documented,
        pg_mw_base_dispatch=base_pg,
    )


def per_machine_table(network) -> list[dict]:
    """Pg, Qg (base power flow), Sn and documented Pmax of every machine."""

    from ..models.ieee39_network import solve_power_flow

    pf = solve_power_flow(network)
    documented = _documented(network)
    rows = []
    for bus in network.generator_buses:
        gen = pf.injection(bus, network.ybus) + network.loads.get(bus, 0j)
        rows.append(
            {
                "bus": bus,
                "slack": bus == network.slack_bus,
                "Pg_MW": gen.real * SYSTEM_BASE_MVA,
                "Qg_Mvar": gen.imag * SYSTEM_BASE_MVA,
                "Sn_MVA": float(network.machines[bus]["Sn"])
                if network.converter_loading <= 0.0
                else math.nan,
                "Pmax_MW": _pmax_pu(network, bus) * SYSTEM_BASE_MVA
                if documented
                else math.nan,
                "Pmax_source": "ANDES PV/Slack.pmax (source case)"
                if documented
                else "not documented",
            }
        )
    return rows


def source_note(network) -> str:
    payload = json.loads(Path(network.config_path).read_text(encoding="utf-8"))
    return json.dumps(payload.get("source", {}))
