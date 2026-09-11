"""Unit data contract: MVA ratings and MW dispatch are separate, measured quantities."""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case
from ibr_cycles.models.ieee39_network import load_network
from ibr_cycles.units import measure, per_machine_table

RESEARCH = Path(__file__).resolve().parents[1]
CORE = (30, 33, 35, 37)


@pytest.fixture(scope="module")
def flagship():
    return solve_case(ReplacementPlan.of({b: 1.0 for b in CORE}))


def test_flagship_rating_is_mva_and_dispatch_is_measured(flagship):
    q = measure(flagship)
    assert q.replaced_sn_mva == pytest.approx(4270.7, abs=1e-9)
    assert q.replaced_sn_mva == pytest.approx(
        flagship.replaced_mw, abs=1e-12
    )  # legacy = MVA
    assert q.retired_sg_sn_mva == pytest.approx(4270.7, abs=1e-9)
    assert q.replaced_pg_mw == pytest.approx(2096.6077, abs=1e-3)
    assert q.replaced_pg_mw == pytest.approx(q.pg_mw_base_dispatch, abs=1e-6)
    assert q.replaced_q_mvar == pytest.approx(420.649, abs=1e-2)
    assert q.replaced_pmax_mw == pytest.approx(1040 + 682 + 697 + 564, abs=1e-9)
    assert q.sg_sn_mva_left_at_replaced_buses == 0.0 and q.condenser_sn_mva == 0.0


def test_partial_replacement_carries_rho_times_dispatch():
    q = measure(solve_case(ReplacementPlan.of({30: 0.5})))
    assert q.replaced_sn_mva == pytest.approx(520.0)
    assert q.replaced_pg_mw == pytest.approx(0.5 * 436.086385, abs=1e-5)
    assert q.sg_sn_mva_left_at_replaced_buses == pytest.approx(520.0)


def test_condenser_rating_is_a_separate_axis():
    case = solve_case(ReplacementPlan.of({30: 1.0}, condenser={30: 0.25}))
    q = measure(case)
    assert q.condenser_sn_mva == pytest.approx(260.0)
    assert q.replaced_pg_mw == pytest.approx(
        436.086385, abs=1e-5
    )  # converter carries all P
    assert abs(q.replaced_q_mvar) < 1e-6  # q_share = 1: the condenser carries Q


def test_pmax_only_where_documented():
    kundur = load_network(RESEARCH / "configs" / "kundur" / "kundur_network.json")
    assert all(math.isnan(r["Pmax_MW"]) for r in per_machine_table(kundur))
    ieee68 = load_network(RESEARCH / "configs" / "ieee68" / "ieee68_network.json")
    assert all(math.isnan(r["Sn_MVA"]) for r in per_machine_table(ieee68))
    rows = {r["bus"]: r for r in per_machine_table(load_network())}
    assert rows[30]["Pmax_MW"] == pytest.approx(1040.0)
    assert rows[37]["Sn_MVA"] == pytest.approx(970.2)
    assert rows[37]["Pg_MW"] == pytest.approx(321.521338, abs=1e-5)
