from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
AUDIT = ROOT / "experiments/tx3/critical_mode_baseline_audit"
sys.path[:0] = [str(ROOT), str(AUDIT)]

from run_audit import state_group  # noqa: E402


def test_state_taxonomy_is_transparent_and_exhaustive_for_known_models() -> None:
    assert state_group("delta GENROU 30") == "synchronous_rotor"
    assert state_group("e2q GENROU 30") == "synchronous_electromagnetic"
    assert state_group("LG2_y GAST 32") == "turbine_governor"
    assert state_group("LAW_y SEXS 34") == "excitation"
    assert state_group("PI_xi PLL2 37") == "pll"
    assert state_group("S0_y REGCP1 38") == "converter_current_command"
    assert state_group("PIQ_xi REECB1 36") == "converter_electrical_control"
    assert state_group("s5_xi REPCA1 37") == "plant_control"
    assert state_group("x UNKNOWN 1") == "other"
