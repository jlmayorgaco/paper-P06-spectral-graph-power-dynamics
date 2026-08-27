from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
VALIDATOR = ROOT / "experiments" / "tx3" / "E05C_margin_conditioned" / "validate_e05c.py"


def test_e05c_preserves_negative_margin_result_and_e06_stop() -> None:
    specification = importlib.util.spec_from_file_location("validate_e05c", VALIDATOR)
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    assert module.validate(ROOT) == []
