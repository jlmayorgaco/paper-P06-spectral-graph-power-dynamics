from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
VALIDATOR = ROOT / "experiments" / "tx3" / "E05B_dynamic_spectral_shift" / "validate_e05b.py"


def test_e05b_artifacts_preserve_frozen_decisions_and_stop_boundary() -> None:
    specification = importlib.util.spec_from_file_location("validate_e05b", VALIDATOR)
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    assert module.validate(ROOT) == []
