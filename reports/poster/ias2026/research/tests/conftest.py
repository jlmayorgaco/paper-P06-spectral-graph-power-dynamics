from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ibr_cycles.models.toy5_case import ToyCase, build_case  # noqa: E402


@pytest.fixture(scope="session")
def toy() -> ToyCase:
    return build_case()
