"""Make the package importable when experiments are run as plain scripts."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"

if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
