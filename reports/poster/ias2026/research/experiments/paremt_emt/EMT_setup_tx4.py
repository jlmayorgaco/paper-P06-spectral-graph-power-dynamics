# ruff: noqa: E501  -- parameter tables and kernel argument lists kept on one line
"""Prepare external/ParaEMT_tx4: apply the documented network patch and copy the TX4 overlay files.

Idempotent. The overlay sources are tracked in experiments/paremt_emt/overlay/; the working copy
gets exact copies (checked byte for byte). Run with any Python.
"""

from __future__ import annotations

import filecmp
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[5]
TX4 = REPO / "external" / "ParaEMT_tx4"
OVERLAY = HERE / "overlay"


def main() -> int:
    subprocess.run([sys.executable, str(HERE / "tx4_patch_paremt.py")], check=True)
    for src in sorted(OVERLAY.glob("*.py")):
        dst = TX4 / src.name
        shutil.copyfile(src, dst)
        assert filecmp.cmp(src, dst, shallow=False)
        print("overlay ->", dst.name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
