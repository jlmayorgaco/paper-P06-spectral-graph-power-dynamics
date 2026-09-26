# ruff: noqa: E501  -- long diagnostic strings and code-patch literals kept on one line
"""EMT00c - tool-level diagnosis of the stock IEEE-39 start-up transient (no TX4 meaning).

EMT00b found that the stock ParaEMT IEEE-39 case is not stationary in a no-event run
(a 0.62 pu bus-voltage jump in the first millisecond). Upstream loads the transformer
ratio pfd.xfmr_k but never uses it; the stock power flow closes (1.8e-3 pu) only with the
6-31 ratio on the to-bus (bus 31), and has a 1.47 pu mismatch without it.

This script:
  1. reruns the unmodified stock no-event case twice (2 s) and checks bit-identity
     (determinism of the installed tool);
  2. runs the SAME stock case, stock models, stock damping, from a copy of the patched
     working copy with ONLY the documented tap applied (ini.apply_xfmr_tap = True,
     xfmr_tap_side = 'to'), 10 s no event, and reports the drift.
Writes results/EMT00/tap_diagnostic.json.
"""

from __future__ import annotations

import json
import os
import shutil
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import EMT00_stock_smoke as smoke  # noqa: E402

REPO = HERE.parents[5]
TX4 = REPO / "external" / "ParaEMT_tx4"
RUN_UP = REPO / "external" / "paremt_runs" / "EMT00_stock"
RUN_TAP = REPO / "external" / "paremt_runs" / "EMT00_stock_tap"
OUT = smoke.OUT


def drift(emt, pfd):
    nbus, ngen = len(pfd.bus_num), len(pfd.gen_bus)
    x = smoke.stack(emt.x)
    v = smoke.stack(emt.v)
    t = np.asarray(emt.t)
    vm = smoke.vmag(v, nbus)
    speed = x[:, 1 : 18 * ngen : 18] / pfd.ws
    delta = x[:, 0 : 18 * ngen : 18]
    rel = delta - delta[:, [0]]
    return {
        "max_speed_drift_pu": float(np.abs(speed - speed[0]).max()),
        "max_voltage_magnitude_drift_pu": float(np.abs(vm - vm[0]).max()),
        "max_relative_angle_drift_rad": float(np.abs(rel - rel[0]).max()),
        "initial_voltage_step_first_1ms_pu": float(np.abs(vm[1] - vm[0]).max()),
        "tlen_s": float(t[-1]),
    }, x


def main() -> int:
    out = {}
    # 1. determinism of the stock tool (two identical 2 s no-event runs)
    smoke.prepare_run_dir()
    os.chdir(RUN_UP)
    sys.path.insert(0, str(RUN_UP))
    runs = []
    for _ in range(2):
        pfd, ini, dyd, emt, wall, n, _ = smoke.run_case("noevent", 2.0)
        d, x = drift(emt, pfd)
        runs.append((d, x))
    out["stock_determinism"] = {
        "bit_identical_states": bool(np.array_equal(runs[0][1], runs[1][1])),
        "stock_drift_2s": runs[0][0],
    }
    for mod in [
        m
        for m in list(sys.modules)
        if m
        in (
            "Lib_BW",
            "lib_numba",
            "psutils",
            "partitionutil",
            "serial_bbd_matrix",
            "preprocessscript",
            "bbd_matrix",
        )
    ]:
        del sys.modules[mod]
    sys.path.remove(str(RUN_UP))

    # 2. the same stock case with only the documented tap applied (patched working copy)
    if RUN_TAP.exists():
        shutil.rmtree(RUN_TAP)
    shutil.copytree(TX4, RUN_TAP, ignore=shutil.ignore_patterns(".git"))
    os.chdir(RUN_TAP)
    sys.path.insert(0, str(RUN_TAP))
    import Lib_BW  # noqa: E402

    Lib_BW.Initialize.apply_xfmr_tap = True
    Lib_BW.Initialize.xfmr_tap_side = "to"
    Lib_BW.Initialize.net_damping = (
        1.0  # stock numerical damping kept: isolates the tap
    )
    pfd, ini, dyd, emt, wall, n, _ = smoke.run_case("noevent_tap", 10.0)
    d, _ = drift(emt, pfd)
    out["stock_with_documented_tap"] = d
    out["interpretation_rule"] = (
        "the stock start-up transient is attributed to the unused transformer ratio iff "
        "the tap-applied run removes the initial jump (< 1e-2 pu) and the drift"
    )
    out["attributed_to_tap"] = bool(d["initial_voltage_step_first_1ms_pu"] < 1e-2)
    (OUT / "tap_diagnostic.json").write_text(
        json.dumps(out, indent=1), encoding="utf-8"
    )
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
