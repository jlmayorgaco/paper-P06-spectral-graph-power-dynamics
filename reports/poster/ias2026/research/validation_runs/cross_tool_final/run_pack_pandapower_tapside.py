"""Translation-corrected Phase B run (label: B-corrected).

The as-supplied run (``run_pack_pandapower_compat.py``) FAILS static parity. The cause
is localised by ``diag_pandapower_ybus.py``. For a ppc transformer whose FROM bus is
the LV bus, pandapower ``from_ppc`` swaps hv/lv and ends up with the off-nominal tap on
the HV winding. The canonical JSON, the project Ybus, ANDES and MATPOWER all put the
tap at the branch FROM bus (bus1). This affects three branches only: 35 (31->6,
t=0.9), 37 (12->11, t=1.006) and 38 (12->13, t=1.006). The other 43 branch
admittances agree to 1e-13.

Passing ``tap_side="lv"`` to from_ppc does not help: from_ppc re-refers the impedance
and produces the identical Ybus (max diff 9.383 pu at (6,6) and (31,31) in both
variants; recorded in docs/FINAL_CROSS_TOOL_VALIDATION.md).

The correction reorients those three ppc rows into the exactly equivalent two-port
with the HV bus as FROM bus. A from-bus tap t with series admittance y gives
    [[y/t^2, -y/t], [-y/t, y]]          (order: LV, HV)
which is identical to a tap tau = 1/t at the HV bus with y' = y/t^2, i.e.
    r' = r t^2,  x' = x t^2,  tap' = 1/t,  from/to swapped.
All three affected rows have b = g = 0 and phi = 0, so no charging or phase shift is
involved; the script asserts that. The canonical pi-model admittance of every
reoriented row is checked against the original to 1e-12 before use. The canonical
JSON, the tolerances and the solver options are not changed.

Outputs go to ``pack_tapside/results`` so the as-supplied results stay untouched.
Usage: python run_pack_pandapower_tapside.py <research_root>
"""

from __future__ import annotations

import runpy
import shutil
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_pack_pandapower_compat as compat  # noqa: E402  (applies the three compat fixes)

cp = compat.canonical_ppc
_rated_canonical_to_ppc = cp.canonical_to_ppc
REORIENTED: list[dict] = []


def _two_port(r: float, x: float, b: float, tap: float) -> np.ndarray:
    y = 1.0 / complex(r, x)
    yc = 1j * b / 2.0
    return np.array([[(y + yc) / tap**2, -y / tap], [-y / tap, y + yc]])


def canonical_to_ppc_hv_from(repo: Path) -> dict:
    ppc = _rated_canonical_to_ppc(repo)
    kv = {int(row[cp.BUS_I]): row[cp.BASE_KV] for row in ppc["bus"]}
    br = ppc["branch"]
    for k in range(br.shape[0]):
        f, t = int(br[k, cp.F_BUS]), int(br[k, cp.T_BUS])
        tap = br[k, cp.TAP]
        if kv[f] >= kv[t] or tap in (0.0, 1.0):
            continue
        r, x, b, shift = br[k, cp.BR_R], br[k, cp.BR_X], br[k, cp.BR_B], br[k, cp.SHIFT]
        assert b == 0.0 and shift == 0.0, f"branch {k}: charging/shift not handled"
        new = (r * tap**2, x * tap**2, 1.0 / tap)
        before = _two_port(r, x, b, tap)  # order (f, t)
        after = _two_port(new[0], new[1], b, new[2])[
            ::-1, ::-1
        ]  # order (t, f) -> (f, t)
        err = float(np.max(np.abs(before - after)))
        assert err < 1e-12, f"branch {k}: reoriented two-port differs by {err}"
        br[k, cp.F_BUS], br[k, cp.T_BUS] = t, f
        br[k, cp.BR_R], br[k, cp.BR_X], br[k, cp.TAP] = new
        REORIENTED.append(
            {
                "branch": k,
                "from": f,
                "to": t,
                "tap": tap,
                "new_tap": new[2],
                "two_port_err": err,
            }
        )
    return ppc


cp.canonical_to_ppc = canonical_to_ppc_hv_from
PACK_TAPSIDE = HERE / "pack_tapside"

if __name__ == "__main__":
    PACK_TAPSIDE.mkdir(exist_ok=True)
    src = compat.PACK / "run_pandapower_parity.py"
    dst = PACK_TAPSIDE / "run_pandapower_parity.py"
    shutil.copyfile(src, dst)  # byte-identical copy, so __file__/results resolve here
    runpy.run_path(str(dst), run_name="__main__")
    for row in REORIENTED:
        print("REORIENTED", row)
