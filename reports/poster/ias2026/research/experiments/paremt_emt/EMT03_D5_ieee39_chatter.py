# ruff: noqa: E501  -- diagnostic tables kept on one line
"""EMT03 diagnostic D5: does the GFL current-source / trapezoidal-inductor alternation seen in the unit
test (D4) also arise in the IEEE-39 realization?  DIAGNOSTIC ONLY: no-event runs, no estimator, no
comparison with any phasor prediction (G3/G4 failed, so EMT04-EMT18 remain BLOCKED).

Cases (dt 50 us, 30 s, no disturbance, preregistered realization): P4 base, P4 H4, H4 at G_S.
Stored every 21 steps (odd, so a step-to-step (-1)^n component survives as a sample-to-sample one).
Reported per case and per window: max over GFL buses of |d|V|| between consecutive stored samples
(alternation proxy) and max over all buses of ||V| - |V(0)||.
Writes results/EMT03/EMT03_D5_ieee39_chatter.csv.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _emt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import tx4_case as C  # noqa: E402
import tx4_emt as K  # noqa: E402

DT, DS, TLEN = 50e-6, 21, 30.0
CASES = {"P4_BASE": ("EMT05_P4_BASE", None), "P4_H4": ("EMT05_P4_30+33+35+37", None),
         "GS_H4": ("EMT05_P4_30+33+35+37", [0.25, 1.425, 1.5, 1.0])}
WINDOWS = ((0.0, 0.1), (0.1, 1.0), (1.0, 10.0), (10.0, 30.0))


def main() -> int:
    data = C.Data(_emt.RESEARCH)
    rows = []
    for name, (cid, theta) in CASES.items():
        case = dict(_emt.CASE[cid])
        if theta is not None:
            case["theta"] = theta
        built = C.build_case(data, case, DT, amplitude=0.0)
        t, sg, gf, v, fin = C.run(built, DT, TLEN, DS)
        vm = np.abs(v[:, :, 0] + 1j * v[:, :, 1])
        gbus = [int(b) for b in built["gfl_par"][:, K.GF["bus"]]] if len(built["gfl_par"]) else []
        for a, b in WINDOWS:
            m = (t >= a) & (t < b)
            if m.sum() < 2:
                continue
            alt = float(np.abs(np.diff(vm[m][:, gbus], axis=0)).max()) if gbus else float("nan")
            alt_all = float(np.abs(np.diff(vm[m], axis=0)).max())
            dev = float(np.abs(vm[m] - vm[0]).max())
            rows.append({"case": name, "t_from": a, "t_to": b, "finite": bool(fin), "t_end": float(t[-1]),
                         "n_gfl": len(gbus), "max_alt_gfl_buses": alt, "max_alt_all_buses": alt_all,
                         "max_dev_vmag_from_t0": dev})
        print(name, "finite", fin, "t_end", t[-1])
    df = pd.DataFrame(rows)
    df.to_csv(_emt.RESULTS / "EMT03" / "EMT03_D5_ieee39_chatter.csv", index=False)
    print(df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
