"""FC04 (amendment I): E14 N6 redone on the Pg-matched controls.

The modal-family definition is the frozen one.

N6 asked whether the matched stable controls hide the same inter-area branch that
goes unstable in the failing quadruples. The controls were matched on the MVA
rating; UC02 re-selected them on the measured active dispatch Pg (failing Pg
range +-5 %, top 25). Here each control, each failing size-4 portfolio and the
flagship is read with the FROZEN v2C modal-family observable (_v2c_common:
anchor = base-case inter-area mode by the frozen nominal shape, family = every
0.3-1.5 Hz mode with overlap >= 0.80 on the common machine coordinates,
endpoint = the family envelope). Never argmax-MAC; a tracking failure is
reported as such.
"""

from __future__ import annotations

import sys
import time

import numpy as np
import pandas as pd
from _fc import RESULTS, FCExperiment, out_dir, write_json

from _v2c_common import base_anchor, nominal_reference, read_family  # noqa: E402
from ibr_cycles.dynamics.modes import eigen_analysis  # noqa: E402
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case  # noqa: E402

OUT = out_dir("FC04_e14_n6_pg")


def parse(m):
    return tuple(int(b) for b in m.split("+"))


def main(argv) -> int:
    exp = FCExperiment(
        name="FC04_e14_n6_pg",
        question="Do the Pg-matched stable controls carry the failing branch damped?",
        config={
            "definition": "frozen v2C family (_v2c_common)",
            "match": "UC02 Pg rule",
        },
    )
    started = time.time()
    uc = pd.read_csv(RESULTS / "UC" / "UC01" / "UC01_census_quantities.csv")
    uc = uc[uc["size"] == 4].copy()
    uc["unstable"] = uc.unstable.map(
        {True: True, False: False, "True": True, "False": False}
    )
    fail = uc[uc.unstable]
    low, high = fail.replaced_pg_mw.min(), fail.replaced_pg_mw.max()
    stable = uc[~uc.unstable]
    controls = stable[
        (stable.replaced_pg_mw >= 0.95 * low) & (stable.replaced_pg_mw <= 1.05 * high)
    ]
    controls = controls.nlargest(25, "replaced_pg_mw")
    nominal = nominal_reference()
    base = solve_case(ReplacementPlan.of({}))
    anchor = base_anchor(eigen_analysis(base.system.A), base, nominal)
    rows = []
    for role, members in [("failing", m) for m in fail.members] + [
        ("control", m) for m in controls.members
    ]:
        case = solve_case(ReplacementPlan.of({b: 1.0 for b in parse(members)}))
        spec = eigen_analysis(case.system.A)
        fr = read_family(anchor, base, case, spec)
        rows.append(
            {
                "role": role,
                "portfolio": members,
                "replaced_pg_mw": float(
                    uc.loc[uc.members == members, "replaced_pg_mw"].iloc[0]
                ),
                "tracked": fr is not None,
                "family_size": None
                if fr is None
                else len(fr.family.members)
                if hasattr(fr.family, "members")
                else None,
                "family_alpha": None if fr is None else fr.alpha,
                "family_worst_hz": None if fr is None else fr.frequency_worst_hz,
                "family_worst_damping": None
                if fr is None
                else float(
                    -fr.alpha / np.hypot(fr.alpha, 2 * np.pi * fr.frequency_worst_hz)
                ),
                "rightmost_nonzero_real": float(
                    max(m.real for m in spec.modes if abs(m.value) > 1e-3)
                ),
            }
        )
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "FC04_n6_pg.csv", index=False)
    ctrl = frame[frame.role == "control"]
    fl = frame[frame.role == "failing"]
    summary = {
        "controls": int(len(ctrl)),
        "controls_tracked": int(ctrl.tracked.sum()),
        "controls_family_alpha_max": float(ctrl.family_alpha.max()),
        "controls_family_damping_min": float(ctrl.family_worst_damping.min()),
        "failing": int(len(fl)),
        "failing_tracked": int(fl.tracked.sum()),
        "failing_family_alpha_range": [
            float(fl.family_alpha.min()),
            float(fl.family_alpha.max()),
        ],
        "failing_family_positive": int((fl.family_alpha > 0).sum()),
        "elapsed_s": round(time.time() - started, 1),
    }
    write_json(OUT / "FC04_summary.json", summary)
    exp.finish("COMPUTED", **summary)
    print(frame.to_string())
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
