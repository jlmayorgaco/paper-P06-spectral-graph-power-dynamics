"""PD0 - timing and reproduction check before the preregistration is frozen.

Not a campaign result. It (a) regenerates one E35 unstable operating point and
checks that the frozen E35 flagship label is reproduced, (b) times one E34-M1
constraint evaluation, and (c) times one census-scale portfolio evaluation and
one gradient of the E9 design engine. The numbers set the compute budget of the
preregistration.
"""

from __future__ import annotations

import os

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_k, "1")

import sys  # noqa: E402
import time  # noqa: E402

import _pd_common as P  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import pd1_adaptive_retune as PD1  # noqa: E402


def main() -> int:
    out = {}
    table = pd.read_parquet(P.E35_TABLE)
    unstable = table[(table.status == "ACCEPTED") & (table.portfolio_unstable == True)]  # noqa: E712
    row = unstable.iloc[0]
    PD1.initialise()
    t0 = time.time()
    ctx = PD1.operating_point(row)
    out["pd1_point_s"] = time.time() - t0
    t0 = time.time()
    flag = PD1.flagship_labels(ctx)
    out["pd1_flagship_s"] = time.time() - t0
    out["pd1_reproduction_alpha_IA"] = [float(row.alpha_IA_exact), flag["alpha_IA"]]
    out["pd1_transverse"] = [flag["status_perp"], flag["alpha_perp"]]
    t0 = time.time()
    v = PD1.abscissa_perp(ctx, np.full(len(PD1.PHYSICAL), -0.15))
    out["pd1_constraint_eval_s"] = time.time() - t0
    out["pd1_constraint_value"] = v

    sys.path.insert(0, str(P.CDW_EXPERIMENTS))
    import _cdw as C  # noqa: E402
    import E09_design as E9  # noqa: E402

    theta = C.policy("D01")
    t0 = time.time()
    r = C.eval_portfolio(C.V9, theta, modes=False)
    out["pd2_portfolio_eval_s"] = time.time() - t0
    out["pd2_alpha_V9_D01"] = r["alpha"]
    plist = E9.params(C.V9, "joint")
    vals = {p: (theta[0] if p[0] == "g" else 1.0) for p in plist}
    t0 = time.time()
    E9.gradient(tuple(C.V9), theta, vals, plist)
    out["pd2_gradient_joint_s"] = time.time() - t0
    for k, val in out.items():
        print(f"{k:32s} {val}")
    P.write_json(P.out_dir("PD0") / "PD0_timing.json", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
