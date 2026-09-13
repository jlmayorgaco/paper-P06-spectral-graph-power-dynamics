# ruff: noqa: E501
"""R1 model audit (BASE CASE ONLY, S = empty, documented excitation k = 1): no converter is present.

Checks, per variant:
- Ybus rebuilt from the JSON equals the frozen network's Ybus;
- power flow against Table 1 (Gate 3 reproduction);
- equilibrium residual; governor/damping wiring (equilibrium identical across variants);
- structural centre: number of transverse-removed zeros; frozen G3 base-mode check (SP33 vs Table 4);
- documented limit margins; state counts; timing.
No portfolio with a converter is evaluated here (preregistration rule).
"""

from __future__ import annotations

import _r68 as R  # noqa: I001

import json
import math
import time

import numpy as np
import pandas as pd

import G3_ieee68 as G3


def main():
    R.env_threads()
    net = R.base_network()
    y2 = R.build_ybus68()
    ybus_err = float(np.abs(y2 - net.ybus).max())
    pf_chk = G3.eligibility()  # frozen G3 code: power flow and base modes of the published model
    rows, cases = [], {}
    for v in R.VARIANTS:
        t0 = time.time()
        case = R.build68((), variant=v, g=0.0, k=1.0)
        cases[v] = case
        rec = R.evaluate(case, modes=True)
        vals = np.linalg.eigvals(R.matrices(case)[0])
        n_zero = int((np.abs(vals) < R.STRUCT_ZERO).sum())
        rows.append({"variant": v, "n_x": rec["n_x"], "n_dead": rec["n_dead"], "eq_residual": rec["eq_residual"],
                     "n_zero_full_A": n_zero, "status": rec["status"], "alpha_perp": rec["alpha"], "lam_hz": rec["lam_hz"],
                     "em_top_re": rec["em_top_re"], "em_top_hz": rec["em_top_hz"], "gz_cond": rec["gz_cond"], "wall_s": time.time() - t0})
    # equilibria: network voltages identical across variants (governors/damping do not move the operating point)
    zr = cases["SP33"].equilibrium.z
    dz = {v: float(np.abs(cases[v].equilibrium.z - zr).max()) for v in R.VARIANTS}
    margins = {v: float(min(getattr(s.device, "limit_margin", math.inf) if not hasattr(s.device, "base") else s.device.base.limit_margin for s in cases[v].dae.slots)) for v in R.VARIANTS}
    gov = []
    for s in cases["REAL"].dae.slots:
        if hasattr(s.device, "tg"):
            gov.append({"bus": s.bus, "kg_device_base": s.device.tg.kg, "pref": s.device.pref, "damping_device_base": s.device.parameters.d,
                        "rating_mva": R.rating(s.bus), "sn_sp_mva": R.sn_sp(s.bus)})
    out = {"ybus_rebuild_max_abs_err": ybus_err, "gate3_power_flow": pf_chk["power_flow"], "gate3_interarea_verdict": pf_chk["interarea_verdict"],
           "gate3_max_abs_df_all15": pf_chk["max_abs_d_f_hz_all15"], "gate3_max_abs_dzeta_all15": pf_chk["max_abs_d_zeta_pct_all15"],
           "variants": rows, "equilibrium_z_diff_vs_SP33": dz, "min_limit_margin": margins, "governors": gov}
    R.write_json(R.RESULTS / "CDW68_R01_audit.json", out)
    print(json.dumps({k: v for k, v in out.items() if k != "governors"}, indent=1, default=str))
    pd.DataFrame(gov).to_csv(R.RESULTS / "CDW68_R01_governors.csv", index=False)


if __name__ == "__main__":
    main()
