# ruff: noqa: E501  -- table labels kept on one line
"""CC06 - Phase 6: stability sensitivity ds*/da vs connected-mechanism sensitivity d chi_H/da.

Preregistered (connected_cumulants_prereg_v1.yaml, de09db3b). Frozen parameters only:
- policy coordinates g, k, t, h at the ten FC18 boundaries (box-normalized by SCALE);
- the 12 frozen holdout lines (seed 20260911) at the FLAG boundary, gamma = 1 -/+ 0.002
  scaling the complete branch two-port (series + charging, tap unchanged; ybus_scaled
  copied verbatim from the frozen F2c script), re-equilibrated C2 realization.

dF(B)/da = tr(adj(I + Q_BB) dQ_BB/da) (valid at the singular boundary); d chi via the
partition product rule; ds*/da = -dF(H)/da / dF(H)/ds. Stability sensitivity for the
policy coordinates is also taken from the FROZEN FC18 grad_port and compared.

H6: argmin Re(ds*/da) == argmin d|chi_H|/da.
    lines: SUPPORTED iff top-1 agrees AND Spearman >= 0.8 over the 12 lines;
    policy: SUPPORTED iff top-1 agrees at >= 8 of 10 boundaries.
"""

from __future__ import annotations

import json
import sys
import time
from dataclasses import replace
from multiprocessing import Pool

import numpy as np
import pandas as pd
from _cc import (
    CORE,
    HOLDOUT,
    NAMES,
    ROOT,
    SCALE,
    STEPS,
    Realization,
    frozen_events,
    out_dir,
    write_json,
)
from scipy.stats import spearmanr

from _overnight import pin_blas_threads
from ibr_cycles.cycles.connected import (
    characteristic_derivatives,
    characteristic_values,
    cumulant_derivatives,
    cumulants,
)
from ibr_cycles.models.ieee39_network import load_network

OUT = out_dir("CC06")
H_S = 1e-5
EPS = 0.002
F2C = (
    ROOT
    / "validation_inputs/cross_tool_handoff/claude_cross_tool_handoff"
    / "dynamic_forest_ieee39_validation/dynamic_forest_ieee39_validation"
    / "data/F2c_holdout_reequilibrated_port_line_sensitivity.csv"
)
payload = json.loads((ROOT / "configs/ias2026/ieee39_network.json").read_text())


# ---- verbatim from dynamic_forest_line_port_holdout.py (frozen F2c) ----
def ybus_scaled(line_idx, gamma):
    order = [int(b["idx"]) for b in payload["buses"]]
    index = {b: i for i, b in enumerate(order)}
    n = len(order)
    y = np.zeros((n, n), complex)
    for ell, line in enumerate(payload["lines"]):
        if float(line["u"]) == 0:
            continue
        f, t = index[int(line["bus1"])], index[int(line["bus2"])]
        scale = gamma if ell == line_idx else 1.0
        series = scale / complex(line["r"], line["x"])
        charging = scale * complex(line["g"], line["b"]) / 2
        m = float(line["tap"]) * np.exp(1j * float(line["phi"]))
        m2 = abs(m) ** 2
        y[f, f] += (series + charging) / m2
        y[t, t] += series + charging
        y[f, t] += -series / np.conj(m)
        y[t, f] += -series / m
    for sh in payload["shunts"]:
        y[index[int(sh["bus"])], index[int(sh["bus"])]] += complex(sh["g"], sh["b"])
    return y


# ---- end verbatim ----


def q_task(task):
    """Return Q (8x8) at the task's s (and at s -/+ H_S for 'center')."""

    pin_blas_threads()
    kind, theta, s, extra = task
    kw = {}
    if kind == "line":
        li, gamma = extra
        net = load_network(ROOT / "configs/ias2026/ieee39_network.json")
        kw["network"] = replace(net, ybus=ybus_scaled(li, gamma))
    r = Realization(CORE, theta, **kw)
    if kind == "center":
        return [r.q_matrix(s + d)[0] for d in (0.0, H_S, -H_S)]
    return r.q_matrix(s)[0]


def h_cols(h):
    return [2 * CORE.index(b) + c for b in h for c in (0, 1)]


def sensitivities(q0, dq, q_sp, q_sm, h):
    cols = h_cols(h)
    blocks = [2] * len(h)
    sub = np.ix_(cols, cols)
    f = characteristic_values(q0[sub], blocks)
    df = characteristic_derivatives(q0[sub], dq[sub], blocks)
    chi = cumulants(f, range(len(h)))
    dchi = cumulant_derivatives(f, df, range(len(h)))
    full = frozenset(range(len(h)))
    fp = characteristic_values(q_sp[sub], blocks)[full]
    fm = characteristic_values(q_sm[sub], blocks)[full]
    dfs = (fp - fm) / (2 * H_S)
    ds = -df[full] / dfs
    dabs = float((np.conj(chi[full]) * dchi[full]).real / abs(chi[full]))
    sigma = dchi[full] / df[full] if abs(df[full]) > 0 else np.nan
    return {
        "ds": complex(ds),
        "dchi": complex(dchi[full]),
        "d_abs_chi": dabs,
        "sigma": complex(sigma),
        "dF": complex(df[full]),
        "abs_chi": abs(chi[full]),
    }


def main(argv) -> int:
    started = time.time()
    events = frozen_events()
    holdout = json.loads(HOLDOUT.read_text())
    assert holdout["seed"] == 20260911
    lines = [int(x) for x in holdout["lines"]]
    flag = next(e for e in events if e["event"].startswith("FLAG"))

    tasks, keys = [], []
    for i, e in enumerate(events):
        tasks.append(("center", e["theta"], e["s_star"], None))
        keys.append(("center", i))
        for j, name in enumerate(NAMES):
            for sign in (+1, -1):
                th = list(e["theta"])
                th[j] += sign * STEPS[name]
                tasks.append(("policy", tuple(th), e["s_star"], None))
                keys.append(("policy", i, name, sign))
    for li in lines:
        for sign, gamma in ((+1, 1 + EPS), (-1, 1 - EPS)):
            tasks.append(("line", flag["theta"], flag["s_star"], (li, gamma)))
            keys.append(("line", li, sign))
    with Pool(16) as pool:
        results = pool.map(q_task, tasks, chunksize=1)
    got = dict(zip(keys, results, strict=True))

    # ---- policy coordinates at the ten boundaries --------------------------------
    pol_rows, agree = [], 0
    for i, e in enumerate(events):
        q0, q_sp, q_sm = got[("center", i)]
        h = tuple(sorted(e["subset"]))
        per = {}
        for j, name in enumerate(NAMES):
            dq = (got[("policy", i, name, +1)] - got[("policy", i, name, -1)]) / (
                2 * STEPS[name]
            )
            sens = sensitivities(q0, dq, q_sp, q_sm, h)
            frozen = e["grad_port"][name]
            per[name] = {
                "Re_ds_scaled": SCALE[j] * sens["ds"].real,
                "Re_ds_frozen_scaled": SCALE[j] * frozen.real,
                "d_abs_chi_scaled": SCALE[j] * sens["d_abs_chi"],
                "abs_sigma": abs(sens["sigma"]),
                "Re_sigma": sens["sigma"].real,
            }
        stab = min(per, key=lambda n: per[n]["Re_ds_frozen_scaled"])
        supp = min(per, key=lambda n: per[n]["d_abs_chi_scaled"])
        agree += stab == supp
        for name, v in per.items():
            pol_rows.append(
                {
                    "event": e["event"],
                    "H": "+".join(map(str, h)),
                    "coordinate": name,
                    **v,
                    "most_stabilizing": stab,
                    "most_suppressing": supp,
                }
            )
    policy = pd.DataFrame(pol_rows)
    policy.to_csv(OUT / "CC06_policy.csv", index=False)
    port_vs_frozen = float(
        np.max(
            np.abs(policy.Re_ds_scaled - policy.Re_ds_frozen_scaled)
            / np.maximum(np.abs(policy.Re_ds_frozen_scaled), 1e-12)
        )
    )

    # ---- lines at FLAG -----------------------------------------------------------------
    i_flag = events.index(flag)
    q0, q_sp, q_sm = got[("center", i_flag)]
    f2c = pd.read_csv(F2C).set_index("line_index")
    line_rows = []
    for li in lines:
        dq = (got[("line", li, +1)] - got[("line", li, -1)]) / (2 * EPS)
        sens = sensitivities(q0, dq, q_sp, q_sm, CORE)
        line_rows.append(
            {
                "line_index": li,
                "line": f"L{li:02d}:{int(payload['lines'][li]['bus1'])}-{int(payload['lines'][li]['bus2'])}",
                "Re_ds_dgamma": sens["ds"].real,
                "Im_ds_dgamma": sens["ds"].imag,
                "d_abs_chi_dgamma": sens["d_abs_chi"],
                "Re_dchi": sens["dchi"].real,
                "Im_dchi": sens["dchi"].imag,
                "abs_sigma": abs(sens["sigma"]),
                "Re_sigma": sens["sigma"].real,
                "frozen_F2c_port_Re_ds": float(f2c.loc[li, "port_reeq_ds_real"]),
                "frozen_F2c_dae_Re_dlambda": float(f2c.loc[li, "dae_dlambda_real"]),
            }
        )
    lines_df = pd.DataFrame(line_rows)
    lines_df.to_csv(OUT / "CC06_lines.csv", index=False)
    top_stab = int(lines_df.loc[lines_df.Re_ds_dgamma.idxmin(), "line_index"])
    top_supp = int(lines_df.loc[lines_df.d_abs_chi_dgamma.idxmin(), "line_index"])
    rho = float(spearmanr(lines_df.Re_ds_dgamma, lines_df.d_abs_chi_dgamma).statistic)
    rho_f2c = float(
        spearmanr(lines_df.Re_ds_dgamma, lines_df.frozen_F2c_port_Re_ds).statistic
    )
    summary = {
        "policy_top1_agreement": f"{agree}/10",
        "policy_H6": "SUPPORTED" if agree >= 8 else "REJECTED",
        "policy_ds_recomputed_vs_frozen_max_rel": port_vs_frozen,
        "lines_top_stabilizing": top_stab,
        "lines_top_suppressing": top_supp,
        "lines_spearman_Re_ds_vs_d_abs_chi": rho,
        "lines_H6": "SUPPORTED"
        if (top_stab == top_supp and rho >= 0.8)
        else "REJECTED",
        "lines_spearman_vs_frozen_F2c_port": rho_f2c,
        "lines_sign_agreement_vs_frozen_F2c_port": int(
            np.sum(
                np.sign(lines_df.Re_ds_dgamma)
                == np.sign(lines_df.frozen_F2c_port_Re_ds)
            )
        ),
        "median_abs_sigma_lines": float(lines_df.abs_sigma.median()),
        "median_abs_sigma_policy": float(policy.abs_sigma.median()),
        "elapsed_s": round(time.time() - started, 1),
    }
    write_json(OUT / "CC06_summary.json", summary)
    pd.set_option("display.width", 250)
    print(policy.to_string(index=False))
    print(lines_df.to_string(index=False))
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
