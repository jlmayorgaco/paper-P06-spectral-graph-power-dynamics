# ruff: noqa: E501  -- table labels kept on one line
"""PCV07 - assemble the deliverable tables from the PCV02-PCV05 outputs (no computation).

- results/20260911_BASELINE_COMPARISON.csv   (Phases 2-3, Task E from PCV04)
- results/20260911_POLICY_COUNTERFACTUAL.csv (Phase 4)
- results/20260911_ROBUSTNESS.csv            (Phase 7)
Deterministic: reads CSV/JSON only.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
PCV = RESEARCH / "results" / "PCV"
OUT = RESEARCH / "results"
H4 = "30+33+35+37"


def read(stage, name):
    return pd.read_csv(PCV / stage / name)


def fmt(x, nd=4):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "N/A"
    return f"{x:.{nd}f}" if isinstance(x, float) else str(x)


# ---------------------------------------------------------------- baselines --
METHODS = [
    # id, column, label, inputs, needs_freq, needs_full_eval, output kind, orientation (+1 higher=riskier)
    (
        "B0",
        "alpha",
        "full transverse eigenanalysis (oracle)",
        "full DAE model of every evaluated portfolio",
        "no",
        "yes (each portfolio solved)",
        "verdict + margin",
        +1,
    ),
    (
        "B1-Pg",
        "pg_mw",
        "displaced active dispatch Pg [MW]",
        "power-flow dispatch",
        "no",
        "no",
        "score (no threshold)",
        +1,
    ),
    (
        "B1-Sn",
        "sn_mva",
        "replaced rating Sn [MVA]",
        "machine ratings",
        "no",
        "no",
        "score (no threshold)",
        +1,
    ),
    (
        "B1-pen",
        "penetration",
        "penetration Pg / total generation",
        "power-flow dispatch",
        "no",
        "no",
        "score (no threshold)",
        +1,
    ),
    (
        "B2-gSCR",
        "gscr",
        "generalized SCR of the replaced set (repository implementation)",
        "network Ybus, ratings",
        "no",
        "no",
        "score (no canonical threshold for this model)",
        -1,
    ),
    (
        "B2-SCR",
        "min_scr",
        "minimum nodal SCR over the replaced buses",
        "network Ybus, ratings",
        "no",
        "no",
        "score",
        -1,
    ),
    (
        "B2-MIIF",
        "max_miif",
        "maximum multi-infeed interaction factor",
        "network Ybus",
        "no",
        "no",
        "score (|S| >= 2)",
        +1,
    ),
    (
        "B3",
        "B3_alpha",
        "first-order modal sensitivity (F10 B5, full theta)",
        "base model + 1 % partial replacement per candidate",
        "no",
        "no",
        "verdict + margin",
        +1,
    ),
    (
        "B4",
        "B4_alpha",
        "additive singles (order-1 Moebius truncation of exact alpha)",
        "exact alpha of base and singles",
        "no",
        "no",
        "verdict + margin",
        +1,
    ),
    (
        "B5",
        "B5_alpha",
        "pairwise truncation (order 2)",
        "exact alpha of all subsets of size <= 2",
        "no",
        "no",
        "verdict + margin",
        +1,
    ),
    (
        "B5b",
        "B5b_alpha",
        "third-order truncation (order 3)",
        "exact alpha of all subsets of size <= 3",
        "no",
        "no",
        "verdict + margin",
        +1,
    ),
    (
        "B6",
        None,
        "participation factors of the base least-damped band mode",
        "base model eigenvectors",
        "no",
        "no",
        "localization only",
        0,
    ),
    (
        "B7",
        "compactness_score",
        "electrical compactness (-mean pairwise Thevenin distance, S excluded)",
        "network Ybus",
        "no",
        "no",
        "score (|S| >= 2)",
        +1,
    ),
    (
        "B8a",
        "closure_nonoracle",
        "closure distance, non-oracle (min over 121 band frequencies)",
        "port operator Q of the candidates (C2 realization)",
        "no",
        "no",
        "score (no threshold)",
        -1,
    ),
    (
        "B8b",
        "closure_oracle",
        "closure distance, oracle-assisted (at the true critical j omega*)",
        "port operator + the true critical frequency",
        "yes",
        "yes (needs omega* of the portfolio)",
        "score (no threshold)",
        -1,
    ),
    (
        "B9a",
        "abs_chi_nonoracle",
        "|chi_S| at the B8a frequency (secondary diagnostic)",
        "port operator",
        "no",
        "no",
        "score (secondary)",
        +1,
    ),
    (
        "B9b",
        "abs_chi_oracle",
        "|chi_S| at the true critical j omega* (secondary diagnostic)",
        "port operator + the true critical frequency",
        "yes",
        "yes",
        "score (secondary)",
        +1,
    ),
    (
        "FW",
        "alpha",
        "exact closure / hypergraph framework (C1-C3)",
        "C2 common realization of the candidates (one equilibrium, device blocks)",
        "no (zero of the exact characteristic function)",
        "no new equilibria; each subset verdict still needs a zero/eigen count",
        "verdict + H + kappa + port derivatives",
        +1,
    ),
]
LINEARIZATIONS = {  # full-model linearizations needed for the P4 16-subset decision (fair count)
    "B0": "16 (one per subset; 1 if only H4 is asked)",
    "B1-Pg": "0",
    "B1-Sn": "0",
    "B1-pen": "0",
    "B2-gSCR": "0 (network only)",
    "B2-SCR": "0 (network only)",
    "B2-MIIF": "0 (network only)",
    "B3": "5 (base + 4 partial replacements)",
    "B4": "5 (base + singles)",
    "B5": "11 (<= pairs)",
    "B5b": "15 (<= triples)",
    "B6": "1 (base eigenvectors)",
    "B7": "0 (network only)",
    "B8a": "5 device linearizations (C2) + 121 frequency evaluations",
    "B8b": "5 + the exact critical frequency (requires B0 of the portfolio)",
    "B9a": "as B8a",
    "B9b": "as B8b",
    "FW": "5 device linearizations at one common equilibrium (C2); exact characteristic function of every subset",
}


def baseline_table() -> pd.DataFrame:
    core = read("PCV02", "PCV02_core_portfolios.csv")
    hyper = read("PCV02", "PCV02_core_hypergraphs.csv")
    ta = read("PCV02", "PCV02_taskA.csv")
    tb = read("PCV02", "PCV02_taskB.csv")
    tc = read("PCV02", "PCV02_taskC.csv")
    td = read("PCV02", "PCV02_taskD.csv")
    te = read("PCV04", "PCV04_taskE.csv")
    summ = json.loads(
        (PCV / "PCV02" / "PCV02_summary.json").read_text(encoding="utf-8")
    )
    p4 = core[(core.point == "P4") & (core.subset == H4)].iloc[0]
    gs = core[(core.point == "G_S") & (core.subset == H4)].iloc[0]
    h_p4 = hyper[hyper.point == "P4"].iloc[0]
    h_gs = hyper[hyper.point == "G_S"].iloc[0]

    def a_metric(ds, col, key):
        r = ta[(ta.dataset == ds) & (ta.method == col)]
        return float(r[key].iloc[0]) if len(r) else float("nan")

    def b_metric(ds, mid):
        r = tb[(tb.dataset == ds) & (tb.method == mid)]
        return (
            f"{int(r.H_exact.iloc[0])}/{int(r.n_points.iloc[0])}" if len(r) else "N/A"
        )

    def c_metric(ds, mid, key):
        r = tc[(tc.dataset == ds) & (tc.method == mid)]
        return float(r[key].iloc[0]) if len(r) else float("nan")

    def e_metric(method, key):
        r = te[
            (te.intervention_set.str.startswith("12 frozen")) & (te.method == method)
        ]
        return float(r[key].iloc[0]) if len(r) else float("nan")

    te_map = {
        "B2-gSCR": "static: delta gSCR(H4) at 1.5x",
        "B7": "static: x-weighted betweenness (frozen F2)",
        "B3": "conventional full-DAE eigen-sensitivity (B3-type)",
        "FW": "port derivative (framework)",
    }
    rows = []
    for mid, col, label, inputs, freq, full, kind, orient in METHODS:
        rec = {
            "method_id": mid,
            "method": label,
            "inputs_required": inputs,
            "needs_critical_frequency": freq,
            "needs_complete_portfolio_evaluated": full,
            "output": kind,
            "orientation": {1: "higher = riskier", -1: "lower = riskier", 0: "N/A"}[
                orient
            ],
            "full_model_linearizations_P4_lattice": LINEARIZATIONS[mid],
        }
        # P4 decision on the complete portfolio
        if mid in ("B0", "FW"):
            rec.update(
                P4_H4_predicted_verdict=p4.status,
                P4_predicted_witness=h_p4.H_true,
                P4_H4_value=p4.alpha,
                G_S_H4_predicted_verdict=gs.status,
                G_S_predicted_witness=h_gs.H_true,
                G_S_H4_value=gs.alpha,
            )
        elif mid in ("B3", "B4", "B5", "B5b"):
            unst = (
                (lambda r: bool(r.B3_verdict_unstable))
                if mid == "B3"
                else (lambda r, c=col: bool(r[c] > 0))
            )
            rec.update(
                P4_H4_predicted_verdict="UNSTABLE" if unst(p4) else "STABLE",
                P4_predicted_witness=h_p4[f"{mid}_H"],
                P4_H4_value=p4[col],
                G_S_H4_predicted_verdict="UNSTABLE" if unst(gs) else "STABLE",
                G_S_predicted_witness=h_gs[f"{mid}_H"],
                G_S_H4_value=gs[col],
            )
        elif mid == "B6":
            rec.update(
                P4_H4_predicted_verdict="N/A (localization only)",
                P4_predicted_witness="N/A",
                P4_H4_value=float("nan"),
                G_S_H4_predicted_verdict="N/A",
                G_S_predicted_witness="N/A",
                G_S_H4_value=float("nan"),
            )
        else:
            rec.update(
                P4_H4_predicted_verdict="N/A (no stability threshold)",
                P4_predicted_witness="N/A",
                P4_H4_value=p4[col],
                G_S_H4_predicted_verdict="N/A",
                G_S_predicted_witness="N/A",
                G_S_H4_value=gs[col],
            )
        rec["value_changes_P4_to_G_S"] = (
            "N/A"
            if not np.isfinite(rec["P4_H4_value"])
            else bool(abs(rec["P4_H4_value"] - rec["G_S_H4_value"]) > 1e-12)
        )
        # Task A
        if col and mid not in ("B0", "FW"):
            rec.update(
                taskA_census4_roc_auc=a_metric("census_size4", col, "roc_auc"),
                taskA_census4_avg_precision=a_metric(
                    "census_size4", col, "average_precision"
                ),
                taskA_census4_spearman_alpha=a_metric(
                    "census_size4", col, "spearman_with_alpha"
                ),
                taskA_census4_perm_p=a_metric("census_size4", col, "permutation_p"),
                taskA_F10pooled_roc_auc=a_metric("F10_points_pooled", col, "roc_auc"),
            )
        else:
            rec.update(
                taskA_census4_roc_auc=float("nan"),
                taskA_census4_avg_precision=float("nan"),
                taskA_census4_spearman_alpha=float("nan"),
                taskA_census4_perm_p=float("nan"),
                taskA_F10pooled_roc_auc=float("nan"),
            )
        # Task B
        if mid in ("B3", "B4", "B5", "B5b"):
            cm = summ["census_minimum_set"][mid]
            rec.update(
                taskB_F10_52_H_exact=b_metric("F10_52", mid),
                taskB_F10_line36_H_exact=b_metric("F10_line36", mid),
                taskB_prereg4_H_exact=b_metric("prereg_points", mid),
                taskB_census_H_exact=str(cm["H_exact"]),
                taskB_census_missed_units=cm["missed_units"],
            )
        elif mid in ("B0", "FW"):
            rec.update(
                taskB_F10_52_H_exact="52/52 (defines truth)"
                if mid == "B0"
                else "= B0 by the C2 identity (not re-run; identity residual <= 3e-10 at P4 and G_S)",
                taskB_F10_line36_H_exact="36/36",
                taskB_prereg4_H_exact="4/4",
                taskB_census_H_exact="True",
                taskB_census_missed_units=0,
            )
        elif mid == "B6":
            rec.update(
                taskB_F10_52_H_exact="N/A",
                taskB_F10_line36_H_exact="N/A",
                taskB_prereg4_H_exact="N/A",
                taskB_census_H_exact=f"N/A (localization: {summ['census_localization']['hyperedges_equal_to_top_k_participation']}/"
                f"{summ['census_localization']['n_hyperedges']} hyperedges equal a top-k participation set)",
                taskB_census_missed_units="N/A",
            )
        else:
            rec.update(
                taskB_F10_52_H_exact="N/A",
                taskB_F10_line36_H_exact="N/A",
                taskB_prereg4_H_exact="N/A",
                taskB_census_H_exact="N/A",
                taskB_census_missed_units="N/A",
            )
        # Task C
        if mid in ("B3", "B4", "B5", "B5b"):
            rec.update(
                taskC_F10_52_sign_accuracy=c_metric("F10_52", mid, "sign_accuracy"),
                taskC_F10_52_mae=c_metric("F10_52", mid, "mae"),
                taskC_F10_52_rmse=c_metric("F10_52", mid, "rmse"),
            )
        else:
            rec.update(
                taskC_F10_52_sign_accuracy=1.0 if mid in ("B0", "FW") else float("nan"),
                taskC_F10_52_mae=0.0 if mid in ("B0", "FW") else float("nan"),
                taskC_F10_52_rmse=0.0 if mid in ("B0", "FW") else float("nan"),
            )
        # Task D
        r = td[td.method == col] if col and mid not in ("B0", "FW") else td.iloc[0:0]
        rec["taskD_within_portfolio_auc"] = (
            float(r.mean_within_portfolio_auc.iloc[0])
            if len(r)
            else (1.0 if mid in ("B0", "FW") else float("nan"))
        )
        # Task E
        m = te_map.get(mid)
        rec["taskE_lines_spearman"] = e_metric(m, "spearman") if m else float("nan")
        rec["taskE_lines_kendall"] = e_metric(m, "kendall") if m else float("nan")
        rec["taskE_method"] = m or "N/A"
        rows.append(rec)
    return pd.DataFrame(rows)


# ---------------------------------------------------------- counterfactual --
def counterfactual_table() -> pd.DataFrame:
    cf = read("PCV03", "PCV03_counterfactual.csv")
    summ = json.loads(
        (PCV / "PCV03" / "PCV03_summary.json").read_text(encoding="utf-8")
    )
    vmax = max(summ["identity"][p]["max_voltage_diff_vs_P4"] for p in ("G_S", "G_S2"))
    rows = []
    for r in cf.itertuples():
        same_clean = bool(r.identical_on_clean_line)
        note = r.note if isinstance(r.note, str) else ""
        if r.quantity == "H4 equilibrium voltages (hash)":
            same_clean = bool(vmax <= 1e-6)
            note = (
                f"identical within power-flow tolerance: max abs difference {vmax:.1e} (hash differs at 1e-10 rounding); "
                + note
            )
        same_pinf = str(r.P4) == str(r.P_inf)
        if r.quantity == "H4 equilibrium voltages (hash)":
            same_pinf = bool(
                summ["identity"]["P_inf"]["max_voltage_diff_vs_P4"] <= 1e-6
            )
        rows.append(
            {
                "quantity": r.quantity,
                "kind": r.kind,
                "P4": r.P4,
                "G_S": r.G_S,
                "G_S2": r.G_S2,
                "P_inf": r.P_inf,
                "unchanged_on_clean_g_line_P4_G_S_G_S2": same_clean,
                "unchanged_P4_vs_P_inf_multi_coordinate": same_pinf,
                "note": note,
            }
        )
    df = pd.DataFrame(rows)
    df.insert(
        0,
        "comparison",
        "clean: only g changes (k = 1.425, t = 1.5, h = 1); P_inf = multi-coordinate (g and k change)",
    )
    return df


# -------------------------------------------------------------- robustness --
SOURCES = {
    "NOMINAL": "nominal P4",
    "EM-f": "E37 machine envelope, fleet-wide (project source)",
    "EM-u": "E37 bounds, per-unit independent factors",
    "EC": "converter stress test (no project source)",
    "EMC": "EM-u + EC",
    "OAT-M": "E37 group extremes, one at a time",
    "OAT-C": "converter group extremes, one at a time",
}


def decision(x):
    return "PASS" if x >= 0.8 else ("PARTIAL" if x >= 0.5 else "FAIL")


def robustness_table() -> pd.DataFrame:
    d = read("PCV05", "PCV05_draws.csv")
    rows = []
    for env in ("NOMINAL", "EM-f", "EM-u", "EC", "EMC", "OAT-M", "OAT-C"):
        g = d[d.envelope == env]
        prereg = env in ("EM-f", "EM-u", "EC", "EMC")
        kc = g.kappa_class.astype(str).value_counts()
        wit = g.H.value_counts()
        a = g.alpha_H4
        f = g.crit_hz_H4
        gs = g.g_star.dropna()
        unst = g[g.alpha_H4 > 0]
        r1, r2 = float(g.non_composable.mean()), float(g.H_is_H4.mean())
        k4 = float((g.kappa_class.astype(str) == "4").mean())
        direction = float((g.alpha_H4_g1 < g.alpha_H4).mean())
        rows.append(
            {
                "envelope": env,
                "source": SOURCES[env],
                "type": "deterministic stress envelope (coverage fraction, not a probability)",
                "n_draws": len(g),
                "n_unresolved_or_failed": int((~g.exact.astype(bool)).sum()),
                "R1_non_composable_fraction": r1,
                "R1_decision": decision(r1) if prereg else "descriptive",
                "R2_witness_H4_fraction": r2,
                "R2_decision": decision(r2) if prereg else "descriptive",
                "R3_kappa_4_fraction": k4,
                **{
                    f"R3_kappa_{k}": int(kc.get(k, 0))
                    for k in ("1", "2", "3", "4", "inf", "base-unstable")
                },
                "witness_frequencies": json.dumps(
                    {k: int(v) for k, v in wit.sort_index().items()}
                ),
                "alpha_H4_min": float(a.min()),
                "alpha_H4_q25": float(a.quantile(0.25)),
                "alpha_H4_median": float(a.median()),
                "alpha_H4_q75": float(a.quantile(0.75)),
                "alpha_H4_max": float(a.max()),
                "crit_hz_H4_min": float(f.min()),
                "crit_hz_H4_median": float(f.median()),
                "crit_hz_H4_max": float(f.max()),
                "R4_n_crossing": len(gs),
                "R4_g_star_min": float(gs.min()) if len(gs) else float("nan"),
                "R4_g_star_q25": float(gs.quantile(0.25)) if len(gs) else float("nan"),
                "R4_g_star_median": float(gs.median()) if len(gs) else float("nan"),
                "R4_g_star_q75": float(gs.quantile(0.75)) if len(gs) else float("nan"),
                "R4_g_star_max": float(gs.max()) if len(gs) else float("nan"),
                "R4_direction_fraction_prereg": direction,
                "R4_direction_robust_prereg": bool(direction >= 0.9),
                "posthoc_restored_at_g1_given_H4_unstable": f"{int((unst.alpha_H4_g1 < 0).sum())}/{len(unst)}",
                "robust_phenomenon": (
                    "yes" if r1 >= 0.8 else "partial" if r1 >= 0.5 else "no"
                )
                if prereg
                else "descriptive",
                "robust_exact_witness": (
                    "yes" if r2 >= 0.8 else "partial" if r2 >= 0.5 else "no"
                )
                if prereg
                else "descriptive",
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    baseline_table().to_csv(OUT / "20260911_BASELINE_COMPARISON.csv", index=False)
    counterfactual_table().to_csv(
        OUT / "20260911_POLICY_COUNTERFACTUAL.csv", index=False
    )
    robustness_table().to_csv(OUT / "20260911_ROBUSTNESS.csv", index=False)
    print("wrote 3 deliverables")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
