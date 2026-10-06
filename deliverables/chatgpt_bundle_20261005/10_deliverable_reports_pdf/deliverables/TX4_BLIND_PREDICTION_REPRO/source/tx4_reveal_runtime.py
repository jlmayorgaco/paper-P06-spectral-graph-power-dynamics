"""Post-reveal comparison and full-order runtime benchmark for TX4.

This script is intentionally separate from tx4_blind_predict.py. It is the
first command allowed to read the frozen full-order answer tables.
"""

from __future__ import annotations

import json
import math
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
EXP = HERE.parent
SRC = EXP.parent / "src"
ROOT = HERE.parents[5]
sys.path.insert(0, str(EXP))
sys.path.insert(0, str(SRC))

from _f7_common import LEAK, Theta  # noqa: E402
from ibr_cycles.certification.symmetry import frequency_partner, rotation_generator  # noqa: E402
from ibr_cycles.certification.transverse import transverse_operator  # noqa: E402
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case  # noqa: E402
from ibr_cycles.models.ieee39_devices import ConverterParameters  # noqa: E402

OUT = ROOT / "results"
V4 = (30, 33, 35, 37)
V9 = tuple(range(30, 39))
P4 = Theta(g=0.03625, k=1.425, t=1.5, h=1.0)


def direct_case(members):
    return solve_case(
        ReplacementPlan.of({b: 1.0 for b in members}),
        converter=ConverterParameters(
            voltage_control=True, voltage_gain=P4.g, voltage_leak=LEAK
        ),
        machine_scaling={"ka": P4.k, "ta": P4.t},
    )


def full_mode(case):
    rx, _ = rotation_generator(case.dae, case.equilibrium.z)
    fp = frequency_partner(case.dae)
    a = transverse_operator(case.system.A, rx, fp.w).a_perp
    vals = np.linalg.eigvals(a)
    freq = abs(vals.imag) / (2 * math.pi)
    keep = (vals.imag >= 0) & (abs(vals) > 1e-3) & (freq >= 0.3) & (freq <= 1.5)
    if not keep.any():
        return float("nan"), float("nan")
    v = vals[np.flatnonzero(keep)]
    x = v[np.argmax(v.real)]
    return float(x.real), float(abs(x.imag) / (2 * math.pi))


def label_map(path: Path, subset_col: str = "members"):
    frame = pd.read_csv(path)
    return {str(r[subset_col]): str(r["status"]) for _, r in frame.iterrows()}


def minimal_edges(labels):
    sets = {frozenset(int(x) for x in k.split("+")) if k != "BASE" else frozenset() for k, v in labels.items() if v == "UNSTABLE"}
    return sorted(
        [s for s in sets if not any(t < s for t in sets)],
        key=lambda s: (len(s), tuple(sorted(s))),
    )


def compare(label: str, pred_path: Path, truth_path: Path, truth_col: str = "members"):
    pred = pd.read_csv(pred_path)
    truth = pd.read_csv(truth_path)
    if "point" in truth.columns:
        truth = truth[truth.point == "P4"].copy()
    truth["portfolio"] = truth[truth_col]
    if "alpha_perp" in truth.columns:
        truth["alpha_truth"] = truth["alpha_perp"]
    else:
        truth["alpha_truth"] = truth["alpha"]
    if "crit_hz" in truth.columns:
        truth["frequency_truth_hz"] = truth["crit_hz"]
    else:
        truth["frequency_truth_hz"] = np.nan
    truth["truth_verdict"] = truth["status"].map({"STABLE": "STABLE_PREDICTED", "UNSTABLE": "UNSTABLE_PREDICTED"})
    # PCV02 is a core-portfolio table and intentionally omits the all-SG BASE
    # row.  Add that frozen reference row from FC01 rather than treating the
    # absent historical row as a prediction miss.
    if label == "V4" and "BASE" not in set(truth["portfolio"].astype(str)):
        base_path = ROOT / "reports/poster/ias2026/research/results/FINAL_CLOSURE/FC01_structure.csv"
        base = pd.read_csv(base_path)
        base = base[base.subset.astype(str) == "BASE"].iloc[0]
        truth = pd.concat(
            [
                truth,
                pd.DataFrame(
                    [
                        {
                            "portfolio": "BASE",
                            "status": "STABLE",
                            "alpha_truth": float(base.alpha_perp),
                            "frequency_truth_hz": np.nan,
                            "truth_verdict": "STABLE_PREDICTED",
                        }
                    ]
                ),
            ],
            ignore_index=True,
        )
    keep = ["portfolio", "status", "truth_verdict", "alpha_truth", "frequency_truth_hz"]
    joined = pred.merge(truth[keep], on="portfolio", how="left")
    joined["verdict_correct"] = joined.predicted_verdict == joined.truth_verdict
    joined["false_safe"] = (joined.predicted_verdict == "STABLE_PREDICTED") & (joined.status == "UNSTABLE")
    joined["false_unstable"] = (joined.predicted_verdict == "UNSTABLE_PREDICTED") & (joined.status == "STABLE")
    joined.to_csv(OUT / f"TX4_{label}_BLIND_VS_FULL.csv", index=False)
    true_labels = {str(r[truth_col]): str(r.status) for _, r in truth.iterrows()}
    pred_labels = {str(r.portfolio): str(r.predicted_verdict).replace("_PREDICTED", "") for _, r in pred.iterrows()}
    h_true = ["+".join(map(str, sorted(s))) or "BASE" for s in minimal_edges(true_labels)]
    h_pred = ["+".join(map(str, sorted(s))) or "BASE" for s in minimal_edges(pred_labels)]
    return {
        "label": label,
        "n": len(joined),
        "correct": int(joined.verdict_correct.sum()),
        "accuracy": float(joined.verdict_correct.mean()),
        "false_safe": int(joined.false_safe.sum()),
        "false_unstable": int(joined.false_unstable.sum()),
        "h_true": h_true,
        "h_pred": h_pred,
        "h_exact": h_true == h_pred,
        "kappa_true": min((len(x) for x in minimal_edges(true_labels)), default=math.inf),
        "kappa_pred": min((len(x) for x in minimal_edges(pred_labels)), default=math.inf),
        "kappa_exact": min((len(x) for x in minimal_edges(true_labels)), default=math.inf) == min((len(x) for x in minimal_edges(pred_labels)), default=math.inf),
        "root_alpha_error_H4": float(joined.loc[joined.portfolio == "30+33+35+37", "predicted_root_real"].iloc[0] - joined.loc[joined.portfolio == "30+33+35+37", "alpha_truth"].iloc[0]) if label == "V4" else float("nan"),
        "root_frequency_error_H4": float(joined.loc[joined.portfolio == "30+33+35+37", "predicted_root_frequency_hz"].iloc[0] - joined.loc[joined.portfolio == "30+33+35+37", "frequency_truth_hz"].iloc[0]) if label == "V4" else float("nan"),
    }


def runtime_full(label, candidates):
    started = time.perf_counter()
    rows = []
    for index, subset in enumerate(
        [tuple(s) for n in range(len(candidates) + 1) for s in combinations(candidates, n)], start=1
    ):
        t0 = time.perf_counter()
        case = direct_case(subset)
        alpha, freq = full_mode(case)
        rows.append({"portfolio": "+".join(map(str, subset)) or "BASE", "alpha": alpha, "frequency_hz": freq, "full_case_s": time.perf_counter() - t0})
        if index == 1 or index % 64 == 0 or index == 2 ** len(candidates):
            print(f"{label} full evaluation {index}/{2 ** len(candidates)}", flush=True)
    total = time.perf_counter() - started
    pd.DataFrame(rows).to_csv(OUT / f"TX4_{label}_FULL_RUNTIME_CASES.csv", index=False)
    return {"label": label, "n_portfolios": len(rows), "full_total_s": total, "full_setup_s": 0.0, "full_median_per_portfolio_s": float(np.median([r["full_case_s"] for r in rows]))}


def main():
    v4 = compare("V4", OUT / "TX4_V4_BLIND_PREDICTIONS.csv", ROOT / "reports/poster/ias2026/research/results/PCV/PCV02/PCV02_core_portfolios.csv", "subset")
    v9 = compare("V9", OUT / "TX4_V9_BLIND_PREDICTIONS.csv", ROOT / "reports/poster/ias2026/research/results/FINAL_CLOSURE/FC10_census_lattice.csv")
    # Full runtime is measured after reveal under the same one-thread launcher.
    rt4 = runtime_full("V4", V4)
    rt9 = runtime_full("V9", V9)
    reduced = pd.read_csv(OUT / "TX4_BLIND_PREDICTION_RUNTIME.csv")
    rt_rows = []
    for row in (rt4, rt9):
        rr = reduced.loc[reduced.label == row["label"]].iloc[0]
        rt_rows.append({**row, "reduced_total_s": float(rr.total_reduced_s), "reduced_setup_s": float(rr.setup_s), "speedup_full_over_reduced": float(row["full_total_s"] / rr.total_reduced_s)})
    pd.DataFrame(rt_rows).to_csv(OUT / "TX4_RUNTIME_SUMMARY.csv", index=False)
    summary = {"v4": v4, "v9": v9, "runtime": rt_rows}
    (OUT / "TX4_REVEAL_SUMMARY.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
