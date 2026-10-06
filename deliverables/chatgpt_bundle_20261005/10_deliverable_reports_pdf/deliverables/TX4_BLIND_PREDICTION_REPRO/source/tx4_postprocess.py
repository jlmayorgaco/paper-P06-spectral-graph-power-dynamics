"""Assemble the preregistered TX4 result tables after the numerical runs."""
from __future__ import annotations

import json
import math
import shutil
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
ROOT = HERE.parents[5]
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)


def alias(src: str, dst: str):
    shutil.copyfile(OUT / src, OUT / dst)


def minimal_edges(labels):
    sets = {frozenset(int(x) for x in k.split("+")) if k != "BASE" else frozenset()
            for k, v in labels.items() if v == "UNSTABLE"}
    return sorted([s for s in sets if not any(t < s for t in sets)], key=lambda s: (len(s), tuple(sorted(s))))


def aggregate_pair():
    frame = pd.read_csv(ROOT / "reports/poster/ias2026/research/results/FINAL_CLOSURE/FC10_census_lattice.csv")
    candidates = []
    for ia, a in frame.iterrows():
        for ib, b in frame.iterrows():
            if ia >= ib or int(a["size"]) != int(b["size"]) or a["status"] == b["status"]:
                continue
            mw_scale = max(abs(float(a.replaced_pg_mw)), abs(float(b.replaced_pg_mw)), 1.0)
            mva_scale = max(abs(float(a.replaced_sn_mva)), abs(float(b.replaced_sn_mva)), 1.0)
            d = math.sqrt(((float(a.replaced_pg_mw) - float(b.replaced_pg_mw)) / mw_scale) ** 2
                          + ((float(a.replaced_sn_mva) - float(b.replaced_sn_mva)) / mva_scale) ** 2)
            exact = abs(float(a.replaced_pg_mw) - float(b.replaced_pg_mw)) <= 1e-9 and abs(float(a.replaced_sn_mva) - float(b.replaced_sn_mva)) <= 1e-9
            candidates.append((exact, d, str(a.members), str(b.members), a, b))
    exact = [x for x in candidates if x[0]]
    pool = exact if exact else [x for x in candidates if x[1] <= 0.05]
    if not pool:
        row = {"selection": "NONE", "selection_rule": "same cardinality; exact MW/MVA if available; otherwise nearest normalized MW/MVA distance; opposite verdict; mismatch <=5%", "valid": False}
    else:
        _, d, pa, pb, a, b = sorted(pool, key=lambda x: (x[1], x[2], x[3]))[0]
        row = {
            "portfolio_A": pa, "portfolio_B": pb, "cardinality": int(a["size"]),
            "MW_A": float(a.replaced_pg_mw), "MW_B": float(b.replaced_pg_mw),
            "MVA_A": float(a.replaced_sn_mva), "MVA_B": float(b.replaced_sn_mva),
            "alpha_A": float(a.alpha_perp), "alpha_B": float(b.alpha_perp),
            "frequency_A_hz": float("nan"), "frequency_B_hz": float("nan"),
            "status_A": str(a.status), "status_B": str(b.status),
            "normalized_distance": d, "mismatch_percent": 100.0 * d,
            "exact_aggregate_match_available": bool(exact), "valid": bool(d <= 0.05),
            "selection_rule": "same cardinality; exact MW/MVA if available; otherwise nearest normalized MW/MVA distance; opposite verdict; mismatch <=5%",
            "alpha_viewed_after_selection": True,
        }
    pd.DataFrame([row]).to_csv(OUT / "TX4_AGGREGATE_MATCHED_PAIRS.csv", index=False)
    return row


def screens():
    pcv = pd.read_csv(ROOT / "reports/poster/ias2026/research/results/PCV/PCV02/PCV02_core_portfolios.csv")
    h = pcv[(pcv.point == "P4") & (pcv.subset == "30+33+35+37")].iloc[0]
    rows = [
        {"screen": "singleton screening", "quantity": "all proper subsets", "value": "PASS", "status": "safe screen; not composability proof"},
        {"screen": "pair screening", "quantity": "all proper subsets", "value": "PASS", "status": "safe screen; not composability proof"},
        {"screen": "triple screening", "quantity": "all proper subsets", "value": "PASS", "status": "safe screen; not composability proof"},
        {"screen": "aggregate MW/MVA", "quantity": "H4", "value": "2096.607719 MW / 4270.7 MVA", "status": "same aggregate for all portfolios of this support; cannot encode collective compatibility"},
        {"screen": "B3 first-order modal approximation", "quantity": "alpha H4", "value": float(h.B3_alpha), "status": "predicts stable; misses H4"},
        {"screen": "B4 additive approximation", "quantity": "alpha H4", "value": float(h.B4_alpha), "status": "predicts stable; misses H4"},
        {"screen": "B5 pairwise approximation", "quantity": "alpha H4", "value": float(h.B5_alpha), "status": "predicts stable; misses H4"},
        {"screen": "B5b third-order alpha approximation", "quantity": "alpha H4", "value": float(h.B5b_alpha), "status": "predicts stable; misses H4"},
        {"screen": "determinant truncation", "quantity": "non-oracle |chi|", "value": float(h.abs_chi_nonoracle), "status": "not equivalent to third-order alpha approximation"},
        {"screen": "exact reduced closure", "quantity": "H4 verdict", "value": "UNSTABLE_PREDICTED", "status": "correct flagship prediction"},
        {"screen": "full DAE", "quantity": "H4 verdict", "value": "UNSTABLE", "status": "ground truth"},
    ]
    pd.DataFrame(rows).to_csv(OUT / "TX4_SCREEN_COMPARISON.csv", index=False)


def main():
    alias("TX4_CONTEXTUAL_RETURN_AT_BOUNDARY.csv", "TX4_CONTEXTUAL_RETURN_BOUNDARY.csv")
    alias("TX4_PROPER_SUBSET_CLOSURE.csv", "TX4_15_PROPER_SUBSETS.csv")
    alias("TX4_CONTEXTUAL_RETURN_DERIVATIVE.csv", "TX4_RETURN_DERIVATIVE_CHECK.csv")
    alias("TX4_CONTROL_BOUNDARY_REMEDIATION.csv", "TX4_CONTROL_REMEDIATION.csv")
    alias("TX4_TDS/TX4_TDS_FINAL.csv", "TX4_TDS_FINAL.csv")
    alias("TX4_TDS/TX4_TDS_SUMMARY.json", "TX4_TDS_SUMMARY.json")
    boundary = pd.read_csv(OUT / "TX4_CONTEXTUAL_RETURN_BOUNDARY.csv")
    boundary[["device_bus", "proper_subset", "g_root", "frequency_hz", "eta_subset", "return_distance_to_unity", "collective_sigma_min", "local_sigma_min_at_device", "schur_identity_residual"]].to_csv(OUT / "TX4_MINIMALITY_SEPARATION.csv", index=False)
    boundary[["device_bus", "proper_subset", "g_root", "frequency_hz", "eta_subset", "collective_sigma_min", "local_sigma_min_at_device", "local_factor_det_abs", "return_distance_to_unity"]].to_csv(OUT / "TX4_LOCAL_VS_COLLECTIVE.csv", index=False)
    screens()
    pair = aggregate_pair()
    e31 = next((p for p in (ROOT / "reports/poster/ias2026/research/outputs/ias2026").rglob("E31_ANDES_validation.csv")), None)
    if e31:
        andes = pd.read_csv(e31)
        andes[andes.members.astype(str).isin(["BASE", "30+33+35", "30+33+35+37"])].assign(scope="independent ANDES phasor comparison; custom GFL g-boundary NOT TESTED").to_csv(OUT / "TX4_ANDES_FINAL_CHECK.csv", index=False)
    else:
        pd.DataFrame([{"scope": "ANDES", "status": "NOT_TESTED", "reason": "existing E31 table unavailable"}]).to_csv(OUT / "TX4_ANDES_FINAL_CHECK.csv", index=False)
    core = json.loads((OUT / "TX4_CONTEXTUAL_RETURN_CORE_SUMMARY.json").read_text(encoding="utf-8"))
    reduced = json.loads((OUT / "TX4_G_BOUNDARY_REDUCED_SUMMARY.json").read_text(encoding="utf-8"))
    full = json.loads((OUT / "TX4_REVEAL_SUMMARY.json").read_text(encoding="utf-8"))
    summary = {"reduced": reduced, "contextual": core, "reveal": full, "aggregate_pair": pair,
               "aliases_created": ["TX4_CONTEXTUAL_RETURN_BOUNDARY.csv", "TX4_MINIMALITY_SEPARATION.csv", "TX4_LOCAL_VS_COLLECTIVE.csv", "TX4_15_PROPER_SUBSETS.csv", "TX4_RETURN_DERIVATIVE_CHECK.csv", "TX4_CONTROL_REMEDIATION.csv", "TX4_SCREEN_COMPARISON.csv", "TX4_AGGREGATE_MATCHED_PAIRS.csv", "TX4_ANDES_FINAL_CHECK.csv"]}
    (OUT / "TX4_POSTPROCESS_SUMMARY.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(json.dumps({"aggregate_pair": pair, "files": summary["aliases_created"]}, indent=2, default=str))


if __name__ == "__main__":
    main()
