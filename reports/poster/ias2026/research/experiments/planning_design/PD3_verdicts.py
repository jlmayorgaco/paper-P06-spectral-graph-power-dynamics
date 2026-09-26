# ruff: noqa: E402, E501
"""PD3 - frozen verdict rules for PD1 and PD2 (written before any campaign run)."""

from __future__ import annotations

import json
import sys

import _pd_common as P
import numpy as np
import pandas as pd
from scipy.stats import binomtest, spearmanr

import pd2_census_design as PD2

TERCILE_EDGES = (0.12215459, 0.24183127)  # frozen E35 edges of alpha_IA over the 240 unstable samples
TARGET = -0.02
REPRO_TOL = 1e-6
REPRO_MAX_FAIL = 0.05


def tercile(a: float) -> str:
    return "mild" if a <= TERCILE_EDGES[0] else ("medium" if a <= TERCILE_EDGES[1] else "severe")


def pd1() -> dict:
    path = P.RESULTS / "PD1" / "PD1_records.jsonl"
    recs = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    d = pd.DataFrame(recs)
    n = len(d)
    ok = d[d.ok == True].copy()  # noqa: E712
    ok["repro_err"] = (ok.flag_alpha_IA - ok.alpha_IA_e35).abs()
    repro_fail = int((ok.repro_err > REPRO_TOL).sum()) + int((d.ok != True).sum())  # noqa: E712
    ok["tercile"] = ok.alpha_IA_e35.map(tercile)
    for x in ("rc", "sc", "ad"):
        st, al, br = f"{x}_status_perp", f"{x}_alpha_perp", f"{x}_band_rhp"
        ok[f"{x}_S0"] = ok[st] == "STABLE"
        ok[f"{x}_S1"] = ok[f"{x}_S0"] & (ok[al] <= TARGET + 1e-6)
        ok[f"{x}_SB"] = ok[br] == 0
    ok.to_csv(P.out_dir("PD1") / "PD1_table.csv", index=False)
    rates = {}
    for x in ("rc", "sc", "ad"):
        for c in ("S0", "S1", "SB"):
            rates[f"{x}_{c}"] = {"all": [int(ok[f"{x}_{c}"].sum()), int(len(ok))]}
            for t in ("mild", "medium", "severe"):
                g = ok[ok.tercile == t]
                rates[f"{x}_{c}"][t] = [int(g[f"{x}_{c}"].sum()), int(len(g))]
    a = rates["ad_S0"]["all"]
    rate_all = a[0] / a[1] if a[1] else float("nan")
    sev = rates["ad_S0"]["severe"]
    rate_sev = sev[0] / sev[1] if sev[1] else float("nan")
    b = int((ok.ad_S0 & ~ok.rc_S0).sum())
    c = int((~ok.ad_S0 & ok.rc_S0).sum())
    p_mcnemar = float(binomtest(b, b + c, 0.5, alternative="greater").pvalue) if b + c else 1.0
    succ = ok[ok.ad_S0]
    rho = spearmanr(succ.alpha_IA_e35, succ.change_norm).statistic if len(succ) > 3 else float("nan")
    out = {
        "n_records": n, "n_ok": int(len(ok)), "reproduction_failures": repro_fail,
        "reproduction_gate": "PASS" if repro_fail <= REPRO_MAX_FAIL * 240 else "FAIL_STOP",
        "max_repro_err": float(ok.repro_err.max()) if len(ok) else float("nan"),
        "flagship_status_perp": ok.flag_status_perp.value_counts().to_dict(),
        "rates": rates,
        "H_PD1a": {"adaptive_S0_rate": rate_all,
                   "verdict": "SUPPORTED" if rate_all >= 0.95 else ("PARTIAL" if rate_all >= 0.80 else "NOT_SUPPORTED")},
        "H_PD1b": {"adaptive_S0_severe_rate": rate_sev, "verdict": "SUPPORTED" if rate_sev >= 0.80 else "NOT_SUPPORTED"},
        "H_PD1c": {"adaptive_only": b, "rc_only": c, "mcnemar_p_one_sided": p_mcnemar,
                   "verdict": "SUPPORTED" if (p_mcnemar <= 0.01 and b > c) else "NOT_SUPPORTED"},
        "descriptive": {
            "spearman_change_norm_vs_alpha_IA_successes": float(rho),
            "bound_active_fraction": float(ok.bound_active.mean()),
            "bound_active_fraction_failures": float(ok[~ok.ad_S0].bound_active.mean()) if (~ok.ad_S0).any() else float("nan"),
            "slsqp_success_fraction": float(ok.slsqp_success.mean()),
            "median_change_norm_by_tercile": {t: float(ok[(ok.tercile == t) & ok.ad_S0].change_norm.median()) for t in ("mild", "medium", "severe")},
            "adaptive_status_counts": ok.ad_status_perp.value_counts().to_dict(),
            "median_wall_s": float(ok.sample_wall_s.median()),
        },
    }
    return out


def pd2() -> dict:
    recs = []
    for task in PD2.tasks():
        path = PD2.TASK_DIR / f"{PD2.task_name(task)}.json"
        if path.exists():
            recs.append(json.loads(path.read_text(encoding="utf-8")))
    g = PD2.gold_d(recs)
    census = [c for c in g["cases"] if c["target"] in ("T2", "T4")]
    controls = [c for c in g["cases"] if c["target"] in ("T1", "T3")]
    any_plan_safe = any(c["plan_safe"] for c in census)
    any_single_safe = any((c["single_T_stable"] and c["single_phi"] < 0) for c in census)
    if any(c["single_leaves_unsafe_subset"] and c["plan_safe"] for c in census):
        verdict = "SUPPORTED"
    elif any_plan_safe:
        verdict = "NOT_SUPPORTED (plan-level safe design exists but single-boundary tuning did not leave an unsafe subset)"
    else:
        verdict = "INCONCLUSIVE (no census-scale plan-level design reached a safe lattice within the budget)"
    rows = []
    for r in recs:
        t = r.get("task", {})
        rows.append({"target": t.get("target"), "family": t.get("family"), "method": t.get("method"), "ok": r.get("ok"),
                     "outcome": r.get("outcome"), "stop": r.get("stop"), "phi_nominal": r.get("phi_nominal"),
                     "phi_final": r.get("phi_final"), "alpha_T_final": r.get("alpha_T_final"), "status_T_final": r.get("status_T_final"),
                     "n_unstable_nominal": r.get("n_unstable_nominal"), "n_unstable_final": len(r.get("unstable_subsets_final") or []),
                     "n_unresolved_final": r.get("n_unresolved_final"), "lattice_evals": r.get("lattice_evals"),
                     "n_iters": len(r.get("history") or []), "wall_s": r.get("wall_s"), "error": r.get("error")})
    pd.DataFrame(rows).to_csv(P.out_dir("PD2") / "PD2_table.csv", index=False)
    return {"n_tasks": len(recs), "GOLD_D_census": verdict, "GOLD_D_pass_any": g["GOLD_D_pass"],
            "any_census_plan_safe": any_plan_safe, "any_census_single_safe": any_single_safe,
            "census_cases": census, "controls": controls}


def main() -> int:
    which = sys.argv[1:] or ["pd1", "pd2"]
    out = {}
    if "pd1" in which:
        out["PD1"] = pd1()
    if "pd2" in which:
        out["PD2"] = pd2()
    P.write_json(P.RESULTS / "PD_VERDICTS.json", out)
    print(json.dumps(out, indent=1, default=str)[:6000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
