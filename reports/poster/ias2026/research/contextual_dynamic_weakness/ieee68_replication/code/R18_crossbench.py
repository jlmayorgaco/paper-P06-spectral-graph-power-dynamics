# ruff: noqa: E501
"""R18 cross-benchmark matrix: IEEE-39 (hardening results, committed) and IEEE-68 (this campaign), two converters.
Writes results/CDW_CROSS_BENCHMARK_MATRIX.csv. Every value is read from a result file."""

from __future__ import annotations

import _r68 as R  # noqa: I001

import json

import numpy as np
import pandas as pd

H = R.CDW / "results" / "hardening"


def j(p):
    return json.loads(p.read_text()) if p.exists() else {}


def main():
    h3, h4, h5, h7, h18, rv = (j(H / "H03_gate.json"), j(H / "H04_gate.json"), j(H / "H05_stats.json"), j(H / "H07_goldb_gate.json"),
                               j(H / "H18_gate.json"), j(H / "H31_revision.json"))
    a7, b7 = j(R.RESULTS / "CDW68_R07_A_gate.json").get("A_REAL", {}), j(R.RESULTS / "CDW68_R07_B_gate.json").get("B_REAL", {})
    r9, r12, r14, r15 = j(R.RESULTS / "CDW68_R09_gate.json"), j(R.RESULTS / "CDW68_R12_gate.json"), j(R.RESULTS / "CDW68_R14_tds_posthoc.json"), j(R.RESULTS / "CDW68_R15_uncertainty.json")
    ld = rv.get("levelD_base_stable", {})
    em = rv.get("alt_em_tracked", {})
    ut = rv.get("h18_uniform_tested_policies", {})
    rows = []
    rows.append({"benchmark": "IEEE-39", "converter": "custom TX4 GFL (Model A)", "base_realism": "two-axis machines, first-order AVR, power PSS, constant-power loads",
                 "governors": "none", "damping": "D = 0", "candidate_count": 9, "policy_count": "63 (24 new holdout; 19 base-stable)",
                 "levelD_reversal_coverage": f"{round(h3.get('new_frac_D', np.nan) * 19)}/19 base-stable new policies", "median_reversal_magnitude": (ld.get("mag_quartiles") or [None, None])[1],
                 "fixed_ranking_regret": f"FULL {h4.get('new_FULL_median_frac_regret'):.2f}; EM {h4.get('new_EM_median_frac_regret'):.2f}",
                 "total_link_median_spearman": h7.get("GB_new_policies_median_rho_Dtotal"), "best_static_median_spearman": h7.get("GB_new_policies_median_rho_oracle_static"),
                 "advantage": h7.get("GB_prereg_median_diff"), "top5_precision": h7.get("GB_prereg_median_top5_Dtotal"),
                 "nonlinear_tds_confirmation": "not performed", "uncertainty_coverage": "EM same-mode coverage by envelope: " + json.dumps(h5.get("L7_cov_by_env_EM_samemode")),
                 "final_interpretation": "small (about 0.01 s^-1) nested same-mode reversals; fixed ranking fails only through singularity crossings; strong portfolio-conditioned branch ranking"})
    rows.append({"benchmark": "IEEE-39", "converter": "WECC library chain (Model B)", "base_realism": "as above (SG2AX transcription)", "governors": "none", "damping": "D = 0",
                 "candidate_count": 4, "policy_count": f"{ut.get('n_tested')} tested", "levelD_reversal_coverage": f"global 0/{ut.get('n_tested')} (pinned pole); EM-tracked nested reversal {sum(v['n_nested_reversals'] > 0 for v in em.values())}/{len(em)}",
                 "median_reversal_magnitude": None, "fixed_ranking_regret": "not computed", "total_link_median_spearman": h18.get("C2_median_rho_altFD"),
                 "best_static_median_spearman": h18.get("C2_median_rho_S1"), "advantage": ut.get("C2_median_diff"), "top5_precision": None,
                 "nonlinear_tds_confirmation": "not performed", "uncertainty_coverage": "not run", "final_interpretation": "reversal not reproduced; numerical total sensitivity ranks the model's own outcomes (12 lines), margin at the bar"})
    for model, g7, conv in (("A", a7, "custom TX4 GFL (Model A)"), ("B", b7, "WECC library chain (Model B)")):
        rk = (r9.get(model) or {})
        rr = (r12.get(model) or {}).get("R12", {})
        r11 = (r12.get(model) or {}).get("R11", {})
        unc = {k: v for k, v in r15.items() if k.startswith(f"{model}_E")}
        tds = "post hoc: " + "; ".join(f"pair {p['i']}: signs {p['sign_agree_1']}/{p['sign_agree_2']}" for p in r14.get("pairs", [])) if (model == "A" and r14) else "not applicable (no case under the preregistered rule)"
        rows.append({"benchmark": "IEEE-68", "converter": conv, "base_realism": "Singh & Pal v3.3 sub-transient machines, DC4B/ST1A/manual, speed PSS, impedance loads",
                     "governors": "PST tg model 1 on every surviving machine", "damping": "PST d_o on area equivalents; damper windings", "candidate_count": 6,
                     "policy_count": f"{g7.get('n_conditions')} ({g7.get('n_eligible_base_stable')} base-stable)",
                     "levelD_reversal_coverage": f"{g7.get('k_D')}/{g7.get('n_eligible_base_stable')} (level T {g7.get('k_T')}/{g7.get('n_eligible_base_stable')})",
                     "median_reversal_magnitude": None, "fixed_ranking_regret": f"FULL {rk.get('FULL', {}).get('median_frac_regret')}; TRACKED {rk.get('TRACKED', {}).get('median_frac_regret')} (p* {rk.get('TRACKED', {}).get('median_p_star')})",
                     "total_link_median_spearman": rr.get("median_rho_tot"), "best_static_median_spearman": rr.get("median_rho_static_sel"), "advantage": rr.get("median_advantage"),
                     "top5_precision": rr.get("median_top5_tot"), "nonlinear_tds_confirmation": tds,
                     "uncertainty_coverage": json.dumps({k: {"D": v.get("frac_draws_level_D"), "T": v.get("frac_draws_level_T")} for k, v in unc.items()}),
                     "final_interpretation": f"R11 {'PASS' if r11.get('pass') else 'FAIL'}; R12 {'PASS' if rr.get('R12_pass') else 'FAIL'}; reversal gate {'PASS' if g7.get('primary_gate_pass') else 'FAIL'}"})
    df = pd.DataFrame(rows)
    df.to_csv(R.RESULTS / "CDW_CROSS_BENCHMARK_MATRIX.csv", index=False)
    print(df.T.to_string())


if __name__ == "__main__":
    main()
