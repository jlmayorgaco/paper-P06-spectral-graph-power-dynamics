# ruff: noqa: E501
"""H19 master evidence table: does a context-independent node ranking suffice? (prereg H19)
Also collects the F_conf Holm family (T1, T2, T3)."""

from __future__ import annotations

import _hinfra as HI  # noqa: I001

import json

import pandas as pd

import _hstats as ST


def g(name):
    p = HI.RESULTS / name
    return json.loads(p.read_text()) if p.exists() else {}


def run():
    h3, h4, h5, h18, h7, h13 = (g("H03_gate.json"), g("H04_gate.json"), g("H05_stats.json"), g("H18_gate.json"),
                                g("H07_goldb_gate.json"), g("H13_mixing.json"))
    rows = [
        ("L1", "arbitrary-context reversal", h3.get("new_frac_arbitrary_global"), ">=0.75", h3.get("L1_pass")),
        ("L2", "nested-context reversal (level A)", h3.get("new_frac_A"), ">=0.75", h3.get("L2_pass")),
        ("L3", "same-mode nested reversal (level B)", h3.get("new_frac_B"), ">=0.50", h3.get("L3_pass")),
        ("L4", "EM same-mode nested reversal (level C)", h3.get("new_frac_C"), ">=0.75", h3.get("L4_pass")),
        ("L5", "optimal fixed ranking failure (FULL)", f"p*={h4.get('new_FULL_median_pstar'):.3f}; regret={h4.get('new_FULL_median_frac_regret'):.3f}" if h4 else None,
         "p*<=0.90 & regret>=0.10", h4.get("L5_pass")),
        ("L6", "ranking transfer failure", f"offdiag regret={h4.get('transfer_new_median_offdiag_regret_by_test'):.3f}; diag={h4.get('transfer_new_median_diag_regret'):.3f}" if h4 else None,
         ">=0.10 and > diagonal", h4.get("L6_pass")),
        ("L7", "envelope robustness (fresh draws)", json.dumps(h5.get("L7_cov_by_env_global")), ">=0.50 in >=3/4 envelopes", h5.get("L7_pass")),
        ("L8", "dynamic-converter holdout (ALT-WECC reversal)", f"{h18.get('A_frac')} ({h18.get('A_verdict')})" if h18 else "BLOCKED", "TRANSFERS", h18.get("L8_pass", False)),
    ]
    df = pd.DataFrame(rows, columns=["layer", "question", "value_new_holdout", "rule", "positive"])
    df["positive"] = df.positive.fillna(False).astype(bool)
    npos = int(df.positive.sum())
    strong = bool(npos >= 5 and df.set_index("layer").loc["L2", "positive"] and df.set_index("layer").loc["L5", "positive"])
    # EM stratification qualifiers (reported next to the layers; not part of the frozen rule)
    em = {"L5_EM_median_pstar": h4.get("new_EM_median_pstar"), "L5_EM_median_regret": h4.get("new_EM_median_frac_regret"),
          "L5_SAME_median_pstar": h4.get("new_SAME_median_pstar"), "L5_SAME_median_regret": h4.get("new_SAME_median_frac_regret"),
          "L5_EM_would_pass_A3_rule": bool(h4.get("new_EM_median_pstar", 1) <= 0.90 and h4.get("new_EM_median_frac_regret", 0) >= 0.10),
          "L7_EM_same_mode_cov": h5.get("L7_cov_by_env_EM_samemode"), "L7_EM_n_envs_ge_0.5": h5.get("L7_EM_n_envs_ge_0.5")}
    df.to_csv(HI.RESULTS / "H19_evidence.csv", index=False)
    pv = {"T1_goldA_sign": h5.get("T1_sign_test", {}).get("p_one_sided"), "T2_goldB_signflip": h7.get("GB_prereg_T2_signflip_p"),
          "T3_mixing": h13.get("T3_primary", {}).get("p")}
    pv = {k: (float(v) if v is not None else float("nan")) for k, v in pv.items()}
    hol = ST.holm(pv)
    t3 = h13.get("T3_primary", {})
    out = {"n_positive": npos, "strong_claim": strong, "layers": df.to_dict("records"), "em_qualifiers": em,
           "F_conf_raw_p": pv, "F_conf_holm_p": hol,
           "T3_correlates": bool(abs(t3.get("rho", 0)) >= 0.5 and hol.get("T3_mixing", 1) <= 0.01)}
    HI.write_json("H19_evidence.json", out)
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=1, default=str))
