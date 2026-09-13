# ruff: noqa: E501
"""H30: supplement.tex for the CDW paper, generated from the result files."""

from __future__ import annotations

import _hinfra as HI  # noqa: I001

import json

import numpy as np
import pandas as pd

import _infra as I

R = HI.RESULTS
PAPER = HI.I.REPO / "reports" / "papers" / "cdw_contextual_dynamic_weakness"


def j(n):
    p = R / n
    return json.loads(p.read_text()) if p.exists() else {}


def f(x, d=3):
    try:
        x = float(x)
        return "--" if not np.isfinite(x) else f"{x:.{d}f}"
    except (TypeError, ValueError):
        return str(x).replace("_", r"\_")


def longtable(df, spec, caption, label, size=r"\scriptsize"):
    head = " & ".join(c if "$" in c or "\\" in c else c.replace("_", r"\_") for c in df.columns) + r" \\"
    body = "\n".join(" & ".join(str(v) for v in r) + r" \\" for r in df.itertuples(index=False))
    return rf"""{{{size}
\begin{{longtable}}{{{spec}}}
\caption{{{caption}}}\label{{{label}}}\\
\toprule
{head}
\midrule
\endfirsthead
\toprule
{head}
\midrule
\endhead
{body}
\bottomrule
\end{{longtable}}}}
"""


def sec_policies():
    rows = json.loads((HI.INPUTS / "hardening_policies.json").read_text())
    pol = pd.read_csv(R / "H03_policy_summary.csv").set_index("pid")
    df = pd.DataFrame([{"id": r["id"].replace("HARDENING_H", "N"), "$g$": f(r["g"], 4), "$k$": f(r["k"], 3), "$t$": f(r["t"], 3), "$h$": f(r["h"], 3),
                        r"$\aperp(\emptyset)$": f(pol.loc[r["id"], "base_alpha"], 3), "base": "stable" if pol.loc[r["id"], "base_stable"] else "unstable"} for r in rows])
    man = json.loads((HI.INPUTS / "manifest.json").read_text())
    txt = (r"\section{Policy holdout}" + "\n" + "The new holdout HARDENING\\_H01--H24 (N01--N24) is the maximin Latin hypercube of "
           f"{man['seeds']['policy_candidates']} candidates (seed {man['seeds']['policy']}; minimum pairwise distance "
           f"{man['policy_design']['min_pairwise_distance_unit_cube']:.4f} in the unit cube). "
           "The file \\path{results/hardening/prereg_inputs/hardening_policies.json} has sha256 "
           f"\\texttt{{{man['sha256']['hardening_policies.json'][:16]}\\ldots}}.\n\n")
    return txt + longtable(df, "lcccccc", "New holdout policies ($g=u^2$).", "tab:s-pol")


def sec_census():
    pol = pd.read_csv(R / "H03_policy_summary.csv")
    rk = pd.read_csv(R / "H04_ranking.csv")
    rk = rk[rk.stratum == "FULL"].set_index("pid")
    rows = []
    for r in pol.itertuples():
        rows.append({"policy": r.pid.replace("HARDENING_H", "N"), "set": r.split, "stable base": "y" if r.base_stable else "n",
                     "rev.": r.n_rev_arbitrary, "A": r.A_n_units, "B": r.B_n_units, "C": r.C_n_units, "C pairs": r.C_n_pairs,
                     "C s$\\to$d/d$\\to$s": f"{r.C_n_s2d}/{r.C_n_d2s}", "C max mag": f(r.C_max_mag),
                     "$p^\\star$": f(rk.loc[r.pid, "p_star"]), "regret": f(rk.loc[r.pid, "frac_regret"])})
    df = pd.DataFrame(rows)
    txt = (r"\section{Per-policy census summary}" + "\nColumns rev./A/B/C: number of units (of nine) with a stable-context reversal "
           "(any pair) or a nested reversal at levels A/B/C. The full marginal table (2304 marginals per policy, 63 policies) is "
           "\\path{results/hardening/H03_marginals.parquet}; all nested pairs with their curvature witnesses are in "
           "\\path{H03_nested_pairs.parquet}.\n\n")
    return txt + longtable(df, "lllrrrrrrrrr", "Per-policy reversal counts and oracle fixed-ranking results (FULL stratum).", "tab:s-census", r"\tiny")


def sec_witness():
    p = pd.read_parquet(R / "H03_nested_pairs.parquet")
    p = p[(p.level == "C") & p.pid.str.startswith("HARDENING")]
    rows = []
    for pid, g in p.groupby("pid"):
        for d in ("s2d", "d2s"):
            q = g[g.dir == d]
            if not len(q):
                continue
            r = q.loc[q.m.idxmin()] if q.m.min() == 1 else q.loc[q.mag.idxmax()]
            rows.append({"policy": pid.replace("HARDENING_H", "N"), "dir": d.replace("s2d", "s$\\to$d").replace("d2s", "d$\\to$s"), "$i$": r.i,
                         "$S_1$": r.S1.replace("+", ","), "$S_2$": r.S2.replace("+", ","), r"$\Delta(S_1)$": f(r.d1), r"$\Delta(S_2)$": f(r.d2),
                         "$j$": r.wit_j, "$d_{ij}$": f(r.wit_d), "EM-clean": "y" if r.em_clean else "n"})
    df = pd.DataFrame(rows)
    txt = (r"\section{Minimal curvature witnesses}" + "\nFor each new-holdout policy and direction, the level-C nested reversal with the "
           "smallest chain length (or, if none has $m=1$, the largest magnitude), and the certifying second difference "
           "$d_{ij}$ on the canonical chain (Proposition~2 of the paper). s$\\to$d needs $d>0$ (a submodularity violation), "
           "d$\\to$s needs $d<0$ (a supermodularity violation).\n\n")
    return txt + longtable(df, "llcllrrcrc", "Minimal nested EM same-mode reversals and certifying curvature terms (s$^{-1}$).", "tab:s-wit", r"\tiny")


def sec_stats():
    return r"""\section{Statistical methodology}
Policies and envelope draws are designed points (maximin Latin hypercube; declared bounds), not samples of a physical
population. Every fraction is a coverage over the tested conditions. Intervals are percentile bootstraps ($10^4$
resamples, seed 20260932) with clusters: one per policy; for the 40 fresh draws, one per envelope. Three secondary
inferential summaries form the Holm family F$_{\rm conf}$: T1, an exact one-sided sign test of $p^\star<0.90$ over new
base-stable policies; T2, a one-sided cluster sign-flip test of the median paired Spearman advantage minus 0.20;
T3, a two-sided permutation test of Spearman($\mu(\theta,\omega_{\rm ref})$, reversal count). All gates are
effect-size thresholds frozen in \path{docs/CDW_HARDENING_PREREG_V1.md}. Deviations (execution fixes, the ALT-WECC
structural-zero and limiter refinements recorded before the cross-model matrix, and the literal application of the
GOLD-B eligibility rule) are in \path{docs/CDW_HARDENING_DEVIATIONS.md}.
"""


def sec_baselines():
    m = pd.read_csv(R / "H06_metrics.csv")
    out = [r"\section{All branch-ranking baselines}"]
    for (tgt, gam, tr), lab in (((("H4", 1.5, "FULL"), "target $V_4$, $\\gamma=1.5$, truth $\\aperp$")), ((("H4", 1.1, "FULL"), "target $V_4$, $\\gamma=1.10$")),
                                ((("H4", 1.25, "FULL"), "target $V_4$, $\\gamma=1.25$")), ((("H4", 1.5, "SAME"), "target $V_4$, $\\gamma=1.5$, tracked-mode truth")),
                                ((("V9", 1.5, "FULL"), "target $V_9$, $\\gamma=1.5$"))):
        q = m[(m.target == tgt) & (np.isclose(m.gamma, gam)) & (m.truth == tr)]
        g = q.groupby("predictor")[["spearman", "kendall", "concord", "top3", "top5", "ndcg5", "sign_acc"]].median()
        g = g.reindex([p for p in ["S1_absP", "S2_absS", "S3_absz", "S4_elecdist", "S5_reff", "S6_fiedler", "S7_betweenness", "S8_dvdq", "S9_dgscr", "Dconv", "Dfrozen", "Dtotal", "Dport"] if p in g.index])
        df = g.reset_index()
        df["predictor"] = df.predictor.str.replace("_", r"\_")
        df = df.map(lambda v: f(v) if not isinstance(v, str) else v)
        df.columns = ["predictor", r"$\rho$", r"$\tau$", "concord.", "top-3", "top-5", "NDCG@5", "sign acc."]
        out.append(longtable(df, "lrrrrrrr", f"Median ranking metrics over {q.cond.nunique()} conditions: {lab}.", f"tab:s-bl-{tgt}-{gam}-{tr}"))
    return "\n".join(out)


def sec_cross():
    d = pd.read_csv(R / "H18_crossmodel.csv")
    cols = {"pid": "policy", "alt_base_status": "ALT base", "alt_H4_alpha": r"ALT $\aperp(V_4)$", "alt_H4_hz": "Hz", "alt_rev_A": "ALT rev.", "gfl_rev_A": "GFL rev.",
            "C_kendall": r"$\tau$ lines", "C2_rho_altFD": r"ALT $\rho_{\rm tot}$", "C2_rho_S1": r"ALT $\rho_{|P|}$", "D_gfl_top": "GFL top", "D_alt_top": "ALT top", "E_agree": "topo agree"}
    df = d[list(cols)].rename(columns=cols)
    df["policy"] = df.policy.str.replace("HARDENING_H", "N")
    df = df.map(lambda v: f(v) if isinstance(v, float) else str(v).replace("_", r"\_"))
    q = j("alt/H17_qualification_andes.json")
    ref = j("alt/H17_internal_reference.json")
    q1 = max(abs(q["Q1"][p]["alpha"] - ref[p]["alpha"]) for p in ref) if q else float("nan")
    txt = (r"\section{Cross-model holdout detail}" + "\nALT-WECC: TX3 WECC library GFL chain (PLL2, REGCP1, REECB1, REPCA1, BusFreq; frozen "
           "parameter set TX3-GFL-0.1, sha256 \\texttt{f07a6a40\\ldots}) in ANDES 2.0.0 with the TX4 machine transcription. Qualification: "
           f"all-SG base $\\aperp$ reproduced within {q1:.1e}~s$^{{-1}}$ at all eight policies; descriptor vs ANDES EIG spectra "
           f"{q.get('Q3_P4', {}).get('eig_cross_rel_err', float('nan')):.1e}; 488/488 cases initialized. Structural treatment: four "
           "decoupled states per converter (REPCA1 s2 integrator, REECB1 Q-PI and V-PI integrators unused under QFLAG=0) are removed "
           "exactly before the structural pair (deviation log).\n\n")
    return txt + longtable(df, "llrrccrrrllr", "Per-policy cross-model results.", "tab:s-cross", r"\tiny")


def sec_mixing():
    d = pd.read_csv(R / "H13_mixing_decomposition.csv")
    h = j("H13_mixing.json")
    rows = [{"test": k.replace("_", r"\_").replace("|", " / "), r"$\rho$": f(v.get("rho")), "$p$": f(v.get("p")), "Holm $p$": f(v.get("p_holm_Fmix")), "$n$": v.get("n")}
            for k, v in h.get("tests", {}).items()]
    txt = (r"\section{Modal-mixing decomposition}" + f"\nTwo-factor decomposition over {len(d)} policies: controller term variance share "
           f"{h.get('var_share_C', float('nan')):.2f}; median $|C|/(|C|+|F|)$ {h.get('share_abs_C_of_abs_C_plus_abs_F_median', float('nan')):.2f}; "
           f"decomposition check {h.get('max_abs_decomposition_check', 0):.1e}. $\\mu_{{10}}=\\mu(\\theta,\\omega_{{\\rm ref}})$, "
           "$\\mu_{11}=\\mu(\\theta,\\omega_c(\\theta))$ (the earlier campaign's quantity).\n\n")
    return txt + longtable(pd.DataFrame(rows), "lrrrr", "Mixing correlation tests (primary: new / mu10 / n\\_rev).", "tab:s-mix")


def sec_corridors():
    p = R / "H11_corridor_table.csv"
    if not p.exists():
        return r"\section{Equal-budget corridors}" + "\nPending.\n"
    t = pd.read_csv(p)
    t = t[t.budget == "L1_0.50"]
    t["corridor_set"] = t.corridor_set.str.replace("_", r"\_")
    t["set"] = t["set"].str.replace("_", r"\_")
    cols = ["set", "corridor_set", "size", "n_cond", "stab_frac", "top3_frac", "median_eff", "median_rank"]
    df = t[cols].map(lambda v: f(v) if isinstance(v, float) else v)
    df.columns = ["set", "corridor", "size", "$n$", "stab. frac.", "top-3 frac.", "median effect", "median rank"]
    g = j("H10_gate.json")
    txt = (r"\section{Equal-budget corridors and size-matched nulls}" + "\nPrimary budget $\\sum_{e\\in C}|\\Delta\\gamma_e|=0.5$, uniform. Null family A: "
           "500 arbitrary branch groups per size ($k\\ge2$), all 46 branches for $k=1$; family B: connected groups (complete for $k\\le4$). "
           f"Robust weak corridors (family-A percentile $\\ge0.95$ in $\\ge75\\%$ of new-holdout conditions): {', '.join(g.get('H10_robust_weak_corridors', [])) or 'none'}.\n\n")
    return txt + longtable(df, "llrrrrrr", "Equal-budget corridor cross-validation (H11), primary budget.", "tab:s-cor", r"\tiny")


def sec_repro():
    man = pd.read_csv(I.RESULTS / "CDW_HARDENING_RUN_MANIFEST.csv")
    df = man[["phase", "state", "n_tasks", "n_errors", "raw_files", "wall_s"]].map(lambda v: f(v, 0) if isinstance(v, float) else str(v).replace("_", r"\_"))
    det = j("CDWH_DETERMINISM.json")
    txt = (r"\section{Reproducibility}" + "\nEnvironment: Python 3.13.14, numpy 2.5.2, scipy 1.18.1, pandas 3.0.5 (\\texttt{.venv/tx3-analysis}); ANDES 2.0.0 "
           "(\\texttt{.venv/xtool-andes-gfl}). Preregistration commit \\texttt{05b507e3}. Determinism reruns: "
           f"{'identical' if det.get('all_identical') else 'NOT identical'} on every summary field except wall-clock. Orchestrator: "
           "\\path{experiments/cdw_hardening/CDWH_MASTER_RUN.py}; analyses: \\texttt{H03}--\\texttt{H19}; figures, tables and numbers: "
           "\\path{CDWH_FIGURES.py}, \\path{CDWH_TABLES.py}, \\path{CDWH_NUMBERS.py}.\n\n")
    return txt + longtable(df, "llrrrr", "Run manifest.", "tab:s-man")


def sec_figs():
    figs = [("CDWH_S1_corridors", "Equal-budget corridor effects by set."), ("CDWH_S2_null", "Size-matched null percentiles of the frozen corridors."),
            ("CDWH_S3_mixing", "Two-factor modal-mixing decomposition and fixed-frequency correlation."), ("CDWH_S4_crossmodel", "Cross-model holdout (ALT-WECC)."),
            ("CDWH_S5_uncertainty", "Fresh envelope draws: reversal coverage and branch ranking."), ("CDWH_S6_stratification", "FULL / EM / same-mode stratification.")]
    out = [r"\section{Supplementary figures}"]
    for name, cap in figs:
        if (HI.FIGS / f"{name}.pdf").exists():
            out.append(rf"\begin{{figure}}[!h]\centering\includegraphics[width=\textwidth]{{{name}.pdf}}\caption{{{cap}}}\end{{figure}}")
    return "\n".join(out)


def main():
    body = "\n\n".join([sec_policies(), sec_census(), sec_witness(), sec_stats(), sec_baselines(), sec_cross(), sec_corridors(), sec_mixing(), sec_repro(), sec_figs()])
    tex = r"""\documentclass[10pt,a4paper]{article}
\usepackage[margin=2cm]{geometry}
\usepackage[T1]{fontenc}
\usepackage{amsmath,amssymb,booktabs,longtable,graphicx,float}
\usepackage[hyphens]{url}
\usepackage[hidelinks]{hyperref}
\graphicspath{{figures/}}
\newcommand{\aperp}{\alpha_{\perp}}
\renewcommand{\thetable}{S\arabic{table}}
\renewcommand{\thefigure}{S\arabic{figure}}
\title{Supplementary Material:\\ Portfolio-Conditioned Dynamic Weakness and Reinforcement Ranking in Inverter-Rich Power Networks}
\author{Jorge L.~Mayorga Taborda}
\date{}
\begin{document}
\maketitle
""" + body + "\n\\end{document}\n"
    (PAPER / "supplement.tex").write_text(tex, encoding="utf-8")
    print("supplement written")


if __name__ == "__main__":
    main()
