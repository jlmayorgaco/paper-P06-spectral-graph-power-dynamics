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
           "\\path{results/hardening/H03_marginals.parquet}; all nested pairs with their located interaction terms are in "
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
    txt = (r"\section{Located interaction terms}" + "\nFor each new-holdout policy and direction, the level-C nested reversal with the "
           "smallest chain length (or, if none has $m=1$, the largest magnitude), and the second difference "
           "$d_{ij}$ located on the canonical chain by the chain identity (Proposition~1 of the paper). The bound holds for every "
           "nested pair by construction; the table shows where the interaction sits. An s$\\to$d reversal locates $d>0$, "
           "a d$\\to$s reversal $d<0$.\n\n")
    return txt + longtable(df, "llcllrrcrc", "Minimal level-C nested reversals and located interaction terms (s$^{-1}$).", "tab:s-wit", r"\tiny")


def sec_stats():
    return r"""\section{Statistical methodology}
Policies and envelope draws are designed points (maximin Latin hypercube; declared bounds), not samples of a physical
population. Every fraction is a coverage over the tested conditions. Intervals are percentile bootstraps ($10^4$
resamples, seed 20260932) with clusters: one per policy; for the 40 fresh draws, one per envelope. Three inferential
summaries form the Holm family F$_{\rm conf}$ (called confirmatory in the plan; all three are reported in the paper): T1, an exact one-sided sign test of $p^\star<0.90$ over new
base-stable policies; T2, a one-sided cluster sign-flip test of the median paired Spearman advantage minus 0.20;
T3, a two-sided permutation test of Spearman($\mu(\theta,\omega_{\rm ref})$, reversal count). All gates are
effect-size thresholds frozen in \path{docs/CDW_HARDENING_PREREG_V1.md}. Deviations are in
\path{docs/CDW_HARDENING_DEVIATIONS.md}:
execution fixes, including a condition-key collision in the corridor analysis that was found by its built-in consistency check;
the ALT-WECC structural-zero and limiter refinements;
the literal application of the GOLD-B eligibility rule;
a correction of hand-written clock times;
the list of post-hoc analyses.
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
           f"{q.get('Q3_P4', {}).get('eig_cross_rel_err', float('nan')):.1e}; 488/488 cases initialized. Structural treatment: three "
           "decoupled states per converter (REPCA1 s2 integrator with $K_p=K_i=0$; REECB1 Q-PI and V-PI integrators, unused under QFLAG=0; "
           "the V-PI state is isolated with diagonal $1.2\\times10^{-6}$~s$^{-1}$) are removed exactly before the structural pair "
           "(deviation log, refinements 1 and 1b). The isolated anti-windup tracking state of the disabled REPCA1 plant active-power "
           "integrator (s5\\_xi, $-0.100$~s$^{-1}$) exceeds the $10^{-3}$ removal bound and is kept; it sets the global $\\aperp$ of most "
           "stable ALT portfolios at H02, H07 and N04, which makes the preregistered global-$\\aperp$ reversal test uninformative. "
           "Exploratory EM-tracked analysis (critical EM mode of $S$ followed by MAC into $S\\cup\\{i\\}$) and a post-hoc QFLAG=1 "
           "variant (library voltage control on, gains unchanged, outside the validated freeze) are in Table~\\ref{tab:s-altem}.\n\n")
    rv = j("H31_revision.json")
    em, q1f = rv.get("alt_em_tracked", {}), rv.get("alt_qflag1", {}).get("em_tracked", {})
    rows = []
    for p, v in em.items():
        w = q1f.get(p, {})
        rows.append({"policy": p.replace("HARDENING_H", "N"), "stable": v["n_stable_portfolios"], "pinned": f(v["frac_stable_alpha_pinned_at_minus_0.1"], 2),
                     "EM marg.": v["n_em_tracked_marginals"], "stab.": v["n_stabilizing"], "destab.": v["n_destabilizing"], "nested rev.": v["n_nested_reversals"],
                     "Q1 stab.": w.get("n_stabilizing", "--"), "Q1 rev.": w.get("n_nested_reversals", "--")})
    t2 = longtable(pd.DataFrame(rows), "lrrrrrrrr", "ALT-WECC on the $V_4$ lattice, exploratory EM-tracked analysis (post hoc). stable: stable portfolios; "
                   "pinned: fraction of them whose $\\aperp$ is the isolated $-0.100$~s$^{-1}$ pole; EM marg.: EM-tracked marginals; Q1: the QFLAG=1 variant.",
                   "tab:s-altem")
    return txt + longtable(df, "llrrccrrrllr", "Per-policy cross-model results.", "tab:s-cross", r"\tiny") + "\n" + t2


def sec_revision():
    rv, ex = j("H31_revision.json"), j("H31_explore.json")
    if not rv:
        return ""
    lc, ld = rv["levelC_base_stable"], rv["levelD_base_stable"]
    out = [r"\section{Analyses requested by internal review (post hoc)}",
           "These analyses were added after the preregistered results had been computed and read (deviation log). None changes a preregistered verdict.\n",
           r"\subsection{Reversal magnitudes and levels}",
           f"Base-stable new-holdout policies. Level C: {lc['pairs']} nested pairs (s$\\to$d {lc['s2d']}, d$\\to$s {lc['d2s']}), "
           f"{lc['distinct_policy_unit']} distinct policy--unit combinations, {lc['em_clean']:.0%} with every critical mode in the EM band, "
           f"{lc['m1_pairs']} one-step pairs; largest magnitude {lc['max_mag']:.3f}~s$^{{-1}}$. Level D: {ld['pairs']} pairs (s$\\to$d {ld['s2d']}, "
           f"d$\\to$s {ld['d2s']}) in {ld['policies']} policies, both directions in {ld['policies_both_dirs']}; magnitude quartiles "
           f"{', '.join(f'{x:.4f}' for x in ld['mag_quartiles'])}~s$^{{-1}}$. Pairs whose rightmost eigenvalue lies within $\\tau_{{\\rm mat}}$ of the next at one "
           f"of the four portfolios: {rv['gap2_sensitivity']['frac_pairs_gap_below_tau']:.0%}; coverage without them: level C "
           f"{rv['gap2_sensitivity']['coverage_C_gap_ge_tau']}, level D {rv['gap2_sensitivity']['coverage_D_gap_ge_tau']}. Continuation of the tracked "
           f"mode under fractional replacement ends on the MAC-matched mode in {ex.get('continuation_match_frac', float('nan')):.0%} of "
           f"{ex.get('continuation_n', 0)} endpoint marginals.\n"]
    tc = pd.read_csv(R / "H31_tau_curve.csv")
    out.append(longtable(tc.astype({c: int for c in ["n", "C", "C_both", "D", "D_both"]}).map(lambda v: f(v, 4) if isinstance(v, float) else v),
                         "rlrrrrr", "Number of base-stable policies with a nested reversal against the materiality threshold.", "tab:s-tau"))
    out.append(r"\subsection{Singularity-type events behind fixed-ranking regret}")
    fm = rv["fast_modes_new"]
    out.append(f"On the new holdout, {fm['n_fast']} of {fm['n_portfolios']} portfolios have $\\aperp\\ge1$~s$^{{-1}}$; {fm['frac_fast_real']:.1%} of these "
               f"critical eigenvalues are real, with 5/50/95\\,\\% quantiles {', '.join(f'{x:.0f}' for x in fm['fast_alpha_quantiles'])}~s$^{{-1}}$. "
               "Converted-unit counts: " + ", ".join(f"{k}: {v}" for k, v in fm["fast_size_counts"].items()) + ". Mean participation of the leading "
               "state groups: " + ", ".join(f"{k.replace('_', chr(92) + '_')} {v:.2f}" for k, v in list(ex["fast_participation_mean_top"].items())[:5]) + ".\n")
    sw = pd.DataFrame(ex["sweeps_recomputed"])
    sw["pid"] = sw.pid.str.replace("HARDENING_H", "N")
    sw = sw[["pid", "i", "rho_onset_alpha_ge_1", "alpha_at_onset", "rho_at_min_sigma_gz", "min_sigma_gz"]]
    sw.columns = ["policy", "unit", r"$\rho$ at onset", r"$\aperp$ at onset", r"$\rho$ at min $\sigma(g_z)$", r"min $\sigma(g_z)$"]
    sw = sw.map(lambda v: (f"{v:.0f}" if abs(v) > 100 else f(v, 4)) if isinstance(v, float) else v)
    out.append(longtable(sw, "lrrrrr", "Fractional-replacement sweeps (replaced fraction $\\rho$ of the added unit, step 0.025): onset of "
                         "$\\aperp\\ge1$~s$^{-1}$ and the minimum singular value of the network Jacobian $g_z$ along the sweep.", "tab:s-sweep"))
    sc, mr, ss = rv["screened_ranking_median"], rv["min_regret_ranking_median"], rv["static_sequencing_median_regret"]
    rows = []
    for sp in ("discovery", "old", "new"):
        rows.append({"split": sp, "no stable option": f(sc[sp]["frac_no_stable_option"], 3), "plain (feasible)": f(sc[sp]["regret_plain_on_feasible"], 3),
                     "screened (feasible)": f(sc[sp]["regret_screened_on_feasible"], 3), "regret-opt. FULL": f(mr[f"{sp}|FULL"], 3),
                     "regret-opt. EM": f(mr[f"{sp}|EM"], 3), "gSCR FULL": f(ss[f"{sp}|FULL|gscr"], 3), "min-SCR FULL": f(ss[f"{sp}|FULL|minscr"], 3)})
    out.append(longtable(pd.DataFrame(rows), "lrrrrrrr", "Median material-regret rates of alternative fixed and static sequencing rules.", "tab:s-seq", r"\tiny"))
    out.append(r"\subsection{Branch ranking against a learned fixed list and large actions}")
    lo = rv["loco_fixed_branch_list"]
    rows = [{"set": k.replace("_", " "), r"$\rho$ list": f(v["rho_fixed_list"]), r"$\rho$ $D^{\rm tot}$": f(v["rho_Dtotal"]), r"$\rho$ $D^{\rm fro}$": f(v["rho_Dfrozen"]),
             "top-5 list": f(v["top5_fixed_list"], 2), "top-5 $D^{\\rm tot}$": f(v["top5_Dtotal"], 2), "$D^{\\rm tot}$ beats list": f(v["frac_Dtotal_beats_fixed_list"], 2)}
            for k, v in lo.items()]
    out.append(longtable(pd.DataFrame(rows), "lrrrrrr", "Leave-one-cluster-out fixed branch list (mean finite-effect rank over the other clusters), medians.", "tab:s-loco"))
    dt = rv["dtot_vs_topology_median"]
    tm = ex.get("timing", {})
    out.append(f"$D^{{\\rm tot}}$ against finite branch actions at six new-holdout policies: doubling, median Spearman {dt['dbl']['rho_full']:.3f} "
               f"(top-5 {dt['dbl']['top5']:.2f}); outage, {dt['out']['rho_full']:.3f} (EM-tracked {dt['out']['rho_em']:.3f}, top-5 {dt['out']['top5']:.2f}). "
               f"Timing of this implementation: {tm.get('sec_per_branch_total_derivative_fd_impl', float('nan')):.2f}~s per branch for the total derivative "
               f"(finite-difference Jacobians) and {tm.get('sec_per_branch_finite_resolve_eig', float('nan')):.2f}~s for a finite re-solve with eigenanalysis.\n")
    tw = rv["topology_within_type"]
    rows = [{"action": k.split("|")[0], "score": k.split("|")[1].replace("_", r"\_"), r"median $\rho$": f(v["median_rho"]), r"frac. $|\rho|\ge0.6$": f(v["frac_abs_ge_0.6"], 2)}
            for k, v in tw.items()]
    out.append(r"\subsection{Topology provenance and within-type static scores}")
    cs = rv["topology_creation_by_split"]
    out.append(f"EM-incompatibility creations by policy set: discovery {cs['discovery']}, old holdout {cs['old']}, new holdout {cs['new']}; all are outages.\n")
    out.append(longtable(pd.DataFrame(rows), "llrr", "Static topology scores against the tracked EM effect, within action type (16 policies).", "tab:s-topo"))
    mv = rv["mixing_variance_terms"]
    out.append(r"\subsection{Mixing variance terms}" + f"\nWith $\\Delta\\mu=C+F$ (controller and frequency terms), "
               f"${{\\rm var}}(C)/{{\\rm var}}(\\Delta\\mu)={mv['var_C_over_var_total']:.2f}$, ${{\\rm var}}(F)/{{\\rm var}}(\\Delta\\mu)={mv['var_F_over_var_total']:.2f}$ "
               f"and $2\\,{{\\rm cov}}(C,F)/{{\\rm var}}(\\Delta\\mu)={mv['2cov_over_var_total']:.2f}$.\n")
    return "\n".join(out)


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
    figs = [("CDWH_S1_corridors", r"Equal-budget corridor effects on $V_4$ by condition set (identical cutsets K3\_01=K4\_03, K3\_02=K4\_02, K3\_12=K4\_23 shown once; band: $\pm\tau_{\rm mat}$)."), ("CDWH_S2_null", "Size-matched null percentiles of the frozen corridors (dashed: 0.95)."),
            ("CDWH_S3_mixing", "Two-factor modal-mixing decomposition and fixed-frequency correlation."), ("CDWH_S4_crossmodel", "Cross-model holdout (ALT-WECC)."),
            ("CDWH_S5_uncertainty", "Fresh envelope draws: reversal coverage and branch ranking."), ("CDWH_S6_stratification", "FULL / EM / tracked-mode stratification."),
            ("CDWH_S7_singularity", "Fractional-replacement sweeps: $\\aperp$ and the minimum singular value of $g_z$ against the replaced fraction.")]
    out = [r"\section{Supplementary figures}"]
    for name, cap in figs:
        if (HI.FIGS / f"{name}.pdf").exists():
            out.append(rf"\begin{{figure}}[!h]\centering\includegraphics[width=\textwidth]{{{name}.pdf}}\caption{{{cap}}}\end{{figure}}")
    return "\n".join(out)


def main():
    body = "\n\n".join([sec_policies(), sec_census(), sec_witness(), sec_stats(), sec_baselines(), sec_cross(), sec_corridors(), sec_mixing(), sec_revision(), sec_repro(), sec_figs()])
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
