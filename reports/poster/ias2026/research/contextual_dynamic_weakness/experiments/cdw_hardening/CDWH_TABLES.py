# ruff: noqa: E501
"""Generate the manuscript result tables (LaTeX) directly from the result files."""

from __future__ import annotations

import _hinfra as HI  # noqa: I001

import json

import numpy as np
import pandas as pd

R = HI.RESULTS
PAPER = HI.I.REPO / "reports" / "papers" / "cdw_contextual_dynamic_weakness"
BS = "\\"


def j(n):
    p = R / n
    return json.loads(p.read_text()) if p.exists() else {}


def f(x, d=2):
    return "--" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{d}f}"


def table(caption, label, spec, head, rows, size="footnotesize", sep="3pt", foot=""):
    L = [BS + "begin{table}[!t]", BS + "centering", BS + "caption{" + caption + "}" + BS + "label{" + label + "}", BS + size,
         BS + "setlength{" + BS + "tabcolsep}{" + sep + "}", BS + "begin{tabular}{" + spec + "}", BS + "toprule", head + " " + BS + BS, BS + "midrule"]
    L += rows
    L += [BS + "bottomrule", BS + "end{tabular}"]
    if foot:
        L.append(BS + BS + "[2pt]" + BS + "parbox{" + BS + "columnwidth}{" + BS + "scriptsize " + foot + "}")
    L.append(BS + "end{table}")
    return "\n".join(L) + "\n"


def tab_nested():
    pol = pd.read_csv(R / "H03_policy_summary.csv")
    rev = j("H31_revision.json")
    dc = rev.get("direction_counts_new", {})
    g2 = rev.get("gap2_sensitivity", {})
    rows = []
    for lab, key in (("any stable-context pair", "arb"), ("A: nested, global $" + BS + "aperp$", "A"), ("B: nested, tracked mode", "B"),
                     ("C: B, EM-band critical modes", "C"), ("D: C, same mode in both contexts", "D")):
        cells = []
        for s in ("discovery", "old", "new"):
            q = pol[(pol.split == s) & pol.base_stable]
            k = int((q.n_rev_arbitrary > 0).sum()) if key == "arb" else int(q[f"{key}_has"].sum())
            cells.append(f"{k}/{len(q)}")
        rows.append(f"{lab} & " + " & ".join(cells) + " " + BS + BS)
    foot = ("New holdout by direction (s$" + BS + "to$d / d$" + BS + "to$s / both): level C "
            f"{dc['C']['s2d_policies']}/{dc['C']['d2s_policies']}/{dc['C']['both_policies']}, level D "
            f"{dc['D']['s2d_policies']}/{dc['D']['d2s_policies']}/{dc['D']['both_policies']} of 19. "
            f"Excluding pairs whose rightmost-eigenvalue gap is below $" + BS + "tau_{" + BS + "rm mat}$ at any of the four portfolios (post hoc): "
            f"level C {g2.get('coverage_C_gap_ge_tau')}, level D {g2.get('coverage_D_gap_ge_tau')}. Base-unstable policies (kept, excluded): "
            "discovery 0, old 6, new 5.")
    return table("Coverage of reversal across the tested base-stable policies at $" + BS + "tau_{" + BS + "rm mat}=0.01$\\,s$^{-1}$. Level C is the preregistered primary gate ($" + BS + "ge0.75$ on the new holdout).",
                 "tab:nested", "@{}lccc@{}", "reversal type & discovery & old holdout & new holdout", rows, foot=foot)


def tab_ranking():
    rk = pd.read_csv(R / "H04_ranking.csv")
    rk = rk[rk.base_stable]
    rows = []
    for st, lab in (("FULL", "all transitions"), ("EM", "EM $" + BS + "to$ EM transitions"), ("SAME", "tracked EM mode (level C)")):
        c = []
        for s in ("old", "new"):
            q = rk[(rk.stratum == st) & (rk.split == s)]
            c += [f(q.p_star.median()), f(q.frac_regret.median())]
        rows.append(f"{lab} & " + " & ".join(c) + " " + BS + BS)
    rows.append(BS + "midrule")
    mr = pd.read_csv(R / "H31_min_regret_ranking.csv")
    c = []
    for s in ("old", "new"):
        c += ["--", f(mr[(mr.split == s) & (mr.stratum == "FULL")].min_regret_rate.median())]
    rows.append("regret-optimal order, all & " + " & ".join(c) + " " + BS + BS)
    sc = pd.read_csv(R / "H31_screened_ranking.csv")
    c = []
    for s in ("old", "new"):
        c += ["--", f(sc[sc.split == s].regret_screened_on_feasible.median())]
    rows.append("stability-screened order$^" + BS + "dagger$ & " + " & ".join(c) + " " + BS + BS)
    ss = pd.read_csv(R / "H31_static_sequencing.csv")
    c = []
    for s in ("old", "new"):
        c += ["--", f(ss[(ss.split == s) & (ss.stratum == "FULL") & (ss.screen == "gscr")].frac_regret.median())]
    rows.append("max gSCR of $S" + BS + "cup" + BS + "{i" + BS + "}$ (static) & " + " & ".join(c) + " " + BS + BS)
    cap = ("Fixed and screened next-replacement rules (medians over base-stable policies). Rows 1--3: agreement-optimal fixed order; "
           "preregistered insufficiency bar $p^" + BS + "star" + BS + "le0.90$ and regret rate $" + BS + "ge0.10$. Rows 4--6 are post hoc. "
           "$^" + BS + "dagger$First ranked unit whose next portfolio is stable, on contexts with at least one stable option.")
    return table(cap, "tab:ranking", "@{}lcccc@{}", "& " + BS + "multicolumn{2}{c}{old holdout} & " + BS + "multicolumn{2}{c}{new holdout}" + BS + BS + "\nrule & $p^" + BS + "star$ & regret & $p^" + BS + "star$ & regret", rows)


def tab_goldb():
    m = pd.read_csv(R / "H06_metrics.csv")
    raw = pd.read_parquet(R / "H06_links_raw.parquet")
    r0 = raw[raw.family == "link"].groupby("cond").R0.max()
    elig = set(r0[r0 <= 1e-8].index)
    base = m[(m.target == "H4") & (m.gamma == 1.5) & (m.truth == "FULL")]
    names = [("S1_absP", "$|P_e|$ (selected on discovery)"), ("S2_absS", "$|S_e|$"), ("S3_absz", "$|z_e|$"), ("S4_elecdist", "electrical distance"),
             ("S5_reff", "effective resistance"), ("S6_fiedler", "Fiedler edge score"), ("S7_betweenness", "edge betweenness"),
             ("S8_dvdq", "endpoint $dV/dQ$"), ("S9_dgscr", "$" + BS + "Delta$gSCR($V_4$)"), ("Dconv", "base-portfolio eig." + BS + " sens."),
             ("Dfrozen", "frozen (portfolio)"), ("Dtotal", "total (portfolio)"), ("Dport", "port form (verification)")]
    rows = []
    for key, lab in names:
        a = base[(base.predictor == key) & base.cond.isin(elig)]
        b = base[base.predictor == key]
        cells = [f(a.spearman.median()), f(b.spearman.median()), f(b.kendall.median()), f(b.concord.median()), f(b.top5.median()), f(b.ndcg5.median())]
        if key == "Dtotal":
            lab = BS + "textbf{" + lab + "}"
            cells = [BS + "textbf{" + c + "}" for c in cells]
        rows.append(f"{lab} & " + " & ".join(cells) + " " + BS + BS)
        if key == "S9_dgscr":
            rows.append(BS + "midrule")
    lc = pd.read_csv(R / "H31_loco_branch_list.csv")
    a = lc[lc.cond.isin(elig)]
    rows.append(BS + "midrule")
    rows.append("fixed list from other clusters$^" + BS + "ddagger$ & " + " & ".join([f(a.rho_fixed_list.median()), f(lc.rho_fixed_list.median()), "--", "--", f(lc.top5_fixed_list.median()), "--"]) + " " + BS + BS)
    n_el, n_all = len(elig & set(base.cond)), base.cond.nunique()
    cap = ("Branch-reinforcement ranking for $V_4$ on the 24 new-holdout policies and the 40 fresh draws at P4 (truth: finite re-equilibrated $" + BS + "times1.5$). Medians over conditions; "
           f"``elig.'': the {n_el} preregistration-eligible conditions ($R_0" + BS + f"le10^{{-8}}$); other columns: all {n_all} solved conditions. "
           "The port form re-derives the total derivative from re-solved Jacobians (its ``elig.'' entry uses 2 conditions). $^" + BS + "ddagger$Post hoc: "
           "mean finite-effect rank over the conditions of the other clusters.")
    return table(cap, "tab:goldb", "@{}lcccccc@{}", "predictor & $" + BS + "rho$ elig. & $" + BS + "rho$ & $" + BS + "tau$ & concord. & top-5 & NDCG@5", rows, size="scriptsize", sep="2.2pt")


def tab_cross():
    g = j("H18_gate.json")
    rev = j("H31_revision.json")
    em = rev.get("alt_em_tracked", {})
    em_rev = sum(1 for v in em.values() if v["n_nested_reversals"] > 0)
    em_stab = sum(1 for v in em.values() if v["n_stabilizing"] > 0)
    q1 = rev.get("alt_qflag1", {}).get("lattice", {})
    rows = [("A", "stable-context reversal, global $" + BS + "aperp$", f"{g['A_frac']:.2f} of {g['n_tested_policies']} (custom {g['A_gfl_frac_same_policies']:.2f})", g["A_verdict"]),
            ("A$'$", "same, tracked EM mode (post hoc)", f"stabilizing at {em_stab}/7, nested reversal at {em_rev}/7", "partial"),
            ("A$''$", "voltage control on (post hoc)", f"{q1.get('global_rev', 0)}/{q1.get('n', 7)} global, {q1.get('tracked_rev', 0)}/{q1.get('n', 7)} tracked", "none"),
            ("C", "12-line finite ranking, Kendall", f"median {g['C_median_kendall']:.2f}", g["C_verdict"]),
            ("C2", "ALT total vs $|P|$ on ALT outcomes", f"{g['C2_median_rho_altFD']:.2f} vs {g['C2_median_rho_S1']:.2f} ({g['C2_n_below_0.20']} of 7 below bar)", g["C2_verdict"]),
            ("D", "corridor top-1 agreement", f"{g['D_frac_top1_agree']:.2f}", g["D_verdict"]),
            ("E", "topology sign agreement", f"{g['E_pooled_sign_agree']:.2f} ({g['E_n_pairs']} pairs)", g["E_verdict"])]
    body = [f"{a} & {b} & {c} & {d.replace('MODEL-SPECIFIC', 'model-specific').replace('TRANSFERS', 'transfers').replace('PARTIAL', 'partial')} " + BS + BS for a, b, c, d in rows]
    cap = ("Cross-model holdout: frozen WECC library GFL chain (ALT-WECC) vs the custom GFL on the 7 of 8 preregistered policies with a stable ALT base. "
           "In ALT the global $" + BS + "aperp$ of most stable portfolios is pinned at $-0.100$" + BS + ",s$^{-1}$ by a decoupled integrator, which makes test A uninformative.")
    return table(cap, "tab:cross", "@{}cL{0.34" + BS + "columnwidth}L{0.34" + BS + "columnwidth}L{0.17" + BS + "columnwidth}@{}", "& test & value & verdict", body, sep="2.5pt")


def tab_layers():
    h = j("H19_evidence.json")
    rows = []
    h5 = j("H05_stats.json")
    qlab = {"L3": "tracked-mode nested reversal (level B)", "L4": "EM-band tracked nested reversal (level C)"}
    for L in h.get("layers", []):
        v = L["value_new_holdout"]
        v = f(v) if isinstance(v, float) else str(v).replace("_", BS + "_").replace('"', "").replace("{", "").replace("}", "")
        if L["layer"] == "L7" and h5.get("L7_cov_by_env_EM_samemode"):
            v += "; EM mode: " + ", ".join(f"{k}: {x:.1f}" for k, x in h5["L7_cov_by_env_EM_samemode"].items())
        rows.append(f"{L['layer']} & {qlab.get(L['layer'], L['question'])} & {v} & {'yes' if L['positive'] else 'no'} " + BS + BS)
    cap = ("Preregistered evidence layers for ``a context-independent node ranking does not suffice'' (new holdout). Frozen rule for a strong claim: "
           "$" + BS + "ge5$ of 8 including L2 and L5. L1--L4 are nested, L6 is positive whenever L5 is, L5 fails in the EM stratum, and L7 uses the "
           "global $" + BS + "aperp$, so the layers are not independent.")
    return table(cap, "tab:layers", "@{}cL{0.44" + BS + "columnwidth}L{0.36" + BS + "columnwidth}c@{}", "& layer & value & pos.", rows, size="scriptsize", sep="2pt")


def tab_corridors():
    g = j("H10_gate.json")
    if not g:
        return ""
    t = pd.read_csv(R / "H11_corridor_table.csv")
    t = t[(t.budget == "L1_0.50") & t["set"].isin(["new", "fresh_draws"])]
    rows = []
    for us in ("TXother", "TXall", "TXcore", "K2_01", "K3_01=K4_03", "K3_02=K4_02", "K3_12=K4_23", "K4_13"):
        q = t[t.corridor_set == us]
        if not len(q):
            continue
        rows.append(f"{COR_LABEL.get(us, us)} & {int(q['size'].iloc[0])} & {f(q.stab_frac.mean())} & {f(q.top3_frac.mean())} & "
                    f"{f(q.median_eff.median(), 3)} & {f(g.get(us + '|A|new_frac_ge95'))} & {f(g.get(us + '|B|new_frac_ge95'))} " + BS + BS)
    cap = ("Frozen corridors under an equal budget ($" + BS + "sum_{e" + BS + "in C}|" + BS + "Delta" + BS + "gamma_e|=0.5$) on the new holdout (24 policies, 40 draws): "
           "fraction of conditions stabilized by $" + BS + "ge" + BS + "tau_{" + BS + "rm mat}$, fraction in the top 3, median effect, and fraction of conditions "
           "in the top 5\\,\\% of the size-matched arbitrary (A) and connected (B) nulls. Gate: A $" + BS + "ge0.75$.")
    return table(cap, "tab:corridors", "@{}L{0.30" + BS + "columnwidth}cccccc@{}", "corridor & $|C|$ & stab. & top-3 & effect & null A & null B", rows, size="scriptsize", sep="2.2pt")


COR_LABEL = {"TXother": "other 8 transformers", "TXall": "all 12 transformers", "TXcore": "step-up transformers of $V_4$", "K2_01": "cutset 1--2, 3--4, 14--15",
             "K3_01=K4_03": "cutset 17--18, 17--27", "K3_02=K4_02": "cutset 3--4, 9--39", "K3_12=K4_23": "line 14--15", "K4_13": "line 16--19"}


def sec_corridors():
    g = j("H10_gate.json")
    if not g:
        return "% corridor subsection: generated after H10\n"
    rob = g.get("H10_robust_weak_corridors", [])
    t = pd.read_csv(R / "H11_corridor_table.csv")
    t = t[(t.budget == "L1_0.50") & t["set"].isin(["new", "fresh_draws"])]
    best = t.groupby("corridor_set").top3_frac.mean().sort_values(ascending=False)
    top = best.index[0]
    fr = {us: g.get(f"{us}|A|new_frac_ge95") for us in COR_LABEL}
    fr_old = {us: g.get(f"{us}|A|old_frac_ge95") for us in COR_LABEL}
    lead = max((k for k in fr if fr[k] is not None and np.isfinite(fr[k])), key=lambda k: fr[k])
    s = [BS + "subsection{Equal-budget corridors}" + BS + "label{sec:corr}"]
    s.append("A first campaign ranked groups of branches (corridors) by the effect of scaling every member by 1.5, which favors large groups. "
             "The preregistered test gives every frozen corridor the same budget $" + BS + "sum_{e" + BS + "in C}|" + BS + "Delta" + BS + "gamma_e|=0.5$ and compares "
             "its stabilizing effect on $V_4$ with 500 arbitrary branch groups of the same size (family A; connected groups, family B, as a secondary null). A "
             "corridor is called robust if it lies in the top 5" + BS + ",\\% of family A in at least 75" + BS + ",\\% of the 64 new-holdout conditions.")
    if rob:
        s.append(f"The rule is met by {', '.join(COR_LABEL.get(r, r) for r in rob)} (Table~" + BS + "ref{tab:corridors}); old-holdout replication: " +
                 ", ".join(f"{COR_LABEL.get(r, r)} {f(fr_old.get(r))}" for r in rob) + ".")
    else:
        s.append(f"No corridor meets the rule (Table~" + BS + f"ref{{tab:corridors}}). The closest is the {COR_LABEL.get(lead, lead)} group, in the top 5" + BS + ",\\% "
                 f"in {f(fr[lead])} of conditions (old holdout {f(fr_old.get(lead))}). The {COR_LABEL.get(top, top)} group is most often among the three "
                 f"strongest frozen corridors ({f(best.iloc[0])} of conditions), so the earlier corridor ranking partly reflects corridor size and "
                 "we make no weak-corridor claim.")
    s.append("")
    s.append(BS + "input{tab_corridors.tex}")
    return "\n".join(s) + "\n"


def main():
    PAPER.mkdir(parents=True, exist_ok=True)
    for name, fn in (("tab_nested", tab_nested), ("tab_ranking", tab_ranking), ("tab_goldb", tab_goldb), ("tab_cross", tab_cross),
                     ("tab_layers", tab_layers), ("tab_corridors", tab_corridors), ("sec_corridors", sec_corridors)):
        (PAPER / f"{name}.tex").write_text(fn(), encoding="utf-8")
        print("wrote", name)


if __name__ == "__main__":
    main()
