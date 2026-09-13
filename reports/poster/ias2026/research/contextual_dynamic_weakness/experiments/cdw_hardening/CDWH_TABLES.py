# ruff: noqa: E501
"""Generate the manuscript result tables (LaTeX) directly from the result files."""

from __future__ import annotations

import _hinfra as HI  # noqa: I001

import json

import numpy as np
import pandas as pd

R = HI.RESULTS
PAPER = HI.I.REPO / "reports" / "papers" / "cdw_contextual_dynamic_weakness"


def j(n):
    p = R / n
    return json.loads(p.read_text()) if p.exists() else {}


def f(x, d=2):
    return "--" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{d}f}"


def tab_nested():
    pol = pd.read_csv(R / "H03_policy_summary.csv")
    rows = []
    for lab, key in (("any stable-context pair", "arb"), ("A: nested, global $\\aperp$", "A"), ("B: nested, tracked mode", "B"),
                     ("C: nested, EM same mode", "C"), ("D: C + same family", "D")):
        cells = []
        for s in ("discovery", "old", "new"):
            q = pol[(pol.split == s) & pol.base_stable]
            k = int((q.n_rev_arbitrary > 0).sum()) if key == "arb" else int(q[f"{key}_has"].sum())
            cells.append(f"{k}/{len(q)}")
        rows.append(f"{lab} & " + " & ".join(cells) + r" \\")
    q = pol[(pol.split == "new") & pol.base_stable]
    both = int(q.C_both_dirs.sum())
    body = "\n".join(rows)
    return rf"""\begin{{table}}[!t]
\centering
\caption{{Coverage of reversal across the tested base-stable policies. Level C is the preregistered primary gate ($\ge 0.75$ on the new holdout).}}\label{{tab:nested}}
\footnotesize
\begin{{tabular}}{{@{{}}lccc@{{}}}}
\toprule
reversal type & discovery & old holdout & new holdout \\
\midrule
{body}
\bottomrule
\end{{tabular}}\\[2pt]
\parbox{{\columnwidth}}{{\scriptsize Both directions (s$\to$d and d$\to$s) at level C occur in {both}/{len(q)} new-holdout policies. Base-unstable policies (kept, excluded from these denominators): discovery 0, old 6, new 5.}}
\end{{table}}
"""


def tab_ranking():
    rk = pd.read_csv(R / "H04_ranking.csv")
    rk = rk[rk.base_stable]
    lines = []
    for st, lab in (("FULL", "all transitions"), ("EM", "EM $\\to$ EM transitions"), ("SAME", "EM same-mode (level C)")):
        c = []
        for s in ("old", "new"):
            q = rk[(rk.stratum == st) & (rk.split == s)]
            c += [f(q.p_star.median()), f(q.frac_regret.median())]
        lines.append(f"{lab} & " + " & ".join(c) + r" \\")
    body = "\n".join(lines)
    return rf"""\begin{{table}}[!t]
\centering
\caption{{Oracle-optimal fixed node ranking (medians over base-stable policies). A3 bar: $p^\star\le0.90$ and regret rate $\ge0.10$.}}\label{{tab:ranking}}
\footnotesize
\begin{{tabular}}{{@{{}}lcccc@{{}}}}
\toprule
& \multicolumn{{2}}{{c}}{{old holdout}} & \multicolumn{{2}}{{c}}{{new holdout}}\\
stratum & $p^\star$ & regret & $p^\star$ & regret \\
\midrule
{body}
\bottomrule
\end{{tabular}}
\end{{table}}
"""


def tab_goldb():
    m = pd.read_csv(R / "H06_metrics.csv")
    raw = pd.read_parquet(R / "H06_links_raw.parquet")
    r0 = raw[raw.family == "link"].groupby("cond").R0.max()
    elig = set(r0[r0 <= 1e-8].index)
    base = m[(m.target == "H4") & (m.gamma == 1.5) & (m.truth == "FULL")]
    names = [("S1_absP", "$|P_e|$ (S1, selected)"), ("S2_absS", "$|S_e|$"), ("S3_absz", "$|z_e|$"), ("S4_elecdist", "electrical distance"),
             ("S5_reff", "effective resistance"), ("S6_fiedler", "Fiedler edge score"), ("S7_betweenness", "edge betweenness"),
             ("S8_dvdq", "endpoint $dV/dQ$"), ("S9_dgscr", "$\\Delta$gSCR(H4)"), ("Dconv", "base-portfolio eig.\\ sens."),
             ("Dfrozen", "frozen (portfolio)"), ("Dtotal", "total (portfolio)"), ("Dport", "port form (policies only)")]
    lines = []
    for key, lab in names:
        a = base[(base.predictor == key) & base.cond.isin(elig)]
        b = base[base.predictor == key]
        cells = [f(a.spearman.median()), f(b.spearman.median()), f(b.kendall.median()), f(b.concord.median()), f(b.top5.median()), f(b.ndcg5.median())]
        if key == "Dtotal":
            lab = r"\textbf{" + lab + "}"
            cells = [r"\textbf{" + c + "}" for c in cells]
        lines.append(f"{lab} & " + " & ".join(cells) + r" \\")
        if key == "S9_dgscr":
            lines.append(r"\midrule")
    body = "\n".join(lines)
    return rf"""\begin{{table}}[!t]
\centering
\caption{{Branch-reinforcement ranking on the new holdout (target H4, truth: finite re-equilibrated $\times1.5$). Medians over conditions; ``elig.'': the {len(elig & set(base.cond))} preregistration-eligible conditions ($R_0\le10^{{-8}}$); other columns: all {base.cond.nunique()} solved conditions.}}\label{{tab:goldb}}
\scriptsize
\setlength{{\tabcolsep}}{{2.2pt}}
\begin{{tabular}}{{@{{}}lcccccc@{{}}}}
\toprule
predictor & $\rho$ elig. & $\rho$ & $\tau$ & concord. & top-5 & NDCG@5 \\
\midrule
{body}
\bottomrule
\end{{tabular}}
\end{{table}}
"""


def tab_cross():
    g = j("H18_gate.json")
    rows = [("A", "stable-context reversal, core lattice", f"{g.get('A_frac'):.2f} of {g.get('n_tested_policies')} (custom {g.get('A_gfl_frac_same_policies'):.2f})", g.get("A_verdict")),
            ("B", "same-mode reversal", f"{g.get('B_frac'):.2f}", g.get("B_verdict")),
            ("C", "12-line finite ranking, Kendall custom vs ALT", f"median {g.get('C_median_kendall'):.2f}", g.get("C_verdict")),
            ("C2", "ALT total vs ALT finite, minus $|P|$", f"{g.get('C2_median_rho_altFD'):.2f} vs {g.get('C2_median_rho_S1'):.2f}", g.get("C2_verdict")),
            ("D", "corridor top-1 agreement", f"{g.get('D_frac_top1_agree'):.2f}", g.get("D_verdict")),
            ("E", "topology sign agreement", f"{g.get('E_pooled_sign_agree'):.2f} ({g.get('E_n_pairs')} pairs)", g.get("E_verdict"))]
    body = "\n".join(rf"{a} & {b} & {c} & {d.replace('MODEL-SPECIFIC', 'model-specific').replace('TRANSFERS', 'transfers')} \\" for a, b, c, d in rows)
    return rf"""\begin{{table}}[!t]
\centering
\caption{{Cross-model holdout: frozen TX3 WECC library GFL chain (ALT-WECC, ANDES) vs the custom GFL, 8 preregistered policies.}}\label{{tab:cross}}
\footnotesize
\setlength{{\tabcolsep}}{{2.5pt}}
\begin{{tabular}}{{@{{}}cL{{0.37\columnwidth}}L{{0.28\columnwidth}}L{{0.20\columnwidth}}@{{}}}}
\toprule
& test & value & verdict \\
\midrule
{body}
\bottomrule
\end{{tabular}}
\end{{table}}
"""


def tab_layers():
    h = j("H19_evidence.json")
    rows = []
    for L in h.get("layers", []):
        v = L["value_new_holdout"]
        v = f(v) if isinstance(v, float) else str(v).replace("_", "\\_").replace('"', "").replace("{", "").replace("}", "")
        rows.append(rf"{L['layer']} & {L['question']} & {v} & {'yes' if L['positive'] else 'no'} \\")
    body = "\n".join(rows)
    return rf"""\begin{{table}}[!t]
\centering
\caption{{Preregistered evidence layers for ``a context-independent node ranking does not suffice'' (new holdout). Strong claim: $\ge5$ of 8 including L2 and L5.}}\label{{tab:layers}}
\scriptsize
\setlength{{\tabcolsep}}{{2pt}}
\begin{{tabular}}{{@{{}}cL{{0.44\columnwidth}}L{{0.36\columnwidth}}c@{{}}}}
\toprule
& layer & value & pos. \\
\midrule
{body}
\bottomrule
\end{{tabular}}
\end{{table}}
"""


def main():
    PAPER.mkdir(parents=True, exist_ok=True)
    for name, fn in (("tab_nested", tab_nested), ("tab_ranking", tab_ranking), ("tab_goldb", tab_goldb), ("tab_cross", tab_cross), ("tab_layers", tab_layers)):
        (PAPER / f"{name}.tex").write_text(fn(), encoding="utf-8")
        print("wrote", name)


if __name__ == "__main__":
    main()
