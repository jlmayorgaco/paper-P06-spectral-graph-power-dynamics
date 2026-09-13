# ruff: noqa: E501
"""Publication figures for the CDW hardening campaign and paper (prereg figure definitions).

Palette (dataviz reference instance, light surface): categorical slots blue #2a78d6, orange
#eb6834, aqua #1baf7a (all-pairs validated); static baselines in neutral grays; sequential blue
ramp; diverging blue<->red with a gray midpoint. Every figure: PDF + SVG + PNG (300 dpi)."""

from __future__ import annotations

import _hinfra as HI  # noqa: I001

import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402

import _cdw as C  # noqa: E402

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
GRAY1, GRAY2 = "#9a9893", "#5f5d58"
SEQ = LinearSegmentedColormap.from_list("seqblue", ["#f4f8fd", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
DIV = LinearSegmentedColormap.from_list("divbr", ["#1c5cab", "#86b6ef", "#f0efec", "#f19a9a", "#c0392f"])
R = HI.RESULTS
OUT = HI.FIGS
plt.rcParams.update({"font.family": "serif", "font.size": 8, "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "pdf.fonttype": 42,
                     "svg.fonttype": "none", "axes.titlesize": 8, "axes.titleweight": "bold"})
COL1, COL2 = 3.45, 7.16  # IEEE column / double-column widths (in)


def save(fig, name):
    for ext in ("pdf", "svg", "png"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


def jload(n):
    return json.loads((R / n).read_text())


# ----------------------------------------------------------------------- Fig 2 --
def fig2():
    mg = pd.read_parquet(R / "H03_marginals.parquet", columns=["pid", "i", "cls", "lvA"])
    pol = pd.read_csv(R / "H03_policy_summary.csv")
    pairs = pd.read_parquet(R / "H03_nested_pairs.parquet", columns=["pid", "level", "i"])
    stab = mg[mg.lvA]
    frac = stab.assign(d=(stab.cls == 1)).groupby(["pid", "i"]).d.mean().unstack()
    order = [p for p in HI.DISC + HI.OLD_HOLD + HI.HPOL if p in frac.index]
    frac = frac.reindex(order)
    cset = set(map(tuple, pairs[pairs.level == "C"][["pid", "i"]].drop_duplicates().to_numpy()))
    fig, ax = plt.subplots(figsize=(COL1, 6.4))
    im = ax.imshow(frac.to_numpy(float), aspect="auto", cmap=SEQ, vmin=0, vmax=1)
    for r, p in enumerate(order):
        for c, u in enumerate(frac.columns):
            if (p, u) in cset:
                ax.plot(c, r, marker="x", ms=3.2, mew=0.8, color="#f28a3a" if frac.iloc[r, c] > 0.55 else INK)
    unst = set(pol[~pol.base_stable].pid)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([(p.replace("HARDENING_H", "N") + (" u" if p in unst else "")) for p in order], fontsize=5)
    ax.set_xticks(range(len(frac.columns)))
    ax.set_xticklabels(frac.columns)
    ax.set_xlabel("replaced unit i (bus)")
    ax.grid(False)
    for y in (len(HI.DISC) - 0.5, len(HI.DISC) + len(HI.OLD_HOLD) - 0.5):
        ax.axhline(y, color=INK, lw=0.8)
    g = jload("H03_gate.json")
    ax.set_title(f"Nested EM same-mode reversal coverage: new {g['new_frac_C']:.2f}, old {g['old_frac_C']:.2f}", fontsize=7)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("fraction of stable contexts where replacing i destabilizes", fontsize=6)
    ax.text(1.02, 1.0, "D: discovery\nH: old holdout\nN: new holdout\nu: unstable base\nx: nested EM\n   same-mode reversal",
            transform=ax.transAxes, fontsize=5, va="top", color=INK2)
    save(fig, "CDWH_F2_contextual_reversal")


# ----------------------------------------------------------------------- Fig 3 --
def strongest_new_C(k=3):
    p = pd.read_parquet(R / "H03_nested_pairs.parquet")
    pol = pd.read_csv(R / "H03_policy_summary.csv")
    ok = set(pol[pol.base_stable].pid)
    p = p[(p.level == "C") & p.pid.str.startswith("HARDENING") & p.pid.isin(ok)].sort_values("mag", ascending=False)
    out, seen = [], set()
    for r in p.itertuples():
        if r.i in seen:
            continue
        seen.add(r.i)
        out.append(r)
        if len(out) == k:
            break
    return out


def fig3():
    mg = pd.read_parquet(R / "H03_marginals.parquet")
    ex = strongest_new_C(3)
    fig, axs = plt.subplots(1, 3, figsize=(COL2, 2.2), sharey=False)
    for ax, r in zip(axs, ex, strict=False):
        g = mg[(mg.pid == r.pid) & (mg.i == r.i) & mg.lvC]
        size = g["mask"].map(lambda m: bin(m).count("1"))
        jit = (np.random.default_rng(1).random(len(g)) - 0.5) * 0.3
        ax.axhspan(-C.TAU_MAT, C.TAU_MAT, color="#f0efec", zorder=0)
        ax.scatter(size + jit, g.delta, s=6, color=MUTED, alpha=0.6, lw=0, label="EM same-mode contexts")
        for lab, d, col, mk in ((r.S1, r.d1, BLUE if r.d1 < 0 else ORANGE, "v" if r.d1 < 0 else "^"),
                                (r.S2, r.d2, BLUE if r.d2 < 0 else ORANGE, "v" if r.d2 < 0 else "^")):
            sz = 0 if lab == "BASE" else lab.count("+") + 1
            ax.scatter([sz], [d], s=40, color=col, marker=mk, edgecolor="white", lw=0.8, zorder=3)
            ax.annotate(("{" + lab.replace("+", ",") + "}") if lab != "BASE" else r"$\emptyset$", (sz, d), xytext=(4, 3), textcoords="offset points", fontsize=6, color=INK)
        ax.axhline(0, color=INK2, lw=0.6)
        ax.set_title(f"unit {r.i}, {r.pid.replace('HARDENING_H', 'N')} ({r.dir.replace('s2d', 'stab→destab').replace('d2s', 'destab→stab')})", fontsize=7)
        ax.set_xlabel("context size |S|")
        ax.set_xticks(range(0, 9))
    axs[0].set_ylabel(r"$\Delta_i\alpha(S)$ on the tracked EM mode (s$^{-1}$)")
    save(fig, "CDWH_F3_same_mode_examples")


# ----------------------------------------------------------------------- Fig 4 --
def fig4():
    rk = pd.read_csv(R / "H04_ranking.csv")
    rk = rk[rk.base_stable]
    reg = pd.read_csv(R / "H04_transfer_regret.csv", index_col=0)
    fig = plt.figure(figsize=(COL2, 2.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.25])
    for k, (col, lab) in enumerate((("p_star", "optimal fixed-ranking agreement p*"), ("frac_regret", "material regret rate"))):
        ax = fig.add_subplot(gs[0, k])
        pos = 0
        ticks, tl = [], []
        for s, sc in (("FULL", BLUE), ("EM", AQUA), ("SAME", ORANGE)):
            for split in ("old", "new"):
                v = rk[(rk.stratum == s) & (rk.split == split)][col].to_numpy(float)
                ax.boxplot([v], positions=[pos], widths=0.6, showfliers=False, patch_artist=True,
                           boxprops={"facecolor": "white", "edgecolor": sc}, medianprops={"color": sc, "lw": 1.5},
                           whiskerprops={"color": sc}, capprops={"color": sc})
                ax.scatter(pos + (np.random.default_rng(2).random(v.size) - 0.5) * 0.3, v, s=5, color=sc, lw=0, alpha=0.7)
                ticks.append(pos)
                tl.append(f"{s}\n{split}")
                pos += 1
            pos += 0.5
        if col == "p_star":
            ax.axhline(0.90, color=INK2, ls="--", lw=0.7)
            ax.text(pos - 0.6, 0.905, "0.90 ceiling", fontsize=5.5, ha="right", color=INK2)
        else:
            ax.axhline(0.10, color=INK2, ls="--", lw=0.7)
        ax.set_xticks(ticks)
        ax.set_xticklabels(tl, fontsize=5)
        ax.set_title(lab, fontsize=7)
    ax = fig.add_subplot(gs[0, 2])
    im = ax.imshow(reg.to_numpy(float), cmap=SEQ, vmin=0, vmax=0.6)
    ax.set_title("regret of ranking fitted at row policy, used at column policy", fontsize=6)
    ax.set_xlabel("test policy")
    ax.set_ylabel("train policy")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02).set_label("material regret rate", fontsize=6)
    fig.tight_layout()
    save(fig, "CDWH_F4_fixed_ranking")


# ----------------------------------------------------------------------- Fig 5 --
def fig5():
    import _hdata as HD

    cen = HD.load_census("H_H01")
    import _analysis as AN

    ds = []
    for pid, d in cen.items():
        v = {tuple(C.V9[b] for b in range(9) if m >> b & 1): d["alpha"][m] for m in range(512)}
        ds += [x[3] for x in AN.second_differences(v, C.V9)]
    ds = np.array(ds)
    p = pd.read_parquet(R / "H03_nested_pairs.parquet")
    pc = p[(p.level == "C") & p.pid.str.startswith("HARDENING")]
    fig, axs = plt.subplots(1, 2, figsize=(COL2, 2.3))
    ax = axs[0]
    x = np.sign(ds) * np.log10(1 + np.abs(ds) / C.TAU_RES)
    ax.hist(x, bins=120, color=GRAY1)
    ax.axvline(np.log10(2), color=INK2, lw=0.6, ls="--")
    ax.axvline(-np.log10(2), color=INK2, lw=0.6, ls="--")
    ax.set_xlabel(r"sign$(d)\,\log_{10}(1+|d_{ij}(S)|/\tau_{res})$")
    ax.set_ylabel("one-step squares (new holdout)")
    ax.set_yscale("log")
    frac_pos = float((ds > C.TAU_RES).mean())
    frac_neg = float((ds < -C.TAU_RES).mean())
    ax.set_title(f"discrete curvature: {frac_pos:.0%} > +τ_res, {frac_neg:.0%} < −τ_res", fontsize=7)
    ax = axs[1]
    for d_, col, lab in (("s2d", ORANGE, "stab→destab (needs d>0)"), ("d2s", BLUE, "destab→stab (needs d<0)")):
        q = pc[pc.dir == d_]
        ax.scatter(q.D.abs() / q.m, q.wit_d.abs(), s=4, color=col, alpha=0.5, lw=0, label=f"{lab}, n={len(q)}")
    lim = [1e-3, max(1.0, float(pc.wit_d.abs().max()) * 1.2)]
    ax.plot(lim, lim, color=INK2, lw=0.6)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"bound $|\Delta_i\alpha(S_2)-\Delta_i\alpha(S_1)|/m$ (s$^{-1}$)")
    ax.set_ylabel(r"largest certifying $|d_{ij_r}(S_{r-1})|$ (s$^{-1}$)")
    ax.set_title("every nested EM reversal has a curvature witness above the bound", fontsize=7)
    ax.legend(fontsize=5.5, loc="upper left")
    fig.tight_layout()
    save(fig, "CDWH_F5_curvature")


# ----------------------------------------------------------------------- Fig 6 --
def fig6():
    m = pd.read_csv(R / "H06_metrics.csv")
    prim = m[(m.target == "H4") & (m.gamma == 1.5) & (m.truth == "FULL")]
    piv = prim.pivot_table(index="cond", columns="predictor", values="spearman")
    statics = [c for c in piv.columns if c.startswith("S")]
    piv["oracle static"] = piv[statics].max(axis=1)
    cols = [("S1_absP", GRAY1, "|P| (discovery-selected static)"), ("oracle static", GRAY2, "best static per condition"),
            ("Dfrozen", ORANGE, "frozen dynamic"), ("Dtotal", BLUE, "total (re-equilibrated) dynamic")]
    fig, axs = plt.subplots(1, 2, figsize=(COL2, 2.5))
    ax = axs[0]
    for k in range(len(piv)):
        ax.plot(range(4), [piv.iloc[k][c] for c, _, _ in cols], color=GRID, lw=0.5, zorder=1)
    raw0 = pd.read_parquet(R / "H06_links_raw.parquet")
    r0 = raw0[raw0.family == "link"].groupby("cond").R0.max()
    elig = piv.index.isin(r0[r0 <= 1e-8].index)
    for j, (c, col, lab) in enumerate(cols):
        v = piv[c].to_numpy(float)
        xx = j + (np.random.default_rng(3).random(v.size) - 0.5) * 0.18
        ax.scatter(xx, v, s=6, color=col, lw=0, zorder=2)
        ax.scatter(xx[elig], v[elig], s=16, facecolor="none", edgecolor=INK, lw=0.6, zorder=4)
        ax.plot([j - 0.25, j + 0.25], [np.nanmedian(v)] * 2, color=INK, lw=1.5, zorder=3)
        ax.text(j, 1.06, f"{np.nanmedian(v):.2f}", ha="center", fontsize=6, color=INK)
    ax.set_xticks(range(4))
    ax.set_xticklabels([lab for _, _, lab in cols], fontsize=5.5, rotation=12)
    ax.set_ylabel("Spearman vs finite x1.5 reinforcement")
    ax.set_ylim(-1.0, 1.12)
    ax.set_title(f"{len(piv)} solved new-holdout conditions ({int(elig.sum())} prereg-eligible, ringed)", fontsize=6.5)
    ax = axs[1]
    raw = pd.read_parquet(R / "H06_links_raw.parquet")
    lk = raw[(raw.family == "link") & (raw.target == "H4")]
    pred = -0.5 * lk.d_total
    true = -(lk["g1.50_alpha"] - lk.alpha0)
    ok = np.isfinite(pred) & np.isfinite(true) & (true.abs() < 1.0)
    ax.scatter(pred[ok], true[ok], s=2, color=BLUE, alpha=0.35, lw=0)
    lim = [min(pred[ok].min(), true[ok].min()), max(pred[ok].max(), true[ok].max())]
    ax.plot(lim, lim, color=INK2, lw=0.6)
    ax.set_xlabel(r"first-order prediction $-0.5\,d\alpha/d\gamma_e$ (s$^{-1}$)")
    ax.set_ylabel(r"finite effect $-(\alpha(\gamma_e{=}1.5)-\alpha_0)$ (s$^{-1}$)")
    ax.set_title("all 46 branches x all new H4 conditions", fontsize=7)
    fig.tight_layout()
    save(fig, "CDWH_F6_goldb")


# ----------------------------------------------------------------------- Fig 7 --
def fig7():
    dec = pd.read_parquet(R / "H08_link_decomposition.parquet")
    dec = dec[dec.target == "H4"]
    ex = jload("H07_goldb_gate.json")["H8_example_modes"]
    nd = pd.read_parquet(R / "H08_node_decomposition.parquet")
    fig, axs = plt.subplots(1, 2, figsize=(COL2, 2.4), gridspec_kw={"width_ratios": [1.3, 1]})
    ax = axs[0]
    pos = 0
    ticks, tl = [], []
    sel = [("agree", ex["agree"]), ("differ", ex["differ"])]
    if ex["flip"] >= 0 and ex["flip"] != ex["differ"]:
        sel.append(("flip", ex["flip"]))
    elif ex["flip"] == ex["differ"]:
        sel[1] = ("differ = flip", ex["differ"])
    for lab, e in sel:
        if e < 0:
            continue
        g = dec[dec.idx == e]
        for j, (c, col) in enumerate((("d_frozen", ORANGE), ("indirect", AQUA), ("d_total", BLUE))):
            v = -g[c].to_numpy(float)
            ax.boxplot([v], positions=[pos + j * 0.8], widths=0.6, showfliers=False, patch_artist=True,
                       boxprops={"facecolor": "white", "edgecolor": col}, medianprops={"color": col, "lw": 1.5},
                       whiskerprops={"color": col}, capprops={"color": col})
        ticks.append(pos + 0.8)
        br = C.branches()[e]
        tl.append(f"{lab}: branch {e}\n({br['f']}-{br['t']})")
        pos += 3.2
    ax.axhline(0, color=INK2, lw=0.6)
    ax.set_xticks(ticks)
    ax.set_xticklabels(tl, fontsize=6)
    ax.set_ylabel(r"$-d\alpha/d\gamma_e$ (s$^{-1}$), new H4 conditions")
    for c, col, lab in (("d_frozen", ORANGE, "frozen"), ("indirect", AQUA, "re-equilibration"), ("d_total", BLUE, "total")):
        ax.plot([], [], color=col, lw=2, label=lab)
    ax.legend(fontsize=6, loc="best")
    ax.set_title("rule-selected branches: frozen + re-equilibration = total", fontsize=7)
    ax = axs[1]
    kinds = [("g", "Q/V gain g"), ("pll", "PLL scale"), ("ka", "AVR K_A"), ("vset", "V_set"), ("load", "load scale"), ("line", "branch γ")]
    vals = []
    for k, _ in kinds:
        if k == "line":
            r = (dec.indirect.abs() / dec.d_frozen.abs().clip(lower=1e-12)).to_numpy(float)
        else:
            r = nd[nd.kind == k].ratio.to_numpy(float)
        vals.append(np.log10(np.clip(r, 1e-12, 1e6)))
    ax.boxplot(vals, showfliers=False, patch_artist=True, boxprops={"facecolor": "white", "edgecolor": INK2}, medianprops={"color": BLUE, "lw": 1.5})
    ax.set_xticks(range(1, len(kinds) + 1))
    ax.set_xticklabels([lab for _, lab in kinds], fontsize=6, rotation=20)
    ax.axhline(0, color=INK2, lw=0.6, ls="--")
    ax.set_ylabel(r"$\log_{10}$ |re-equilibration| / |frozen|")
    ax.annotate("frozen ≡ 0", (4, 5.6), fontsize=6, ha="center", va="top", color=INK2)
    ax.set_title("controller coordinates: frozen = total (SPR)", fontsize=7)
    fig.tight_layout()
    save(fig, "CDWH_F7_why_total")


# ----------------------------------------------------------------------- Fig 8 --
def fig8():
    cl = pd.read_parquet(R / "H14_topology_em.parquet")
    h4 = cl[(cl.S == "30+33+35+37") & cl.em_class.isin(["EM-STABILIZING", "EM-DESTABILIZING", "EM-NEUTRAL"]) & ~cl.fast_flag]
    st = pd.read_csv(R / "H16_static_scores.csv", index_col=0)
    fig, axs = plt.subplots(1, 2, figsize=(COL2, 2.4), gridspec_kw={"width_ratios": [1.6, 1]})
    ax = axs[0]
    pids = [p for p in HI.DISC + HI.OLD_HOLD + HI.HPOL if p in set(h4.pid)]
    for k, p in enumerate(pids):
        g = h4[h4.pid == p]
        out = g[g.action.str.startswith("out")]
        dbl = g[g.action.str.startswith("dbl")]
        ax.scatter(k - 0.12 + (np.random.default_rng(k).random(len(out)) - 0.5) * 0.15, out.d_em, s=5, color=ORANGE, lw=0, alpha=0.8)
        ax.scatter(k + 0.12 + (np.random.default_rng(k + 9).random(len(dbl)) - 0.5) * 0.15, dbl.d_em, s=5, color=BLUE, lw=0, alpha=0.8)
    ax.axhspan(-C.TAU_MAT, C.TAU_MAT, color="#f0efec", zorder=0)
    ax.set_xticks(range(len(pids)))
    ax.set_xticklabels([p.replace("HARDENING_H", "N") for p in pids], rotation=60, fontsize=5.5)
    ax.set_yscale("symlog", linthresh=0.01)
    ax.set_ylabel(r"$\Delta\alpha_{EM}$ of H4 tracked EM mode (s$^{-1}$)")
    ax.scatter([], [], color=ORANGE, s=8, label="single-branch outage")
    ax.scatter([], [], color=BLUE, s=8, label="single-branch doubling")
    ax.legend(fontsize=6, loc="upper left")
    ax.set_title("admissible topology actions, fast-mode-dominated events removed", fontsize=7)
    ax = axs[1]
    g = h4[h4.pid == "D01"]
    x = st.loc[g.action, "fiedler"].to_numpy(float)
    ax.scatter(x, g.d_em, s=7, color=np.where(g.action.str.startswith("out"), ORANGE, BLUE), lw=0)
    ax.axhline(0, color=INK2, lw=0.6)
    ax.set_xlabel(r"$\Delta$ Fiedler value of $L_B$")
    ax.set_ylabel(r"$\Delta\alpha_{EM}$ (s$^{-1}$), P4")
    gt = jload("H15_topology_gate.json")["H16"]["fiedler"]
    ax.set_title(f"P4: static score vs EM effect (median ρ over policies {gt['median_rho_em']:.2f})", fontsize=7)
    fig.tight_layout()
    save(fig, "CDWH_F8_topology_em")


# --------------------------------------------------------------------- Fig S3 --
def figS3():
    d = pd.read_csv(R / "H13_mixing_decomposition.csv")
    fig, axs = plt.subplots(1, 2, figsize=(COL2, 2.3))
    ax = axs[0]
    for s, col in (("discovery", GRAY1), ("old", ORANGE), ("new", BLUE)):
        q = d[d.split == s]
        ax.scatter(q.F, q.C, s=10, color=col, lw=0, label=s)
    ax.axhline(0, color=INK2, lw=0.6)
    ax.axvline(0, color=INK2, lw=0.6)
    ax.set_xlabel("frequency contribution F")
    ax.set_ylabel("controller contribution C")
    ax.legend(fontsize=6)
    j = jload("H13_mixing.json")
    ax.set_title(f"two-factor decomposition of μ(θ,ω_c(θ)) − μ(θ_ref,ω_ref); var share C = {j['var_share_C']:.2f}", fontsize=6)
    ax = axs[1]
    q = d[(d.split == "new") & d.base_stable]
    ax.scatter(q.mu10, q.n_rev, s=12, color=BLUE, lw=0, label=f"fixed ω_ref (new): ρ={j['tests']['new|mu10|n_rev']['rho']:.2f}")
    q = d[(d.split == "old") & d.base_stable]
    ax.scatter(q.mu10, q.n_rev + 0.15, s=12, color=ORANGE, lw=0, marker="s", label=f"fixed ω_ref (old): ρ={j['tests']['old|mu10|n_rev']['rho']:.2f}")
    ax.set_xlabel(r"$\mu(\theta,\omega_{ref})$")
    ax.set_ylabel("# units with a stable-context reversal")
    ax.legend(fontsize=6)
    ax.set_title(f"old E12 at ω_c(θ), old holdout: ρ={j['tests']['old|mu11|n_rev']['rho']:.2f}", fontsize=7)
    fig.tight_layout()
    save(fig, "CDWH_S3_mixing")


# ----------------------------------------------------------------------- Fig 1 --
def fig1():
    import networkx as nx

    r = strongest_new_C(1)[0]
    G = nx.Graph()
    for b in C.branches():
        G.add_edge(b["f"], b["t"])
    pos = nx.kamada_kawai_layout(G)
    S1 = set() if r.S1 == "BASE" else {int(u) for u in r.S1.split("+")}
    S2 = set() if r.S2 == "BASE" else {int(u) for u in r.S2.split("+")}
    fig, axs = plt.subplots(1, 2, figsize=(COL2, 2.9))
    for ax, S, d, name in ((axs[0], S1, r.d1, "S_1"), (axs[1], S2, r.d2, "S_2")):
        ax.set_axis_off()
        nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#c3c2b7", width=0.8)
        load = [n for n in G.nodes if n not in C.SG_BUSES]
        nx.draw_networkx_nodes(G, pos, nodelist=load, ax=ax, node_size=10, node_color="#c3c2b7")
        sg = [n for n in C.SG_BUSES if n not in S and n != r.i]
        nx.draw_networkx_nodes(G, pos, nodelist=sg, ax=ax, node_size=55, node_color="white", edgecolors=INK2, linewidths=0.8)
        nx.draw_networkx_nodes(G, pos, nodelist=sorted(S), ax=ax, node_size=55, node_color=AQUA, edgecolors=INK2, linewidths=0.8)
        col = BLUE if d < 0 else ORANGE
        nx.draw_networkx_nodes(G, pos, nodelist=[r.i], ax=ax, node_size=120, node_color=col, edgecolors=INK, linewidths=1.2)
        nx.draw_networkx_labels(G, pos, labels={n: str(n) for n in C.SG_BUSES}, ax=ax, font_size=5.5)
        verb = "stabilizes" if d < 0 else "destabilizes"
        ax.set_title(rf"${name}=\{{{','.join(map(str, sorted(S))) or ''}\}}$: replacing {r.i} {verb}" + "\n"
                     + rf"$\Delta_{{{r.i}}}\alpha({name})={d:+.3f}\ \mathrm{{s}}^{{-1}}$ (same tracked EM mode)", fontsize=7)
    fig.text(0.5, 0.01, rf"nested contexts $S_1\subset S_2$ at policy {r.pid.replace('HARDENING_H', 'N')}; aqua: already replaced; white: synchronous; bus 39: interconnection",
             ha="center", fontsize=6, color=INK2)
    save(fig, "CDWH_F1_concept")


# --------------------------------------------------------------------- Fig S4 --
def figS4():
    d = pd.read_csv(R / "H18_crossmodel.csv")
    gt = jload("H18_gate.json")
    fig, axs = plt.subplots(1, 3, figsize=(COL2, 2.3))
    ax = axs[0]
    x = np.arange(len(d))
    ax.bar(x - 0.2, d.gfl_rev_A.astype(float), 0.38, color=BLUE, label="custom GFL")
    ax.bar(x + 0.2, d.alt_rev_A.astype(float) + 0.02, 0.38, color=ORANGE, label="ALT-WECC")
    ax.set_xticks(x)
    ax.set_xticklabels([p.replace("HARDENING_H", "N") for p in d.pid], fontsize=5.5, rotation=45)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["no", "yes"])
    ax.set_title(f"A reversal: {gt['A_verdict']}", fontsize=6.5)
    ax.legend(fontsize=5.5, loc="center right")
    ax = axs[1]
    ax.scatter(x, d.C_kendall, color=GRAY2, s=14, label="custom vs ALT finite ranking (C)")
    ax.scatter(x, d.C2_rho_altFD, color=BLUE, s=14, marker="s", label="ALT total vs ALT finite (C2)")
    ax.scatter(x, d.C2_rho_S1, color=GRAY1, s=14, marker="^", label="|P| vs ALT finite (C2)")
    ax.axhline(0, color=INK2, lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels([p.replace("HARDENING_H", "N") for p in d.pid], fontsize=5.5, rotation=45)
    ax.set_ylabel("rank correlation (12 lines)")
    ax.legend(fontsize=5, loc="lower left")
    ax.set_title(f"C: {gt['C_verdict']} / C2: {gt['C2_verdict']}", fontsize=6.5)
    ax = axs[2]
    for a in ("dbl45", "dbl41", "out2", "out33", "out26", "out20"):
        ax.scatter(d[f"E_{a}_gfl"], d[f"E_{a}_alt"], s=10, color=BLUE if a.startswith("dbl") else ORANGE, lw=0)
    ax.axhline(0, color=INK2, lw=0.6)
    ax.axvline(0, color=INK2, lw=0.6)
    ax.set_xlabel(r"custom GFL $\Delta\alpha$(H4)")
    ax.set_ylabel(r"ALT-WECC $\Delta\alpha$(H4)")
    ax.set_title(f"E topology: {gt['E_pooled_sign_agree']:.2f} agree ({gt['E_verdict']})", fontsize=6.5)
    fig.tight_layout()
    save(fig, "CDWH_S4_crossmodel")


# --------------------------------------------------------------------- Fig S5 --
def figS5():
    l7 = pd.read_csv(R / "H05_L7_fresh_draws.csv")
    m = pd.read_csv(R / "H06_metrics.csv")
    fig, axs = plt.subplots(1, 2, figsize=(COL2, 2.2))
    ax = axs[0]
    cov = l7.groupby("env")[["rev_A", "rev_C"]].mean().reindex(list(HI.ENVS))
    x = np.arange(len(cov))
    ax.bar(x - 0.2, cov.rev_A, 0.38, color=BLUE, label="stable-context reversal (global)")
    ax.bar(x + 0.2, cov.rev_C, 0.38, color=AQUA, label="EM same-mode reversal")
    ax.axhline(0.5, color=INK2, ls="--", lw=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(cov.index)
    ax.set_ylabel("coverage over 10 fresh draws")
    ax.legend(fontsize=5.5, loc="lower right")
    ax.set_title("P4 core lattice, fresh envelope draws", fontsize=7)
    ax = axs[1]
    q = m[(m.target == "H4") & (m.gamma == 1.5) & (m.truth == "FULL") & m.cond.str.startswith("D01|")]
    q = q.assign(env=q.cond.str.split("|").str[1])
    for k, (pr, col) in enumerate((("S1_absP", GRAY1), ("Dfrozen", ORANGE), ("Dtotal", BLUE))):
        for j, e in enumerate(HI.ENVS):
            v = q[(q.predictor == pr) & (q.env == e)].spearman
            ax.scatter(np.full(len(v), j + (k - 1) * 0.25), v, s=6, color=col, lw=0, label=pr if j == 0 else None)
    ax.set_xticks(range(4))
    ax.set_xticklabels(HI.ENVS)
    ax.set_ylabel("Spearman vs finite x1.5")
    ax.legend(fontsize=5.5, loc="lower right")
    ax.set_title("GOLD-B per fresh draw", fontsize=7)
    fig.tight_layout()
    save(fig, "CDWH_S5_uncertainty")


# --------------------------------------------------------------------- Fig S6 --
def figS6():
    g = jload("H03_gate.json")
    h7 = jload("H07_goldb_gate.json")
    fig, axs = plt.subplots(1, 2, figsize=(COL2, 2.2))
    ax = axs[0]
    lv = ["A", "B", "C", "D"]
    for k, (s, col) in enumerate((("discovery", GRAY1), ("old", ORANGE), ("new", BLUE))):
        ax.bar(np.arange(4) + (k - 1) * 0.27, [g[f"{s}_frac_{x}"] for x in lv], 0.25, color=col, label=s)
    ax.axhline(0.75, color=INK2, ls="--", lw=0.7)
    ax.set_xticks(range(4))
    ax.set_xticklabels(["A global", "B tracked", "C EM same-mode", "D same family"], fontsize=6)
    ax.set_ylabel("fraction of base-stable policies")
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=6, loc="lower left")
    ax.set_title("nested reversal by stratum", fontsize=7)
    ax = axs[1]
    labs = [("GB_prereg", "prereg-eligible"), ("GB_primary", "all solved"), ("GB_EM", "EM conditions"), ("GB_SAME", "tracked truth")]
    for k, (key, lab) in enumerate(labs):
        ax.bar(k - 0.2, h7[f"{key}_median_rho_Dtotal"], 0.38, color=BLUE, label="total" if k == 0 else None)
        ax.bar(k + 0.2, h7[f"{key}_median_rho_S1"], 0.38, color=GRAY1, label="|P|" if k == 0 else None)
    ax.set_xticks(range(4))
    ax.set_xticklabels([lab for _, lab in labs], fontsize=6)
    ax.set_ylabel("median Spearman")
    ax.legend(fontsize=6)
    ax.set_title("GOLD-B by stratum", fontsize=7)
    fig.tight_layout()
    save(fig, "CDWH_S6_stratification")


# --------------------------------------------------------------------- Fig S1 --
def figS1():
    h9 = pd.read_csv(R / "H09_corridors.csv")
    h9 = h9[h9.budget == "L1_0.50"]
    uniq = ["TXother", "TXall", "TXcore", "K2_01", "K3_01", "K3_02", "K3_12", "K4_13"]
    sets = [("old", ORANGE), ("old_draws", "#f4a582"), ("new", BLUE), ("fresh_draws", "#86b6ef")]
    fig, ax = plt.subplots(figsize=(COL2, 2.4))
    for k, c in enumerate(uniq):
        for s_i, (st, col) in enumerate(sets):
            v = h9[(h9.corridor == c) & (h9.set == st)].eff.to_numpy(float)
            x = k + (s_i - 1.5) * 0.18
            ax.scatter(np.full(v.size, x) + (np.random.default_rng(k * 7 + s_i).random(v.size) - 0.5) * 0.1, v, s=4, color=col, lw=0,
                       label=st.replace("_", " ") if k == 0 else None)
    ax.axhspan(-C.TAU_MAT, C.TAU_MAT, color="#f0efec", zorder=0)
    ax.set_xticks(range(len(uniq)))
    ax.set_xticklabels(uniq)
    ax.set_yscale("symlog", linthresh=0.01)
    ax.set_ylabel(r"$-\Delta\alpha(V_4)$ at $\sum|\Delta\gamma_e|=0.5$ (s$^{-1}$)")
    ax.legend(fontsize=6, ncol=4, loc="upper right")
    ax.set_title("equal-budget corridor effects (identical cutsets K3_01=K4_03, K3_02=K4_02, K3_12=K4_23 shown once)", fontsize=7)
    save(fig, "CDWH_S1_corridors")


# --------------------------------------------------------------------- Fig S2 --
def figS2():
    nl = pd.read_csv(R / "H10_null_percentiles.csv")
    order = [c for c in ("TXother", "TXall", "TXcore", "K2_01", "K3_01=K4_03", "K3_02=K4_02", "K3_12=K4_23", "K4_13") if c in set(nl.corridor_set)]
    fig, axs = plt.subplots(1, 2, figsize=(COL2, 2.3), sharey=True)
    for ax, fam in zip(axs, ("A", "B"), strict=True):
        for k, c in enumerate(order):
            for s_i, (sts, col) in enumerate(((("old", "old_draws"), ORANGE), (("new", "fresh_draws"), BLUE))):
                v = nl[(nl.corridor_set == c) & (nl.family == fam) & nl.set.isin(sts)].percentile.to_numpy(float)
                x = k + (s_i - 0.5) * 0.3
                ax.scatter(np.full(v.size, x) + (np.random.default_rng(k + 3 * s_i).random(v.size) - 0.5) * 0.15, v, s=4, color=col, lw=0,
                           label=("old holdout" if s_i == 0 else "new holdout") if k == 0 else None)
        ax.axhline(0.95, color=INK2, ls="--", lw=0.7)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([o.split("=")[0] for o in order], rotation=30, fontsize=6)
        ax.set_title(f"null family {fam} ({'arbitrary' if fam == 'A' else 'connected'} size-matched groups)", fontsize=7)
    axs[0].set_ylabel("percentile of corridor effect in null")
    axs[0].legend(fontsize=6, loc="lower left")
    save(fig, "CDWH_S2_null")


FIGS = {"S1": figS1, "S2": figS2, "F1": fig1, "S4": figS4, "S5": figS5, "S6": figS6, "F2": fig2, "F3": fig3, "F4": fig4, "F5": fig5, "F6": fig6, "F7": fig7, "F8": fig8, "S3": figS3}

if __name__ == "__main__":
    for k in (sys.argv[1:] or FIGS):
        FIGS[k]()
