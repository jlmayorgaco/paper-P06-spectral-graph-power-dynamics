# ruff: noqa: E501
"""R21 figures F1-F12 (PDF/SVG/PNG) for the CDW68 replication. Palette: blue #2a78d6, orange #eb6834, aqua #1baf7a,
grays (dataviz reference instance). Every value is read from results/ (and raw/ for the census)."""

from __future__ import annotations

import _r68 as R  # noqa: I001

import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import networkx as nx  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import _rev68 as V  # noqa: E402

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
BAND = "#f0efec"
plt.rcParams.update({"font.family": "serif", "font.size": 8, "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5, "axes.spines.top": False,
                     "axes.spines.right": False, "legend.frameon": False, "pdf.fonttype": 42, "svg.fonttype": "none",
                     "axes.titlesize": 8, "axes.titleweight": "bold"})
COL1, COL2 = 3.45, 7.16
RES = R.RESULTS


def save(fig, name):
    for ext in ("pdf", "svg", "png"):
        fig.savefig(R.FIGS / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


def j(n):
    p = RES / n
    return json.loads(p.read_text()) if p.exists() else {}


def design():
    return json.loads((R.INPUTS / "cdw68_design_v1.json").read_text())


def f1():
    net = R.base_network()
    G = nx.Graph()
    for b in R.branches():
        G.add_edge(b["f"], b["t"])
    pos = nx.kamada_kawai_layout(G)
    fig, ax = plt.subplots(figsize=(COL2, 4.0))
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#c8c6bd", width=0.8)
    loads = set(net.loads)
    others = [n for n in G if n > 16]
    nx.draw_networkx_nodes(G, pos, nodelist=others, node_size=[14 if n in loads else 7 for n in others], node_color="#c8c6bd", ax=ax)
    cand = [b for b in R.V68]
    phys = [b for b in R.PHYSICAL if b not in cand]
    eqv = [13, 14, 15, 16]
    nx.draw_networkx_nodes(G, pos, nodelist=phys, node_size=120, node_color="white", edgecolors=INK2, linewidths=1.0, ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=cand, node_size=150, node_color=AQUA, edgecolors=INK, linewidths=1.0, ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=eqv, node_size=150, node_color=MUTED, edgecolors=INK, linewidths=1.0, node_shape="s", ax=ax)
    nx.draw_networkx_labels(G, pos, labels={b: f"G{b}" for b in range(1, 17)}, font_size=6, ax=ax)
    ax.set_axis_off()
    ax.text(0.0, -0.04, "aqua: candidate set V68 (largest scheduled output among G1-G12); white: other physical plants (NETS G1-G9, NYPS G10-G12);\n"
            "gray squares: area equivalents G13-G16 (G16 slack). Small dots: load buses.", transform=ax.transAxes, fontsize=6, color=INK2, va="top")
    save(fig, "CDW68_F1_network")


def census_marginals(model):
    p = RES / f"CDW68_R07_{model}_REAL_marginals.parquet"
    return pd.read_parquet(p)


def f2():
    fig, axs = plt.subplots(1, 2, figsize=(COL2, 3.3), sharey=False)
    for ax, model in zip(axs, ("A", "B"), strict=True):
        mg = census_marginals(model)
        t = mg[mg.lvT & np.isfinite(mg.t_delta)]
        frac = t.assign(d=t.t_delta >= R.TAU_MAT).groupby(["cond", "i"]).d.mean().unstack()
        stab = t.assign(s=t.t_delta <= -R.TAU_MAT).groupby(["cond", "i"]).s.sum().unstack()
        im = ax.imshow(frac.to_numpy(float), aspect="auto", cmap="Blues", vmin=0, vmax=1)
        ax.set_xticks(range(len(frac.columns)))
        ax.set_xticklabels([f"G{c}" for c in frac.columns])
        ax.set_yticks(range(len(frac.index)))
        ax.set_yticklabels([c.split("|")[1] for c in frac.index], fontsize=6)
        for r_, c_ in zip(*np.nonzero(stab.to_numpy() > 0), strict=True):
            ax.plot(c_, r_, marker="o", ms=3, color=ORANGE)
        ax.grid(False)
        ax.set_title(f"Model {model}: share of contexts where the replacement\nmaterially destabilizes the tracked EM mode", fontsize=7)
    fig.colorbar(im, ax=axs, fraction=0.03, pad=0.02).set_label("fraction of stable contexts", fontsize=6)
    fig.text(0.01, -0.02, "No level-D or level-T nested reversal exists (no crosses). Orange dots: cells with any materially stabilizing EM-tracked marginal.", fontsize=6, color=INK2)
    save(fig, "CDW68_F2_reversal_heatmap")


def f3():
    ref = "P68_10"
    mg = census_marginals("A")
    g = mg[(mg.cond == f"REAL|{ref}") & mg.lvT]
    fig, axs = plt.subplots(1, 6, figsize=(COL2, 2.0), sharey=True)
    for ax, i in zip(axs, R.V68, strict=True):
        q = g[g.i == i]
        size = q["mask"].map(lambda m: bin(m).count("1"))
        ax.axhspan(-R.TAU_MAT, R.TAU_MAT, color=BAND, zorder=0)
        ax.scatter(size + (np.random.default_rng(i).random(len(q)) - 0.5) * 0.3, q.t_delta, s=6, color=BLUE, lw=0, alpha=0.8)
        ax.axhline(0, color=INK2, lw=0.6)
        ax.set_title(f"replace G{i}", fontsize=7)
        ax.set_xticks(range(0, 6))
        ax.set_xlabel("|S|")
    axs[0].set_ylabel(r"$\Delta_i$ Re(tracked EM mode) (s$^{-1}$)")
    fig.suptitle(f"No nested reversal: EM-tracked marginals by context (Model A, REAL, {ref}); band = $\\pm\\tau$", fontsize=7, y=1.04)
    save(fig, "CDW68_F3_tracked_marginals")


def f4():
    ref = "P68_10"
    recs = [(r["task"]["pid"], r["task"]["S"], r["record"]) for r in R.Store("R06A").all() if r["task"]["variant"] == "REAL" and r["task"]["pid"] == ref]
    cen = V.census_from_records(recs)
    d = cen[ref]
    fT, fa = d["em_re"], d["alpha"]
    dT, da = [], []
    for S in range(V.N):
        for a in range(V.NB):
            for b in range(a + 1, V.NB):
                if S >> a & 1 or S >> b & 1:
                    continue
                i_, j_, ij = S | 1 << a, S | 1 << b, S | 1 << a | 1 << b
                dT.append(fT[ij] - fT[i_] - fT[j_] + fT[S])
                da.append(fa[ij] - fa[i_] - fa[j_] + fa[S])
    tds = j("CDW68_R14_tds_posthoc.json")
    fig, axs = plt.subplots(1, 2, figsize=(COL2, 2.3))
    ax = axs[0]
    ax.hist(np.array(dT), bins=60, color=BLUE, alpha=0.85, label="rightmost EM mode")
    ax.hist(np.array(da), bins=60, color=ORANGE, alpha=0.7, label=rf"global $\alpha_\perp$ (all $|d| \leq$ {np.nanmax(np.abs(da)):.0e})")
    ax.set_yscale("log")
    ax.set_xlabel(r"second difference $d_{ij}(S)$ (s$^{-1}$)")
    ax.set_ylabel("one-step squares")
    ax.legend(fontsize=6)
    ax.set_title(f"Interaction terms exist but stay small ({ref})", fontsize=7)
    ax = axs[1]
    for k, p in enumerate(tds.get("pairs", [])):
        i = p["i"]
        S1 = R.mask_of(() if p["S1"] == "BASE" else tuple(int(b) for b in p["S1"].split("+")))
        S2 = R.mask_of(tuple(int(b) for b in p["S2"].split("+")))
        ch = V.chain(fT, S1, S2, R.V68.index(i))
        add = [bb for bb in range(V.NB) if (S2 >> bb & 1) and not (S1 >> bb & 1)]
        vals = [fT[S1 | 1 << R.V68.index(i)] - fT[S1]]
        Sr = S1
        for bb in add:
            vals.append(vals[-1] + (fT[Sr | 1 << R.V68.index(i) | 1 << bb] - fT[Sr | 1 << R.V68.index(i)] - fT[Sr | 1 << bb] + fT[Sr]))
            Sr |= 1 << bb
        ax.plot(range(len(vals)), vals, "-o", color=(BLUE, AQUA)[k % 2], ms=3, label=f"G{i}: {p['S1']} -> {p['S2']} (resid {ch['resid']:.0e})")
    ax.axhspan(-R.TAU_MAT, R.TAU_MAT, color=BAND, zorder=0)
    ax.axhline(0, color=INK2, lw=0.6)
    ax.set_xlabel("chain step r")
    ax.set_ylabel(r"$\Delta_i f(S_0)+\sum_{k\leq r} d_{ij_k}$ (s$^{-1}$)")
    ax.legend(fontsize=5.5)
    ax.set_title("Chain identity: no sign crossing along nested chains", fontsize=7)
    fig.tight_layout()
    save(fig, "CDW68_F4_curvature")


def f5():
    h4 = json.loads((R.CDW / "results" / "hardening" / "H04_gate.json").read_text())
    rv = json.loads((R.CDW / "results" / "hardening" / "H31_revision.json").read_text())
    r9 = j("CDW68_R09_gate.json")
    labels = ["IEEE-39 A\nFULL", "IEEE-39 A\nEM", "IEEE-39 A\nscreened", "IEEE-68 A\nFULL", "IEEE-68 A\nTRACKED", "IEEE-68 B\nFULL", "IEEE-68 B\nTRACKED"]
    vals = [h4["new_FULL_median_frac_regret"], h4["new_EM_median_frac_regret"], rv["screened_ranking_median"]["new"]["regret_screened_on_feasible"],
            r9["A"]["FULL"]["median_frac_regret"], r9["A"]["TRACKED"]["median_frac_regret"], r9["B"]["FULL"]["median_frac_regret"], r9["B"]["TRACKED"]["median_frac_regret"]]
    cols = [ORANGE, ORANGE, ORANGE, BLUE, BLUE, AQUA, AQUA]
    fig, ax = plt.subplots(figsize=(COL2, 2.2))
    ax.bar(range(len(vals)), vals, color=cols, width=0.6)
    for k, v in enumerate(vals):
        ax.text(k, v + 0.01, f"{v:.3f}", ha="center", fontsize=6, color=INK2)
    ax.axhline(0.10, color=INK2, ls="--", lw=0.7)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(labels, fontsize=6)
    ax.set_ylabel("material next-action regret (median)")
    ax.set_title("Fixed node ranking: the IEEE-39 failure does not replicate on IEEE-68 (dashed: insufficiency bar 0.10)", fontsize=7)
    save(fig, "CDW68_F5_fixed_ranking")


def f6():
    mt = pd.read_csv(RES / "CDW68_R12_metrics.csv")
    g12 = j("CDW68_R12_gate.json")
    fig, axs = plt.subplots(1, 2, figsize=(COL2, 2.6), sharey=True, gridspec_kw={"width_ratios": [5, 2.4]})
    for ax, model in zip(axs, ("A", "B"), strict=True):
        q = mt[(mt.model == model) & (mt.gamma == 1.5) & (mt.split == "holdout")]
        sel = g12[model]["R12"]["selected_static"]
        preds = [("P_tot", "total (IFT)" if model == "A" else "total (numerical)", BLUE)]
        if model == "A":
            preds += [("P_fro", "frozen", ORANGE), ("P_conv", "base-\ncase", MUTED), ("P_em_tot", "EM-tracked\ntotal\n(post hoc)", AQUA)]
        preds += [(sel, f"best static\n({sel.split('_', 1)[1]})", INK2)]
        data = [q[q.predictor == p].spearman.dropna().to_numpy() for p, _, _ in preds]
        bp = ax.boxplot(data, widths=0.5, showfliers=False, patch_artist=True, medianprops={"color": INK, "lw": 1.0})
        for patch, (_, _, c) in zip(bp["boxes"], preds, strict=True):
            patch.set(facecolor="white", edgecolor=c)
        for k, (dd, (_, _, c)) in enumerate(zip(data, preds, strict=True), 1):
            ax.scatter(k + (np.random.default_rng(k).random(dd.size) - 0.5) * 0.25, dd, s=6, color=c, lw=0)
        ax.set_xticks(range(1, len(preds) + 1))
        ax.set_xticklabels([lab for _, lab, _ in preds], fontsize=6)
        ax.axhline(0.70, color=INK2, ls="--", lw=0.7)
        ax.set_title(f"Model {model}", fontsize=7)
    axs[0].set_ylabel(r"Spearman $\rho$")
    fig.suptitle("Spearman correlation with the finite x1.5 branch effect (12 holdout conditions per model; dashed: gate 0.70)", fontsize=7, y=1.02)
    save(fig, "CDW68_F6_ranking_methods")


def f7():
    mt = pd.read_csv(RES / "CDW68_R12_metrics.csv")
    g12 = j("CDW68_R12_gate.json")
    fig, ax = plt.subplots(figsize=(COL1, 2.4))
    for k, model in enumerate(("A", "B")):
        q = mt[(mt.model == model) & (mt.gamma == 1.5) & (mt.split == "holdout")]
        sel = g12[model]["R12"]["selected_static"]
        a = q[q.predictor == "P_tot"].set_index("cond").spearman - q[q.predictor == sel].set_index("cond").spearman
        ax.scatter(np.full(a.size, k) + (np.random.default_rng(k).random(a.size) - 0.5) * 0.2, a, s=10, color=(BLUE, AQUA)[k], lw=0)
        ax.hlines(a.median(), k - 0.25, k + 0.25, color=INK, lw=1.2)
    ax.axhline(0.20, color=INK2, ls="--", lw=0.7)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Model A", "Model B"])
    ax.set_ylabel("paired advantage over selected static")
    ax.set_title("Paired Spearman advantage (dashed: gate 0.20)", fontsize=7)
    save(fig, "CDW68_F7_advantage")


def f8():
    df = pd.read_csv(RES / "CDW68_R10_A_branch_table.csv")
    fig, axs = plt.subplots(1, 3, figsize=(COL2, 2.3), sharey=False)
    for ax, gm in zip(axs, (1.10, 1.25, 1.50), strict=True):
        xs, ys = df.P_tot * (gm - 1) * 1e3, df[f"T_{gm}"] * 1e3
        ax.scatter(xs, ys, s=3, color=BLUE, lw=0, alpha=0.6)
        lim = np.nanmax(np.abs(np.concatenate([xs, ys])))
        ax.plot([-lim, lim], [-lim, lim], color=INK2, lw=0.6, ls="--")
        ax.set_xlabel(r"first-order $-D_{tot}(\gamma-1)$ ($10^{-3}$ s$^{-1}$)", fontsize=6.5)
        ax.set_title(f"x{gm:.2f}", fontsize=7)
    axs[0].set_ylabel(r"finite $-\Delta\alpha_\perp$ ($10^{-3}$ s$^{-1}$)")
    mx = np.nanmax(np.abs(df[[f"T_{g}" for g in (1.10, 1.25, 1.50)]].to_numpy()))
    fig.suptitle(f"Model A: first-order prediction vs finite effect (all branches, 16 policies); largest finite effect {mx:.1e} s$^{{-1}}$, against $\\tau = 10^{{-2}}$", fontsize=7, y=1.05)
    fig.tight_layout()
    save(fig, "CDW68_F8_prediction_vs_finite")


def f9():
    p = pd.read_parquet(RES / "CDW68_PORTFOLIOS_A.parquet")
    d = design()
    abl = set(d["ablation_policies"])
    fig, axs = plt.subplots(1, 2, figsize=(COL2, 2.4))
    for k, v in enumerate(("REAL", "NOGOV", "SP33")):
        q = p[(p.variant == v) & p.pid.isin(abl)]
        axs[0].scatter(np.full(len(q), k) + (np.random.default_rng(k).random(len(q)) - 0.5) * 0.3, q.alpha, s=3, color=(BLUE, ORANGE, AQUA)[k], lw=0)
        mg = pd.read_parquet(RES / f"CDW68_R13_{v}_marginals.parquet")
        t = mg[mg.lvT].t_delta.dropna()
        axs[1].hist(t, bins=50, histtype="step", color=(BLUE, ORANGE, AQUA)[k], label=f"{v} (stabilizing >= tau: {(t <= -R.TAU_MAT).sum()})")
    axs[0].set_xticks([0, 1, 2])
    axs[0].set_xticklabels(["REAL", "NOGOV (gov off)", "SP33 (published)"], fontsize=6)
    axs[0].set_ylabel(r"$\alpha_\perp$ (s$^{-1}$)")
    axs[0].set_title("Rightmost transverse eigenvalue: slow real poles in every variant", fontsize=7)
    axs[1].axvspan(-R.TAU_MAT, R.TAU_MAT, color=BAND, zorder=0)
    axs[1].set_xlabel(r"EM-tracked marginal (s$^{-1}$)")
    axs[1].set_ylabel("marginals")
    axs[1].legend(fontsize=5.5)
    axs[1].set_title("EM-tracked marginals: almost only destabilizing", fontsize=7)
    fig.tight_layout()
    save(fig, "CDW68_F9_ablation")


def f10():
    tds = j("CDW68_R14_tds_posthoc.json")
    pairs = tds.get("pairs", [])
    if not pairs or "trace_absq" not in next(iter(pairs[0]["portfolios"].values())):
        print("F10 skipped: no traces")
        return
    fig, axs = plt.subplots(1, len(pairs), figsize=(COL2, 2.4), sharey=False)
    axs = np.atleast_1d(axs)
    for ax, p in zip(axs, pairs, strict=True):
        for tag, c, ls in (("S1", BLUE, "-"), ("S1i", BLUE, "--"), ("S2", ORANGE, "-"), ("S2i", ORANGE, "--")):
            q = p["portfolios"][tag]
            lab = "{" + q["S"].replace("+", ",") + "}" if q["S"] != "BASE" else "BASE"
            ax.semilogy(q["trace_t"], np.maximum(q["trace_absq"], 1e-12), color=c, ls=ls, lw=0.9, label=f"{lab}: decay {q['tds_decay']:+.4f} s$^{{-1}}$")
        ax.axvspan(0, 10, color=BAND, zorder=0)
        ax.set_xlabel("t (s)")
        ax.set_title(f"replace G{p['i']}: TDS change {p['tds_d1']:+.4f} / {p['tds_d2']:+.4f}\nlinear change {p['lin_d1']:+.4f} / {p['lin_d2']:+.4f} (s$^{{-1}}$)", fontsize=6.5)
        ax.set_ylim(1e-9, 1e-4)
        ax.legend(fontsize=5.5, loc="lower left")
    axs[0].set_ylabel("|modal coordinate| $|q(t)|$")
    fig.tight_layout()
    fig.suptitle("Post-hoc nonlinear phasor TDS (0.5 pu reactor at bus 3 for 10 s): same-sign effects reproduced", fontsize=7, y=1.04)
    save(fig, "CDW68_F10_tds")


def f11():
    mt = pd.read_csv(RES / "CDW68_R12_metrics.csv")
    fig, axs = plt.subplots(1, 2, figsize=(COL2, 2.4))
    for k, model in enumerate(("A", "B")):
        mg = census_marginals(model)
        t = mg[mg.lvT].t_delta.dropna()
        axs[0].hist(t, bins=60, histtype="step", color=(BLUE, AQUA)[k], label=f"Model {model}")
    axs[0].axvspan(-R.TAU_MAT, R.TAU_MAT, color=BAND, zorder=0)
    axs[0].set_xlabel(r"EM-tracked marginal (s$^{-1}$)")
    axs[0].legend(fontsize=6)
    axs[0].set_title("Replacement effects on the tracked EM mode", fontsize=7)
    a = mt[(mt.model == "A") & (mt.gamma == 1.5) & (mt.predictor == "P_tot")].set_index("cond").spearman
    b = mt[(mt.model == "B") & (mt.gamma == 1.5) & (mt.predictor == "P_tot")].set_index("cond").spearman
    b.index = [c.replace("B68", "P68") for c in b.index]
    both = pd.concat([a.rename("A"), b.rename("B")], axis=1).dropna()
    axs[1].scatter(both.A, both.B, s=12, color=BLUE, lw=0)
    axs[1].set_xlabel(r"$\rho$ total vs finite, Model A")
    axs[1].set_ylabel(r"$\rho$ total vs finite, Model B")
    axs[1].set_title("Ranking method per condition, A vs B", fontsize=7)
    fig.tight_layout()
    save(fig, "CDW68_F11_converters")


def f12():
    cb = pd.read_csv(RES / "CDW_CROSS_BENCHMARK_MATRIX.csv")
    fig, ax = plt.subplots(figsize=(COL2, 2.4))
    labs = [f"{r.benchmark}\n{r.converter.split(' (')[0]}" for r in cb.itertuples()]
    x = np.arange(len(cb))
    ax.bar(x - 0.2, pd.to_numeric(cb.total_link_median_spearman, errors="coerce"), width=0.35, color=BLUE, label="total-sensitivity Spearman")
    ax.bar(x + 0.2, pd.to_numeric(cb.best_static_median_spearman, errors="coerce"), width=0.35, color=MUTED, label="static comparator")
    for k, r in enumerate(cb.itertuples()):
        s = str(r.levelD_reversal_coverage).replace(" base-stable new policies", " policies").replace(" (pinned pole)", "").replace("EM-tracked nested reversal", "EM-tracked").replace("; ", "\n").replace(" (level T", "\n(level T")
        ax.text(k, 1.03, s, ha="center", va="bottom", fontsize=5.5, color=INK2)
    ax.set_xticks(x)
    ax.set_xticklabels(labs, fontsize=6)
    ax.set_ylim(0, 1.3)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_ylabel(r"median Spearman $\rho$")
    ax.legend(fontsize=6, loc="center left", bbox_to_anchor=(1.0, 0.5))
    ax.set_title("Cross-benchmark summary (text above bars: level-D reversal coverage)", fontsize=7)
    save(fig, "CDW68_F12_crossbench")


FIGS = {"F1": f1, "F2": f2, "F3": f3, "F4": f4, "F5": f5, "F6": f6, "F7": f7, "F8": f8, "F9": f9, "F10": f10, "F11": f11, "F12": f12}

if __name__ == "__main__":
    for k in (sys.argv[1:] or FIGS):
        FIGS[k]()
