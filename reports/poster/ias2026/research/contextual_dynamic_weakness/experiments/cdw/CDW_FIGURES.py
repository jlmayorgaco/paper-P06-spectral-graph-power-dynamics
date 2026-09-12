# ruff: noqa: E501
"""CDW figures F1-F13 from the result tables (no new model evaluation). PDF + SVG + PNG."""

from __future__ import annotations

import json

import _infra as I
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import _cdw as C  # noqa: E402

R = I.RESULTS
plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8, "legend.fontsize": 7,
                     "figure.dpi": 150, "savefig.bbox": "tight", "axes.spines.top": False, "axes.spines.right": False})
BLUE, ORANGE, GREEN, RED, GREY = "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#7f7f7f"


def save(fig, name):
    for ext in ("pdf", "svg", "png"):
        fig.savefig(I.FIGURES / f"{name}.{ext}")
    plt.close(fig)


def f1():
    s = pd.read_csv(R / "CDW_E1_summary.csv")
    pol = pd.read_csv(R / "CDW_E1_policy_summary.csv")
    order = list(C.DISCOVERY) + [f"H{i:02d}" for i in range(1, 25)]
    piv = s.pivot(index="i", columns="pid", values="frac_destab")[order]
    rev = s.pivot(index="i", columns="pid", values="rev_stable")[order]
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.8), gridspec_kw={"width_ratios": [3, 1.2]})
    im = ax[0].imshow(piv.to_numpy(float), aspect="auto", cmap="RdBu_r", vmin=0, vmax=1)
    yy, xx = np.where(rev.to_numpy(bool))
    ax[0].scatter(xx, yy, marker="x", s=10, c="k", lw=0.7, label="stable-context reversal")
    ax[0].set_yticks(range(len(piv.index)), piv.index)
    ax[0].set_xticks(range(len(order)), order, rotation=90, fontsize=5)
    ax[0].axvline(14.5, color="k", lw=0.8)
    ax[0].set_ylabel("replaced unit i")
    ax[0].set_title("fraction of contexts where replacing i destabilizes (|Δα| ≥ 0.01)")
    ax[0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), fontsize=6, frameon=False)
    fig.colorbar(im, ax=ax[0], fraction=0.03)
    unst = ~pol.base_stable
    ax[1].bar(range(len(pol)), pol.n_rev_interventions, color=np.where(pol.split == "holdout", ORANGE, BLUE))
    ax[1].set_title("reversing interventions per policy")
    ax[1].set_xlabel("policy (blue discovery, orange holdout)")
    ax[1].set_xticks([])
    for k in np.where(unst)[0]:
        ax[1].text(k, 0, "u", ha="center", va="bottom", fontsize=5)
    save(fig, "CDW_F1_contextual_reversal_map")


def f2():
    """Same-mode (tracked) stable-context reversals: the same intervention stabilizing in one
    context and destabilizing in another, both on the same modal family."""

    tr = pd.read_csv(R / "CDW_E1_tracked_reversals_holdout.csv").sort_values("strength", ascending=False)
    mg = pd.read_parquet(R / "CDW_E1_contextual_marginals.parquet")
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.5))
    for ax, (_, r) in zip(axs, tr.head(3).iterrows(), strict=False):
        g = mg[(mg.pid == r.pid) & (mg.i == r.i) & mg.stable_S & mg.same_mode].sort_values("size")
        jit = np.random.default_rng(0).uniform(-0.15, 0.15, len(g))
        col = np.where(g.delta >= C.TAU_MAT, RED, np.where(g.delta <= -C.TAU_MAT, BLUE, GREY))
        ax.scatter(g["size"] + jit, g.delta, s=5, c=col, alpha=0.8)
        for lab, d, mk in ((r.S_stab, r.d_stab, "v"), (r.S_destab, r.d_destab, "^")):
            sz = 0 if lab == "BASE" else lab.count("+") + 1
            ax.scatter([sz], [d], marker=mk, s=40, facecolors="none", edgecolors="k", lw=0.8)
            ax.annotate(f"S={lab}", (sz, d), fontsize=5, xytext=(3, 3), textcoords="offset points")
        ax.axhline(0, color="k", lw=0.5)
        ax.axhspan(-C.TAU_MAT, C.TAU_MAT, color=GREY, alpha=0.15)
        ax.set_title(f"{r.pid}: unit {int(r.i)}, same mode", fontsize=7)
        ax.set_xlabel("context size |S| (stable S)")
    axs[0].set_ylabel("Δ_i α(S) [1/s]")
    save(fig, "CDW_F2_same_intervention_two_signs")


def f3():
    nd = pd.read_csv(R / "CDW_E5_static_vs_dynamic_nodes.csv")
    fig, ax = plt.subplots(1, 2, figsize=(6.4, 2.6))
    ax[0].scatter(nd.SCR, nd.NG_norm_abs_median, c=BLUE)
    for r in nd.itertuples():
        ax[0].annotate(str(r.bus), (r.SCR, r.NG_norm_abs_median), fontsize=6)
    ax[0].set_xlabel("SCR (static)")
    ax[0].set_ylabel("|dα/dg_i| (holdout median, V9)")
    ax[1].scatter(nd.SCR, nd.frac_destab_median, c=ORANGE)
    for r in nd.itertuples():
        ax[1].annotate(str(r.bus), (r.SCR, r.frac_destab_median), fontsize=6)
    ax[1].set_xlabel("SCR (static)")
    ax[1].set_ylabel("destabilizing context fraction (median)")
    save(fig, "CDW_F3_static_vs_dynamic_nodes")


def f4():
    b = pd.read_csv(R / "CDW_E4_baselines.csv")
    h = b[b.split.str.startswith("holdout") & (b.target == "H4")]
    med = h.groupby("predictor").spearman.median().sort_values()
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.6))
    cols = [RED if p.startswith("D") else GREY for p in med.index]
    ax[0].barh(range(len(med)), med.to_numpy(), color=cols)
    ax[0].set_yticks(range(len(med)), med.index)
    ax[0].axvline(0, color="k", lw=0.5)
    ax[0].set_xlabel("median Spearman vs finite ×1.5 effect (holdout, H4)")
    lk = pd.read_parquet(R / "CDW_E4_link_sensitivity.parquet")
    d = lk[(lk.pid == "D01") & lk.env.isna() & (lk.target == "H4") & (lk.semantics == "SPR")]
    eff = -(d["large_x1.5_alpha"] - d.alpha0)
    ax[1].scatter(-d.d_total * 0.5, eff, s=8, c=RED, label="total")
    ax[1].scatter(-d.d_frozen * 0.5, eff, s=8, c=BLUE, marker="x", label="frozen")
    lim = [min(eff.min(), (-d.d_total * 0.5).min()), max(eff.max(), (-d.d_total * 0.5).max())]
    ax[1].plot(lim, lim, color="k", lw=0.5)
    ax[1].set_xlabel("first-order prediction −dα/dγ·0.5")
    ax[1].set_ylabel("finite −Δα (×1.5)")
    ax[1].set_title("P4 (D01), H4, 46 branches")
    ax[1].legend(frameon=False)
    save(fig, "CDW_F4_static_vs_total_line_ranking")


def f5():
    e2 = pd.read_csv(R / "CDW_E2_robust_summary.csv")
    e2 = e2[e2.pid == "D01"]
    mets = ["cov_rev_stable", "cov_witness_changed", "sign_persistence", "kendall_ctx", "top1_ctx", "cov_rev_given_witness_changed"]
    lab = [f"{r.source}:{r.env}" for r in e2.itertuples()]
    M = e2[mets].to_numpy(float)
    fig, ax = plt.subplots(figsize=(5.2, 2.6))
    im = ax.imshow(M, aspect="auto", cmap="viridis", vmin=0, vmax=1)
    for (i, j), v in np.ndenumerate(M):
        ax.text(j, i, f"{v:.2f}" if np.isfinite(v) else "–", ha="center", va="center", fontsize=6, color="w" if v < 0.6 else "k")
    ax.set_yticks(range(len(lab)), lab)
    ax.set_xticks(range(len(mets)), ["reversal", "witness\nchanged", "sign\npersist.", "Kendall\n(ctx)", "top-1\n(ctx)", "reversal |\nwitness chg"], fontsize=6)
    ax.set_title("P4 core lattice under envelopes (coverage fractions, not probabilities)")
    fig.colorbar(im, ax=ax, fraction=0.03)
    save(fig, "CDW_F5_robustness_heatmap")


def f6():
    import networkx as nx

    cor = json.loads((R / "CDW_E6_corridor_definitions.json").read_text())
    gate = json.loads((R / "CDW_E6_gates.json").read_text())
    G = nx.Graph()
    for b in C.branches():
        G.add_edge(b["f"], b["t"], e=b["e"])
    pos = nx.kamada_kawai_layout(G)
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.0))
    top = sorted(gate["top3_frequency"].items(), key=lambda x: -x[1])
    ax = axs[0]
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color=GREY, width=0.6)
    colors = [RED, ORANGE, GREEN]
    styles = ["solid", "dashed", "dotted"]
    widths = [3.2, 2.4, 1.8]
    # draw the largest-cardinality corridor first so a smaller, overlapping one stays visible on top
    ordered = sorted(zip(top[:3], colors, styles, widths, strict=False), key=lambda x: -len(cor[x[0][0]]))
    for (name, fr), col, ls, lw in ordered:
        es = [(b["f"], b["t"]) for b in C.branches() if b["e"] in cor[name]]
        nx.draw_networkx_edges(G, pos, edgelist=es, ax=ax, edge_color=col, width=lw, style=ls, label=f"{name} ({fr:.2f}, {len(cor[name])} branches)")
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=[40 if n in C.SG_BUSES else 8 for n in G.nodes], node_color=[BLUE if n in C.V4 else "k" for n in G.nodes])
    ax.legend(fontsize=6, frameon=False, loc="lower left")
    ax.set_title("most frequent top-3 corridors (holdout)")
    ax.axis("off")
    df = pd.read_csv(R / "CDW_E6_corridors.csv")
    h = df[df.split.str.startswith("holdout")]
    ax = axs[1]
    ax.scatter(h.sum_single_finite, h.delta_finite, s=6, c=BLUE)
    lim = [min(h.sum_single_finite.min(), h.delta_finite.min()), max(h.sum_single_finite.max(), h.delta_finite.max())]
    ax.plot(lim, lim, color="k", lw=0.5)
    ax.set_xlabel("Σ single-branch finite Δα (×1.5)")
    ax.set_ylabel("corridor finite Δα (×1.5)")
    ax.set_title("non-additivity of corridor effects")
    save(fig, "CDW_F6_corridors")


def f7():
    t = pd.read_parquet(R / "CDW_E7_topology.parquet")
    t = t[t.action != "NOMINAL"]
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.6))
    pids = sorted(t.pid.unique())
    for k, pid in enumerate(pids):
        g = t[t.pid == pid]
        for kind, col, off in (("out", BLUE, -0.15), ("dbl", RED, 0.15)):
            gg = g[g.action.str.startswith(kind)]
            axs[0].scatter(np.full(len(gg), k + off), gg.d_alpha_H4, s=4, c=col, alpha=0.6, label=kind if k == 0 else None)
    axs[0].axhline(0, color="k", lw=0.5)
    axs[0].set_xticks(range(len(pids)), pids, rotation=90, fontsize=6)
    axs[0].set_ylabel("Δα(H4) after single action")
    axs[0].legend(frameon=False, title="outage / doubling", fontsize=6)
    g = t[t.pid == "D01"]
    axs[1].scatter(g.d_fiedler, g.d_alpha_H4, s=6, c=np.where(g.action.str.startswith("out"), BLUE, RED))
    axs[1].set_xlabel("Δ Fiedler value of L_B")
    axs[1].set_ylabel("Δα(H4) at P4")
    save(fig, "CDW_F7_topology")


def f8():
    e8 = pd.read_csv(R / "CDW_E8_exchange.csv")
    fig, ax = plt.subplots(figsize=(4.4, 2.6))
    for k, (key, g) in enumerate(e8.groupby(["pid", "size"])):
        v = g.compensation_ratio.dropna().clip(upper=3)
        ax.scatter(np.full(len(v), k) + np.random.default_rng(k).uniform(-0.2, 0.2, len(v)), v, s=6)
    ax.axhline(0.2, color=RED, lw=0.8, ls="--")
    ax.set_xticks(range(len(e8.groupby(["pid", "size"]))), [f"{a}\n{b}" for a, b in e8.groupby(["pid", "size"]).groups])
    ax.set_ylabel("|Δα pair| / |Δα control only|")
    ax.set_title("local stability-equivalent exchange (paired finite)")
    save(fig, "CDW_F8_exchange_rate")


def f9():
    rows = []
    for r in I.Store("E09").all():
        if r.get("ok") and r.get("history"):
            for h in r["history"]:
                rows.append({**r["task"], **{k: h[k] for k in ("iter", "phi", "alpha_T")}})
            rows.append({**r["task"], "iter": len(r["history"]), "phi": r["phi_final"], "alpha_T": r["alpha_T_final"]})
    df = pd.DataFrame(rows)
    if df.empty:
        return
    tg = sorted(df.target.unique())
    fig, axs = plt.subplots(1, len(tg), figsize=(2.3 * len(tg), 2.4), squeeze=False)
    for ax, t in zip(axs[0], tg, strict=True):
        for (fam, m), g in df[df.target == t].groupby(["family", "method"]):
            ls = "-" if m == "plan" else "--"
            col = {"control": BLUE, "topology": RED, "joint": GREEN}[fam]
            ax.plot(g.iter, g.phi, ls=ls, color=col, marker="o", ms=2, label=f"{fam}/{m}")
        ax.axhline(0, color="k", lw=0.5)
        ax.set_title(f"{t}: Φ_T = max_S α(S)")
        ax.set_xlabel("SQP iteration")
    axs[0][0].set_ylabel("Φ_T [1/s]")
    axs[0][-1].legend(fontsize=5, frameon=False)
    save(fig, "CDW_F9_single_vs_plan_design")


def f10():
    m = pd.read_csv(R / "CDW_E12_mixing.csv")
    h4 = m[m.S == "30+33+35+37"]
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.5))
    axs[0].scatter(range(len(h4)), h4.mu_mix, c=np.where(h4.pid.str.startswith("H"), ORANGE, BLUE), s=10)
    axs[0].set_xticks(range(len(h4)), h4.pid, rotation=90, fontsize=5)
    axs[0].set_ylabel("μ_mix of T(jω_c), H4")
    axs[0].set_title("graph-modal mixing vs policy (fixed network)")
    axs[1].scatter(m.low3_share, m.alpha, s=4, c=GREY)
    axs[1].set_xlabel("share of critical angle-mode energy in L_B modes 1–3")
    axs[1].set_ylabel("α(S)")
    save(fig, "CDW_F10_graph_modal_mixing")


def f11():
    r = pd.read_csv(R / "CDW_E13_reduction.csv")
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.5))
    for fam, col in (("GM", BLUE), ("POD", ORANGE)):
        g = r[r.family.str.startswith(fam)].copy()
        g["r"] = g.family.str.replace(fam, "").astype(int)
        g = g.sort_values("r")
        axs[0].plot(g.r, g.verdict_acc, marker="o", color=col, label=fam)
    ts = r[r.family.str.startswith("TS")]
    for rr in ts.itertuples():
        axs[0].scatter([39], [rr.verdict_acc], marker="s", s=14, label=rr.family)
    axs[0].set_xlabel("retained network modes r")
    axs[0].set_ylabel("verdict accuracy (holdout)")
    axs[0].legend(fontsize=6, frameon=False)
    axs[1].scatter(r.state_ratio_median, r.end_to_end_speedup_median, c=[BLUE if f.startswith("GM") else ORANGE if f.startswith("POD") else GREEN for f in r.family])
    for rr in r.itertuples():
        axs[1].annotate(rr.family, (rr.state_ratio_median, rr.end_to_end_speedup_median), fontsize=5)
    axs[1].axhline(3, color=RED, ls="--", lw=0.8)
    axs[1].axvline(0.7, color=RED, ls="--", lw=0.8)
    axs[1].set_xlabel("state ratio")
    axs[1].set_ylabel("end-to-end speedup")
    save(fig, "CDW_F11_reduction")


def f12():
    s = pd.read_csv(R / "CDW_E14_summary.csv")
    fig, ax = plt.subplots(figsize=(4.4, 2.5))
    x = np.arange(len(s))
    cert = s.coverage * s.n
    ax.bar(x, cert, color=GREEN, label="sampled-certified")
    ax.bar(x, s.abstain_pole, bottom=cert, color=GREY, label="abstain: pole in Γ")
    ax.bar(x, s.abstain_bound, bottom=cert + s.abstain_pole, color=ORANGE, label="abstain: bound ≥ 0.9")
    ax.bar(x, s.abstain_scope, bottom=cert + s.abstain_pole + s.abstain_bound, color=RED, label="abstain: scope")
    for k, r in enumerate(s.itertuples()):
        ax.text(k, r.n * 1.01, f"false={r.false_certifications}", ha="center", fontsize=6)
    ax.set_xticks(x, s.fam)
    ax.set_ylabel("cases (holdout core)")
    ax.legend(fontsize=6, frameon=False)
    save(fig, "CDW_F12_certificate")


def f13():
    e = pd.read_csv(R / "CDW_E16_modal_energy.csv")
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.5))
    d = e[e.differ]
    axs[0].scatter(e.limit_re, e.dom_re, s=4, c=np.where(e.differ, RED, GREY))
    axs[0].plot([e.limit_re.min(), 0], [e.limit_re.min(), 0], color="k", lw=0.5)
    axs[0].set_xlabel("Re λ of stability-limiting mode")
    axs[0].set_ylabel("Re λ of energy-dominant mode")
    per = e.groupby("pid").differ.mean()
    axs[1].bar(range(len(per)), per.to_numpy(), color=np.where(per.index.str.startswith("H"), ORANGE, BLUE))
    axs[1].axhline(0.25, color=RED, ls="--", lw=0.8)
    axs[1].set_xticks([])
    axs[1].set_ylabel("fraction of stable portfolios where they differ")
    axs[1].set_xlabel("policy")
    del d
    save(fig, "CDW_F13_modal_energy_vs_limiting")


if __name__ == "__main__":
    for f in (f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13):
        try:
            f()
            print("ok", f.__name__)
        except Exception as e:  # noqa: BLE001
            print("FAILED", f.__name__, type(e).__name__, e)
