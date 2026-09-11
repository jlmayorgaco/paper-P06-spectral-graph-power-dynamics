"""MC03: figures of the TX4 paper at IEEE column size, each with its source CSV.

Reads only frozen campaign outputs (final-closure run, paper Monte Carlo run) and
writes PDF + CSV to reports/papers/tx4_policy_dependent_incompatibility/figures.

    fig_policy_map   transverse order kappa_perp on the F7A policy plane
    fig_governed     frozen vs governed kappa on the coarse FC03 plane
    fig_closure      network-closure anatomy of the ten frozen boundaries (FC18)
    fig_port         symmetry-deflated zero-frequency port (Kundur holdout line)
    fig_mc           Monte Carlo validation summary (MC01 + MC02)
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from _mc import OUT_ROOT, ROOT  # noqa: E402

PAPER = ROOT.parents[2] / "papers" / "tx4_policy_dependent_incompatibility" / "figures"
PAPER.mkdir(parents=True, exist_ok=True)
FC = Path((OUT_ROOT / "FINAL_CLOSURE_CURRENT_RUN").read_text(encoding="utf-8").strip())
MC = Path((OUT_ROOT / "PAPER_MC_CURRENT_RUN").read_text(encoding="utf-8").strip())

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "legend.fontsize": 6.5,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "axes.linewidth": 0.6,
        "lines.linewidth": 1.1,
        "pdf.fonttype": 42,
    }
)
INK, MUTED, GRID = "#1f2328", "#5c6370", "#d8dde3"
CAT = ["#1f5fa8", "#d95f02", "#1b9e77", "#7a3fb0", "#c2185b", "#6b6b6b"]
KCOL = {1: "#67000d", 2: "#cb181d", 3: "#fb6a4a", 4: "#fcbba1", -1: "#e8f0f8"}
KLAB = {
    1: r"$\kappa=1$",
    2: r"$\kappa=2$",
    3: r"$\kappa=3$",
    4: r"$\kappa=4$",
    -1: "empty",
}
COL1, COL2 = 3.5, 7.16


def style(ax):
    ax.grid(color=GRID, lw=0.4)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def save(fig, name, source):
    fig.savefig(PAPER / f"{name}.pdf", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(PAPER / f"{name}.png", dpi=220, bbox_inches="tight", pad_inches=0.02)
    source.to_csv(PAPER / f"{name}_source.csv", index=False)
    plt.close(fig)


def sqrt_axis(ax, ticks=(0, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0)):
    ax.set_xticks([np.sqrt(t) for t in ticks], [f"{t:g}" for t in ticks])


# ------------------------------------------------------------------- maps --
def fig_policy_map():
    p = pd.read_csv(FC / "FC15_figures" / "FIG1_policy_map_source.csv")
    fig, ax = plt.subplots(figsize=(COL1, 2.55))
    for kap in (-1, 4, 3, 2, 1):
        sel = p[(p.kappa == kap) & p.exact_final.astype(bool)]
        ax.scatter(
            np.sqrt(sel.g),
            sel.k,
            s=1.3,
            c=KCOL[kap],
            marker="s",
            lw=0,
            rasterized=True,
            label=KLAB[kap],
        )
    un = p[~p.exact_final.astype(bool)]
    ax.scatter(
        np.sqrt(un.g),
        un.k,
        s=0.6,
        c=INK,
        marker=".",
        lw=0,
        rasterized=True,
        label="unresolved",
    )
    for name, g, k in (("P4", 0.03625, 1.425), (r"P$_\infty$", 1.0, 0.5)):
        ax.plot(np.sqrt(g), k, marker="o", ms=4, mfc="none", mec=INK, mew=0.9)
        ax.annotate(
            name, (np.sqrt(g), k), xytext=(4, 3), textcoords="offset points", fontsize=7
        )
    sqrt_axis(ax)
    ax.set_xlabel(r"reactive-policy (Q/V) gain $g$ (square-root scale)")
    ax.set_ylabel(r"excitation-gain scale $k$")
    ax.legend(
        markerscale=4,
        ncol=3,
        loc="upper right",
        frameon=True,
        framealpha=0.9,
        edgecolor=GRID,
        handletextpad=0.2,
        columnspacing=0.6,
    )
    save(fig, "fig_policy_map", p[["g", "k", "kappa", "H_perp", "exact_final"]])


def fig_governed():
    p = pd.read_csv(FC / "FC03_governed_replication" / "FC03_points.csv")
    p = p[p.tag == "PLANE"].copy()
    fig, axes = plt.subplots(1, 2, figsize=(COL1, 1.85), sharey=True)
    for ax, col, title in (
        (axes[0], "kappa_frozen", "frozen (no governor)"),
        (axes[1], "kappa_governed", "documented TGOV1N"),
    ):
        hcol = "H_frozen" if col == "kappa_frozen" else "H_governed"
        for kap in (-1, 4, 3, 2, 1):
            sel = p[(p[col] == kap) & ~p[hcol].astype(str).str.startswith("BASE")]
            ax.scatter(
                np.sqrt(sel.g),
                sel.k,
                s=9,
                c=KCOL[kap],
                marker="s",
                lw=0,
                label=KLAB[kap],
            )
        base = p[p[hcol].astype(str).str.startswith("BASE")]
        ax.scatter(
            np.sqrt(base.g),
            base.k,
            s=9,
            c="#9aa3ad",
            marker="x",
            lw=0.7,
            label="base unstable",
        )
        ax.set_title(title, pad=2)
        sqrt_axis(ax, (0, 0.1, 0.5, 1.0))
        ax.set_xlabel(r"$g$")
    axes[0].set_ylabel(r"$k$")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=6,
        bbox_to_anchor=(0.5, -0.2),
        frameon=False,
        handletextpad=0.1,
        columnspacing=0.5,
        markerscale=1.2,
    )
    save(
        fig,
        "fig_governed",
        p[["g", "k", "H_frozen", "kappa_frozen", "H_governed", "kappa_governed"]],
    )


# ---------------------------------------------------------------- closure --
def fig_closure():
    s = json.loads(
        (FC / "FC18_targeted_port_checks" / "FC18_summary.json").read_text(
            encoding="utf-8"
        )
    )
    rows = []
    for e in s["events"]:
        tr = {int(k): v for k, v in e["check1"]["interaction_truncation_abs"].items()}
        n = max(tr)
        rows.append(
            {
                "event": e["event"].split()[0].replace("FLAG", "F"),
                "subset": e["subset"],
                "size": n,
                "local": e["check1"]["local_factor_abs"],
                **{f"trunc_le_{k}": tr[k] for k in tr},
            }
        )
    df = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(COL1, 2.0))
    x = np.arange(len(df))
    w = 0.2
    ax.bar(x - 1.5 * w, df.local, w, color=CAT[5], label="local factor")
    ax.bar(x - 0.5 * w, df.trunc_le_2, w, color=CAT[0], label=r"order $\leq 2$")
    ax.bar(
        x + 0.5 * w,
        df.get("trunc_le_3", np.nan).where(df["size"] > 3),
        w,
        color=CAT[1],
        label=r"order $\leq 3$ ($|S|=4$)",
    )
    full = [r[f"trunc_le_{int(r['size'])}"] for _, r in df.iterrows()]
    ax.bar(x + 1.5 * w, full, w, color=CAT[2], label=r"full order $|S|$")
    ax.set_yscale("log")
    ax.set_ylim(1e-15, 2)
    ax.set_xticks(
        x,
        [f"{a}\n{int(b)}" for a, b in zip(df.event, df["size"], strict=True)],
        fontsize=6,
    )
    ax.set_ylabel(r"$|$truncated $\det(I+Q_{SS})|$")
    ax.set_xlabel(
        r"boundary event and size $|S|$ of the changing coalition", labelpad=1
    )
    ax.legend(
        ncol=4,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.0),
        frameon=False,
        handletextpad=0.3,
        columnspacing=0.8,
    )
    style(ax)
    save(fig, "fig_closure", df)


# ------------------------------------------------------------------- port --
def fig_port():
    p = pd.read_csv(FC / "FC15_figures" / "FIG6_zero_frequency_port_source.csv")
    fig, axes = plt.subplots(2, 1, figsize=(COL1, 2.4), sharex=True)
    for i, (sub, blk) in enumerate(p.groupby("subset")):
        axes[0].plot(
            blk.g * 1e3,
            blk.rightmost_real_transverse,
            color=CAT[i],
            label=f"subset {sub}",
        )
        axes[1].step(
            blk.g * 1e3,
            blk.reloc_sign * (1 - 0.06 * i),
            where="mid",
            color=CAT[i],
            label=rf"relocated port, {sub}",
        )
        axes[1].plot(
            blk.g * 1e3,
            np.sign(blk.raw_sign) * 0.5,
            color=CAT[i],
            ls=":",
            lw=0.7,
            label="plain port (sign/2)" if i == 0 else None,
        )
    axes[0].axhline(0, color=MUTED, lw=0.6)
    axes[0].set_ylabel("rightmost real\n" + r"eig. of $A_\perp$ [s$^{-1}$]")
    axes[1].set_ylabel(r"sign $\det T^{\#\#}(0)$")
    axes[1].set_yticks([-1, 0, 1])
    axes[1].set_ylim(-1.25, 1.25)
    axes[1].set_xlabel(r"$g\times 10^{3}$ (Kundur, $k=1.1125$, holdout line)")
    axes[0].legend(loc="upper left", frameon=False)
    axes[1].legend(
        loc="center right", frameon=True, framealpha=0.95, edgecolor=GRID, fontsize=6
    )
    for ax in axes:
        style(ax)
    save(fig, "fig_port", p)


# --------------------------------------------------------------------- MC --
def fig_mc():
    s1 = pd.read_csv(MC / "MC01_synthetic" / "S1.csv")
    s6 = pd.read_csv(MC / "MC01_synthetic" / "S6.csv")
    p1 = pd.read_csv(MC / "MC02_ieee39" / "P1.csv")
    p3 = pd.read_csv(MC / "MC02_ieee39" / "P3.csv")
    p4 = pd.read_csv(MC / "MC02_ieee39" / "P4.csv")
    p2_path = MC / "MC02_ieee39" / "P2.csv"
    p2 = pd.read_csv(p2_path) if p2_path.exists() else None
    syn = json.loads(
        (MC / "MC01_synthetic" / "MC01_summary.json").read_text(encoding="utf-8")
    )
    fig, axes = plt.subplots(
        1, 3, figsize=(COL2, 2.05), gridspec_kw={"width_ratios": [1.45, 1.0, 1.0]}
    )
    # (a) residuals of the exact identities
    checks = [
        ("S1 det. identity", s1.det_identity_rel, 1e-8),
        ("S6 minors", s6.vertex_identity_rel, 1e-10),
        ("S6 invariance", s6.invariance_rel, 1e-10),
        ("P1 spectrum", p1.spectral_identity, 1e-6),
        ("P1 rotation (NL)", p1.nonlinear_rotation, 1e-8),
        ("P1 drift (NL)", p1.nonlinear_drift, 1e-8),
        ("P3 port ratio", p3.vertex_identity_rel, 1e-8),
        ("P3 Moebius", p3.moebius_minor_rel, 1e-8),
    ]
    ax = axes[0]
    rng = np.random.default_rng(0)
    rows = []
    for i, (name, vals, thr) in enumerate(checks):
        v = np.log10(np.maximum(vals.dropna().to_numpy(), 1e-18))
        ax.scatter(
            v,
            i + 0.25 * (rng.random(v.size) - 0.5),
            s=1.2,
            color=CAT[0],
            lw=0,
            rasterized=True,
        )
        ax.plot([np.log10(thr)] * 2, [i - 0.35, i + 0.35], color=CAT[1], lw=1.2)
        rows.append(
            {
                "check": name,
                "n": v.size,
                "max_log10": float(v.max()),
                "threshold_log10": float(np.log10(thr)),
            }
        )
    ax.set_yticks(range(len(checks)), [c[0] for c in checks], fontsize=6.3)
    ax.set_xlabel(r"$\log_{10}$ residual (bar: frozen pass threshold)")
    ax.set_title("(a) exact identities", loc="left")
    style(ax)
    # (b) TDS growth rate vs alpha_perp
    ax = axes[1]
    for lab, col, mk in (
        ("RECOVERS", CAT[2], "o"),
        ("FAILS_DECLARED_SECURITY", CAT[1], "^"),
        ("OUTSIDE_MODEL_SCOPE", CAT[3], "s"),
    ):
        sel = p4[p4.label == lab]
        if sel.empty:
            continue
        ax.scatter(
            sel.alpha_perp,
            sel.envelope_rate,
            s=7,
            color=col,
            marker=mk,
            lw=0,
            label=lab.split("_")[0].lower(),
        )
    lim = [
        min(p4.alpha_perp.min(), p4.envelope_rate.min()) - 0.02,
        max(p4.alpha_perp.max(), p4.envelope_rate.max()) + 0.02,
    ]
    ax.plot(lim, lim, color=MUTED, lw=0.6, ls="--")
    ax.axvline(0, color=MUTED, lw=0.5)
    ax.axhline(0, color=MUTED, lw=0.5)
    ax.set_xlabel(r"transverse abscissa $\alpha_\perp$ [s$^{-1}$]")
    ax.set_ylabel(r"TDS envelope rate [s$^{-1}$]")
    ax.set_title("(b) nonlinear TDS, 5 MW pulse", loc="left")
    ax.legend(frameon=False, loc="upper left", handletextpad=0.1)
    style(ax)
    # (c) crossing frequencies on random IEEE-39 segments
    ax = axes[2]
    rb = MC / "MC04_p2_rebisect" / "MC04_events.csv"
    if rb.exists():  # exact crossings from the post-hoc re-bisection
        cr = pd.read_csv(rb)
        cr["freq_hz"] = cr.im / (2 * np.pi)
    else:
        cr = (
            p2[p2.subset.fillna("").astype(bool)]
            if p2 is not None
            else pd.DataFrame({"freq_hz": []})
        )
    if p2 is None:
        ax.text(
            0.5, 0.5, "P2 pending", ha="center", va="center", transform=ax.transAxes
        )
    ax.hist(cr.freq_hz, bins=np.linspace(0, 1.6, 33), color=CAT[0])
    ax.set_xlabel("crossing frequency [Hz]")
    ax.set_ylabel("located crossings")
    ax.set_title("(c) crossings closing changes of $\\mathcal{H}$", loc="left")
    style(ax)
    fig.tight_layout(w_pad=0.8)
    src = pd.DataFrame(rows)
    src["S2_S4_summary"] = json.dumps(
        {k: syn[k] for k in ("S2", "S4") if k in syn}, default=str
    )[:30000]
    save(fig, "fig_mc", src)
    p4.to_csv(PAPER / "fig_mc_tds_source.csv", index=False)
    cr.to_csv(PAPER / "fig_mc_crossings_source.csv", index=False)


def fig_path():
    g1 = pd.read_csv(ROOT / "results" / "G1" / "G1_ieee39_points.csv.gz")
    line = g1[(g1["map"] == "F7B") & (np.isclose(g1.t, 0.8515625))].sort_values("g")
    line = line[line.g <= 0.35]
    s = json.loads(
        (FC / "FC18_targeted_port_checks" / "FC18_summary.json").read_text(
            encoding="utf-8"
        )
    )
    ev = pd.DataFrame(
        [
            {
                "event": e["event"].split()[0],
                "g": e["g_star"],
                "sigma1": e["check2"]["sv_Q"][0],
                "sigma8": e["check2"]["sv_Q"][-1],
                "rank_1e-2": e["check2"]["eff_rank_Q_1e-2"],
                "size": len(e["subset"].split("+")),
            }
            for e in s["events"]
            if e["event"] != "FLAG P4-line"
        ]
    )
    fig, axes = plt.subplots(
        2, 1, figsize=(COL1, 2.35), sharex=True, gridspec_kw={"height_ratios": [1, 1.1]}
    )
    kap = line.kappa_RHP.replace(-1, 5)
    axes[0].step(np.sqrt(line.g), kap, where="mid", color=CAT[0])
    axes[0].set_yticks([2, 3, 4, 5], ["2", "3", "4", "none"])
    axes[0].set_ylabel(r"$\kappa_\perp$")
    for _, r in ev.iterrows():
        axes[0].axvline(np.sqrt(r.g), color=GRID, lw=0.6, zorder=0)
    axes[1].semilogy(
        np.sqrt(ev.g), ev.sigma1, "o-", color=CAT[1], ms=3, label=r"$\sigma_1(Q)$"
    )
    axes[1].semilogy(
        np.sqrt(ev.g),
        ev.sigma8,
        "s-",
        color=CAT[2],
        ms=3,
        label=r"$\sigma_8(Q)$ (smallest)",
    )
    for _, r in ev.iterrows():
        axes[1].annotate(
            str(int(r["rank_1e-2"])),
            (np.sqrt(r.g), r.sigma1),
            xytext=(0, 4),
            textcoords="offset points",
            ha="center",
            fontsize=6,
        )
    axes[1].set_ylabel("singular value")
    axes[1].legend(frameon=False, loc="lower left", ncol=2)
    sqrt_axis(axes[1], (0.01, 0.05, 0.1, 0.2, 0.3))
    axes[1].set_xlabel(
        r"$g$ on the F7B line $t=0.852$ (numbers: effective rank of $Q$)"
    )
    for ax in axes:
        style(ax)
    save(
        fig,
        "fig_path",
        pd.concat(
            [
                line[["g", "kappa_RHP", "H_RHP"]].assign(kind="map"),
                ev.assign(kind="event"),
            ],
            ignore_index=True,
        ),
    )


def fig_nonlinear():
    rk = pd.read_csv(FC / "FC06_resilience_complex" / "FC06_R_k.csv")
    s6 = json.loads(
        (FC / "FC06_resilience_complex" / "FC06_summary.json").read_text(
            encoding="utf-8"
        )
    )
    fam = {
        "D1": "D1 load pulse [MW]",
        "D2": "D2 $P_m$ dip [MW]",
        "D3": "D3 GFL dip [p.u.]",
    }
    fig, axes = plt.subplots(1, 3, figsize=(COL2, 1.6))
    rows = []
    for j, f in enumerate(("D1", "D2", "D3")):
        ax = axes[j]
        for i, pt in enumerate(("P4", "P_inf")):
            blk = rk[(rk.point == pt) & (rk.family == f)].sort_values("k")
            scope = s6[f"{pt}|{f}"]["rho_scope"]
            ax.plot(
                blk.k,
                blk.R_k,
                "o-",
                color=CAT[i],
                ms=3,
                label="P4" if pt == "P4" else r"P$_\infty$",
            )
            ax.axhline(scope, color=CAT[i], ls=":", lw=0.9)
            rows.append(
                {
                    "family": f,
                    "point": pt,
                    "rho_scope": scope,
                    **{f"R_{int(k)}": v for k, v in zip(blk.k, blk.R_k, strict=True)},
                }
            )
        ax.set_title(fam[f], pad=2)
        ax.set_xticks([1, 2, 3, 4])
        ax.set_xlabel(r"portfolio size $k$")
        style(ax)
    axes[0].set_ylabel(r"$R_k=\min_{|S|\leq k}\, r_S$")
    axes[0].legend(frameon=False)
    fig.tight_layout(w_pad=0.6)
    save(fig, "fig_nonlinear", pd.DataFrame(rows))


def fig_curvature():
    e = pd.read_csv(FC / "FC07_second_order_curvature" / "FC07_errors.csv")
    case = "critical triple 30+33+35 P4"
    blk = e[e.case == case]
    fig, axes = plt.subplots(1, 2, figsize=(COL1, 1.75))
    names = {
        "M0 linear": "linear",
        "M1 full second order": "full 2nd order",
        "M2 no KCL curvature": r"no $f_zD^2\psi$",
        "M2t no KCL curvature, transverse forcing": r"no $f_zD^2\psi$, transv.",
    }
    for ax, col, lab in (
        (axes[0], "err_freq_hz", "COI frequency [Hz]"),
        (axes[1], "err_v16_pu", r"$|V_{16}|$ [p.u.]"),
    ):
        for i, (mdl, b) in enumerate(blk.groupby("model")):
            ax.loglog(
                b.eps_mw, b[col], marker="o", ms=2.5, color=CAT[i], label=names[mdl]
            )
        ax.set_title(lab, pad=2)
        ax.set_xlabel(r"pulse $\varepsilon$ [MW]")
        style(ax)
    axes[0].set_ylabel("max error vs. full DAE")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=4,
        bbox_to_anchor=(0.5, -0.12),
        frameon=False,
        handletextpad=0.3,
        columnspacing=0.7,
    )
    fig.tight_layout(w_pad=0.5)
    save(fig, "fig_curvature", blk)


def fig_dispatch():
    d = pd.read_csv(FC / "FC04_e14_n6_pg" / "FC04_n6_pg.csv")
    fig, ax = plt.subplots(figsize=(COL1, 1.9))
    for role, col, mk, lab in (
        ("failing", CAT[1], "^", "unstable four-unit portfolios"),
        ("control", CAT[0], "o", "stable, matched on dispatch"),
    ):
        sel = d[d.role == role]
        ax.scatter(
            sel.replaced_pg_mw,
            sel.family_alpha,
            s=12,
            color=col,
            marker=mk,
            lw=0,
            label=lab,
        )
    ax.axhline(0, color=MUTED, lw=0.6)
    ax.set_xlabel("displaced active dispatch [MW]")
    ax.set_ylabel(r"inter-area family $\max\mathrm{Re}$ [s$^{-1}$]")
    ax.legend(frameon=False, loc="lower left")
    style(ax)
    save(fig, "fig_dispatch", d)


def fig_retune():
    pts = pd.read_csv(FC / "FC03_governed_replication" / "FC03_points.csv")
    line = pts[pts.tag == "PATH k=1.425"].sort_values("g")
    s = json.loads(
        (FC / "FC18_targeted_port_checks" / "FC18_summary.json").read_text(
            encoding="utf-8"
        )
    )
    it = pd.DataFrame(s["retune_P4"])
    fig, ax = plt.subplots(figsize=(COL1, 1.85))
    ax.plot(
        line.g, line.alpha_flag_frozen, color=CAT[0], label=r"direct $\alpha_\perp(g)$"
    )
    g0, a0, d0 = it.g.iloc[0], it.re_s.iloc[0], it.dRe_dg_port.iloc[0]
    gg = np.linspace(g0, g0 + 0.09, 10)
    ax.plot(
        gg,
        a0 + d0 * (gg - g0),
        color=CAT[1],
        ls="--",
        lw=0.9,
        label="one linear step (port slope)",
    )
    ax.plot(it.g, it.re_s, "o", color=CAT[2], ms=3.5, label="port-Newton iterates")
    for i, r in it.iterrows():
        if i < 4:
            ax.annotate(
                str(i),
                (r.g, r.re_s),
                xytext=(3, 3),
                textcoords="offset points",
                fontsize=6,
            )
    ax.axhline(0, color=MUTED, lw=0.6)
    ax.axvline(s["retune_direct_boundary_g"], color=MUTED, lw=0.6, ls=":")
    ax.set_xlim(0, 0.32)
    ax.set_ylim(-0.12, 0.2)
    ax.set_xlabel(r"Q/V gain $g$ on the P4 line ($k=1.425$, $t=1.5$)")
    ax.set_ylabel(r"$\mathrm{Re}\,s^*$, $\alpha_\perp$ [s$^{-1}$]")
    ax.legend(frameon=False, loc="upper right")
    style(ax)
    save(
        fig,
        "fig_retune",
        pd.concat(
            [
                line[["g", "alpha_flag_frozen"]].assign(kind="direct"),
                it.assign(kind="newton"),
            ],
            ignore_index=True,
        ),
    )


FIGS = {
    "retune": fig_retune,
    "dispatch": fig_dispatch,
    "policy": fig_policy_map,
    "governed": fig_governed,
    "closure": fig_closure,
    "port": fig_port,
    "mc": fig_mc,
    "path": fig_path,
    "nonlinear": fig_nonlinear,
    "curvature": fig_curvature,
}


def main(argv) -> int:
    for n in argv or list(FIGS):
        try:
            FIGS[n]()
            print(n, "ok")
        except FileNotFoundError as err:
            print(n, "skipped:", err)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
