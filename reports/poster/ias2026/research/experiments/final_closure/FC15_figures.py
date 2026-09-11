"""FC15 (step 18): the eight paper-ready figures, each with its source-data CSV.

    FIG1  policy-dependent transverse incompatibility map (F7A slice, kappa_perp)
    FIG2  kappa_NL(rho) overlaid with kappa_RHP_perp
    FIG3  R_k curves per disturbance family
    FIG4  H_RHP_perp vs H_NL at one decisive severity
    FIG5  TDS just below / just above a nonlinear threshold
    FIG6  symmetry-deflated real-axis port closure (Kundur holdout line)
    FIG7  linear vs full second order vs no-KCL curvature error scaling
    FIG8  planning: spectral-only vs nonlinear-safe vs support-assisted

Usage: python FC15_figures.py FIG1 FIG6 ...
(no argument: every figure whose data exist)
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
from _fc import RESULTS, out_dir  # noqa: E402

OUT = out_dir("FC15_figures")
INK, MUTED, GRID = "#1f2328", "#5c6370", "#d8dde3"
# fixed categorical order (never cycled); status colors are not reused for series
CAT = ["#1f5fa8", "#d95f02", "#1b9e77", "#7a3fb0", "#c2185b", "#6b6b6b"]
KAPPA_COLORS = {1: "#7a1c1c", 2: "#c0392b", 3: "#e67e22", 4: "#f5c16c", -1: "#dfe9f3"}


def style(ax):
    ax.grid(color=GRID, lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(colors=INK, labelsize=9)


def save(fig, name, source: pd.DataFrame):
    fig.savefig(OUT / f"{name}.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    source.to_csv(OUT / f"{name}_source.csv", index=False)
    plt.close(fig)


# ------------------------------------------------------------------- FIG1 --
def fig1():
    p = pd.read_csv(
        out_dir("FC01_transverse_quotient") / "FC01_points.csv.gz", low_memory=False
    )
    p = p[p["map"] == "F7A"].copy()
    kd = pd.to_numeric(p.get("kappa_perp_direct"), errors="coerce")
    p["kappa"] = kd.fillna(p.kappa_perp)
    ex = p.get("exact_direct")
    p["exact_final"] = (
        ex.map({True: True, False: False, "True": True, "False": False}).fillna(p.exact)
    ).astype(bool)
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    for kap, col in KAPPA_COLORS.items():
        sel = p[(p.kappa == kap) & p.exact_final]
        lbl = (
            "no hyperedge (transversely composable)"
            if kap == -1
            else f"kappa_perp = {kap}"
        )
        ax.scatter(
            sel.g, sel.k, s=2.2, c=col, marker="s", lw=0, label=lbl, rasterized=True
        )
    un = p[~p.exact_final]
    ax.scatter(
        un.g,
        un.k,
        s=1.5,
        c=INK,
        marker=".",
        lw=0,
        label="within classifier resolution (unresolved)",
        rasterized=True,
    )
    ax.set_xlabel("reactive-policy gain g (leaky Q/V regulator)", color=INK)
    ax.set_ylabel("surviving-AVR gain scale k", color=INK)
    ax.set_title(
        "Transverse incompatibility order kappa_perp over the policy plane\n"
        "IEEE-39, F7A slice (t = 1.5, h = 1), core 30/33/35/37, matched dispatch, "
        "no governor",
        fontsize=10,
        color=INK,
    )
    style(ax)
    leg = ax.legend(markerscale=5, fontsize=8, loc="upper right", frameon=True)
    leg.get_frame().set_edgecolor(GRID)
    save(fig, "FIG1_policy_map", p[["g", "k", "kappa", "H_perp", "exact_final"]])


# ------------------------------------------------------------------- FIG6 --
def fig6():
    sys.path.insert(0, str(Path(__file__).parent))
    from F12_kundur import solve as k_solve
    from ibr_cycles.certification.port_origin import relocated_port
    from ibr_cycles.certification.symmetry import frequency_partner, rotation_generator
    from ibr_cycles.certification.transverse import transverse_operator
    from ibr_cycles.dynamics.linearize import central_difference_jacobians

    rows = []
    for members in ((3, 4), (2, 3, 4)):
        for g in np.linspace(0.0, 0.0012, 49):
            case = k_solve(members, {"g": float(g), "k": 1.1125, "t": 1.0})
            jac = central_difference_jacobians(
                case.dae, case.equilibrium.x, case.equilibrium.z, {}
            )
            a = jac.fx - jac.fz @ np.linalg.solve(jac.gz, jac.gx)
            r_x, _ = rotation_generator(case.dae, case.equilibrium.z)
            w = frequency_partner(case.dae).w
            ev = np.linalg.eigvals(transverse_operator(a, r_x, w).a_perp)
            real = ev[np.abs(ev.imag) < 1e-9]
            rp = relocated_port(jac, r_x, w, 1.0, 1.0)
            n = jac.fx.shape[0]
            # unrelocated port at s = 0 (structural double zero, split by rounding)
            s0 = 1e-9
            t_raw = jac.gz + jac.gx @ np.linalg.solve(s0 * np.eye(n) - jac.fx, jac.fz)
            sr, lr = np.linalg.slogdet(t_raw)
            rows.append(
                {
                    "subset": "+".join(map(str, members)),
                    "g": g,
                    "rightmost_real_transverse": float(real.real.max()),
                    "reloc_sign": rp.t0_sign,
                    "reloc_logabs": rp.t0_logabs,
                    "raw_sign": float(np.real(sr)),
                    "raw_logabs": float(lr),
                }
            )
    df = pd.DataFrame(rows)
    hold = pd.read_csv(RESULTS / "zero_frequency_port_validation.csv")
    fig, axes = plt.subplots(
        1, 3, figsize=(12.5, 3.8), gridspec_kw={"width_ratios": [1.1, 1.1, 0.9]}
    )
    for i, (sub, blk) in enumerate(df.groupby("subset")):
        axes[0].plot(
            blk.g * 1e3,
            blk.rightmost_real_transverse,
            color=CAT[i],
            lw=2,
            label=f"subset {sub}",
        )
        axes[1].plot(
            blk.g * 1e3,
            blk.reloc_sign,
            color=CAT[i],
            lw=2,
            drawstyle="steps-mid",
            label=f"relocated port, {sub}",
        )
        axes[1].plot(
            blk.g * 1e3,
            np.sign(blk.raw_sign) * 0.5,
            color=CAT[i],
            lw=1,
            ls=":",
            label=f"original port (s -> 0), {sub}",
        )
    axes[0].axhline(0, color=MUTED, lw=1)
    axes[0].set_xlabel("g x 1e3", color=INK)
    axes[0].set_ylabel("rightmost real eigenvalue of A_perp [1/s]", color=INK)
    axes[0].set_title(
        "Physical real crossing (Kundur, holdout line k = 1.1125)", fontsize=9
    )
    axes[1].set_xlabel("g x 1e3", color=INK)
    axes[1].set_ylabel("sign det T##(0)  (dotted: original, scaled)", color=INK)
    axes[1].set_title("The relocated port changes sign at the crossing", fontsize=9)
    for ax in axes[:2]:
        style(ax)
        ax.legend(fontsize=7)
    axes[2].axis("off")
    txt = "\n".join(
        f"{r.source}: {int(r.detected)}/{int(r.holdout_crossings)} detected, "
        f"{int(r.false_positives)} false positives / "
        f"{int(r.negative_controls)} controls"
        for r in hold.itertuples()
    )
    axes[2].text(
        0.0,
        0.6,
        "Holdout (frozen port_relocation_v1.1)\n\n"
        + txt
        + "\n\nKundur miss: a real-pair coalescence\n"
        "into a complex pair (not a zero crossing)",
        fontsize=8,
        color=INK,
        va="center",
    )
    save(fig, "FIG6_zero_frequency_port", df)


FAM = {"D1": "D1 load pulse [MW]", "D2": "D2 Pm dip [MW]", "D3": "D3 PV dip [fraction]"}


# ------------------------------------------------------------------- FIG2 --
def fig2():
    k = pd.read_csv(out_dir("FC06_resilience_complex") / "FC06_kappa_NL.csv")
    fig, axes = plt.subplots(2, 3, figsize=(12.5, 6.2), sharey=True)
    for i, pt in enumerate(("P4", "P_inf")):
        for j, fam in enumerate(("D1", "D2", "D3")):
            ax = axes[i, j]
            blk = k[(k.point == pt) & (k.family == fam)].sort_values("rho_lo")
            kr = blk.kappa_RHP_perp.iloc[0] if len(blk) else np.nan
            kr_plot = 5 if kr == -1 else kr
            scope = float(blk.rho_scope.iloc[0]) if len(blk) else np.nan
            xmax = (
                scope
                if np.isfinite(scope)
                else (blk.rho_hi.replace(np.inf, np.nan).max() or 1.0)
            )
            for r in blk.itertuples():
                hi = min(r.rho_hi, xmax) if np.isfinite(r.rho_hi) else xmax
                kn = 5 if r.kappa_NL == -1 else r.kappa_NL
                ax.hlines(
                    kn,
                    r.rho_lo,
                    hi,
                    color=CAT[0],
                    lw=3,
                    label="kappa_NL(rho)" if r.Index == blk.index[0] else None,
                )
            ax.axhline(kr_plot, color=CAT[1], lw=1.5, ls="--", label="kappa_RHP_perp")
            if np.isfinite(scope):
                ax.axvspan(
                    scope,
                    scope * 1.6 + 1e-9,
                    color=GRID,
                    alpha=0.8,
                    label="outside model scope (rho >= rho_scope)",
                )
                ax.set_xlim(0, scope * 1.6)
            ax.set_yticks([1, 2, 3, 4, 5], ["1", "2", "3", "4", "none"])
            ax.set_title(f"{pt}, {FAM[fam]}", fontsize=9)
            ax.set_xlabel("disturbance severity rho", fontsize=8)
            style(ax)
            if i == 0 and j == 0:
                ax.set_ylabel("minimum incompatible order")
                ax.legend(fontsize=7, loc="lower left")
    fig.suptitle(
        "Nonlinear composability order vs the transverse spectral order "
        "(IEEE-39, TDS, not EMT)",
        fontsize=10,
    )
    fig.subplots_adjust(hspace=0.45)
    save(fig, "FIG2_kappa_NL", k)


# ------------------------------------------------------------------- FIG3 --
def fig3():
    rk = pd.read_csv(out_dir("FC06_resilience_complex") / "FC06_R_k.csv")
    thr = pd.read_csv(out_dir("FC05_nonlinear_thresholds") / "FC05_thresholds.csv")
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8))
    for j, fam in enumerate(("D1", "D2", "D3")):
        ax = axes[j]
        for i, pt in enumerate(("P4", "P_inf")):
            blk = rk[(rk.point == pt) & (rk.family == fam)].sort_values("k")
            ax.plot(
                blk.k,
                blk.R_k.replace(np.inf, np.nan),
                marker="o",
                color=CAT[i],
                lw=2,
                label=pt,
            )
            tb = thr[(thr.point == pt) & (thr.family == fam) & (thr["size"] >= 1)]
            for s in tb.itertuples():
                face = CAT[i] if s.type == "FAILS_DECLARED_SECURITY" else "white"
                ax.scatter(
                    s.size + (0.08 if i else -0.08),
                    s.r if np.isfinite(s.r) else np.nan,
                    s=16,
                    facecolors=face,
                    edgecolors=CAT[i],
                    lw=0.8,
                    zorder=3,
                )
        ax.set_title(FAM[fam], fontsize=9)
        ax.set_xlabel("maximum portfolio size k", fontsize=8)
        ax.set_xticks([1, 2, 3, 4])
        style(ax)
        if j == 0:
            ax.set_ylabel("R_k = min_{|S|<=k} r_S")
            ax.legend(fontsize=8)
    fig.suptitle(
        "Resilience radii by portfolio size. Line: R_k. Dots: r_S per portfolio "
        "(open = scope-censored, OUTSIDE_MODEL_SCOPE; filled = in-scope failure; "
        "r_S = 0: transversely unstable)",
        fontsize=9,
    )
    save(fig, "FIG3_R_k", rk)


# ------------------------------------------------------------------- FIG4 --
def fig4():
    import yaml
    from _fc import CONFIGS

    cfg = yaml.safe_load(
        (CONFIGS / "final_nonlinear_composability_v1.yaml").read_text(encoding="utf-8")
    )
    thr = pd.read_csv(out_dir("FC05_nonlinear_thresholds") / "FC05_thresholds.csv")
    fam, pt = "D1", "P4"
    rho = cfg["planning"]["rho_star"][fam]
    blk = thr[(thr.point == pt) & (thr.family == fam)].copy()

    def status(r):
        if r.alpha_perp > 0:
            return "transversely unstable (spectral hyperedge)"
        if np.isfinite(r.r) and r.r <= rho:
            return (
                "fails in scope (nonlinear)"
                if r.type == "FAILS_DECLARED_SECURITY"
                else "not observable: model scope exceeded"
            )
        return "tolerates rho*"

    blk["status_at_rho"] = [status(r) for r in blk.itertuples()]
    colors = {
        "transversely unstable (spectral hyperedge)": "#7a1c1c",
        "fails in scope (nonlinear)": "#e67e22",
        "not observable: model scope exceeded": "#9aa3ad",
        "tolerates rho*": "#dfe9f3",
    }
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    by_size = {}
    for r in blk.itertuples():
        by_size.setdefault(r.size, []).append(r)
    pos = {}
    for sz, items in by_size.items():
        for i, r in enumerate(sorted(items, key=lambda q: q.subset)):
            x = (i - (len(items) - 1) / 2) * 1.6
            pos[r.subset] = (x, sz)
            ax.scatter(
                x,
                sz,
                s=900,
                c=colors[r.status_at_rho],
                edgecolors=INK,
                lw=0.8,
                zorder=3,
            )
            ax.text(
                x,
                sz,
                r.subset.replace("+", "\n"),
                ha="center",
                va="center",
                fontsize=6.5,
                zorder=4,
            )
    for s, (x, y) in pos.items():
        mem = set() if s == "BASE" else set(s.split("+"))
        for t, (x2, y2) in pos.items():
            mt = set() if t == "BASE" else set(t.split("+"))
            if y2 == y + 1 and mem < mt:
                ax.plot([x, x2], [y, y2], color=GRID, lw=0.8, zorder=1)
    for lbl, c in colors.items():
        ax.scatter([], [], s=80, c=c, edgecolors=INK, label=lbl)
    ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    ax.set_yticks([0, 1, 2, 3, 4])
    ax.set_ylim(-0.45, 4.5)
    ax.set_ylabel("portfolio size")
    ax.set_xticks([])
    ax.set_title(
        f"Spectral H_RHP_perp vs nonlinear status at rho* = {rho:g} MW (D1, {pt})",
        fontsize=9,
    )
    for sname in ("top", "right", "bottom"):
        ax.spines[sname].set_visible(False)
    save(
        fig,
        "FIG4_H_vs_HNL",
        blk[["subset", "size", "alpha_perp", "r", "type", "status_at_rho"]],
    )


# ------------------------------------------------------------------- FIG5 --
def fig5():
    thr = pd.read_csv(out_dir("FC05_nonlinear_thresholds") / "FC05_thresholds.csv")
    tr = out_dir("FC05_nonlinear_thresholds/traces")
    cand = thr[(thr.alpha_perp < 0) & np.isfinite(thr.r) & (thr.r > 0)]
    fails = cand[cand.type == "FAILS_DECLARED_SECURITY"]
    pick = (
        fails if len(fails) else cand[(cand.point == "P4") & (cand.family == "D1")]
    ).sort_values("size")
    r = pick.iloc[-1] if len(pick) else cand.iloc[0]
    stem = f"{r.point}_{r.family}_{r.subset}".replace("+", "p")
    below = pd.read_csv(tr / f"{stem}_below.csv.gz")
    above = pd.read_csv(tr / f"{stem}_above.csv.gz")
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.6))
    for data, name, c in (
        (below, "just below r_S", CAT[0]),
        (above, "just above r_S", CAT[1]),
    ):
        axes[0].plot(data.t, data.D, color=c, lw=1.4, label=name)
        axes[1].plot(
            data.t, data.gfl_vmin if "gfl_vmin" in data else data.vmin, color=c, lw=1.4
        )
        axes[1].plot(
            data.t,
            data.gfl_vmax if "gfl_vmax" in data else data.vmax,
            color=c,
            lw=1.0,
            ls="--",
        )
        axes[2].plot(data.t, data.pss_max, color=c, lw=1.4)
    axes[1].axhline(0.9, color=MUTED, ls=":", lw=1)
    axes[1].axhline(1.1, color=MUTED, ls=":", lw=1)
    axes[2].axhline(0.1, color=MUTED, ls=":", lw=1)
    axes[0].set_ylabel("deviation D(t)")
    axes[1].set_ylabel("GFL terminal |V| min / max [pu]")
    axes[2].set_ylabel("max |PSS output| [pu]")
    for ax in axes:
        ax.set_xlabel("time [s]")
        style(ax)
    axes[0].set_yscale("log")
    axes[0].legend(fontsize=8, loc="lower left")
    fig.subplots_adjust(wspace=0.35)
    fig.suptitle(
        f"{r.point}, {r.family}, portfolio {r.subset}: r_S = {r.r:.3g} ({r.type}). "
        "Dotted lines: model-scope guards (documented / declared limits)",
        fontsize=9,
    )
    save(
        fig,
        "FIG5_threshold_traces",
        pd.concat([below.assign(case="below"), above.assign(case="above")]),
    )


# ------------------------------------------------------------------- FIG7 --
def fig7():
    e = pd.read_csv(out_dir("FC07_second_order_curvature") / "FC07_errors.csv")
    case = "critical triple 30+33+35 P4"
    blk = e[e.case == case]
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8))
    obs = [
        ("err_freq_hz", "COI frequency [Hz]"),
        ("err_v16_pu", "|V| bus 16 [pu]"),
        ("err_i30_pu", "current at bus 30 [pu]"),
    ]
    for j, (col, name) in enumerate(obs):
        ax = axes[j]
        for i, (mdl, b) in enumerate(blk.groupby("model")):
            ax.loglog(b.eps_mw, b[col], marker="o", color=CAT[i], lw=1.8, label=mdl)
        x = np.array([blk.eps_mw.min(), blk.eps_mw.max()])
        ref = blk[blk.model.str.startswith("M0")][col].iloc[0]
        ax.loglog(x, ref * (x / x[0]) ** 2, color=MUTED, ls=":", lw=1, label="slope 2")
        ax.loglog(
            x, ref * 1e-2 * (x / x[0]) ** 3, color=MUTED, ls="--", lw=1, label="slope 3"
        )
        ax.set_title(name, fontsize=9)
        ax.set_xlabel("pulse amplitude eps [MW]")
        style(ax)
        if j == 0:
            ax.set_ylabel("max |model − full DAE|")
            handles, names = ax.get_legend_handles_labels()
            fig.legend(
                handles,
                names,
                fontsize=7,
                loc="lower center",
                ncol=6,
                bbox_to_anchor=(0.5, -0.06),
                frameon=False,
            )
    fig.suptitle(
        "Second-order reduced response: linear (M0), full (M1), no KCL curvature "
        f"(M2; M2t with transverse forcing) — {case}",
        fontsize=9,
    )
    save(fig, "FIG7_curvature_scaling", e)


# ------------------------------------------------------------------- FIG8 --
def fig8():
    s = json.loads(
        (out_dir("FC12_planning") / "FC12_summary.json").read_text(encoding="utf-8")
    )
    plans = pd.DataFrame(s["plans"])
    p3 = pd.DataFrame(s["P3"])
    fig, axes = plt.subplots(
        1, 2, figsize=(12.5, 4.0), gridspec_kw={"width_ratios": [1.6, 1]}
    )
    lbls, vals, cols = [], [], []
    for r in plans.itertuples():
        mw = getattr(r, "replaced_pg_mw", np.nan)
        lbls.append(f"{r.point}\n{r.constraint}")
        vals.append(mw if isinstance(mw, float) and np.isfinite(mw) else 0.0)
        cols.append(
            CAT[0]
            if "P1" in r.constraint
            else (CAT[5] if "SCOPE" in str(getattr(r, "status", "")) else CAT[1])
        )
    axes[0].barh(range(len(vals)), vals, color=cols)
    for i, r in enumerate(plans.itertuples()):
        st = str(getattr(r, "status", ""))
        txt = getattr(r, "portfolio", "") if "SCOPE" not in st else "out of model scope"
        axes[0].text(max(vals[i], 5), i, f" {txt}", va="center", fontsize=7)
    axes[0].set_yticks(range(len(lbls)), lbls, fontsize=6.5)
    axes[0].set_xlabel(
        "replaced active dispatch of the optimal any-order-safe portfolio [MW]"
    )
    style(axes[0])
    y = range(len(p3))
    fam = p3.family.fillna("").astype(str)
    p3_labels = [
        "transversely stable\n(small-signal)"
        if not f
        else f"{f}: stays inside the declared\nenvelope at {r.requirement.split()[-1]}"
        for f, r in zip(fam, p3.itertuples(), strict=True)
    ]
    colors = [CAT[2] if not f else CAT[5] for f in fam]
    axes[1].barh(list(y), p3.condenser_mva.fillna(0), color=colors)
    axes[1].set_yticks(list(y), p3_labels, fontsize=7)
    axes[1].set_xlabel("minimum damped-condenser rating [MVA] (P4 flagship)")
    for i, r in enumerate(p3.itertuples()):
        done = np.isfinite(r.condenser_mva)
        axes[1].text(
            (r.condenser_mva if done else 0) + 8,
            i,
            f"{r.condenser_mva:.0f} MVA" if done else "not reached at 25 % (1068 MVA)",
            va="center",
            fontsize=7,
        )
    axes[1].set_xlim(0, 1150)
    axes[1].set_title(
        "grey: envelope (guard) requirement, not a stability requirement",
        fontsize=7,
        color=MUTED,
    )
    style(axes[1])
    fig.subplots_adjust(wspace=0.55)
    fig.suptitle(
        "Planning: spectral-only vs nonlinear-safe vs support-assisted "
        "(MW and MVA on separate axes)",
        fontsize=9,
    )
    save(fig, "FIG8_planning", plans)


FIGS = {
    "FIG1": fig1,
    "FIG2": fig2,
    "FIG3": fig3,
    "FIG4": fig4,
    "FIG5": fig5,
    "FIG6": fig6,
    "FIG7": fig7,
    "FIG8": fig8,
}


def main(argv) -> int:
    names = argv or list(FIGS)
    for n in names:
        try:
            FIGS[n]()
            print(n, "ok")
        except FileNotFoundError as e:  # data not produced yet
            print(n, "skipped:", e)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
