from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
R = ROOT / "results"
OUT = ROOT / "figures" / "PD39_confirmatory"
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.titlesize": 11, "figure.dpi": 130})

def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.png", dpi=220, bbox_inches="tight")
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.svg", bbox_inches="tight")
    plt.close(fig)

def disc():
    return pd.read_csv(R / "PD39_CLASSICAL_BASELINES.csv")

def portfolio_landscape():
    d = disc(); fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.scatter(d.converted_mw, d.m_9, s=13, c=d.cardinality, cmap="viridis", alpha=.75, edgecolor="none")
    v8 = d[d.portfolio == "30;32;33;34;35;36;37;38"]
    p78 = d[d.cardinality == 7]
    ax.scatter(p78.converted_mw, p78.m_9, facecolors="none", edgecolors="black", s=45, label="7-of-8")
    ax.scatter(v8.converted_mw, v8.m_9, marker="*", c="crimson", s=110, label="V8")
    ax.axhline(.05, color="tab:orange", ls="--", label="robust target")
    ax.axhline(0, color="black", lw=.8, label="true-stability boundary")
    ax.set(xlabel="Converted MW", ylabel="Discovery m₉ = min(-α)", title="F1 — Complete portfolio landscape")
    ax.legend(frameon=False, ncol=3, fontsize=8); save(fig, "F1_complete_portfolio_landscape")

def blocker_lattice():
    d = pd.read_csv(R / "PD39_PORTFOLIO_BLOCKER_STRUCTURE.csv")
    g = d.groupby("cardinality").agg(n=("portfolio", "size"), true_blockers=("minimal_true_blocker", "sum"), robust_blockers=("minimal_0_05_blocker", "sum"))
    fig, ax = plt.subplots(figsize=(6.8, 4.0)); x=np.arange(len(g))
    ax.bar(x-.18, g.true_blockers, .36, label="true stability blockers", color="crimson")
    ax.bar(x+.18, g.robust_blockers, .36, label="0.05 robustness blockers", color="tab:orange")
    ax.set_xticks(x, g.index); ax.set_xlabel("Replacement cardinality"); ax.set_ylabel("Minimal blockers")
    ax.set_title("F2 — 255 + 1 structure (exact finite benchmark)"); ax.legend(frameon=False); save(fig, "F2_255_plus_1_structure")

def modal_plane():
    d = pd.read_csv(R / "PD39_7OF8_TO_8OF8_MODAL_ANALYSIS.csv")
    d["im_hz"] = d.critical_frequency_hz
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for pid, q in d.groupby("portfolio_id"):
        if pid == "V8": ax.scatter(q.alpha, q.im_hz, c=q.alpha, cmap="coolwarm", s=48, marker="*", label="V8")
        else: ax.scatter(q.alpha, q.im_hz, s=30, alpha=.75, label=pid.replace("missing_", "missing "))
    ax.axvline(0, color="black", lw=.8); ax.axvline(-.05, color="tab:orange", ls="--")
    ax.set(xlabel="Critical α (s⁻¹)", ylabel="Critical frequency (Hz)", title="F3 — Critical-mode plane")
    ax.legend(frameon=False, ncol=3, fontsize=7); save(fig, "F3_mode_transition_eigenvalue_plane")

def participation():
    d = pd.read_csv(R / "PD39_7OF8_TO_8OF8_MODAL_ANALYSIS.csv")
    picks=[]
    for pid in ["missing_37", "missing_34", "V8"]:
        q=d[(d.portfolio_id==pid) & (d.scenario=="nominal")]
        if len(q): picks.append(q.iloc[0])
    fig, axes=plt.subplots(1, len(picks), figsize=(11,4), sharey=False)
    if len(picks)==1: axes=[axes]
    for ax,row in zip(axes,picks):
        vals=[]; labs=[]
        for token in str(row.critical_participation).split("|")[:8]:
            lab,val=token.rsplit(":",1); labs.append(lab.replace("VIndex(","").replace(")","")); vals.append(float(val))
        order=np.argsort(vals)[::-1]; ax.barh(np.array(labs)[order][::-1], np.array(vals)[order][::-1], color="slateblue")
        ax.set_title(row.portfolio_id); ax.set_xlabel("participation")
    fig.suptitle("F4 — Common state participation (nominal critical mode)"); save(fig, "F4_mode_participation")

def radius_scatter():
    p=R/"PD39_STRUCTURED_RADIUS_SEARCH.csv"
    if not p.exists(): return
    d=pd.read_csv(p); q=d[(d.target=="rho_0") & (d.status=="boundary_found")]
    if q.empty: return
    best=q.loc[q.groupby("selection_id").radius_inf.idxmin()]
    fig,ax=plt.subplots(figsize=(6.8,4.2))
    alphas={"V8":-.001444778, "missing_30":-.130545972, "missing_32":-.1380, "missing_33":-.1383, "missing_34":-.0858, "missing_35":-.1383, "missing_36":-.1381, "missing_37":-.1383, "missing_38":-.1380}
    for _,r in best.iterrows():
        x=-alphas.get(r.selection_id, np.nan); ax.scatter(x,r.radius_inf,s=60,label=r.selection_id)
        ax.annotate(r.selection_id,(x,r.radius_inf),textcoords="offset points",xytext=(4,4),fontsize=8)
    ax.set(xlabel="Nominal mα = −α", ylabel="Smallest boundary found ρ₀ (normalized)", title="F5 — Nominal margin vs structured radius")
    ax.legend(frameon=False, fontsize=7, ncol=2); save(fig,"F5_nominal_margin_vs_structured_radius")
    fig,ax=plt.subplots(figsize=(6.8,4.2));
    q2=d[(d.target=="rho_0.05") & (d.status=="boundary_found")]
    if q2.empty: ax.text(.5,.5,"No positive ρ₀.₀₅ boundary: tested cases were already below target or uncrossed",ha="center",va="center",wrap=True)
    else: ax.scatter(q2.radius_inf, q2.radius_inf, s=25)
    ax.set_title("F6 — ρ₀.₀₅ vs nominal margin (available boundaries)"); save(fig,"F6_rho005_vs_nominal_margin")

def repair_front():
    fig,ax=plt.subplots(figsize=(7,4.2))
    c=pd.read_csv(R/"PD39_V8_CONTROLLER_REPAIR_SUMMARY.csv")
    ax.scatter(c.delta_pll.abs()+c.delta_xf.abs()+c.delta_cc.abs(), c.worst_margin, s=12, alpha=.65, label="controller grid")
    l=R/"PD39_V8_LINE_REPAIR_SUMMARY.csv"
    if l.exists():
        q=pd.read_csv(l); ax.scatter(q.minimum_gamma-1,q.max_gamma_worst_margin,s=35,label="single-line γ=1.25")
        for _,r in q[q.robust_at_1_25].iterrows(): ax.annotate(f"L{int(r.line)}",(r.minimum_gamma-1,r.max_gamma_worst_margin),fontsize=8)
    ax.axhline(.05,color="tab:orange",ls="--"); ax.set(xlabel="Effort proxy",ylabel="Worst discovery margin",title="F9 — V8 repair front")
    ax.legend(frameon=False); save(fig,"F9_v8_repair_front")

def pareto():
    d=disc(); fig,ax=plt.subplots(figsize=(7,4.2)); q=d[d.cardinality>0]
    sc=ax.scatter(q.converted_mw,q.m_9,c=q.nominal_alpha,cmap="coolwarm",s=20)
    v8=q[q.cardinality==8]; ax.scatter(v8.converted_mw,v8.m_9,marker="*",c="black",s=100,label="V8")
    ax.axhline(.05,color="tab:orange",ls="--"); ax.set(xlabel="Replaced MW",ylabel="Discovery robust margin m₉",title="F10 — IBR penetration vs robust margin")
    fig.colorbar(sc,ax=ax,label="nominal α (s⁻¹)"); ax.legend(frameon=False); save(fig,"F10_ibr_penetration_pareto")

def planner():
    p=R/"PD39_PLANNER_SELECTIONS.csv"
    if not p.exists(): return
    d=pd.read_csv(p); q=d[d.status=="selected"]
    fig,ax=plt.subplots(figsize=(8,4.5));
    for strat,g in q.groupby("strategy"):
        ax.plot(g.target, g.nominal_alpha, marker="o", label=strat.replace("_"," "))
    ax.axhline(-.05,color="tab:orange",ls="--"); ax.set(ylabel="Nominal α",title="F11 — Planner selections on discovery screen"); ax.legend(frameon=False,fontsize=7); save(fig,"F11_planner_comparison")

def tds():
    files=list(R.glob("PD39_TDS_TRACE_*.csv"))
    if not files: return
    fig,ax=plt.subplots(figsize=(7.4,4.4))
    for f in files:
        q=pd.read_csv(f); ax.plot(q.time,q.frequency_spread_hz,label=f.stem.replace("PD39_TDS_TRACE_", ""))
    ax.set(xlabel="Time (s)",ylabel="Frequency spread (Hz)",title="F12 — TDS common-disturbance frequency spread"); ax.legend(frameon=False,fontsize=7); save(fig,"F12_tds_hero")

def summary():
    fig,ax=plt.subplots(figsize=(10,6)); ax.axis("off")
    items=[("SG grid","25 / 50 / 75 / 100% SG→GFL"),("Nominal","stability screen"),("Hidden margin","structured perturbation radius"),("Weakness","dynamic node / link"),("Co-design","control + network + topology"),("Validation","24 holdout + common TDS")]
    y=.86
    for i,(a,b) in enumerate(items):
        ax.text(.08,y,a,fontsize=13,weight="bold",bbox=dict(boxstyle="round,pad=.35",fc="#e8f0f7",ec="#527a9e"),va="center")
        ax.text(.42,y,b,fontsize=12,va="center")
        if i<len(items)-1: ax.annotate("",xy=(.25,y-.08),xytext=(.25,y-.03),arrowprops=dict(arrowstyle="->",lw=1.5))
        y-=.13
    ax.set_title("F13 — PD39 confirmatory campaign summary",fontsize=16,pad=18); save(fig,"F13_one_page_summary")

portfolio_landscape(); blocker_lattice(); modal_plane(); participation(); radius_scatter(); repair_front(); pareto(); planner(); tds(); summary()
print(f"figures written to {OUT}")
