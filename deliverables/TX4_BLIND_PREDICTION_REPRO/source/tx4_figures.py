"""Publication figures for the narrow TX4 blind-prediction campaign."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
ROOT = HERE.parents[5]
OUT = ROOT / "figures" / "tx4_blind_prediction"
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.titlesize": 11, "axes.labelsize": 10, "figure.dpi": 140})


def save(fig, stem):
    fig.tight_layout()
    fig.savefig(OUT / f"{stem}.png", dpi=260, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def main():
    res = ROOT / "results"
    v4 = pd.read_csv(res / "TX4_V4_BLIND_VS_FULL.csv")
    v9 = pd.read_csv(res / "TX4_V9_BLIND_VS_FULL.csv")
    sweep = pd.read_csv(res / "TX4_CONTEXTUAL_RETURN_SWEEP.csv")
    boundary = pd.read_csv(res / "TX4_CONTEXTUAL_RETURN_BOUNDARY.csv")
    proper = pd.read_csv(res / "TX4_15_PROPER_SUBSETS.csv")
    pair = pd.read_csv(res / "TX4_AGGREGATE_MATCHED_PAIRS.csv").iloc[0]
    tds = pd.read_csv(res / "TX4_TDS/TX4_TDS_FINAL.csv")

    # P1: candidate support and local-vs-collective decomposition.
    fig, ax = plt.subplots(figsize=(7.0, 3.8))
    xs = np.arange(4)
    ax.scatter(xs, np.zeros(4), s=500, c=["#277da1", "#f9844a", "#90be6d", "#f9c74f"], edgecolor="black", zorder=3)
    for x, b in zip(xs, [30, 33, 35, 37]):
        ax.text(x, 0, str(b), ha="center", va="center", color="white", weight="bold")
    for a in xs:
        for b in xs:
            if a < b:
                ax.plot([a, b], [0, 0], color="#777", lw=0.7, alpha=0.35)
    ax.text(1.5, 0.23, "H4 = {30, 33, 35, 37}", ha="center", weight="bold", fontsize=13)
    ax.text(1.5, -0.20, "four local SG→GFL substitutions; collective closure tested", ha="center")
    ax.set_xlim(-0.5, 3.5); ax.set_ylim(-0.35, 0.4); ax.set_xticks(xs, ["bus 30", "bus 33", "bus 35", "bus 37"]); ax.set_yticks([])
    ax.set_title("P1 — TX4 candidate support")
    for s in ax.spines.values(): s.set_visible(False)
    save(fig, "P1_H4_support")

    # P2: 16-portfolio lattice colored by full-order verdict.
    fig, ax = plt.subplots(figsize=(7.0, 4.3))
    for _, r in v4.iterrows():
        x = int(r.cardinality); y = int(r.portfolio == "30+33+35+37")
        color = "#d95f02" if r.status == "UNSTABLE" else "#1b9e77"
        ax.scatter(x + (np.random.default_rng(0).uniform(-0.23, 0.23) if x not in (0,4) else 0), y + np.random.default_rng(hash(r.portfolio) % (2**32)).uniform(-0.03,0.03), s=80, c=color, alpha=.85, edgecolor="white")
    ax.axhline(0.5, color="#888", lw=.7); ax.set_xticks(range(5)); ax.set_xlabel("replacement cardinality"); ax.set_yticks([0,1], ["proper / stable", "H4 / unstable"]); ax.set_title("P2 — 15 safe proper portfolios, one unsafe full portfolio")
    save(fig, "P2_16_portfolio_lattice")

    # P3: blind versus full verdicts.
    fig, ax = plt.subplots(figsize=(7.0, 4.3))
    for x, (df, label) in enumerate([(v4, "V4 (16)"), (v9, "V9 (512)")]):
        cm = pd.crosstab(df.status, df.predicted_verdict.replace({"STABLE_PREDICTED":"STABLE", "UNSTABLE_PREDICTED":"UNSTABLE"}))
        for j, truth in enumerate(["STABLE", "UNSTABLE"]):
            for i, pred in enumerate(["STABLE", "UNSTABLE"]):
                n = int(cm.get(pred, pd.Series()).get(truth, 0)) if hasattr(cm.get(pred, pd.Series()), 'get') else 0
                ax.text(x + i*.25, j*.28, str(n), ha="center", va="center", fontsize=16, weight="bold", color="#1b9e77" if pred==truth else "#d95f02")
        ax.text(x+.125, -.18, label, ha="center")
    ax.set_xlim(-.35, 1.6); ax.set_ylim(-.35,.7); ax.set_xticks([]); ax.set_yticks([0,.28], ["full stable", "full unstable"]); ax.set_title("P3 — blind reduced closure versus full-order verdicts")
    ax.text(.05,.62,"columns within each panel: predicted stable / unstable",fontsize=8)
    save(fig, "P3_blind_vs_full")

    # P4: full alpha and return distance on g path.
    fig, ax = plt.subplots(figsize=(7.0, 4.3)); ax2 = ax.twinx()
    ax.plot(sweep.g, sweep.alpha, "o-", color="#277da1", label="full-order α")
    ax.axhline(0, color="black", lw=.7); ax.axvline(.207681405, color="#d95f02", ls="--", label="boundary")
    ax2.semilogy(sweep.g, sweep.return_distance_to_unity, "s--", color="#90be6d", label="min |μ−1|")
    ax.set_xlabel("g"); ax.set_ylabel("α⊥ [s⁻¹]"); ax2.set_ylabel("contextual-return distance"); ax.set_title("P4 — controller boundary and contextual return")
    h1,l1=ax.get_legend_handles_labels(); h2,l2=ax2.get_legend_handles_labels(); ax.legend(h1+h2,l1+l2,loc="best")
    save(fig, "P4_alpha_return")

    # P5: collective closure singularity versus local regularity.
    fig, ax = plt.subplots(figsize=(7.0, 4.1)); x=np.arange(len(boundary)); w=.34
    ax.semilogy(x-w/2, boundary.collective_sigma_min, "o", color="#d95f02", label="collective σmin")
    ax.semilogy(x+w/2, boundary.local_sigma_min_at_device, "s", color="#1b9e77", label="local-factor σmin")
    ax.set_xticks(x, boundary.device_bus.astype(str)); ax.set_xlabel("device removed from H4 context"); ax.set_ylabel("singular value"); ax.set_title("P5 — collective closure closes while local factors remain regular"); ax.legend()
    save(fig, "P5_local_collective")

    # P6: all proper subsets remain separated at H4 root frequency.
    fig, ax = plt.subplots(figsize=(7.0, 4.4)); proper2=proper[proper.subset!="30+33+35+37"].copy(); proper2["label"]=proper2.subset
    ax.barh(np.arange(len(proper2)), proper2.closure_sigma_min_at_H4_root, color="#277da1")
    ax.set_yticks(np.arange(len(proper2)), proper2.label, fontsize=7); ax.set_xlabel("σmin(I+Q_S) at H4 root frequency"); ax.set_title("P6 — proper-subset closure separation")
    save(fig, "P6_proper_closure")

    # P7: selected aggregate-matched counterexample.
    fig, ax = plt.subplots(figsize=(7.0, 4.1)); names=[str(pair.portfolio_A),str(pair.portfolio_B)]; vals=[float(pair.alpha_A),float(pair.alpha_B)]; colors=["#d95f02","#1b9e77"]
    ax.bar(names, vals, color=colors); ax.axhline(0,color="black",lw=.7); ax.set_ylabel("α⊥ [s⁻¹]"); ax.tick_params(axis="x",labelrotation=25); ax.set_title(f"P7 — same cardinality, 0.25% normalized MW/MVA mismatch")
    save(fig, "P7_aggregate_pair")

    # P8: TDS envelope proxy from the frozen traces.
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    for _, r in tds.iterrows():
        key=r.case.replace(" ","_").replace("+","p").replace("=","_")
        paths=list((res/"TX4_TDS").glob(f"trace_{key}_D2.csv.gz"))
        if paths:
            tr=pd.read_csv(paths[0]); ax.plot(tr.t, np.abs(tr.signal-tr.signal.mean()), label=f"{r['case']} ({r['outcome']})")
    ax.set_yscale("log"); ax.set_xlabel("time [s]"); ax.set_ylabel("|speed-difference deviation|"); ax.set_title("P8 — common +2% bus-20 load-pulse TDS"); ax.legend(fontsize=7)
    save(fig, "P8_TDS_final")

    # P9: benchmark accuracy and runtime.
    fig, ax = plt.subplots(figsize=(7.0, 4.0)); labels=["V4", "V9"]; acc=[float(v4.verdict_correct.mean()), float(v9.verdict_correct.mean())]; fs=[int(v4.false_safe.sum()), int(v9.false_safe.sum())]
    xx=np.arange(2); ax.bar(xx-.18,acc,width=.36,color="#277da1",label="accuracy"); ax2=ax.twinx(); ax2.bar(xx+.18,fs,width=.36,color="#d95f02",label="false-safe count"); ax.set_xticks(xx,labels); ax.set_ylim(0,1.05); ax.set_ylabel("accuracy"); ax2.set_ylabel("false-safe count"); ax.set_title("P9 — V4 exact flagship prediction; V9 negative scalability result")
    save(fig, "P9_V9_benchmark")

    # P10: one-page graphical abstract.
    fig, ax = plt.subplots(figsize=(8.2, 4.7)); ax.axis("off")
    boxes=[(0.04,.56,.21,.25,"SG grid\nall-SG baseline"),(0.29,.56,.21,.25,"local ΔY₃₀…ΔY₃₇\n4 local models"),(0.54,.56,.21,.25,"exact closure\nI + Q_H"),(0.79,.56,.17,.25,"15 SAFE\n1 UNSAFE")]
    for x,y,w,h,txt in boxes: ax.add_patch(plt.Rectangle((x,y),w,h,facecolor="#e8f1f5",edgecolor="#277da1",lw=1.5)); ax.text(x+w/2,y+h/2,txt,ha="center",va="center",weight="bold")
    for x in [.25,.50,.75]: ax.annotate("",xy=(x,.685),xytext=(x-.025,.685),arrowprops=dict(arrowstyle="->",lw=1.5,color="#555"))
    ax.text(.5,.31,"det(I + Q_H)=0  ⇔  μ(R_i|R)=1",ha="center",fontsize=15,weight="bold",color="#d95f02")
    ax.text(.5,.17,"g = 0.03625: unstable   →   g = 0.25: stable\nnonlinear TDS agrees; V9 benchmark exposes scope limits",ha="center",fontsize=11)
    ax.set_title("P10 — Safe alone, unsafe together: blind closure + contextual return",weight="bold",fontsize=14)
    save(fig, "P10_graphical_abstract")


if __name__ == "__main__":
    main()
