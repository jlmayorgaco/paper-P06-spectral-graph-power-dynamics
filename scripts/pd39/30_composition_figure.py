from pathlib import Path
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

static = pd.read_csv(R / "PD39_CLASSICAL_BASELINES.csv")
discovery = pd.read_csv(R / "pd39" / "portfolio_campaign" / "portfolio_scenario_results.csv")
portfolio_m9 = discovery.groupby("portfolio", as_index=False).agg(
    m9=("dynamic_margin", "min"), alpha_worst=("max_real", "max")
)
d = static.merge(portfolio_m9, on="portfolio", how="inner")

fig, ax = plt.subplots(figsize=(7.4, 4.6))
colors = {6: "#2878b5", 7: "#59a14f", 8: "#e15759"}
for k in (6, 7, 8):
    q = d[d.cardinality == k]
    ax.scatter(q.converted_mw, q.m9, s=32 if k < 8 else 115,
               color=colors[k], marker="*" if k == 8 else "o",
               alpha=.85, edgecolor="black" if k == 8 else "none",
               label=f"{k}/8")

pairs = pd.read_csv(R / "PD39_PENETRATION_COMPOSITION_PAIRS.csv")
pair_ids = ["k6_p2", "k6_p5", "k7_p2", "k7_p3"]
for pid in pair_ids:
    q = pairs[pairs.pair_id == pid].iloc[0]
    aa = d[d.portfolio == q.portfolio_a].iloc[0]
    bb = d[d.portfolio == q.portfolio_b].iloc[0]
    ax.plot([aa.converted_mw, bb.converted_mw], [aa.m9, bb.m9],
            color="0.35", lw=.8, alpha=.8)

ax.axhline(.05, color="tab:orange", ls="--", lw=1.0, label="0.05 target")
ax.axhline(0, color="black", lw=.8, label="true-stability boundary")
ax.set(xlabel="Converted dispatch (MW)", ylabel="m₉ = min over 9 scenarios (−α), s⁻¹",
       title="F14 — Penetration versus composition in the frozen discovery set")
ax.legend(frameon=False, ncol=3, fontsize=8)
fig.tight_layout()
for ext, kwargs in (("png", {"dpi": 220, "bbox_inches": "tight"}),
                    ("pdf", {"bbox_inches": "tight"}),
                    ("svg", {"bbox_inches": "tight"})):
    fig.savefig(OUT / f"F14_penetration_vs_composition.{ext}", **kwargs)
plt.close(fig)
print(f"figure written to {OUT}")
