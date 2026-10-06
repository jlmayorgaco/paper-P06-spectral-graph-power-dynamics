from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"reports"/"experiment_F1"
FIG=OUT/"figures"
FIG.mkdir(parents=True,exist_ok=True)
df=pd.read_csv(OUT/"tables"/"TABLE_F01_modal_critical_damping.csv")
order=["underdamped","modal_critical","overdamped"]
fig,ax=plt.subplots(figsize=(8.2,4.5))
for scenario in order:
    rows=df[df.scenario==scenario]
    ax.plot(rows.nu,rows.analytic_decay_rate,marker="o",linewidth=1.8,label=scenario.replace("_"," "))
ax.set(xlabel=r"Modal stiffness $\nu_k$",ylabel="Asymptotic decay rate [s⁻¹]",
       title="Critical modal damping maximizes each isolated mode's decay")
ax.grid(alpha=.25); ax.legend(frameon=False)
for ext in ("png","pdf","svg"):
    fig.savefig(FIG/f"F1_01_modal_damping.{ext}",dpi=220,bbox_inches="tight")
plt.close(fig)
print(f"F1 figure written to {FIG}")
