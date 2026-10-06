from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"reports"/"experiment_F_program"
FIG=OUT/"figures"
FIG.mkdir(parents=True,exist_ok=True)


def save(fig,stem):
    for ext in ("png","pdf","svg"):
        fig.savefig(FIG/f"{stem}.{ext}",dpi=220,bbox_inches="tight")
    plt.close(fig)


# F00: pre-registered pipeline and the observed gate stop.
fig,ax=plt.subplots(figsize=(12,5.2))
ax.set_xlim(0,12); ax.set_ylim(0,5.2); ax.axis("off")
nodes={
    "F0\nphysical graph\nPASS":(1.1,3.7,"#d8f0e6"),
    "F1\nmodal damping\nPASS":(3.4,3.7,"#d8f0e6"),
    "F2\ndetailed modal test\nNON_MODAL":(5.8,3.7,"#ffe1d6"),
    "F3–F4\ngain law / graph controller\nBLOCKED":(8.5,3.7,"#eceff4"),
    "F6–F7\ncapacity / validation\nBLOCKED":(10.9,3.7,"#eceff4"),
    "F5\nexact $\\Gamma_k$ diagnostic\nDIAGNOSTIC ONLY":(5.8,1.1,"#fff0c2"),
}
for label,(x,y,color) in nodes.items():
    ax.text(x,y,label,ha="center",va="center",fontsize=10,
            bbox={"boxstyle":"round,pad=0.65","facecolor":color,"edgecolor":"#475467","linewidth":1.2})
def arrow(a,b,style="-",color="#475467"):
    ax.annotate("",xy=b,xytext=a,arrowprops={"arrowstyle":"-|>","lw":1.8,"linestyle":style,"color":color})
arrow((2.0,3.7),(2.55,3.7)); arrow((4.35,3.7),(4.85,3.7)); arrow((6.8,3.7),(7.45,3.7),"--")
arrow((9.55,3.7),(10.0,3.7),"--")
arrow((5.8,3.05),(5.8,1.85),"--",color="#a36b00")
ax.text(6.8,4.75,"F0–F7 experimental program",ha="center",fontsize=16,weight="bold",color="#172b4d")
ax.text(8.55,2.6,"No $K_p=h_p(L_c)$ or $K_i=h_i(L_c)$ claim",ha="center",fontsize=9,color="#b42318")
ax.text(5.8,0.35,"Reactive cross-mode coupling = 29.14% (gate: 5%)",ha="center",fontsize=9,color="#7a4a00")
save(fig,"FIG_F00_research_pipeline")

# F25: final gate evidence, with unrun performance kept visibly unreported.
fig,(ax0,ax1)=plt.subplots(1,2,figsize=(12,5.5),gridspec_kw={"width_ratios":[1.1,1]})
labels=["Reactive network", "Homogeneous GFL transfer", "Heterogeneous descriptor", "F2 threshold"]
values=[0.29138391,0.28079173,0.06141531,0.05]
colors=["#c36a34","#276c91","#7356a6","#b42318"]
ax0.barh(labels,values,color=colors)
ax0.set_xlim(0,.34); ax0.set_xlabel("Off-diagonal Frobenius ratio")
ax0.set_title("Modal coupling (threshold = 5%)")
ax0.grid(axis="x",alpha=.25)
ax1.axis("off")
summary=("F0 PASS  ·  Lc SPD: 0.00113–5.013 p.u.\n"
         "Kron port-current residual: 4.1×10⁻¹⁶\n"
         "F1 PASS  ·  modal critical law reproduced\n"
         "F2 STOP  ·  29.14% reactive cross-mode coupling\n\n"
         "No graph gain law was derived.\n"
         "F3–F4 and F6–F7 were not run.\n"
         "Exact Γₖ was retained as a coupled-mode diagnostic.\n"
         "HARD_CURRENT_LIMIT = NOT_MODELED")
ax1.text(.02,.94,"Program result: PARTIAL",transform=ax1.transAxes,fontsize=16,weight="bold",color="#172b4d",va="top")
ax1.text(.02,.75,summary,transform=ax1.transAxes,fontsize=12,va="top",linespacing=1.7,
         bbox={"boxstyle":"round,pad=.7","facecolor":"#f5f7fa","edgecolor":"#cbd5e1"})
fig.suptitle("Explicit graph-structured PLL design: pre-registered gate outcome",fontsize=15,weight="bold")
fig.tight_layout()
save(fig,"FIG_F25_final_poster_summary")
print(f"Final program figures written to {FIG}")
