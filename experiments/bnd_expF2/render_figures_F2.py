from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
F0=ROOT/"reports"/"experiment_F0"
F2=ROOT/"reports"/"experiment_F2"
FIG=F2/"figures"
FIG.mkdir(parents=True,exist_ok=True)
U=pd.read_csv(F0/"matrices"/"Lc_eigenvectors.csv").iloc[:,1:].to_numpy(float)
G=pd.read_csv(F0/"matrices"/"Yport_real.csv").iloc[:,1:].to_numpy(float)
B=pd.read_csv(F0/"matrices"/"Yport_imag.csv").iloc[:,1:].to_numpy(float)
GH=U.T@G@U
BH=U.T@B@U
metrics=pd.read_csv(F2/"tables"/"TABLE_F01_modal_structure.csv").set_index("metric").value
labels=["Reactive network", "Homogeneous GFL ports", "Heterogeneous full pencil"]
values=[float(metrics["susceptance_offdiag_ratio"]),
        float(metrics["homogeneous_closed_loop_port_offdiag_ratio"]),
        float(pd.read_csv(F2/"tables"/"TABLE_F02_heterogeneous_coupling.csv").set_index("metric").loc["heterogeneous_descriptor_offdiag_ratio","value"])]
fig,(ax0,ax1)=plt.subplots(1,2,figsize=(12,5.1),gridspec_kw={"width_ratios":[1.45,1]})
vmax=max(np.max(np.abs(GH)),np.max(np.abs(BH)))
im=ax0.imshow(BH,cmap="RdBu_r",vmin=-vmax,vmax=vmax)
ax0.set(title=r"$U^T\,\mathrm{Im}(Y_{port})\,U$",xlabel="Conductance graph mode",ylabel="Conductance graph mode")
ax0.set_xticks(range(len(BH)),range(1,len(BH)+1)); ax0.set_yticks(range(len(BH)),range(1,len(BH)+1))
fig.colorbar(im,ax=ax0,label="Reactive coupling [p.u.]")
ax1.barh(labels,values,color=["#c36a34","#276c91","#7356a6"])
ax1.axvline(.05,color="#b42318",linestyle="--",linewidth=1.5,label="predeclared 5% threshold")
ax1.set(xlim=(0,max(.35,1.15*max(values))),xlabel="Off-diagonal Frobenius ratio",
        title="Residual coupling in the graph basis")
ax1.grid(axis="x",alpha=.25); ax1.legend(frameon=False,loc="lower right")
fig.suptitle("The passive conductance graph does not diagonalize the detailed GFL network pencil",y=1.02)
fig.tight_layout()
for ext in ("png","pdf","svg"):
    fig.savefig(FIG/f"F2_01_modal_characteristic_structure.{ext}",dpi=220,bbox_inches="tight")
plt.close(fig)
print(f"F2 figure written to {FIG}")
