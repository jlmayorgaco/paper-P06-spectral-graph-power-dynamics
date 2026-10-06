"""Plot measured local poles and SG dispatch from the fixed-rho architecture audit."""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"reports"/"experiment_M"


def rightmost(raw:str)->float:
    if not raw:return np.nan
    vals=[]
    for token in raw.split(";"):
        try:vals.append(complex(token.replace("im","j")).real)
        except ValueError:continue
    return max(vals) if vals else np.nan


def main()->None:
    with (OUT/"tables"/"TABLE_M14_fractional_rho_semantics.csv").open(newline="",encoding="utf-8") as fh:
        rows=list(csv.DictReader(fh))
    fig,axes=plt.subplots(2,2,figsize=(9,7),constrained_layout=True)
    for col,bus in enumerate((36,38)):
        sub=sorted((r for r in rows if int(r["bus"])==bus and r["status"]=="EVALUATED"),
                   key=lambda r:float(r["rho"]))
        x=[float(r["rho"]) for r in sub]
        axes[0,col].plot(x,[rightmost(r["SG_internal_poles"]) for r in sub],"o-",label="SG-dominant local poles")
        axes[0,col].plot(x,[rightmost(r["GFL_internal_poles"]) for r in sub],"s-",label="GFL-dominant local poles")
        axes[0,col].axhline(0,color="black",lw=.7)
        axes[0,col].set_title(f"Bus {bus}")
        axes[0,col].set_ylabel("Local open-loop rightmost Re(λ), s⁻¹")
        axes[0,col].legend(fontsize=7)
        axes[1,col].plot(x,[float(r["SG_P_MW"]) for r in sub],"o-",label="Initialized PD SG MW")
        axes[1,col].plot(x,[float(r["expected_frozen_SG_P_MW"]) for r in sub],"k--",label="Frozen share MW")
        axes[1,col].set_ylabel("SG active power (MW)")
        axes[1,col].set_xlabel("GFL fraction ρ")
        axes[1,col].legend(fontsize=7)
    fig.savefig(OUT/"FIG_M02_internal_poles_vs_rho.png",dpi=170)
    print(OUT/"FIG_M02_internal_poles_vs_rho.png")


if __name__=="__main__":main()
