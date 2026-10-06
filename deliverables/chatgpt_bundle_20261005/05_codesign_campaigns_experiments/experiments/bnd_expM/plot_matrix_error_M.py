"""Visualize physically aligned reduced-Jacobian discrepancies."""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import solve

from compare_models_M import pd_key, an_key


ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"reports"/"experiment_M"


def read(path:Path)->list[dict]:
    with path.open(newline="",encoding="utf-8") as fh:return list(csv.DictReader(fh))


def aligned_error(case:str)->np.ndarray:
    p=OUT/"matrices"/case
    ps=read(p/"PD_state_map.csv");ans=read(p/"AN_state_map.csv")
    keys=[pd_key(r["state_name"]) for r in ps]
    ankeys=[an_key(r) for r in ans]
    a=np.loadtxt(p/"PD_A.csv",delimiter=",")
    m=np.loadtxt(p/"PD_M.csv",delimiter=",")
    d=np.flatnonzero(np.diag(m)==1);z=np.flatnonzero(np.diag(m)==0)
    ar=a[np.ix_(d,d)]-a[np.ix_(d,z)]@solve(a[np.ix_(z,z)],a[np.ix_(z,d)])
    position={keys[i]:j for j,i in enumerate(d)}
    order=[position[ankeys[i]] for i,row in enumerate(ans) if row["differential"]=="true"]
    return ar[np.ix_(order,order)]-np.loadtxt(p/"AN_Ared.csv",delimiter=",")


def main()->None:
    cases=("all_SG","ExpG_candidate","ExpK_nominal")
    fig,axes=plt.subplots(1,3,figsize=(12,3.8),constrained_layout=True)
    for ax,case in zip(axes,cases):
        delta=aligned_error(case)
        im=ax.imshow(np.log10(np.abs(delta)+1e-14),vmin=-14,vmax=4,cmap="magma",origin="lower",aspect="auto")
        ax.set_title(case.replace("_"," "))
        ax.set_xlabel("Physical state column");ax.set_ylabel("Physical state row")
    fig.colorbar(im,ax=axes,label="log10 |A_PD − A_AN| (s⁻¹)",shrink=.82)
    fig.savefig(OUT/"FIG_M01_A_error_heatmap.png",dpi=170)
    print(OUT/"FIG_M01_A_error_heatmap.png")


if __name__=="__main__":main()
