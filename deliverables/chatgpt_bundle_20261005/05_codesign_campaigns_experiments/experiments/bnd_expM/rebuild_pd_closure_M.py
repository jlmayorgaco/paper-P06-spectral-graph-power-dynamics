"""Reclose installed PD bus/line port operators with the library's tested sign."""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy.linalg import eigvals, solve
from scipy.optimize import linear_sum_assignment


ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"reports"/"experiment_M"
TABLES=OUT/"tables"


def mat(path:Path)->np.ndarray:
    return np.atleast_2d(np.loadtxt(path,delimiter=","))


def finite_poles(a:np.ndarray,m:np.ndarray)->np.ndarray:
    d=np.flatnonzero(np.diag(m)==1)
    z=np.flatnonzero(np.diag(m)==0)
    red=a[np.ix_(d,d)]-a[np.ix_(d,z)]@solve(a[np.ix_(z,z)],a[np.ix_(z,d)])
    return eigvals(red)


def main()->None:
    summaries=[]
    for case in ("all_SG","ExpG_candidate","ExpK_nominal"):
        p=OUT/"matrices"/case
        if not (p/"PD_Zbus_A.csv").exists():continue
        a=mat(p/"PD_Zbus_A.csv");m=mat(p/"PD_Zbus_M.csv")
        b=mat(p/"PD_Zbus_B.csv");c=mat(p/"PD_Zbus_C.csv");d=mat(p/"PD_Zbus_D.csv")
        network=mat(p/"PD_Ynw_D.csv")
        if (p/"PD_Yinj_D.csv").exists():
            network+=mat(p/"PD_Yinj_D.csv")
        correction=solve(np.eye(network.shape[0])-network@d,network@c)
        closed=a+b@correction
        poles=finite_poles(closed,m)
        with (p/"PD_closure_positive_poles.csv").open("w",newline="",encoding="utf-8") as fh:
            writer=csv.writer(fh);writer.writerow(("real","imag"))
            writer.writerows(zip(poles.real,poles.imag))
        direct=np.loadtxt(p/"PD_poles.csv",delimiter=",",skiprows=1)
        reference=direct[:,0]+1j*direct[:,1]
        if len(poles)==len(reference):
            row,col=linear_sum_assignment(np.abs(reference[:,None]-poles[None,:]))
            errors=np.abs(reference[row]-poles[col])
            summaries.append(dict(case=case,status="EVALUATED_POSITIVE_FEEDBACK",
                PD_poles=len(reference),closure_poles=len(poles),
                max_pole_error=float(errors.max()),median_pole_error=float(np.median(errors)),
                pass_1e_8=bool(errors.max()<1e-8)))
        else:
            summaries.append(dict(case=case,status="DIMENSION_MISMATCH",
                PD_poles=len(reference),closure_poles=len(poles),
                max_pole_error="",median_pole_error="",pass_1e_8=False))
    if summaries:
        TABLES.mkdir(parents=True,exist_ok=True)
        with (TABLES/"TABLE_M08_PD_open_loop_closure_identity.csv").open("w",newline="",encoding="utf-8") as fh:
            writer=csv.DictWriter(fh,fieldnames=list(summaries[0]));writer.writeheader();writer.writerows(summaries)
    for row in summaries:print(row)


if __name__=="__main__":main()
