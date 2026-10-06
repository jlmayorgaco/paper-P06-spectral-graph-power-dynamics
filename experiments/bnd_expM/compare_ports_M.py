"""Compare initialized PD bus port admittance with the unchanged analytic port."""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import solve


ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"reports"/"experiment_M"
TABLES=OUT/"tables"
FREQUENCIES=np.geomspace(1e-4,1e3,51)


def matrix(path:Path)->np.ndarray:
    return np.atleast_2d(np.loadtxt(path,delimiter=","))


def transfer(m:np.ndarray,a:np.ndarray,b:np.ndarray,c:np.ndarray,d:np.ndarray,s:complex)->np.ndarray:
    return c@solve(s*m-a,b)+d


def main()->None:
    summaries=[]; traces=[]
    with (ROOT/"reports"/"experiment_E"/"tables"/"TABLE_E03_initialized_ZIP_loads.csv").open(newline="",encoding="utf-8") as fh:
        loads={int(r["bus"]):complex(float(r["initialized_admittance_real"]),
                                     float(r["initialized_admittance_imag"])) for r in csv.DictReader(fh)}
    for case in ("all_SG","ExpG_candidate","ExpK_nominal"):
        p=OUT/"matrices"/case
        if not (p/"AN_A_local.csv").exists():
            continue
        with (p/"AN_state_map.csv").open(newline="",encoding="utf-8") as fh:
            states=list(csv.DictReader(fh))
        aan=matrix(p/"AN_A_local.csv")
        ban=matrix(p/"AN_B_port.csv")
        can=matrix(p/"AN_C_port.csv")
        dan=matrix(p/"AN_D_port.csv")
        for bus in (31,36,38,39):
            idx=np.array([int(row["index"])-1 for row in states if int(row["bus"])==bus and row["differential"]=="true"])
            if len(idx)==0:
                continue
            y=np.array([2*bus-2,2*bus-1])
            mp=matrix(p/f"PD_bus{bus}_M.csv")
            ap=matrix(p/f"PD_bus{bus}_A.csv")
            bp=matrix(p/f"PD_bus{bus}_B.csv")
            cp=matrix(p/f"PD_bus{bus}_C.csv")
            dp=matrix(p/f"PD_bus{bus}_D.csv")
            local=[]
            for w in FREQUENCIES:
                s=1j*w
                zan=transfer(np.eye(len(idx)),aan[np.ix_(idx,idx)],ban[np.ix_(idx,y)],
                    can[np.ix_(y,idx)],dan[np.ix_(y,y)],s)
                if bus in loads:
                    adm=loads[bus]
                    zan=zan+np.array([[adm.real,-adm.imag],[adm.imag,adm.real]])
                zpd=transfer(mp,ap,bp,cp,dp,s)
                ypd=np.linalg.inv(zpd)
                sign=1 if np.linalg.norm(ypd-zan)<np.linalg.norm(-ypd-zan) else -1
                err=np.linalg.norm(sign*ypd-zan)/max(np.linalg.norm(zan),1e-15)
                local.append(err)
                traces.append(dict(case=case,bus=bus,omega_rad_s=float(w),
                    relative_port_error=float(err),PD_admittance_sign=sign))
            summaries.append(dict(case=case,bus=bus,component="SG_or_GFL_bus_port",
                samples=len(local),max_relative_error=float(np.max(local)),
                median_relative_error=float(np.median(local)),p95_relative_error=float(np.quantile(local,0.95)),
                omega_max_rad_s=float(FREQUENCIES[int(np.argmax(local))]),
                sign_consistent=len(set(x["PD_admittance_sign"] for x in traces[-len(local):]))==1))
    TABLES.mkdir(parents=True,exist_ok=True)
    if summaries:
        with (TABLES/"TABLE_M15_port_transfer_identity.csv").open("w",newline="",encoding="utf-8") as fh:
            writer=csv.DictWriter(fh,fieldnames=list(summaries[0]));writer.writeheader();writer.writerows(summaries)
        with (TABLES/"TABLE_M15_port_transfer_trace.csv").open("w",newline="",encoding="utf-8") as fh:
            writer=csv.DictWriter(fh,fieldnames=list(traces[0]));writer.writeheader();writer.writerows(traces)
        plt.figure(figsize=(7,4))
        for case,bus in sorted({(r["case"],r["bus"]) for r in traces}):
            rows=[r for r in traces if r["case"]==case and r["bus"]==bus]
            plt.loglog([r["omega_rad_s"] for r in rows],[r["relative_port_error"] for r in rows],
                label=f"{case}, bus {bus}")
        plt.xlabel("Angular frequency (rad/s)");plt.ylabel("Relative port admittance error")
        plt.legend(fontsize=7);plt.grid(True,which="both",alpha=.3);plt.tight_layout()
        plt.savefig(OUT/"FIG_M03_GFL_port_error_vs_frequency.png",dpi=160)
    for row in summaries:print(row)


if __name__=="__main__":main()
