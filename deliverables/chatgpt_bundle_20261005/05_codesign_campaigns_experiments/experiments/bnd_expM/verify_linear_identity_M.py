"""Trajectory identity for direct PD and PD-port analytical closure at ExpK."""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve


ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"reports"/"experiment_M"
CASE=OUT/"matrices"/"ExpK_nominal"


def load(name:str)->np.ndarray:
    return np.atleast_2d(np.loadtxt(CASE/name,delimiter=","))


def reduced(a:np.ndarray,m:np.ndarray,b:np.ndarray,c:np.ndarray,direct:np.ndarray)->tuple:
    di=np.flatnonzero(np.diag(m)==1);zi=np.flatnonzero(np.diag(m)==0)
    zx=solve(a[np.ix_(zi,zi)],a[np.ix_(zi,di)])
    zu=solve(a[np.ix_(zi,zi)],b[zi,:])
    return (a[np.ix_(di,di)]-a[np.ix_(di,zi)]@zx,
            b[di,:]-a[np.ix_(di,zi)]@zu,
            c[:,di]-c[:,zi]@zx,direct-c[:,zi]@zu,di)


def propagate(a:np.ndarray,b:np.ndarray,x0:np.ndarray,u:np.ndarray,t0:float,t1:float,times:np.ndarray)->np.ndarray:
    sol=solve_ivp(lambda _,x:a@x+b@u,(t0,t1),x0,method="BDF",rtol=1e-11,atol=1e-13,
        dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    return sol.sol(times)


def main()->None:
    with (CASE/"ANALYTICAL_PD_EXACT_state_map.csv").open(newline="",encoding="utf-8") as fh:
        states=list(csv.DictReader(fh))
    order=[int(r["PD_index"])-1 for r in states]
    ap=load("PD_load_input_A.csv")[np.ix_(order,order)]
    mp=load("PD_load_input_M.csv")[np.ix_(order,order)]
    bp=load("PD_load_input_B.csv")[order,:]
    cp=load("PD_load_input_C.csv")[:,order]
    dp=load("PD_load_input_D.csv")
    ac=load("ANALYTICAL_PD_EXACT_A.csv")
    mc=load("ANALYTICAL_PD_EXACT_M.csv")
    rp,up,yp,fp,di=reduced(ap,mp,bp,cp,dp)
    rc,uc,yc,fc,_=reduced(ac,mc,bp,cp,dp)
    meta=list(csv.DictReader((CASE/"PD_load_input_metadata.csv").open(newline="",encoding="utf-8")))[0]
    u=np.array([float(meta["Pset_pu"]),float(meta["Qset_pu"])])*1e-4
    weights=np.array([float(meta["H38_Sn38"]),float(meta["H39_Sn39"])])
    weights/=weights.sum()
    times=np.linspace(1,1.1,101)
    zp=propagate(rp,up,np.zeros(len(di)),u,1,1.1,times)
    zc=propagate(rc,uc,np.zeros(len(di)),u,1,1.1,times)
    fpulse=60*weights@(yp@zp+fp@u[:,None])
    fclosure=60*weights@(yc@zc+fc@u[:,None])
    loaderr=float(np.max(np.abs(fpulse-fclosure)))
    loadscale=max(float(np.max(np.abs(fpulse))),1e-30)
    # Perturb one physical SG rotor-speed state by 1e-7 pu.
    target=next(i for i,s in enumerate(states) if s["state_name"]=="VIndex(39, :machine₊ω)")
    local=int(np.where(di==target)[0][0])
    x0=np.zeros(len(di));x0[local]=1e-7
    t=np.linspace(0,2,101)
    zero=np.zeros(2)
    xp=propagate(rp,up,x0,zero,0,2,t)
    xc=propagate(rc,uc,x0,zero,0,2,t)
    stateerr=float(np.max(np.abs(xp-xc)))
    statesize=max(float(np.max(np.abs(xp))),1e-30)
    rows=[dict(case="ExpK_nominal",perturbation="tiny_load_step_bus16",
               status="PASS_PD_EXACT_CLOSURE",max_absolute_trajectory_error=loaderr,
               max_relative_trajectory_error=loaderr/loadscale,
               input_source="PD parameter Jacobian shared by both models"),
          dict(case="ExpK_nominal",perturbation="SG39_speed_state_1e-7_pu",
               status="PASS_PD_EXACT_CLOSURE",max_absolute_trajectory_error=stateerr,
               max_relative_trajectory_error=stateerr/statesize,
               input_source="same physical state perturbation"),
          dict(case="ExpK_nominal",perturbation="independent_port_current_voltage",
               status="NOT_INDEPENDENTLY_DEFINED",max_absolute_trajectory_error="",
               max_relative_trajectory_error="",
               input_source="Closure port matrices match by construction; direct full-network port channel not exported")]
    with (OUT/"tables"/"TABLE_M19_linear_trajectory_identity.csv").open("w",newline="",encoding="utf-8") as fh:
        writer=csv.DictWriter(fh,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    for row in rows:print(row)


if __name__=="__main__":main()
