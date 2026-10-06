"""Compare the corrected PD-exact linear response with frozen ExpK nonlinear traces."""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve


ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"reports"/"experiment_M"
CASE=OUT/"matrices"/"ExpK_nominal"


def mat(name:str)->np.ndarray:
    return np.atleast_2d(np.loadtxt(CASE/f"PD_load_input_{name}.csv",delimiter=","))


def read(path:Path)->list[dict]:
    with path.open(newline="",encoding="utf-8") as fh:return list(csv.DictReader(fh))


def write(path:Path,rows:list[dict])->None:
    with path.open("w",newline="",encoding="utf-8") as fh:
        writer=csv.DictWriter(fh,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def main()->None:
    A=mat("A");M=mat("M");B=mat("B");C=mat("C");D=mat("D")
    d=np.flatnonzero(np.diag(M)==1);z=np.flatnonzero(np.diag(M)==0)
    zz=solve(A[np.ix_(z,z)],np.hstack((A[np.ix_(z,d)],B[z,:])))
    zx=zz[:,:len(d)];zu=zz[:,len(d):]
    R=A[np.ix_(d,d)]-A[np.ix_(d,z)]@zx
    U=B[d,:]-A[np.ix_(d,z)]@zu
    Y=C[:,d]-C[:,z]@zx
    F=D-C[:,z]@zu
    meta=read(CASE/"PD_load_input_metadata.csv")[0]
    pset=float(meta["Pset_pu"]);qset=float(meta["Qset_pu"])
    w=np.array([float(meta["H38_Sn38"]),float(meta["H39_Sn39"])])
    w/=w.sum()
    u=np.array([pset,qset])
    source=U@u
    new_traces=OUT/"tables"/"TABLE_M20_nonlinear_traces.csv"
    if new_traces.exists():
        archive=read(new_traces)
        amplitudes=sorted({float(r["pulse_fraction"]) for r in archive})
        trace_source="ExpM_high_accuracy_nonlinear_trace_vs_PD_exact_linearization"
    else:
        archive=read(ROOT/"reports"/"experiment_K"/"tables"/"TABLE_K16_tds_trajectories.csv")
        amplitudes=(0.0005,0.001)
        trace_source="frozen_ExpK_nonlinear_trace_vs_PD_exact_linearization"
    trace=[];summary=[]
    for amplitude in amplitudes:
        selected=sorted((r for r in archive if r.get("candidate","nominal")=="nominal" and r["event_bus"]=="16" and
            abs(float(r["pulse_fraction"])-amplitude)<1e-12),key=lambda r:float(r["time_s"]))
        times=np.array([float(r["time_s"]) for r in selected])
        within=np.flatnonzero((times>=1.0)&(times<=1.1))
        after=np.flatnonzero(times>1.1)
        def rhs_on(_,x):return R@x+amplitude*source
        initial=np.zeros(len(d))
        on=solve_ivp(rhs_on,(1.0,1.1),initial,method="BDF",rtol=1e-10,atol=1e-12,
                     dense_output=True)
        if not on.success:raise RuntimeError(on.message)
        at_off=on.sol(1.1)
        off=solve_ivp(lambda _,x:R@x,(1.1,max(times)),at_off,method="BDF",rtol=1e-10,
                      atol=1e-12,dense_output=True)
        if not off.success:raise RuntimeError(off.message)
        linear=np.zeros(len(times))
        for index,t in enumerate(times):
            if t<1.0:continue
            if t<=1.1:
                state=on.sol(t);input_=amplitude*u
            else:
                state=off.sol(t);input_=np.zeros(2)
            local=Y@state+F@input_
            linear[index]=60*w@local
        nonlinear=np.array([float(r["COI_frequency_Hz"]) for r in selected])
        err=np.abs(nonlinear-linear)
        scale=max(np.max(np.abs(linear)),1e-30)
        summary.append(dict(case="ExpK_nominal",event_bus=16,pulse_fraction=amplitude,
            trace_points=len(times),linear_peak_Hz=float(np.max(np.abs(linear))),
            nonlinear_peak_Hz=float(np.max(np.abs(nonlinear))),
            max_absolute_error_Hz=float(np.max(err)),max_relative_error=float(np.max(err)/scale),
            source=trace_source))
        trace.extend(dict(case="ExpK_nominal",event_bus=16,pulse_fraction=amplitude,
            time_s=float(t),linear_frequency_Hz=float(l),nonlinear_frequency_Hz=float(n),
            error_Hz=float(abs(l-n))) for t,l,n in zip(times,linear,nonlinear))
    write(OUT/"tables"/"TABLE_M19_linear_nonlinear_trace.csv",trace)
    write(OUT/"tables"/"TABLE_M20_nonlinear_linear_scaling.csv",summary)
    if len(summary)>=3:
        x=np.array([r["pulse_fraction"] for r in summary],dtype=float)
        y=np.array([r["max_absolute_error_Hz"] for r in summary],dtype=float)
        slope,intercept=np.polyfit(np.log(x),np.log(y),1)
        fit=[dict(case="ExpK_nominal",event_bus=16,samples=len(x),
            absolute_error_scaling_exponent=float(slope),
            relative_error_scaling_exponent=float(slope-1),
            minimum_amplitude=float(np.min(x)),maximum_amplitude=float(np.max(x)))]
        write(OUT/"tables"/"TABLE_M20_scaling_fit.csv",fit)
        fig,ax=plt.subplots(figsize=(5.5,4))
        ax.loglog(x,y,"o",label="Measured nonlinear − linear")
        grid=np.geomspace(min(x),max(x),100)
        ax.loglog(grid,np.exp(intercept)*grid**slope,"--",label=f"Fit exponent {slope:.2f}")
        ax.set_xlabel("Pulse fraction δ")
        ax.set_ylabel("Maximum absolute frequency error (Hz)")
        ax.grid(True,which="both",alpha=.3);ax.legend();fig.tight_layout()
        fig.savefig(OUT/"FIG_M04_error_scaling.png",dpi=170)
    for row in summary:print(row)


if __name__=="__main__":main()
