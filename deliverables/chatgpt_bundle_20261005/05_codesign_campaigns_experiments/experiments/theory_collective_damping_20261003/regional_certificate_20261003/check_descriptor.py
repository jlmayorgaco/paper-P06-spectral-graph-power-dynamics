"""Algebraic descriptor refinement; reuses the exact frozen points and box."""
import os
os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["MKL_NUM_THREADS"]="1"
from pathlib import Path
import json
import tomllib
import numpy as np
import pandas as pd
from scipy import linalg as la
from checks import OUT,ROOT,MODEL,save,rel,freeze

def main():
    freeze()
    M={name:pd.read_csv(MODEL/(name+".csv")).to_numpy() for name in
       ["Adev","Bv","Cs","Cf","Ds","Y","Etheta","Hv","Bp","Bi"]}
    before=tomllib.loads((ROOT/"experiments/graph_gsp_codesign_20261003/baseline.toml").read_text())
    after=json.loads((OUT/"FROZEN_ALGEBRA_POINT.json").read_text())
    dr=np.array(after["rho"])-before["rho"]
    dp=np.array(after["Kp"])-before["Kp"]
    di=np.array(after["Ki"])-before["Ki"]
    tau=np.array(after["fixed_tau_seconds"])
    theta=np.diag(np.r_[np.repeat(dr,2),dp,di])
    radii=np.r_[np.full(20,.001),.01*np.array(before["Kp"]),.01*np.array(before["Ki"])]
    U=np.block([[np.zeros((204,20)),M["Bp"],M["Bi"]],
                [np.eye(20),np.zeros((20,10)),np.zeros((20,10))],
                [np.zeros((10,20)),np.zeros((10,10)),np.zeros((10,10))]])
    def matrices(s,d,E):
        w=np.repeat(1-np.array(d["rho"]),2)
        G=M["Y"]+w[:,None]*M["Ds"]
        N=w[:,None]*M["Cs"]+(1-w[:,None])*M["Cf"]
        B=M["Bp"]*np.array(d["Kp"])+M["Bi"]*np.array(d["Ki"])
        F=np.block([[s*np.eye(204)-M["Adev"],-M["Bv"],-B@E],
                    [N,G,np.zeros((20,10))],
                    [-M["Etheta"],-M["Hv"],np.eye(10)]])
        V=-la.solve(G,N)
        Delta=s*np.eye(204)-M["Adev"]-M["Bv"]@V-B@E@(M["Etheta"]+M["Hv"]@V)
        return F,Delta,G
    rows=[]
    for f in [.5,5.,10.]:
        s=-.05+2j*np.pi*f
        E=np.diag(np.exp(-s*tau))
        V=np.block([[-(M["Cs"]-M["Cf"]),-M["Ds"],np.zeros((20,10))],
                    [np.zeros((10,204)),np.zeros((10,20)),-E],
                    [np.zeros((10,204)),np.zeros((10,20)),-E]])
        F0,D0,G0=matrices(s,before,E)
        F1,D1,G1=matrices(s,after,E)
        T=V@la.solve(F0,U)
        H=theta@T
        def dratio(A,B):
            sa,laa=np.linalg.slogdet(A);sb,lbb=np.linalg.slogdet(B)
            return sa/sb*np.exp(laa-lbb)
        ratio=dratio(F1,F0)
        small=la.det(np.eye(40)+H)
        eliminated=dratio(G1,G0)*dratio(D1,D0)
        positive=radii[:,None]*abs(T)
        perron=float(max(abs(la.eigvals(positive))))
        q0=1.1*perron+1e-12
        weights=la.solve(q0*np.eye(40)-positive,np.ones(40))
        assert weights.min()>0
        weights/=weights.max()
        q=float(max((positive@weights)/weights))
        actual=float(max(abs((H*weights[None,:])/weights[:,None]).sum(axis=1)))
        assert actual<=q*(1+1e-10)
        err=rel(F1-F0,U@theta@V)
        errdet=abs(ratio-small)/max(abs(ratio),abs(small))
        errelim=abs(ratio-eliminated)/max(abs(ratio),abs(eliminated))
        rows.append(dict(frequency_Hz=f,affine_factor_error=err,
            determinant_lemma_error=float(errdet),eliminated_equivalence_error=float(errelim),
            unweighted_box_infinity_norm=float(positive.sum(axis=1).max()),
            perron_numeric=perron,weighted_box_bound_at_this_s=q,
            actual_weighted_norm_at_frozen_point=actual,
            small_gain_condition_at_this_s=q<1,
            uniform_contour_verified=False,pass_check=max(err,errdet,errelim)<1e-9))
        assert rows[-1]["pass_check"],rows[-1]
    pd.DataFrame(rows).to_csv(OUT/"TABLE_03_AFFINE_DESCRIPTOR.csv",index=False)
    save("DESCRIPTOR_SUMMARY.json",{
        "status":"NUMERICALLY_VALIDATED","scope":"three preregistered algebraic frequency points",
        "points_passed":int(sum(r["pass_check"] for r in rows)),
        "matrix_parameter_remainder":"EXACTLY_ZERO_IN_AFFINE_DESCRIPTOR",
        "weighted_box_bound_min":min(r["weighted_box_bound_at_this_s"] for r in rows),
        "weighted_box_bound_max":max(r["weighted_box_bound_at_this_s"] for r in rows),
        "pointwise_conditions_under_one":int(sum(r["small_gain_condition_at_this_s"] for r in rows)),
        "maximum_algebra_error":max(max(r["affine_factor_error"],r["determinant_lemma_error"],r["eliminated_equivalence_error"]) for r in rows),
        "IEEE39_certificate":"BLOCKED_CONTINUOUS_CONTOUR_ENCLOSURES",
        "IEEE39_replacement_upper_MW":None})
    print((OUT/"DESCRIPTOR_SUMMARY.json").read_text())

if __name__=="__main__":
    main()
