"""Rebuild the full nonlinear certificate after changing rho and PLL gains.

Parameters modify the physical equations, equilibrium offset, complete Jacobian,
modal coordinates and interval bounds. No event/location-specific objective.
"""
from pathlib import Path
import json
import time
import numpy as np
import mpmath as mp
import coupled_dq_model as md
from modal_comparison_contract import Comparison,block_metric,vector_norm_bound
from coupled_interval_model import mm,magnitude,up

OUT=md.OUT.parent/"robust_parameter_search"


class ParametricComparison(Comparison):
    def __init__(self,rho,kp,ki):
        start=time.monotonic();self.out=OUT;self.coordinate="dq";self.module=md
        self.model=m=md.Model()
        for name,value in (("rho",rho),("Kp",kp),("Ki",ki)):
            value=np.asarray(value,dtype=float);assert value.shape==(10,) and np.isfinite(value).all()
            assert np.all(value>0)
            if name=="rho":assert np.all(value<1)
            m.raw[name]=value.tolist()
        offset=np.array(json.loads((md.OUT/"refined_center.json").read_text())["offset"])
        m.center=offset.copy()
        point=m.intervals(np.zeros(m.nx),np.zeros(m.nw))[0]
        A=np.array([(q.d[0][:m.nx]+q.d[1][:m.nx])/2 for q in point])
        B=np.array([(q.d[0][m.nx:]+q.d[1][m.nx:])/2 for q in point])
        spectral=float(np.linalg.eigvals(A).real.max())
        if spectral>=0:raise ValueError(f"Nominal Jacobian unstable: {spectral}")
        history=[]
        for iteration in range(4):
            f=m.evaluate([mp.iv.mpf(0)]*m.nx,[mp.iv.mpf(0)]*m.nw)[0]
            flo=np.array([md.lower(q) for q in f]);fhi=np.array([md.upper(q) for q in f])
            residual=float(np.maximum(abs(flo),abs(fhi)).max());history.append(residual)
            if residual<1e-20:break
            offset-=np.linalg.solve(A,(flo+fhi)/2);m.center=offset.copy()
        assert residual<1e-18,(residual,history)
        self.center=dict(offset=offset.tolist(),f0lo=flo.tolist(),f0hi=fhi.tolist(),residual_history=history)
        m.raw["A"]=A.tolist();m.raw["B"]=B.tolist()
        self.H,self.Ti,self.groups,self.modes,self.inverse_check=block_metric(A,len(m.keep))
        self.ng=len(self.groups)
        self.U=np.column_stack([up(np.sqrt(up(np.sum(magnitude(self.Ti)[:,g]**2,axis=1)*(1+8*len(g)*np.finfo(float).eps)))) for g in self.groups])
        hi=mm((self.H,self.H),(flo[:,None],fhi[:,None]))
        self.drift=np.array([vector_norm_bound(magnitude(hi)[g]) for g in self.groups])
        initial=mm((self.H,self.H),(-m.center[:,None],-m.center[:,None]))
        self.initial=np.array([vector_norm_bound(magnitude(initial)[g]) for g in self.groups])
        self.nominal_abscissa=spectral;self.build_seconds=time.monotonic()-start
        self.nominal_A=A;self.nominal_B=B


def parameters():
    d=md.Model().raw
    return [np.array(d[name]) for name in ("rho","Kp","Ki")]
