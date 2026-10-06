"""Condition modal coordinates for a decay target, not exact normal form.

Real invariant two-dimensional spaces are retained. Within each complex pair,
use an orthonormal basis and a Lyapunov metric with a fractional decay target.
This avoids forcing almost-real pairs into badly conditioned rotation form.
"""
import numpy as np
import scipy.linalg as la
from modal_comparison_contract import Comparison,vector_norm_bound
from reverse_interval_contract import ReverseComparison
from coupled_interval_model import inverse_enclosure,magnitude,mm,up


def balanced_metric(A,nplant,fraction):
    ap=A[:nplant,:nplant];am=A[nplant:,nplant:];c=A[nplant:,:nplant]
    ev,vec=la.eig(ap);columns=[];groups=[];modes=[]
    for i,e in enumerate(ev):
        if e.imag>1e-7:
            Q,_=la.qr(np.column_stack([vec[:,i].real,vec[:,i].imag]),mode="economic")
            local=Q.T@ap@Q
            target=fraction*(-e.real)
            assert target>0
            P=la.solve_continuous_lyapunov((local+target*np.eye(2)).T,-np.eye(2))
            P=(P+P.T)/2
            root=la.cholesky(P,lower=False)
            pair=Q@la.solve_triangular(root,np.eye(2),lower=False)
            pair/=la.norm(pair)
            groups.append(list(range(len(columns),len(columns)+2)))
            columns.extend([pair[:,0],pair[:,1]]);modes.append([float(e.real),float(e.imag)])
        elif abs(e.imag)<=1e-7:
            groups.append([len(columns)]);columns.append(vec[:,i].real/la.norm(vec[:,i].real))
            modes.append([float(e.real),0.])
    tp=np.column_stack(columns);assert tp.shape==ap.shape
    cross=la.solve_sylvester(am,-ap,-c)
    tm=np.diag(np.r_[np.ones(39),5*np.ones(39)])
    T=np.block([[tp,np.zeros((nplant,78))],[cross@tp,tm]])
    for bus in range(39):groups.append([nplant+bus,nplant+39+bus]);modes.append([-10.,0.])
    H=np.linalg.inv(T);Ti,check=inverse_enclosure((H,H))
    return H,Ti,groups,modes,check


def configure(cc,fraction):
    m=cc.model
    cc.H,cc.Ti,cc.groups,cc.modes,cc.inverse_check=balanced_metric(np.array(m.raw["A"]),len(m.keep),fraction)
    cc.ng=len(cc.groups);cc.decay_fraction=fraction
    cc.U=np.column_stack([up(np.sqrt(up(np.sum(magnitude(cc.Ti)[:,g]**2,axis=1)*(1+8*len(g)*np.finfo(float).eps)))) for g in cc.groups])
    fi=(np.array(cc.center["f0lo"])[:,None],np.array(cc.center["f0hi"])[:,None])
    hi=mm((cc.H,cc.H),fi)
    cc.drift=np.array([vector_norm_bound(magnitude(hi)[g]) for g in cc.groups])
    initial=mm((cc.H,cc.H),(-m.center[:,None],-m.center[:,None]))
    cc.initial=np.array([vector_norm_bound(magnitude(initial)[g]) for g in cc.groups])
    return cc


class BalancedComparison(Comparison):
    def __init__(self,fraction=.5):
        super().__init__("dq");configure(self,fraction)


class BalancedReverseComparison(ReverseComparison):
    def __init__(self,fraction=.5):
        super().__init__();configure(self,fraction)
