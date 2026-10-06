"""Seed interval derivatives in modal coordinates before evaluating the DAE.

This retains cancellations in right multiplication by H^-1 within each nonlinear
expression and the algebraic Schur solve. It changes neither the model nor the
physical box used for value enclosures. The original certificate stays frozen.
"""
import numpy as np
from modal_comparison_contract import Comparison,vector_norm_bound
from coupled_interval_model import positive_product,mm,magnitude,add,up,down
from synthesize_coupled_contract import symmetric_upper


class SeededComparison(Comparison):
    def evaluate(self,r,epsilon):
        m=self.model;Jet=self.module.Jet
        xr=positive_product(self.U,r[:,None]).ravel()
        Jet.n=m.nx+m.nw;zz=[];ww=[]
        for j in range(m.nx):
            lo=np.r_[self.Ti[0][j],np.zeros(m.nw)]
            hi=np.r_[self.Ti[1][j],np.zeros(m.nw)]
            zz.append(Jet(-xr[j],xr[j],lo,hi))
        for j in range(m.nw):
            grad=np.zeros(Jet.n);grad[m.nx+j]=1.
            ww.append(Jet(-epsilon,epsilon,grad,grad.copy()))
        out,v,ff,rr=m.evaluate(zz,ww)
        # These derivatives are already with respect to y=Hz, not physical z.
        J=(np.array([q.d[0][:m.nx] for q in out]),np.array([q.d[1][:m.nx] for q in out]))
        B=(np.array([q.d[0][m.nx:] for q in out]),np.array([q.d[1][m.nx:] for q in out]))
        AJ=mm((self.H,self.H),J);HB=mm((self.H,self.H),B)
        aa=magnitude(AJ);bb=magnitude(HB)
        M=np.zeros((self.ng,self.ng));a=np.zeros(self.ng);b=np.zeros(self.ng)
        for i,g in enumerate(self.groups):
            diag=(AJ[0][np.ix_(g,g)],AJ[1][np.ix_(g,g)])
            sym=add(diag,(diag[0].T,diag[1].T))
            a[i]=down(-symmetric_upper(sym)/2)
            b[i]=vector_norm_bound(bb[g])
            for j,h in enumerate(self.groups):
                if i!=j:M[i,j]=vector_norm_bound(aa[np.ix_(g,h)])
        rhs=up(positive_product(M,r[:,None]).ravel()+up(b*epsilon))
        rhs=up(rhs+self.drift)
        fmax=max(magnitude(q.v) for q in ff);rmax=max(magnitude(q.v) for q in rr)
        vv=[v[2*i]**2+v[2*i+1]**2 for i in range(39)]
        vmin=min(down(np.sqrt(max(0.,q.v[0]))) for q in vv)
        vmax=max(up(np.sqrt(max(0.,q.v[1]))) for q in vv)
        self.last_AJ=AJ
        return dict(a=a,M=M,b=b,rhs=rhs,margins=down(down(a*r)-rhs),x_radius=xr,
            frequency=float(fmax),rocof=float(rmax),voltage_min=float(vmin),voltage_max=float(vmax),
            output_pass=bool(fmax<.5 and rmax<.5 and vmin>.9 and vmax<1.1),
            inverse_network_check=m.last_inverse_residual)
