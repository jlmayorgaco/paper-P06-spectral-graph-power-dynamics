"""Reverse interval differentiation preserves shared-output correlations.

Evaluate the unchanged full dq DAE expression graph. Propagate H from outputs
backward before interval derivatives are multiplied, including one implicit
network-solve node. The result encloses H F_z and H F_w on the physical box.
"""
import numpy as np
import coupled_dq_model as md
from modal_comparison_contract import Comparison,vector_norm_bound
from synthesize_coupled_contract import symmetric_upper
from coupled_interval_model import add,neg,mul,inv,ivfun,mm,magnitude,inverse_enclosure,positive_product,up,down


class ReverseJet:
    tape=[]
    def __init__(self,lo,hi=None,parents=(),variable=False):
        self.v=(float(lo),float(lo if hi is None else hi));self.idx=None
        parents=tuple((p.idx,d) for p,d in parents if p.idx is not None)
        if variable or parents:
            self.idx=len(self.tape);self.tape.append(parents)
    @classmethod
    def make(cls,value,parents=()):return cls(value[0],value[1],parents)
    @staticmethod
    def cast(x):return x if isinstance(x,ReverseJet) else ReverseJet(x)
    def __add__(self,other):
        b=self.cast(other);return self.make(add(self.v,b.v),((self,(1.,1.)),(b,(1.,1.))))
    __radd__=__add__
    def __neg__(self):return self.make(neg(self.v),((self,(-1.,-1.)),))
    def __sub__(self,other):return self+(-self.cast(other))
    def __rsub__(self,other):return self.cast(other)+(-self)
    def __mul__(self,other):
        b=self.cast(other);return self.make(mul(self.v,b.v),((self,b.v),(b,self.v)))
    __rmul__=__mul__
    def reciprocal(self):
        value=inv(self.v);return self.make(value,((self,neg(mul(value,value))),))
    def __truediv__(self,other):return self*self.cast(other).reciprocal()
    def __rtruediv__(self,other):return self.cast(other)*self.reciprocal()
    def __pow__(self,n):
        assert n==2
        value=mul(self.v,self.v)
        if self.v[0]<=0<=self.v[1]:value=(0.,value[1])
        return self.make(value,((self,mul((2.,2.),self.v)),))
    def fun(self,name):
        value=ivfun(name,self.v)
        if name=="sin":derivative=ivfun("cos",self.v)
        elif name=="cos":derivative=neg(ivfun("sin",self.v))
        elif name=="sqrt":derivative=mul((.5,.5),inv(value))
        elif name=="atan":derivative=inv(add((1.,1.),(self**2).v))
        else:raise ValueError(name)
        return self.make(value,((self,derivative),))


class ReverseModel(md.Model):
    def voltage(self,coeff,source):
        gl=self.Y.copy();gh=self.Y.copy()
        for i,k in enumerate(coeff):
            j=2*(i+29)
            gl[j,j+1],gh[j,j+1]=add((gl[j,j+1],gh[j,j+1]),neg(k.v))
            gl[j+1,j],gh[j+1,j]=add((gl[j+1,j],gh[j+1,j]),k.v)
        kin,check=inverse_enclosure((gl,gh));self.last_inverse_residual=check
        hv=(np.array([s.v[0] for s in source])[:,None],np.array([s.v[1] for s in source])[:,None])
        vc=-np.linalg.solve(gl/2+gh/2,(hv[0]+hv[1]).ravel()/2)
        residual=add(hv,mm((gl,gh),(vc[:,None],vc[:,None])))
        vin=add((vc[:,None],vc[:,None]),neg(mm(kin,residual)))
        values=[]
        for row in range(78):
            parents=[(s,(-kin[1][row,j],-kin[0][row,j])) for j,s in enumerate(source)]
            for i,k in enumerate(coeff):
                j=2*(i+29)
                partial=add(mul((kin[0][row,j],kin[1][row,j]),(vin[0][j+1,0],vin[1][j+1,0])),
                            neg(mul((kin[0][row,j+1],kin[1][row,j+1]),(vin[0][j,0],vin[1][j,0]))))
                parents.append((k,partial))
            values.append(ReverseJet(vin[0][row,0],vin[1][row,0],parents))
        return values

    def reverse_enclosure(self,xr,epsilon,H):
        previous=md.Jet;md.Jet=ReverseJet;ReverseJet.tape=[]
        try:
            z=[ReverseJet(-float(r),float(r),variable=True) for r in xr]
            w=[ReverseJet(-epsilon,epsilon,variable=True) for _ in range(self.nw)]
            outputs=self.evaluate(z,w)
            n=H.shape[0];count=len(ReverseJet.tape)
            lo=np.zeros((count,n));hi=np.zeros((count,n));active=np.zeros(count,bool)
            for j,q in enumerate(outputs[0]):
                if q.idx is None:continue
                idx=q.idx;v=H[:,j]
                if active[idx]:lo[idx],hi[idx]=add((lo[idx],hi[idx]),(v,v))
                else:lo[idx]=v;hi[idx]=v;active[idx]=True
            for idx in range(count-1,-1,-1):
                if not active[idx]:continue
                adj=(lo[idx],hi[idx])
                for parent,derivative in ReverseJet.tape[idx]:
                    contribution=mul(adj,derivative)
                    if active[parent]:lo[parent],hi[parent]=add((lo[parent],hi[parent]),contribution)
                    else:lo[parent],hi[parent]=contribution;active[parent]=True
            self.last_tape_size=count
            self.last_HJ=(lo[:self.nx].T,hi[:self.nx].T)
            self.last_HB=(lo[self.nx:self.nx+self.nw].T,hi[self.nx:self.nx+self.nw].T)
            return outputs
        finally:md.Jet=previous


class ReverseComparison(Comparison):
    def __init__(self):
        super().__init__("dq")
        self.model=ReverseModel(center=self.center["offset"])
    def evaluate(self,r,epsilon):
        m=self.model;xr=positive_product(self.U,r[:,None]).ravel()
        out,v,ff,rr=m.reverse_enclosure(xr,epsilon,self.H)
        AJ=mm(m.last_HJ,self.Ti);HB=m.last_HB
        aa=magnitude(AJ);bb=magnitude(HB)
        M=np.zeros((self.ng,self.ng));a=np.zeros(self.ng);b=np.zeros(self.ng)
        for i,g in enumerate(self.groups):
            diag=(AJ[0][np.ix_(g,g)],AJ[1][np.ix_(g,g)])
            sym=add(diag,(diag[0].T,diag[1].T))
            a[i]=down(-symmetric_upper(sym)/2)
            b[i]=vector_norm_bound(bb[g])
            for j,h in enumerate(self.groups):
                if i!=j:M[i,j]=vector_norm_bound(aa[np.ix_(g,h)])
        rhs=up(positive_product(M,r[:,None]).ravel()+up(b*epsilon));rhs=up(rhs+self.drift)
        fmax=max(magnitude(q.v) for q in ff);rmax=max(magnitude(q.v) for q in rr)
        vv=[q.v for q in v];vmin=1e99;vmax=0.
        for i in range(39):
            def square(q):
                val=mul(q,q);return (0.,val[1]) if q[0]<=0<=q[1] else val
            sq=add(square(vv[2*i]),square(vv[2*i+1]))
            vmin=min(vmin,down(np.sqrt(max(0.,sq[0]))));vmax=max(vmax,up(np.sqrt(max(0.,sq[1]))))
        self.last_AJ=AJ
        return dict(a=a,M=M,b=b,rhs=rhs,margins=down(down(a*r)-rhs),x_radius=xr,
            frequency=float(fmax),rocof=float(rmax),voltage_min=float(vmin),voltage_max=float(vmax),
            output_pass=bool(fmax<.5 and rmax<.5 and vmin>.9 and vmax<1.1),
            inverse_network_check=m.last_inverse_residual)
