"""Full-network Laurent moments. Explicitly restores only rotational identities."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
from pathlib import Path
import sys, json, tomllib, hashlib
import numpy as np
import pandas as pd
from scipy import linalg as la

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
OLD=OUT.parent/'interaction_decision_20261004'
sys.path.insert(0,str(OLD))
from model import Model, M, PORTS, ROOTS
DESIGN=tomllib.loads((OLD/'designs/corrected.toml').read_text())
P=np.r_[DESIGN['rho'],DESIGN['Kp'],DESIGN['Ki']]
TF=1/(300*2*np.pi); W=.5

def dump(name,obj):
    (OUT/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf8')

def freeze():
    paths=[OUT/'PROTOCOL.md',OUT/'ADDENDUM_01_ALL_SITES.md',OUT/'moments.py',OUT/'export_inputs.jl',OLD/'designs/corrected.toml']
    paths+=sorted((OUT/'model').glob('*.csv'))
    info={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    f=OUT/'MOMENT_INPUT_LOCK.json'
    if f.exists(): assert json.loads(f.read_text())==info
    else: dump(f.name,info)

class Moments:
    def __init__(self):
        mats={p.stem:pd.read_csv(p).to_numpy() for p in (OUT/'model').glob('*.csv')}
        self.mats=mats
        th=PORTS.pll_angle_index.to_numpy()-1
        hidden=np.setdiff1d(np.arange(204),np.r_[th,th+1,th+2])
        A0=mats['A0']; C=mats['C']; O=mats['phase_state']
        self.A=-A0[np.ix_(hidden,hidden)]
        self.B=A0[np.ix_(hidden,th)];self.Bd=mats['Binput'][hidden]
        self.C=C[:,hidden];self.D=C[:,th]
        self.O=O[:,hidden];self.J=O[:,th]
        self.de=mats['detector_input'];self.df=mats['phase_input']
        self.lu=la.lu_factor(self.A)
        XX=la.lu_solve(self.lu,np.c_[self.B,self.Bd])
        self.G=[];self.g=[];self.T=[];self.Z=[]
        for j in range(4):
            self.G.append(self.C@XX[:,:10]+(self.D if j==0 else 0))
            self.g.append(self.C@XX[:,10:]+(self.de if j==0 else 0))
            self.T.append(self.O@XX[:,:10]+(self.J if j==0 else 0))
            self.Z.append(self.O@XX[:,10:]+(self.df if j==0 else 0))
            XX=-la.lu_solve(self.lu,XX)
        r=self.r=np.ones(10)
        self.Gfix=np.outer(self.G[0]@r,r)/10
        self.Tfix=np.outer(np.ones(39)-self.T[0]@r,r)/10
        self.G[0]=self.G[0]-self.Gfix;self.T[0]=self.T[0]+self.Tfix
        K0=-self.G[0];U,sv,Vh=la.svd(K0)
        self.ell=U[:,-1];self.ell/=self.ell@r
        self.d=float(self.ell@(-self.G[1])@r)
        self.weights=self.ell/self.d
        self.border=np.block([[K0,r[:,None]],[self.ell[None,:],np.zeros((1,1))]])
        self.blu=la.lu_factor(self.border)
        ab,scales=la.matrix_balance(self.A)
        parity=Model(P)
        self.audit={'hidden_dimension':len(hidden),'hidden_condition_raw':float(np.linalg.cond(self.A)),
            'hidden_condition_balanced':float(np.linalg.cond(ab)),
            'hidden_min_sv':float(la.svdvals(self.A)[-1]),'gauge_rank_second_smallest_sv':float(sv[-2]),
            'gauge_rank_smallest_sv_restored':float(sv[-1]),'d':self.d,
            'G0_gauge_correction_norm':float(la.norm(self.Gfix)),
            'T0_gauge_correction_norm':float(la.norm(self.Tfix)),
            'G0_gauge_correction_relative':float(la.norm(self.Gfix)/la.norm(self.G[0])),
            'A0_export_model_relative':float(la.norm(A0-parity.A0)/la.norm(A0)),
            'C_export_model_relative':float(la.norm(C-parity.C)/la.norm(C)),
            'status':'NUMERICALLY_VALIDATED; explicit symmetry-restored germ, not interval certificate'}

    def solve_border(self,b):
        return la.lu_solve(self.blu,np.r_[b,np.zeros((1,b.shape[1]))])[:10]

    def coefficients(self,kp,ki,tau):
        K1=-self.G[1];K2=-self.G[2]+np.diag(1/ki)
        K3=-self.G[3]+np.diag((TF+tau)/ki-kp/ki**2)
        xm=np.outer(self.r,(self.ell@self.g[0])/self.d)
        z0=self.solve_border(self.g[0]-K1@xm)
        a0=self.ell@(self.g[1]-K1@z0-K2@xm)/self.d
        x0=z0+np.outer(self.r,a0)
        z1=self.solve_border(self.g[1]-K1@x0-K2@xm)
        a1=self.ell@(self.g[2]-K1@z1-K2@x0-K3@xm)/self.d
        x1=z1+np.outer(self.r,a1)
        t0=self.T[0]@xm
        t1=self.T[0]@x0+self.T[1]@xm+self.Z[0]
        t2=self.T[0]@x1+self.T[1]@x0+self.T[2]@xm+self.Z[1]
        return np.array([t0,t1-W/2*t0,t2-W/2*t1+W**2/6*t0])/(2*np.pi)

    def maps(self,s,restore=True):
        X=la.solve(self.A+s*np.eye(174),np.c_[self.B,self.Bd])
        G=self.C@X[:,:10]+self.D
        g=self.C@X[:,10:]+self.de
        T=self.O@X[:,:10]+self.J
        Z=self.O@X[:,10:]+self.df
        if restore:G=G-self.Gfix;T=T+self.Tfix
        return G,g,T,Z

    def transfer(self,s,kp,ki,tau,restore=True):
        G,g,T,Z=self.maps(s,restore)
        Q=s*s*(1+TF*s)*np.exp(s*tau)/(s*kp+ki)
        return -np.expm1(-s*W)/(2*np.pi*W)*(T@la.solve(np.diag(Q)-G,g)+Z)

    def difference(self,s,kpa,kia,ta,kpb,kib,tb,restore=True):
        G,g,T,Z=self.maps(s,restore)
        Qa=s*s*(1+TF*s)*np.exp(s*ta)/(s*kpa+kia)
        dQ=Qa*((s*kpa+kia)*np.expm1(s*(tb-ta))-s*(kpb-kpa)-(kib-kia))/(s*kpb+kib)
        Ka=np.diag(Qa)-G;Kb=Ka+np.diag(dQ)
        return np.expm1(-s*W)/(2*np.pi*W)*(T@la.solve(Kb,dQ[:,None]*la.solve(Ka,g)))

def patterns():
    out={f'uniform_{x:g}ms':np.full(10,x/1000) for x in [.1,.5,1.]}
    out.update({f'site{30+i}_1ms':np.eye(10)[i]*.001 for i in range(10)})
    out['ramp_1ms']=np.linspace(0,.001,10)
    return out

def main():
    freeze();m=Moments();dump('STRUCTURAL_AUDIT.json',m.audit)
    kp=P[10:20];ki=P[20:];tau=np.full(10,.04);H=m.coefficients(kp,ki,tau)
    pd.DataFrame({'bus':np.arange(30,40),'weight':m.weights,'moment_Kp_coefficient':m.weights/ki**2,
                  'moment_tau_coefficient':m.weights/ki}).to_csv(OUT/'TABLE_02_NETWORK_WEIGHTS.csv',index=False)
    np.savez(OUT/'LAURENT_BASE.npz',H=H,weights=m.weights)
    rows=[];fd=[]
    cases=[(n,kp.copy(),ki.copy(),tau+dt) for n,dt in patterns().items()]
    cases+=[(f'Ki_site{i+30}',kp.copy(),ki*(1+.1*np.eye(10)[i]),tau.copy()) for i in [0,2,3,9]]
    for n,kpb,kib,tb in cases:
        Hb=m.coefficients(kpb,kib,tb)
        iski=n.startswith('Ki_');power=1 if iski else 2
        kap=float(m.weights@(1/kib-1/ki)) if iski else float(m.weights@((tb-tau)/ki-(kpb-kp)/ki**2))
        pred=-kap*H[0];actual=Hb[power]-H[power]
        rows.append(dict(case=n,kappa=kap,power=power,coefficient_abs_error=float(la.norm(actual-pred)),
            coefficient_relative_error=float(la.norm(actual-pred)/max(la.norm(pred),1e-30)),
            lower_coeff_max_error=float(np.max(abs(Hb[:power]-H[:power])))))
        for s in [.02,.01,.005,.0025]:
            diff=m.difference(s,kp,ki,tau,kpb,kib,tb)/s**power
            raw=m.difference(s,kp,ki,tau,kpb,kib,tb,False)/s**power
            fd.append(dict(case=n,s=s,power=power,relative_asymptotic_error=float(la.norm(diff-pred)/la.norm(pred)),
                raw_vs_restored_relative=float(la.norm(raw-diff)/la.norm(pred))))
    pd.DataFrame(rows).to_csv(OUT/'TABLE_03_MOMENT_IDENTITIES.csv',index=False)
    pd.DataFrame(fd).to_csv(OUT/'TABLE_04_FINITE_FREQUENCY_CHECK.csv',index=False)
    print(json.dumps(m.audit,indent=2));print(pd.DataFrame(rows).to_string(index=False))
    print('DC responses Hz per unit MW parameter:',H[0,0,:])

if __name__=='__main__':main()
