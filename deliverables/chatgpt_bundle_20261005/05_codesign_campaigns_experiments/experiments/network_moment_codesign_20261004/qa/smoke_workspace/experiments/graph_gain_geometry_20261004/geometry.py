"""Gain-box exclusion with free complex PLL patterns; exact exponential delay."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
from pathlib import Path
import sys, json, hashlib, time, platform
sys.dont_write_bytecode=True
import numpy as np
import pandas as pd
from scipy import linalg as la
import cvxpy as cp
OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
OLD=ROOT/'experiments/interaction_decision_20261004'
sys.path.insert(0,str(OLD))
from model import Model, M, PORTS, P0
TH=PORTS.pll_angle_index.to_numpy(int)
OM=PORTS.pll_frequency_index.to_numpy(int)
XI=np.argmax(abs(M['Bi']),axis=0)
HIDDEN=np.setdiff1d(np.arange(204),np.r_[TH,OM,XI])
TF=1/M['Bp'][OM,np.arange(10)]
LO=np.r_[np.full(10,.25*2*np.pi*5),np.full(10,.25*(2*np.pi*5)**2/4)]
HI=np.r_[np.full(10,4*2*np.pi*5),np.full(10,4*(2*np.pi*5)**2/4)]

def save(name,obj):
    (OUT/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def freeze():
    paths=[OUT/'PROTOCOL.md',OUT/'geometry.py',OLD/'model.py']
    paths+=list((ROOT/'experiments/graph_gsp_codesign_20261003/model').glob('*.csv'))
    m={'inputs':{str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
       'python':platform.python_version(),'numpy':np.__version__,'cvxpy':cp.__version__}
    p=OUT/'INPUT_MANIFEST.json'
    if p.exists(): assert json.loads(p.read_text())==m,'Frozen source changed; use an addendum.'
    else: save(p.name,m)

def return_matrix(m,s):
    A=s*np.eye(len(HIDDEN))-m.A0[np.ix_(HIDDEN,HIDDEN)]
    X=la.solve(A,m.A0[np.ix_(HIDDEN,TH)],check_finite=False)
    return m.C[:,TH]+m.C[:,HIDDEN]@X

def herm(x): return (x+x.conj().T)/2

def inequalities(G,s,lo=LO,hi=HI,tau=None):
    assert s.imag>0
    tau=np.full(10,.04) if tau is None else np.asarray(tau)
    a=s*s*(1+s*TF)*np.exp(s*tau)
    Q=[];labels=[]
    for i,g in enumerate(G):
        e=np.eye(10)[i];V=np.outer(g.conj(),g)
        Z=a[i]*np.outer(g.conj(),e)
        Re=herm(Z);Im=(Z-Z.conj().T)/(2j)
        p0,p1=lo[i],hi[i];j0,j1=lo[10+i],hi[10+i]
        raw=[Im-s.imag*p0*V,s.imag*p1*V-Im,
             Re-s.real/s.imag*Im-j0*V,j1*V-Re+s.real/s.imag*Im]
        center=s*(p0+p1)/2+(j0+j1)/2
        radius=max(abs(s*x+y-center) for x in [p0,p1] for y in [j0,j1])
        row=a[i]*e-center*g
        raw.append(radius**2*V-np.outer(row.conj(),row))
        for kind,q in zip(['p_lower','p_upper','i_lower','i_upper','disk'],raw):
            q=herm(q);norm=la.norm(q,'fro')
            if norm>1e-28: Q.append(q/norm);labels.append(f'{30+i}:{kind}')
    return np.array(Q),labels

def exclude(Q):
    n=Q.shape[1];m=Q.shape[0]
    w=cp.Variable(m,nonneg=True);t=cp.Variable()
    W=sum(w[i]*Q[i] for i in range(m))
    prob=cp.Problem(cp.Maximize(t),[cp.sum(w)==1, -W-t*np.eye(n)>>0])
    start=time.perf_counter()
    prob.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_feas=1e-9,tol_gap_rel=1e-9,max_iter=200)
    if w.value is None: raise RuntimeError(prob.status)
    weights=np.maximum(w.value,0);weights/=sum(weights)
    slack=-la.eigvalsh(np.einsum('i,ijk->jk',weights,Q))[-1]
    return dict(slack=float(slack),objective=float(t.value),status=prob.status,
                seconds=time.perf_counter()-start,weights=weights.tolist())

def main():
    freeze();rows=[];witnesses={}
    for rho in [.875,.90,.95,.975,.99]:
        p=P0.copy();p[:10]=rho;m=Model(p)
        for hz in [.25,.5,1,2,3,4,5,6,8,10]:
            s=-.05+2j*np.pi*hz;G=return_matrix(m,s)
            Q,lab=inequalities(G,s);res=exclude(Q)
            key=f'rho{rho}_hz{hz}';witnesses[key]={'s':[s.real,s.imag],'rho':rho,'labels':lab,**res}
            rows.append(dict(id=key,rho=rho,hz=hz,slack=res['slack'],seconds=res['seconds'],status=res['status']))
            print(key,'slack',res['slack'],flush=True)
        pd.DataFrame(rows).to_csv(OUT/'TABLE_01_POINT_SCREEN.csv',index=False)
        save('POINT_WITNESSES.json',witnesses)
    checks=[]
    for name in ['anchor','single30','single37','joint','complex_pair','corrected']:
        d=json.loads((OLD/f'CERT_{name}.json').read_text());m=Model(np.array(d['p']));s=complex(*d['center'])
        G=return_matrix(m,s);Q,labels=inequalities(G,s);res=exclude(Q)
        F=np.diag(s*s*(1+s*TF))-np.diag((s*m.p[10:20]+m.p[20:])*np.exp(-s*.04))@G
        _,sv,vh=la.svd(F);q=vh[-1].conj();vals=np.einsum('i,kij,j->k',q.conj(),Q,q).real
        checks.append(dict(design=name,root_residual=float(sv[-1]/sv[0]),min_constraint=float(min(vals)),slack=res['slack'],pass_check=bool(res['slack']<=1e-7 and min(vals)>=-1e-7)))
    pd.DataFrame(checks).to_csv(OUT/'TABLE_02_TRUE_ROOT_CHECK.csv',index=False)
    assert all(r['pass_check'] for r in checks)
    print('TRUE_ROOT_CHECKS_PASS',flush=True)

if __name__=='__main__':main()
