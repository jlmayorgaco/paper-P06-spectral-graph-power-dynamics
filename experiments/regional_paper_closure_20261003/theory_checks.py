"""Nontrajectory checks of the second-order contour/closed-walk identities."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json,tomllib
import numpy as np
import pandas as pd
from scipy import linalg as la
import mpmath as mp
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
from paths import resolve
MODEL=resolve('experiments/graph_gsp_codesign_20261003/model')
M={k:pd.read_csv(MODEL/(k+'.csv')).to_numpy() for k in ['Adev','Bv','Cs','Cf','Ds','Y','Etheta','Hv','Bp','Bi']}
base=tomllib.loads((MODEL.parent/'baseline.toml').read_text())
frozen=json.loads(resolve('experiments/theory_collective_damping_20261003/regional_certificate_20261003/FROZEN_ALGEBRA_POINT.json').read_text())
radii=np.r_[np.full(20,.001),.01*np.array(base['Kp']),.01*np.array(base['Ki'])]
dr=np.array(frozen['rho'])-base['rho'];dp=np.array(frozen['Kp'])-base['Kp'];di=np.array(frozen['Ki'])-base['Ki']
d=np.r_[np.repeat(dr,2),dp,di]
U=np.block([[np.zeros((204,20)),M['Bp'],M['Bi']],[np.eye(20),np.zeros((20,20))],[np.zeros((10,40))]])
mp.mp.dps=50
def psi(A,order):
    aa=mp.matrix(A.tolist());val=-mp.log(mp.det(mp.eye(A.shape[0])-aa))
    powa=mp.eye(A.shape[0])
    for ell in range(1,order):
        powa=powa*aa;val-=sum(powa[i,i] for i in range(A.shape[0]))/ell
    return float(mp.re(val))
rows=[]
def canonical_log(H):
    eig=la.eigvals(H)
    assert max(abs(eig))<1
    return sum(np.log1p(eig))
for freq in [.5,5.,10.]:
    s=-.05+2j*np.pi*freq;e=np.diag(np.exp(-s*np.array(frozen['fixed_tau_seconds'])))
    w=np.repeat(1-np.array(base['rho']),2);G=M['Y']+w[:,None]*M['Ds'];N=w[:,None]*M['Cs']+(1-w[:,None])*M['Cf']
    B=M['Bp']*base['Kp']+M['Bi']*base['Ki']
    F=np.block([[s*np.eye(204)-M['Adev'],-M['Bv'],-B@e],[N,G,np.zeros((20,10))],[-M['Etheta'],-M['Hv'],np.eye(10)]])
    V=np.block([[-(M['Cs']-M['Cf']),-M['Ds'],np.zeros((20,10))],[np.zeros((10,224)),-e],[np.zeros((10,224)),-e]])
    T=V@la.solve(F,U);H=d[:,None]*T;MA=radii[:,None]*abs(T)
    perron=float(max(abs(la.eigvals(MA))));ww=la.solve(1.1*perron*np.eye(40)-MA,np.ones(40));q=float(max((MA@ww)/ww))
    logdet=canonical_log(H);p1=np.trace(H);p2=p1-.5*np.trace(H@H)
    major2=psi(MA,2);major3=psi(MA,3)
    same=0
    for i in range(10):
        inds=[2*i,2*i+1,20+i,30+i];same+=psi(MA[np.ix_(inds,inds)],2)
    diagonal=np.zeros_like(H)
    for i in range(10):
        inds=[2*i,2*i+1,20+i,30+i];diagonal[np.ix_(inds,inds)]=H[np.ix_(inds,inds)]
    omitted=abs(logdet-canonical_log(diagonal))
    # Independent mixed derivative of log determinant vs tr(Ea T Eb T).
    # Coordinates scaled to their physical half-widths for conditioning.
    aa,bb=0,20
    A=np.zeros_like(H);C=np.zeros_like(H)
    A[0]=radii[0]*T[0];A[1]=radii[1]*T[1];C[bb]=radii[bb]*T[bb]
    h=1e-3
    fun=lambda x,y:canonical_log(x*A+y*C)
    fd=(fun(h,h)-fun(h,-h)-fun(-h,h)+fun(-h,-h))/(4*h*h)
    exact=-np.trace(A@C)
    row=dict(frequency_Hz=freq,q=q,scalar_norm_remainder=40*(-np.log1p(-q)-q),
       walk_remainder_order1=major2,walk_remainder_order2=major3,
       actual_error_order1=abs(logdet-p1),actual_error_order2=abs(logdet-p2),
       cross_bus_walk_remainder=major2-same,actual_bus_block_omission=omitted,
       mixed_Hessian_absolute_error=abs(fd-exact),mixed_Hessian_relative_error=abs(fd-exact)/max(abs(exact),1e-30),
       status='NUMERICALLY_VALIDATED_POINTWISE_NOT_CONTOUR_CERTIFICATE')
    assert row['actual_error_order1']<=major2 and row['actual_error_order2']<=major3
    assert omitted<=major2-same and row['mixed_Hessian_absolute_error']<1e-8
    rows.append(row)
pd.DataFrame(rows).to_csv(OUT/'TABLE_01_WALK_REMAINDER.csv',index=False)

# The previous 2-state pure-DDE fixture; exact implicit root derivatives at p=0.
lam=mp.mpf(-1);ep=mp.exp(mp.mpf('.1'))
J=ep;Q=mp.mpf('1.8')*ep**2
toy=[]
for nn in [256,512]:
    contour=[-1+.25*np.exp(2j*np.pi*k/nn) for k in range(nn)]
    weights=[.25*np.exp(2j*np.pi*k/nn)/nn for k in range(nn)]
    Ts=[-np.exp(-s*.1)*np.diag([1/(s+1),1/(s+2)])@np.array([[1,1],[1,-1]]) for s in contour]
    JJ=-sum(np.trace(T)*w for T,w in zip(Ts,weights))
    QQ=sum(np.trace(T@T)*w for T,w in zip(Ts,weights))
    assert abs(JJ-complex(J))<1e-12 and abs(QQ-complex(Q))<1e-12
    for p in ['-.01','-.005','0','.005','.01']:
        pp=mp.mpf(p)
        det=lambda x:(x+1)*(x+2)-pp*mp.exp(-x/10)-2*pp**2*mp.exp(-x/5)
        rt=mp.findroot(det,(-1,-.99))
        prediction=lam+J*pp+Q*pp**2/2
        kappa=mp.mpf(12)/175
        bound=mp.mpf('.25')*2*kappa**3/(3*(1-kappa))
        assert abs(rt-prediction)<bound
        toy.append(dict(nodes=nn,p=float(pp),root=float(rt),J_error=abs(JJ-complex(J)),Q_error=abs(QQ-complex(Q)),
           quadratic_error=float(abs(rt-prediction)),rigorous_analytic_cubic_remainder_upper=float(bound)))
pd.DataFrame(toy).to_csv(OUT/'TABLE_02_QUADRATIC_TOY.csv',index=False)
summary={'pointwise_checks':len(rows),'toy_checks':len(toy),
         'walk_vs_scalar_improvement':[r['scalar_norm_remainder']/r['walk_remainder_order1'] for r in rows],
         'max_toy_quadratic_error':max(r['quadratic_error'] for r in toy),
         'rigorous_toy_bound':toy[0]['rigorous_analytic_cubic_remainder_upper'],
         'scope':'algebraic/contour identity validation; no new physical trajectory or optimization'}
(OUT/'CHECKS_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
