"""Exploratory two-Kp interpolation check. No all-root/nonlinear certificate."""
import sys,json,itertools
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.linalg import solve,svdvals
sys.path.insert(0,str(Path(__file__).parent/'BND_IEEE9_Python'/'src'))
from model import Grid9, root_catalog
out=Path(__file__).parent
m=Grid9(rho=.75)
kp=np.full(3,2*.7071*2*np.pi*5)
ki=np.full(3,(2*np.pi*5)**2)
t0=np.full(3,.020)
roots,errs=root_catalog(m,kp,ki,t0,order=8)
print('SAFE BASE ROOTS',roots,flush=True)
rows=[]
for lam in roots[roots.imag>1e-3]:
 gp,gi,gt,_,_=m.gain_sensitivity(lam,kp,ki,t0)
 for support in itertools.combinations(range(3),2):
  J=np.array([gp[list(support)].real,gp[list(support)].imag])
  eta=abs(np.linalg.det(J))/(abs(gp[support[0]])*abs(gp[support[1]]))
  for policy in ['all','first_of_pair']:
   for h in [1e-5,1e-4,1e-3]:
    tau=t0.copy()
    if policy=='all':tau+=h
    else:tau[support[0]]+=h
    F=m.char(lam,kp,ki,tau)
    P=np.zeros((24,2),complex)
    for j,site in enumerate(support):P[12+site,j]=-np.exp(-lam*tau[site])/m.tf[site]
    R=m.C[list(support)]@solve(F,P)
    aa,bb,cc=R[0,0],R[1,1],np.linalg.det(R)
    coeff=np.array([np.imag(aa*np.conj(cc)),np.imag(np.conj(cc)+aa*np.conj(bb)),np.imag(np.conj(bb))])
    # numerator is 1+aa*x, denominator bb+cc*x
    cand=np.roots(np.trim_zeros(coeff,'f')) if np.max(np.abs(coeff))>0 else []
    anyreal=False
    for x in cand:
     if abs(x.imag)>1e-7*(1+abs(x.real)):continue
     x=float(x.real); den=bb+cc*x
     if abs(den)<1e-14:continue
     y=-(1+aa*x)/den
     if abs(y.imag)>1e-6*(1+abs(y.real)):continue
     pp=kp.copy();pp[list(support)]+=np.array([x,y.real])
     relative=pp/kp
     admissible=bool(np.all((relative>=.25)&(relative<=2)))
     Fnew=m.char(lam,pp,ki,tau)
     sig=svdvals(Fnew)
     rec=dict(target_real=lam.real,target_imag=lam.imag,pair=''.join(str(i+1) for i in support),policy=policy,h=h,
       phase_diversity=eta,J_cond=np.linalg.cond(J),kp1=pp[0],kp2=pp[1],kp3=pp[2],ki_fixed=ki[0],
       gain_admissible=admissible,relative_smin=sig[-1]/sig[0],absolute_smin=sig[-1],relative_effort=np.linalg.norm(relative-1),
       interpolation_residual=abs(1+aa*x+bb*y.real+cc*x*y.real))
     if admissible:
      rr,ee=root_catalog(m,pp,ki,tau,order=8)
      rec['catalog_alpha']=rr[0].real
      rec['catalog_root_count']=len(rr)
     rows.append(rec);anyreal=True
    if not anyreal:rows.append(dict(target_real=lam.real,target_imag=lam.imag,pair=''.join(str(i+1) for i in support),policy=policy,h=h,phase_diversity=eta,J_cond=np.linalg.cond(J),gain_admissible=False,no_real_solution=True))
df=pd.DataFrame(rows)
df.to_csv(out/'two_kp_finite_probe.csv',index=False)
for pol in ['all','first_of_pair']:
 s=df[(df.policy==pol)&df.gain_admissible.fillna(False)&np.isclose(df.h,.001)]
 print('\n',pol,'h=1ms admissible',len(s))
 if len(s): print(s[['target_real','target_imag','pair','kp1','kp2','kp3','phase_diversity','catalog_alpha','absolute_smin']].sort_values('catalog_alpha').to_string(index=False))
print('Saved',len(df),'rows. No complete-spectrum or nonlinear validations performed.')
