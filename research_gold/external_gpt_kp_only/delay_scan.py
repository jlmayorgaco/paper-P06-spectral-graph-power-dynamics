import sys,json
from pathlib import Path
from itertools import combinations
import numpy as np,pandas as pd
from scipy.linalg import solve,svdvals
sys.path.insert(0,str(Path(__file__).parent/'BND_IEEE9_Python'/'src'))
from model import Grid9,root_catalog
out=Path(__file__).parent;g=Grid9(rho=.75)
kp=np.full(3,2*.7071*2*np.pi*5);ki=np.full(3,(2*np.pi*5)**2);t0=np.full(3,.020)
z,_=root_catalog(g,kp,ki,t0,order=8);lam=z[np.argmax(z.imag)]
rows=[]
for delay_site in range(3):
 for h in [.0001,.00025,.0005,.001,.0015,.002,.003,.004,.005,.008]:
  t=t0.copy();t[delay_site]+=h
  roots,_=root_catalog(g,kp,ki,t,order=8);unch=roots[0].real
  F=g.char(lam,kp,ki,t)
  for support in combinations(range(3),2):
   P=np.zeros((24,2),complex)
   for j,site in enumerate(support): P[12+site,j]=-np.exp(-lam*t[site])/g.tf[site]
   R=g.C[list(support)]@solve(F,P)
   a,b,c=R[0,0],R[1,1],np.linalg.det(R)
   cf=np.array([np.imag(a*np.conj(c)),np.imag(np.conj(c)+a*np.conj(b)),np.imag(np.conj(b))])
   for x in np.roots(np.trim_zeros(cf,'f')):
    if abs(x.imag)>1e-7*(1+abs(x.real)):continue
    x=x.real;den=b+c*x
    if abs(den)<1e-14:continue
    y=-(1+a*x)/den
    if abs(y.imag)>1e-7*(1+abs(y.real)):continue
    p=kp.copy();p[list(support)]+=[x,y.real]
    if not np.all((p/kp>=.25)&(p/kp<=2)):continue
    rr,ee=root_catalog(g,p,ki,t,order=8)
    rows.append(dict(delay_site=delay_site+1,h=h,pair=''.join(str(i+1) for i in support),untuned_alpha=unch,candidate_alpha=rr[0].real,kp1=p[0],kp2=p[1],kp3=p[2],target_real=lam.real,target_imag=lam.imag))
df=pd.DataFrame(rows);df.to_csv(out/'extra_local_delay_scan.csv',index=False)
print(df[(df.untuned_alpha>-.1)&(df.candidate_alpha<=-.1)].to_string(index=False))
print('all admissible rows',len(df))
