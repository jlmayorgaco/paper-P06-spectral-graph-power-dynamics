import sys,json
from pathlib import Path
import numpy as np,pandas as pd
sys.path.insert(0,str(Path(__file__).parent/'BND_IEEE9_Python'/'src'))
from model import Grid9,root_catalog
from validate import contour_count
out=Path(__file__).parent
g=Grid9(rho=.75)
kp=np.full(3,2*.7071*2*np.pi*5);ki=np.full(3,(2*np.pi*5)**2)
t0=np.full(3,.020);t1=t0.copy();t1[0]=.021
D=pd.read_csv(out/'extra_local_delay_scan.csv')
s=D[(D.delay_site==1)&(D.pair.astype(str)=='23')&np.isclose(D.h,.001)].sort_values('candidate_alpha').iloc[0]
p=np.array([s.kp1,s.kp2,s.kp3]);lam=complex(s.target_real,s.target_imag)
a,w=lam.real,lam.imag;dh=.001
pt=kp.copy();it=ki.copy();i=0
pt[i]=np.exp(a*dh)*(kp[i]*np.cos(w*dh)+(ki[i]+a*kp[i])/w*np.sin(w*dh))
it[i]=np.exp(a*dh)*(ki[i]*np.cos(w*dh)-(abs(lam)**2*kp[i]+a*ki[i])/w*np.sin(w*dh))
results={}
for name,P,I,T in [('base20',kp,ki,t0),('delay_bus1_unretuned',kp,ki,t1),('two_Kp_fixed_Ki',p,ki,t1),('standard_transport_bus1',pt,it,t1)]:
 rec=dict(kp=P.tolist(),ki=I.tolist(),tau=T.tolist(),target=[lam.real,lam.imag])
 roots_by_order={}
 for order in [6,8,10]:
  roots,errs=root_catalog(g,P,I,T,order=order)
  roots_by_order[str(order)]=[[float(z.real),float(z.imag)] for z in roots]
 rec['roots_by_seed_order']=roots_by_order
 rec['catalog_alpha']=roots[0].real
 rec['nearest_target_error']=min(abs(roots-lam))
 counts=[]
 for sigma in [1e-8,.1]:
  for density in [1,2]:counts.append(dict(density=density,**contour_count(g,P,I,T,sigma=sigma,density=density)))
 rec['numeric_counts']=counts
 results[name]=rec
 print(name,rec['kp'],rec['ki'],'alpha',rec['catalog_alpha'],'targeterr',rec['nearest_target_error'],'counts',[c['count'] for c in counts],flush=True)
(out/'remote_repair_numeric_checks.json').write_text(json.dumps(results,indent=2))
