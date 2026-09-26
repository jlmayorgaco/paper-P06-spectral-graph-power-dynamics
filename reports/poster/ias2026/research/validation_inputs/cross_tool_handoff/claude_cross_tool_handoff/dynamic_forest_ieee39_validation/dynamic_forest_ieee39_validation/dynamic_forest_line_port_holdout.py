import os,sys,json
from pathlib import Path
from dataclasses import replace
import numpy as np,pandas as pd
from scipy.linalg import eig
REPO=Path('/mnt/data/research_unzip/research');sys.path[:0]=[str(REPO/'src'),str(REPO/'experiments')];os.chdir(REPO)
from _f7_common import Theta, LEAK
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case
from ibr_cycles.models.ieee39_devices import ConverterParameters
from ibr_cycles.models.ieee39_network import load_network
from ibr_cycles.models.port_admittance import build_action_space
from ibr_cycles.certification.symmetry import rotation_generator,frequency_partner
from ibr_cycles.certification.transverse import transverse_operator

OUT=Path('/mnt/data/dynamic_forest_ieee39_validation');DATA=OUT/'data'
payload=json.loads((REPO/'configs/ias2026/ieee39_network.json').read_text());net=load_network()
def ybus_scaled(line_idx,gamma):
 order=[int(b['idx']) for b in payload['buses']];index={b:i for i,b in enumerate(order)};n=len(order);y=np.zeros((n,n),complex)
 for ell,line in enumerate(payload['lines']):
  if float(line['u'])==0:continue
  f,t=index[int(line['bus1'])],index[int(line['bus2'])]; scale=gamma if ell==line_idx else 1.0
  series=scale/complex(line['r'],line['x']);charging=scale*complex(line['g'],line['b'])/2
  m=float(line['tap'])*np.exp(1j*float(line['phi']));m2=abs(m)**2
  y[f,f]+=(series+charging)/m2;y[t,t]+=series+charging;y[f,t]+=-series/np.conj(m);y[t,f]+=-series/m
 for sh in payload['shunts']:y[index[int(sh['bus'])],index[int(sh['bus'])]]+=complex(sh['g'],sh['b'])
 return y

def solve(members,network):
 return solve_case(ReplacementPlan.of({b:1 for b in members}),network=network,
                   converter=ConverterParameters(voltage_control=True,voltage_gain=0.20768,voltage_leak=0.05),
                   machine_scaling={'ka':1.425,'ta':1.5})
def qs(space,s):
 m=space.m(s);I=np.eye(m.shape[0],dtype=complex);total=I+m;sb=np.zeros_like(total)
 for k in range(space.order):
  sl=slice(2*k,2*k+2);sb[sl,sl]=total[sl,sl]
 return np.linalg.solve(sb,total)-I
def mu(space,s,target=-1):
 vals=np.linalg.eigvals(qs(space,s));return vals[np.argmin(np.abs(vals-target))]

def trans_eig(case,target):
 rx,_=rotation_generator(case.dae,case.equilibrium.z);w=frequency_partner(case.dae).w
 ev=np.linalg.eigvals(transverse_operator(case.system.A,rx,w).a_perp);c=ev[ev.imag>=0]
 return c[np.argmin(np.abs(c-target))]

sel=[3,9,15,17,22,26,31,32,35,39,42,44]
base0=solve((),net);full0=solve((30,33,35,37),net);sp0=build_action_space(base0,full0,(30,33,35,37))
sstar=1j*2*np.pi*0.706;target=trans_eig(full0,sstar);hs=1e-5;mu_s=(mu(sp0,sstar+hs)-mu(sp0,sstar-hs))/(2*hs)
eps=0.002;rows=[]
for li in sel:
 try:
  spaces=[];lams=[]
  for gam in (1-eps,1+eps):
   nn=replace(net,ybus=ybus_scaled(li,gam));b=solve((),nn);f=solve((30,33,35,37),nn);spaces.append(build_action_space(b,f,(30,33,35,37)));lams.append(trans_eig(f,target))
  mug=(mu(spaces[1],sstar)-mu(spaces[0],sstar))/(2*eps);ds=-mug/mu_s;dl=(lams[1]-lams[0])/(2*eps)
  line=payload['lines'][li]
  rows.append({'line_index':li,'line':f"L{li:02d}:{line['bus1']}-{line['bus2']}",'port_reeq_ds_real':ds.real,'port_reeq_ds_imag':ds.imag,'dae_dlambda_real':dl.real,'dae_dlambda_imag':dl.imag,'error_abs':abs(ds-dl)})
  print(rows[-1],flush=True)
 except Exception as e: print('ERR',li,e,flush=True)
df=pd.DataFrame(rows);df.to_csv(DATA/'F2c_holdout_reequilibrated_port_line_sensitivity.csv',index=False)
from scipy.stats import spearmanr
print('rho',spearmanr(df.port_reeq_ds_real,df.dae_dlambda_real).statistic,'maxerr',df.error_abs.max(),'relmed',np.median(df.error_abs/np.maximum(np.abs(df.dae_dlambda_real+1j*df.dae_dlambda_imag),1e-12)))
