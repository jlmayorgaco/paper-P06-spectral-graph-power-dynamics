import os,sys,itertools,json
from pathlib import Path
import numpy as np,pandas as pd
REPO=Path('/mnt/data/research_unzip/research'); sys.path[:0]=[str(REPO/'src'),str(REPO/'experiments')]; os.chdir(REPO)
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case
from ibr_cycles.models.port_admittance import build_action_space

def auc_rank(y, score):
 y=np.asarray(y,bool); s=np.asarray(score,float); p=s[y]; n=s[~y]
 if len(p)==0 or len(n)==0:return float('nan')
 return float((p[:,None]>n[None,:]).mean()+0.5*(p[:,None]==n[None,:]).mean())

OUT=Path('/mnt/data/dynamic_forest_ieee39_validation'); DATA=OUT/'data'; DATA.mkdir(parents=True,exist_ok=True)
CANDS=tuple(range(30,39))
base=solve_case(ReplacementPlan.of({}))
all9=solve_case(ReplacementPlan.of({b:1 for b in CANDS}))
sp=build_action_space(base,all9,CANDS)
census=pd.read_csv(REPO/'results/tables/E12_compatibility_census_portfolios.csv')
uc=pd.read_csv(REPO/'results/UC/UC01/UC01_census_quantities.csv')[['members','replaced_pg_mw']]
census=census.merge(uc,on='members',how='left')
rows=[]
for r in census[census['size']==4].itertuples():
 S=tuple(map(int,str(r.members).split('+')))
 f=float(r.critical_frequency_hz)
 if not np.isfinite(f) or f<1e-6: f=0.65
 ss=1j*2*np.pi*f
 M=sp.m(ss); K=sp.k(ss); dY=sp.update(ss)
 ids=[]
 for b in S:
  k=CANDS.index(b); ids += [2*k,2*k+1]
 ids=np.array(ids)
 loc=np.linalg.det(dY[np.ix_(ids,ids)]); net=np.linalg.det(K[np.ix_(ids,ids)]); high=np.linalg.det(M[np.ix_(ids,ids)])
 rows.append({'members':r.members,'unstable':bool(r.unstable),'alpha':float(r.spectral_abscissa),'freq_hz':f,
              'local_logabs':np.log10(max(abs(loc),1e-300)),'network_logabs':np.log10(max(abs(net),1e-300)),
              'fullorder_logabs':np.log10(max(abs(high),1e-300)),'replaced_pg_mw':r.replaced_pg_mw,
              'factor_error':abs(high-loc*net)/max(abs(high),1e-300)})
df=pd.DataFrame(rows); df.to_csv(DATA/'F3_size4_selfenergy_vs_network_factor.csv',index=False)
aucs={c:auc_rank(df.unstable,df[c]) for c in ['local_logabs','network_logabs','fullorder_logabs','replaced_pg_mw']}
print('counts',df.unstable.value_counts().to_dict())
print('AUC raw',aucs)
print('AUC best orientation',{k:max(v,1-v) for k,v in aucs.items()})
print('factor error',df.factor_error.max())
print(df.sort_values('network_logabs',ascending=False).head(15).to_string(index=False))
(Path('/mnt/data/dynamic_forest_ieee39_validation/F3_summary.json')).write_text(json.dumps({'auc_raw':aucs,'auc_best':{k:max(v,1-v) for k,v in aucs.items()},'unstable':int(df.unstable.sum()),'n':len(df),'max_factor_error':float(df.factor_error.max())},indent=2))
