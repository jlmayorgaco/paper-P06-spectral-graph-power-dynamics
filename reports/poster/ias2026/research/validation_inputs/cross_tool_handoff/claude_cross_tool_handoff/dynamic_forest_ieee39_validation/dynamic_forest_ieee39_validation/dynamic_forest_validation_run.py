from __future__ import annotations
import copy, itertools, json, math, os, sys, time
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.linalg import eig
from scipy.stats import spearmanr

REPO = Path('/mnt/data/research_unzip/research')
sys.path[:0] = [str(REPO/'src'), str(REPO/'experiments')]
os.chdir(REPO)

from _f7_common import Theta, solve_subset, LEAK
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case
from ibr_cycles.models.ieee39_devices import ConverterParameters
from ibr_cycles.models.ieee39_network import load_network
from ibr_cycles.models.port_admittance import build_action_space, PortActionSpace, _real_ybus
from ibr_cycles.certification.symmetry import rotation_generator, frequency_partner
from ibr_cycles.certification.transverse import transverse_operator
from ibr_cycles.units.quantities import per_machine_table

OUT = Path('/mnt/data/dynamic_forest_ieee39_validation')
DATA=OUT/'data'; FIG=OUT/'figures'
DATA.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------ utilities

def transverse_spectrum_case(case):
    rx,_=rotation_generator(case.dae, case.equilibrium.z)
    w=frequency_partner(case.dae).w
    tr=transverse_operator(case.system.A, rx, w)
    return np.linalg.eigvals(tr.a_perp)

def tracked_transverse_eig(case, target):
    ev=transverse_spectrum_case(case)
    # positive imag if target positive imag
    if np.imag(target) >= 0:
        cand=ev[ev.imag>=0]
    else:
        cand=ev[ev.imag<=0]
    return complex(cand[np.argmin(np.abs(cand-target))])

def generic_case(members, theta:Theta, network=None):
    return solve_case(
        ReplacementPlan.of({b:1.0 for b in members}),
        network=network,
        converter=ConverterParameters(voltage_control=True, voltage_gain=theta.g, voltage_leak=LEAK),
        machine_scaling={'ka':theta.k,'ta':theta.t},
    )

def q_matrix(space, s):
    m=space.m(s)
    n=m.shape[0]
    total=np.eye(n,dtype=complex)+m
    selfb=np.zeros_like(total)
    for k in range(space.order):
        sl=slice(2*k,2*k+2)
        selfb[sl,sl]=total[sl,sl]
    return np.linalg.solve(selfb,total)-np.eye(n,dtype=complex)

def eig_near_minus1(space,s):
    q=q_matrix(space,s)
    vals,vl,vr=eig(q,left=True,right=True)
    i=int(np.argmin(np.abs(vals+1)))
    return vals[i], vl[:,i], vr[:,i]

def slogdet_ratio(num, den):
    sn,ln=np.linalg.slogdet(num); sd,ld=np.linalg.slogdet(den)
    return complex(sn/sd*np.exp(ln-ld))

def auc_rank(y, score):
    y=np.asarray(y,dtype=bool); score=np.asarray(score,float)
    pos=score[y]; neg=score[~y]
    if len(pos)==0 or len(neg)==0: return float('nan')
    comp=(pos[:,None]>neg[None,:]).mean()
    ties=(pos[:,None]==neg[None,:]).mean()
    return float(comp+0.5*ties)

# ------------------------------------------------------------ line builder
payload=json.loads((REPO/'configs/ias2026/ieee39_network.json').read_text())
base_network=load_network()

def ybus_scaled(line_idx:int, gamma:float):
    order=[int(b['idx']) for b in payload['buses']]; index={b:i for i,b in enumerate(order)}
    n=len(order); y=np.zeros((n,n),complex)
    for ell,line in enumerate(payload['lines']):
        if float(line['u'])==0: continue
        f,t=index[int(line['bus1'])],index[int(line['bus2'])]
        series=1.0/complex(line['r'],line['x']); charging=complex(line['g'],line['b'])/2
        scale=gamma if ell==line_idx else 1.0
        series*=scale; charging*=scale
        m=float(line['tap'])*np.exp(1j*float(line['phi'])); m2=abs(m)**2
        y[f,f]+=(series+charging)/m2; y[t,t]+=series+charging
        y[f,t]+=-series/np.conj(m); y[t,f]+=-series/m
    for sh in payload['shunts']:
        y[index[int(sh['bus'])],index[int(sh['bus'])]]+=complex(sh['g'],sh['b'])
    return y

line_labels=[]
for i,l in enumerate(payload['lines']):
    line_labels.append(f"L{i:02d}:{int(l['bus1'])}-{int(l['bus2'])}")

# ------------------------------------------------------------ F0 matrix-tree / forest enumeration on complex graphs

def forest_sum(weights, roots):
    # nodes 0..n, 0 ground. roots have no outgoing parent. Each nonroot chooses parent != itself.
    n=weights.shape[0]-1
    roots=set(roots)
    non=[i for i in range(n+1) if i not in roots]
    choices={i:[j for j in range(n+1) if j!=i and abs(weights[i,j])>0] for i in non}
    total=0j
    for vals in itertools.product(*(choices[i] for i in non)):
        parent=dict(zip(non,vals))
        ok=True
        for i in non:
            seen=set(); cur=i
            while cur not in roots:
                if cur in seen or cur not in parent:
                    ok=False; break
                seen.add(cur); cur=parent[cur]
            if not ok: break
        if ok:
            w=1+0j
            for i,j in parent.items(): w*=weights[i,j]
            total+=w
    return total

rng=np.random.default_rng(20260911)
rows=[]
for rep in range(60):
    n=4
    W=np.zeros((n+1,n+1),complex)
    for i in range(1,n+1):
        for j in range(0,n+1):
            if i==j: continue
            if rng.random()<0.75:
                W[i,j]=(0.2+rng.random())+1j*rng.normal(0,0.3)
        if not np.any(W[i]): W[i,0]=1+0.1j
    A=np.zeros((n,n),complex)
    for i in range(1,n+1):
        A[i-1,i-1]=sum(W[i,j] for j in range(n+1) if j!=i)
        for j in range(1,n+1):
            if i!=j: A[i-1,j-1]=-W[i,j]
    tree=forest_sum(W,{0})
    e0=abs(np.linalg.det(A)-tree)/max(1,abs(np.linalg.det(A)))
    for rsize in (1,2):
        for I0 in itertools.combinations(range(1,n+1),rsize):
            keep=[i-1 for i in range(1,n+1) if i not in I0]
            minor=A[np.ix_(keep,keep)]
            fs=forest_sum(W,{0,*I0})
            err=abs(np.linalg.det(minor)-fs)/max(1,abs(np.linalg.det(minor)))
            rows.append({'rep':rep,'roots':'+'.join(map(str,I0)),'tree_error':e0,'forest_error':err})
F0=pd.DataFrame(rows); F0.to_csv(DATA/'F0_complex_matrix_tree.csv',index=False)

# ------------------------------------------------------------ event list and Jacobi/forest-ratio identities
EVENTS=[
 ('E1',(30,33,35),0.0102,1.0,0.852,1.0,0.628),
 ('E2',(30,33,37),0.0188,1.0,0.852,1.0,0.633),
 ('E3',(33,35,37),0.0609,1.0,0.852,1.0,0.671),
 ('E4',(30,33),0.1034,1.0,0.852,1.0,0.656),
 ('E5',(33,35,37),0.1870,1.0,0.852,1.0,0.695),
 ('E6',(30,33),0.2491,1.0,0.852,1.0,0.672),
 ('E7',(30,33,35),0.2778,1.0,0.852,1.0,0.706),
 ('E8',(30,33,37),0.2931,1.0,0.852,1.0,0.685),
 ('E9',(30,33,35,37),0.2973,1.0,0.852,1.0,0.718),
 ('F',(30,33,35,37),0.20768,1.425,1.5,1.0,0.706),
]

jac_rows=[]; spaces={}
for name,S,g,k,t,h,f in EVENTS:
    th=Theta(g,k,t,h)
    base=solve_subset((),th); full=solve_subset(S,th)
    sp=build_action_space(base,full,S); spaces[name]=(sp,base,full)
    ss=1j*2*np.pi*f; T=sp.t0.evaluate(ss); K=sp.k(ss); dY=sp.update(ss); M=sp.m(ss)
    sel=sp.selector(); inds=np.flatnonzero(np.any(np.abs(sel)>0,axis=1)); comp=np.array([i for i in range(T.shape[0]) if i not in set(inds)])
    ratio=slogdet_ratio(T[np.ix_(comp,comp)],T)
    dk=np.linalg.det(K)
    rel=abs(dk-ratio)/max(1e-300,abs(dk),abs(ratio))
    local=np.linalg.det(dY); network=dk; high=np.linalg.det(M)
    jac_rows.append({'event':name,'size':len(S),'freq_hz':f,'jacobi_relerr':rel,
                     'local_abs':abs(local),'network_forest_ratio_abs':abs(network),'full_order_abs':abs(high),
                     'factor_relerr':abs(high-local*network)/max(1e-300,abs(high))})
JAC=pd.DataFrame(jac_rows); JAC.to_csv(DATA/'F1_ieee39_jacobi_forest_ratio.csv',index=False)

# complex-linearity of offdiagonal bus blocks at flagship
sp,base,full=spaces['F']; ss=1j*2*np.pi*0.706; T=sp.t0.evaluate(ss)
cl=[]
for i in range(base.dae.network.n_bus):
  for j in range(base.dae.network.n_bus):
    B=T[2*i:2*i+2,2*j:2*j+2]
    if np.linalg.norm(B)==0: continue
    a=0.5*(B[0,0]+B[1,1]); b=0.5*(B[1,0]-B[0,1])
    P=np.array([[a,-b],[b,a]],complex)
    res=np.linalg.norm(B-P)/np.linalg.norm(B)
    cl.append({'i_bus':base.dae.network.bus_idx[i],'j_bus':base.dae.network.bus_idx[j],
               'diagonal':i==j,'norm':np.linalg.norm(B),'complex_linear_residual':res})
CL=pd.DataFrame(cl); CL.to_csv(DATA/'F1_complex_linearity_blocks.csv',index=False)

# ------------------------------------------------------------ flagship forest leverage and line sensitivity
sp,base,full=spaces['F']; sstar=1j*2*np.pi*0.706; T=sp.t0.evaluate(sstar); K=sp.k(sstar); sel=sp.selector()
inds=np.flatnonzero(np.any(np.abs(sel)>0,axis=1)); iset=set(inds); comp=np.array([i for i in range(T.shape[0]) if i not in iset])
Tcomp=T[np.ix_(comp,comp)]
# eigenvalue nearest -1 and mu_s
mu0,pl,pr=eig_near_minus1(sp,sstar); normlr=np.vdot(pl,pr); pl=pl/np.conj(normlr)
hs=1e-5
mup=eig_near_minus1(sp,sstar+hs)[0]; mum=eig_near_minus1(sp,sstar-hs)[0]
mu_s=(mup-mum)/(2*hs)

lev_rows=[]
eps=1e-4
for li,label in enumerate(line_labels):
    yp=ybus_scaled(li,1+eps); ym=ybus_scaled(li,1-eps)
    dYreal=(_real_ybus(yp)-_real_ybus(ym))/(2*eps)
    dT=dYreal
    dTc=dT[np.ix_(comp,comp)]
    lev=np.trace(np.linalg.solve(Tcomp,dTc))-np.trace(np.linalg.solve(T,dT))
    # FD check of log det K using fixed ports
    def sp_y(y):
        return PortActionSpace(replace(sp.t0,ybus_real=_real_ybus(y)), replace(sp.ts,ybus_real=_real_ybus(y)), sp.buses)
    spp=sp_y(yp); spm=sp_y(ym)
    kp=np.linalg.det(spp.k(sstar)); km=np.linalg.det(spm.k(sstar))
    # branch-safe complex log ratio
    lev_fd=(np.log(kp/K.det() if False else kp/np.linalg.det(K))-np.log(km/np.linalg.det(K)))/(2*eps)
    # mu derivative and implicit s derivative at frozen equilibrium
    mu_p=eig_near_minus1(spp,sstar)[0]; mu_m=eig_near_minus1(spm,sstar)[0]
    mu_g=(mu_p-mu_m)/(2*eps); ds=-mu_g/mu_s
    lev_rows.append({'line_index':li,'line':label,'bus1':int(payload['lines'][li]['bus1']),'bus2':int(payload['lines'][li]['bus2']),
                     'forest_leverage_real':lev.real,'forest_leverage_imag':lev.imag,'forest_leverage_abs':abs(lev),
                     'forest_fd_error':abs(lev-lev_fd),'port_ds_dgamma_real':ds.real,'port_ds_dgamma_imag':ds.imag,
                     'port_alpha_sensitivity':ds.real})
LEV=pd.DataFrame(lev_rows)

# Re-equilibrated full DAE line sensitivities for all lines; physical branch scale, same policy/dispatch.
thetaF=Theta(0.20768,1.425,1.5,1.0)
lam_target=tracked_transverse_eig(full,1j*2*np.pi*0.706)
fd_eps=0.002
phys=[]
for li,label in enumerate(line_labels):
    try:
        npw=replace(base_network,ybus=ybus_scaled(li,1+fd_eps))
        nmw=replace(base_network,ybus=ybus_scaled(li,1-fd_eps))
        cp=generic_case((30,33,35,37),thetaF,npw); cm=generic_case((30,33,35,37),thetaF,nmw)
        lp=tracked_transverse_eig(cp,lam_target); lm=tracked_transverse_eig(cm,lam_target)
        der=(lp-lm)/(2*fd_eps)
        phys.append({'line_index':li,'dae_dlambda_real':der.real,'dae_dlambda_imag':der.imag,'dae_alpha_sensitivity':der.real,
                     'fd_plus_alpha':lp.real,'fd_minus_alpha':lm.real})
    except Exception as e:
        phys.append({'line_index':li,'dae_dlambda_real':np.nan,'dae_dlambda_imag':np.nan,'dae_alpha_sensitivity':np.nan,
                     'error':str(e)[:120]})
PHYS=pd.DataFrame(phys)
LEV=LEV.merge(PHYS,on='line_index',how='left')
LEV.to_csv(DATA/'F2_flagship_dynamic_forest_line_sensitivity.csv',index=False)
valid=LEV.dropna(subset=['dae_alpha_sensitivity'])
r_port=spearmanr(valid.port_alpha_sensitivity,valid.dae_alpha_sensitivity).statistic
r_lev=spearmanr(valid.forest_leverage_abs,np.abs(valid.dae_alpha_sensitivity)).statistic
# Top stabilizing overlap (more negative derivative means strengthening stabilizes)
top_port=set(valid.nsmallest(10,'port_alpha_sensitivity').line_index)
top_dae=set(valid.nsmallest(10,'dae_alpha_sensitivity').line_index)
top_lev=set(valid.nlargest(10,'forest_leverage_abs').line_index)

# ------------------------------------------------------------ 9-candidate size-4 census local self energy vs network factor
CANDS=tuple(range(30,39))
thetaP4=Theta(0.03625,1.425,1.5,1.0)
base9=generic_case((),thetaP4); all9=generic_case(CANDS,thetaP4)
space9=build_action_space(base9,all9,CANDS)
# metadata Pg per bus
mt=pd.DataFrame(per_machine_table(base9.dae.network)); pgb={int(r.bus):float(r.pg_mw) for r in mt.itertuples()} if 'pg_mw' in mt.columns else {}

census=[]
for idx,S in enumerate(itertools.combinations(CANDS,4)):
    c=generic_case(S,thetaP4)
    ev=transverse_spectrum_case(c)
    # rightmost oscillatory in 0.3-1.5Hz if exists; otherwise rightmost nonzero
    f=np.abs(ev.imag)/(2*np.pi); band=ev[(ev.imag>=0)&(f>=0.3)&(f<=1.5)]
    if len(band): lam=band[np.argmax(band.real)]
    else: lam=ev[np.argmax(ev.real)]
    ss=1j*abs(lam.imag) if abs(lam.imag)>1e-6 else 1j*2*np.pi*0.706
    m=space9.m(ss); K9=space9.k(ss); dY9=space9.update(ss)
    rows=[]
    for b in S:
        k=CANDS.index(b); rows.extend([2*k,2*k+1])
    rows=np.array(rows)
    loc=np.linalg.det(dY9[np.ix_(rows,rows)])
    net=np.linalg.det(K9[np.ix_(rows,rows)])
    high=np.linalg.det(m[np.ix_(rows,rows)])
    census.append({'subset':'+'.join(map(str,S)),'alpha':float(lam.real),'freq_hz':float(abs(lam.imag)/(2*np.pi)),
                   'unstable':bool(lam.real>0),'local_logabs':float(np.log10(max(abs(loc),1e-300))),
                   'network_logabs':float(np.log10(max(abs(net),1e-300))),
                   'fullorder_logabs':float(np.log10(max(abs(high),1e-300))),
                   'factor_error':float(abs(high-loc*net)/max(abs(high),1e-300)),
                   'dispatch_mw':sum(pgb.get(b,0.0) for b in S)})
CEN=pd.DataFrame(census); CEN.to_csv(DATA/'F3_size4_selfenergy_vs_network_factor.csv',index=False)
aucs={col:auc_rank(CEN.unstable,CEN[col]) for col in ['local_logabs','network_logabs','fullorder_logabs','dispatch_mw']}
# Also inverted AUC if direction unknown
aucs_best={k:max(v,1-v) if np.isfinite(v) else v for k,v in aucs.items()}

# Node scores from line sensitivities
nodes=[]
for b in base_network.bus_idx:
    inc=valid[(valid.bus1==b)|(valid.bus2==b)]
    nodes.append({'bus':b,'incident_abs_dae_sensitivity':float(np.nansum(np.abs(inc.dae_alpha_sensitivity))),
                  'incident_abs_forest_leverage':float(np.nansum(inc.forest_leverage_abs)),
                  'incident_stabilizing_strength':float(np.nansum(np.maximum(-inc.dae_alpha_sensitivity,0)))})
NODE=pd.DataFrame(nodes).sort_values('incident_abs_dae_sensitivity',ascending=False); NODE.to_csv(DATA/'F2_node_dynamic_risk_scores.csv',index=False)

summary={
 'F0_matrix_tree':{'cases':len(F0),'max_tree_error':float(F0.tree_error.max()),'max_forest_error':float(F0.forest_error.max())},
 'F1_ieee39_forest_ratio':{'events':len(JAC),'max_jacobi_relerr':float(JAC.jacobi_relerr.max()),'max_factor_relerr':float(JAC.factor_relerr.max()),
                           'offdiag_complex_linearity_max':float(CL.loc[~CL.diagonal,'complex_linear_residual'].max()),
                           'diag_complex_linearity_median':float(CL.loc[CL.diagonal,'complex_linear_residual'].median())},
 'F2_line_diagnosis':{'n_lines':len(valid),'spearman_port_vs_rebalanced_dae':float(r_port),'spearman_forestabs_vs_abs_dae':float(r_lev),
                      'top10_port_dae_overlap':len(top_port&top_dae),'top10_forest_dae_overlap':len(top_lev&top_dae),
                      'max_forest_fd_error':float(valid.forest_fd_error.max())},
 'F3_census':{'size4_portfolios':len(CEN),'unstable':int(CEN.unstable.sum()),'auc_raw':aucs,'auc_orientation_free':aucs_best,
              'max_factor_error':float(CEN.factor_error.max())},
}
(OUT/'SUMMARY.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
print('\nTop line sensitivities (DAE stabilizing by strengthening):')
print(valid.nsmallest(12,'dae_alpha_sensitivity')[['line','dae_alpha_sensitivity','port_alpha_sensitivity','forest_leverage_abs']].to_string(index=False))
print('\nTop node scores:')
print(NODE.head(12).to_string(index=False))
