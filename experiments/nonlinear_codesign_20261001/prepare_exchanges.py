"""Finite experiments for the derived, local gain-compensation formula."""
import csv,json,tomllib
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/'reports/nonlinear_codesign_20261001'
audit=BASE/'stationarity';out=BASE/'exchange_checkpoint';out.mkdir(exist_ok=True)
d=tomllib.loads((audit/'0001_grad_constraints.toml').read_text())
x=np.array(d['parameters']);g=np.array(d['constraints'])
A=np.genfromtxt(audit/'0001_grad_jacobian.csv',delimiter=',',skip_header=1)
active=np.flatnonzero(g>=-.03);Ak=A[active,10:]
weights=np.array([float(r['P_gen_MW']) for r in csv.DictReader((ROOT/'reports/experiment_N/TABLE_N01_original_operating_point.csv').open())])
kp0=10*np.pi;ki0=kp0**2/4
lo=np.r_[np.full(10,np.log(kp0/4)),np.full(10,np.log(ki0/4))]
hi=np.r_[np.full(10,np.log(kp0*4)),np.full(10,np.log(ki0*4))]
lower=np.maximum(lo-x[10:],-.25);upper=np.minimum(hi-x[10:],.25)
rows=[{'id':'anchor','increase_bus':0,'decrease_bus':0,'MW':0.,'tuned':False,**{f'p{i+1}':v for i,v in enumerate(x)}}]
for inc,dec in [(31,32),(35,36)]:
    for mw in [.1,1.]:
        dr=np.zeros(10);dr[inc-30]=mw/weights[inc-30];dr[dec-30]=-mw/weights[dec-30]
        target=A[active,:10]@dr
        opt=minimize(lambda z:float(np.sum((target+Ak@z)**2)+1e-12*(z@z)),np.zeros(20),
            jac=lambda z:2*Ak.T@(target+Ak@z)+2e-12*z,bounds=list(zip(lower,upper)),
            method='L-BFGS-B',options={'ftol':1e-15,'gtol':1e-12,'maxiter':400})
        for tuned in [False,True]:
            p=x+np.r_[dr,opt.x if tuned else np.zeros(20)]
            rows.append({'id':f'{inc}_{dec}_{mw}_{"tuned" if tuned else "fixed"}',
                'increase_bus':inc,'decrease_bus':dec,'MW':mw,'tuned':tuned,
                **{f'p{i+1}':v for i,v in enumerate(p)}})
with (out/'parameters.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
(out/'scope.json').write_text(json.dumps({'gradient_parameters':x.tolist(),'active_rows_zero_based':active.tolist(),
    'claim':'Finite validation of a local constraint-signature compensation, not complete trajectory equivalence.'},indent=2))
