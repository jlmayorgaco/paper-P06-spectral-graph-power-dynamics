"""Report stationarity and tuning-aware equal-MW exchange directions honestly."""
import csv,json,tomllib,argparse
from pathlib import Path
import numpy as np
from scipy.optimize import nnls,minimize
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'reports/nonlinear_codesign_20261001'
pa=argparse.ArgumentParser();pa.add_argument('candidate');pa.add_argument('--label',default='stationarity');args=pa.parse_args()
d=tomllib.loads(Path(args.candidate).read_text());out=BASE/args.label
cdata=tomllib.loads((out/'0001_grad_constraints.toml').read_text())
g=np.array(cdata['constraints']);x=np.array(cdata['parameters'])
A=np.genfromtxt(out/'0001_grad_jacobian.csv',delimiter=',',skip_header=1)
expected=np.r_[np.array(d['rho']),np.log(d['Kp']),np.log(d['Ki'])]
if not np.allclose(expected,x,rtol=0,atol=1e-12):raise RuntimeError('Candidate and gradient parameters differ')
weights=np.array([float(r['P_gen_MW']) for r in csv.DictReader((ROOT/'reports/experiment_N/TABLE_N01_original_operating_point.csv').open())])
cost=np.r_[-weights/weights.sum(),np.zeros(20)];scale=np.r_[np.ones(10),np.full(20,6.)]
kp0=10*np.pi;ki0=kp0**2/4
lo=np.r_[np.full(10,.001),np.full(10,np.log(kp0/4)),np.full(10,np.log(ki0/4))]
hi=np.r_[np.full(10,.999),np.full(10,np.log(kp0*4)),np.full(10,np.log(ki0*4))]
records=[]
for active_tol in (.001,.003,.01):
    ids=np.flatnonzero(g>=-active_tol);rows=[A[i]*scale for i in ids];names=[f'constraint_{i}' for i in ids]
    for j in range(30):
        if x[j]-lo[j]<1e-6:rows.append(-np.eye(30)[j]);names.append(f'lower_{j}')
        if hi[j]-x[j]<1e-6:rows.append(np.eye(30)[j]);names.append(f'upper_{j}')
    if rows:
        M=np.array(rows).T;mu=nnls(M,-cost*scale,maxiter=10000)[0]
        residual=cost*scale+M@mu
    else:mu=np.array([]);residual=cost*scale
    records.append({'active_tolerance':active_tol,'active_constraints':names,'multipliers':mu.tolist(),
        'stationarity_inf':float(max(abs(residual))),'relative_stationarity_inf':float(max(abs(residual))/max(abs(cost*scale))),
        'primal_violation':float(max(0.,max(g))),'constraint_complementarity_inf':float(max(abs(mu[:len(ids)]*g[ids]),default=0.))})
# Select near-active constraints to define a local security signature. This is
# not an equivalence claim for complete trajectories or finite exchanges.
active=np.flatnonzero(g>=-.03)
if len(active)==0:active=np.array([np.argmax(g)])
Ar=A[active,:10];Ak=A[active,10:];pairs=[]
lower=np.maximum(lo[10:]-x[10:],-.25);upper=np.minimum(hi[10:]-x[10:],.25)
for i in range(10):
    for j in range(10):
        if i==j:continue
        dr=np.zeros(10);dr[i]=1/weights[i];dr[j]=-1/weights[j]
        if not np.all((x[:10]+dr>=lo[:10])&(x[:10]+dr<=hi[:10])):continue
        target=Ar@dr
        uncon=-np.linalg.pinv(Ak,rcond=1e-9)@target
        opt=minimize(lambda z:float(np.sum((target+Ak@z)**2)+1e-12*(z@z)),np.zeros(20),
            jac=lambda z:2*Ak.T@(target+Ak@z)+2e-12*z,bounds=list(zip(lower,upper)),
            method='L-BFGS-B',options={'ftol':1e-15,'gtol':1e-12,'maxiter':200})
        pairs.append({'increase_GFL_bus':30+i,'decrease_GFL_bus':30+j,'exchanged_MW':1.,
            'raw_signature_norm':float(np.linalg.norm(target)),
            'unbounded_gain_change_norm':float(np.linalg.norm(uncon)),
            'bounded_signature_residual':float(np.linalg.norm(target+Ak@opt.x)),
            'gain_change_norm':float(np.linalg.norm(opt.x)),
            'finite_nonlinear_exchange_validated':False})
if pairs:
    with (out/'local_exchangeability.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=pairs[0]);writer.writeheader();writer.writerows(pairs)
report={'status':'LOCAL_DIAGNOSTIC_NOT_GLOBAL_OPTIMALITY_CERTIFICATE',
    'gradient_horizon_s':cdata['gradient_horizon_s'],'gradient_step_s':cdata['gradient_step_s'],
    'KKT_audits':records,'exchange_active_rows':active.tolist(),
    'cautions':['Peak ties and active-set changes require gradient bundles.',
                'KKT with a loose active tolerance is only a diagnostic.',
                'Exchange signatures require finite nonlinear validation.',
                'No global lower bound on retained SG has been computed.']}
(out/'stationarity_report.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
