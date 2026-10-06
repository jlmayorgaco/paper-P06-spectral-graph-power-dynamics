"""Direct full-matrix falsification of the finite latency Schur identity."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json,tomllib
import numpy as np
import pandas as pd
from scipy import linalg as la
OUT=Path(__file__).resolve().parent; ROOT=OUT.parents[1]
MODEL=ROOT/'experiments/graph_gsp_codesign_20261003/model'
M={k:pd.read_csv(MODEL/(k+'.csv')).to_numpy() for k in ['Adev','Bv','Cs','Cf','Ds','Y','Etheta','Hv','Bp','Bi']}
ports=pd.read_csv(MODEL/'ports.csv'); ii=ports.state_index.to_numpy(dtype=int)
zz=np.setdiff1d(np.arange(204),ii)
rows=[]
for name in ['baseline','analytic']:
    d=tomllib.loads((ROOT/'experiments/all_pll_gain_map_validation_20261003/designs'/f'{name}.toml').read_text())
    rho=np.array(d['rho']); w=np.repeat(1-rho,2)
    V=-la.solve(M['Y']+w[:,None]*M['Ds'],w[:,None]*M['Cs']+(1-w[:,None])*M['Cf'])
    C=M['Etheta']+M['Hv']@V; B=M['Bp']*d['Kp']+M['Bi']*d['Ki']
    A0=M['Adev']+M['Bv']@V
    mass=ports.M.to_numpy()*(1-rho)/(1-ports.rho.to_numpy())
    for freq in [.5,5.,10.]:
        s=2j*np.pi*freq; e=np.exp(-s*.04)
        D=s*np.eye(204)-A0-e*B@C
        nn,nz,zn,cc=D[np.ix_(ii,ii)],D[np.ix_(ii,zz)],D[np.ix_(zz,ii)],D[np.ix_(zz,zz)]
        invzn=la.solve(cc,zn); Z=mass[:,None]*(nn-nz@invzn)
        for j in range(10):
            b=B[:,j];c=C[j];k=e-np.exp(-s*.041)
            invb=la.solve(cc,b[zz]); h=c[zz]@invb
            u=mass*(b[ii]-nz@invb); vrow=c[ii]-c[zz]@invzn
            gamma=k/(1+k*h); update=gamma*np.outer(u,vrow)
            DD=D+k*np.outer(b,c)
            cc2=DD[np.ix_(zz,zz)]
            Z2=mass[:,None]*(DD[np.ix_(ii,ii)]-DD[np.ix_(ii,zz)]@la.solve(cc2,DD[np.ix_(zz,ii)]))
            direct=Z2-Z; DH=(direct+direct.conj().T)/2
            exact=(update+update.conj().T)/2
            eig=la.eigvalsh(DH); scale=max(la.norm(exact,2),1e-30)
            threshold=1e-7*scale
            err=la.norm(direct-update)/max(la.norm(update),1e-30)
            ranktail=max(abs(eig[1:-1]))/scale
            # v is the conjugate transpose of the analytic row vrow.
            v=vrow.conj();a=gamma*u;inner=np.vdot(v,a)
            rad=max(0,la.norm(a)**2*la.norm(v)**2-inner.imag**2)
            minus=(inner.real-np.sqrt(rad))/2;plus=(inner.real+np.sqrt(rad))/2
            angle=max(0,1-abs(np.vdot(u,v))**2/(la.norm(u)**2*la.norm(v)**2))
            rows.append(dict(design=name,frequency_Hz=freq,bus=30+j,
                tau_base_ms=40,delta_tau_ms=1,update_relative_error=err,
                hidden_condition_before=float(np.linalg.cond(cc)),hidden_condition_after=float(np.linalg.cond(cc2)),
                lambda_min=float(eig[0]),lambda_max=float(eig[-1]),predicted_min=float(minus),predicted_max=float(plus),
                rank_tail_relative=float(ranktail),path_sine_squared=float(angle),
                positive_count=int((eig>threshold).sum()),negative_count=int((eig < -threshold).sum()),
                pass_identity=bool(err<1e-7 and ranktail<1e-7),
                pass_sign=bool(eig[0]<-threshold and eig[-1]>threshold),
                status='NUMERICALLY_VALIDATED_POINTWISE'))
df=pd.DataFrame(rows);df.to_csv(OUT/'TABLE_02_FINITE_DELAY_RANK_LAW.csv',index=False)
summary={'cases':len(df),'identity_pass':bool(df.pass_identity.all()),'all_indefinite':bool(df.pass_sign.all()),
    'max_relative_update_error':float(df.update_relative_error.max()),'max_rank_tail':float(df.rank_tail_relative.max()),
    'min_path_sine_squared':float(df.path_sine_squared.min()),
    'scope':'60 frozen harmonic algebra checks, not a nonlinear or full-spectrum stability certificate'}
(OUT/'RANK_LAW_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
