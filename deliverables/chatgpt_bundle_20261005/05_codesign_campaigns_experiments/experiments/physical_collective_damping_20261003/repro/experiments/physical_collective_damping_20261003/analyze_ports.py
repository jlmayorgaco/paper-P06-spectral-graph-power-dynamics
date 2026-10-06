"""Exact exponential local roots and physical finite-frequency port diagnostics."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
import json, hashlib, subprocess, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import linalg as la, optimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent
ROOT=OUT.parent.parent

class Model:
    def __init__(self,key):
        self.key=key
        p=OUT/'models'/key
        for name in ['A','A0','B','C']:
            setattr(self,name,pd.read_csv(p/(name+'.csv')).to_numpy())
        self.ports=pd.read_csv(p/'ports.csv')
        self.n=self.ports.state_index.to_numpy(int)
        self.h=np.setdiff1d(np.arange(len(self.A)),self.n)
        self.M=self.ports.M.to_numpy()
        self.I=np.eye(len(self.A))
    def delta(self,s,tau):
        return s*self.I-self.A0-np.exp(-s*tau)*(self.B@self.C)
    def refine(self,s,tau):
        for iteration in range(50):
            rb=la.solve(s*self.I-self.A0,self.B,check_finite=False)
            e=np.exp(-s*tau)
            K=np.eye(10)-e*(self.C@rb)
            u,z,vh=la.svd(K,check_finite=False)
            right=vh[-1].conj();left=u[:,-1]
            if z[-1]<3e-11: return s,float(z[-1])
            Ks=e*(tau*self.C@rb+self.C@la.solve(s*self.I-self.A0,rb,check_finite=False))
            ds=-np.vdot(left,K@right)/np.vdot(left,Ks@right)
            if abs(ds)>15: ds*=15/abs(ds)
            s+=ds
            if not np.isfinite(s) or abs(s)>1e5: break
        raise RuntimeError(f'root failed {self.key} {tau} {s} {z[-1]}')
    def impedance(self,s,tau):
        D=self.delta(s,tau)
        cc=D[np.ix_(self.h,self.h)]
        cn=D[np.ix_(self.h,self.n)]
        S=D[np.ix_(self.n,self.n)]-D[np.ix_(self.n,self.h)]@la.solve(cc,cn,check_finite=False)
        return self.M[:,None]*S
    def follow(self,s,t0,t1):
        t=t0;direction=np.sign(t1-t0);step=.0001
        for _ in range(3000):
            if abs(t1-t)<1e-12: return s
            rb=la.solve(s*self.I-self.A0,self.B,check_finite=False)
            e=np.exp(-s*t);G=self.C@rb
            u,sv,vh=la.svd(np.eye(10)-e*G,check_finite=False)
            z=vh[-1].conj();left=u[:,-1];v0=rb@z;v0/=la.norm(v0)
            Ks=e*(t*G+self.C@la.solve(s*self.I-self.A0,rb,check_finite=False))
            slope=-np.vdot(left,(s*e*G)@z)/np.vdot(left,Ks@z)
            dt=direction*min(step,abs(t1-t));pred=s+dt*slope
            nxt,_=self.refine(pred,t+dt)
            rb1=la.solve(nxt*self.I-self.A0,self.B,check_finite=False)
            _,_,vh1=la.svd(np.eye(10)-np.exp(-nxt*(t+dt))*(self.C@rb1),check_finite=False)
            v1=rb1@vh1[-1].conj();v1/=la.norm(v1)
            if abs(np.vdot(v0,v1))<.95 or abs(nxt-pred)>.1:
                step/=2
                if step<1e-8: raise RuntimeError('branch tracking unresolved')
                continue
            t+=dt;s=nxt;step=min(.0001,step*1.5)
        raise RuntimeError('continuation iteration budget')
        return s
    def mode(self,s,tau):
        D=self.delta(s,tau)
        u,z,vh=la.svd(D,check_finite=False)
        rb=la.solve(s*self.I-self.A0,self.B,check_finite=False)
        _,_,kh=la.svd(np.eye(10)-np.exp(-s*tau)*(self.C@rb),check_finite=False)
        v=rb@kh[-1].conj() if tau>0 else vh[-1].conj()
        w=u[:,-1]
        # Fixed PLL angle scaling makes physical perturbations comparable.
        v/=max(abs(v[self.ports.pll_angle_index.to_numpy(int)]))
        return v,w,float(la.norm(D@v)/(la.norm(D)*la.norm(v)))

def decomposition(m,s,tau,v=None):
    Z=m.impedance(1j*s.imag,tau)
    H=(Z+Z.conj().T)/2
    if v is None: v=m.mode(s,tau)[0]
    q=v[m.n];q=q/la.norm(q)
    total=float(np.vdot(q,H@q).real)
    nodal=float(np.sum(abs(q)**2*np.diag(H).real))
    cross=total-nodal
    material=abs(cross)/max(1e-20,abs(nodal)+abs(cross))
    # Baseline decomposition evaluated in exactly the same physical direction.
    Z0=m.impedance(1j*s.imag,0)
    H0=(Z0+Z0.conj().T)/2
    delayloss=float(np.vdot(q,(H0-H)@q).real)
    D=m.delta(s,tau);cc=D[np.ix_(m.h,m.h)]
    return dict(mode_damping_total=total,mode_damping_nodal=nodal,
        mode_damping_cross=cross,cross_absolute_fraction=material,
        diagonal_sign_wrong=bool(total*nodal<0),delay_damping_loss=delayloss,
        lambda_min_H=float(la.eigvalsh(H)[0]),
        offdiag_fraction=float(la.norm(H-np.diag(np.diag(H)))/la.norm(H)),
        speed_port_norm=float(la.norm(v[m.n])),hidden_condition=float(np.linalg.cond(cc)),
        schur_root_residual=float(la.norm(m.impedance(s,tau)@q)/max(1,la.norm(m.impedance(s,tau)))))

def run():
    gates=pd.read_csv(OUT/'TABLE_01_PARITY.csv')
    if not ((gates.equilibrium_inf<1e-7)&(gates.jacobian_relative<1e-8)&(gates.eigenvalue_max_set_difference<1e-5)).all():
        raise RuntimeError('BLOCKED_BASELINE_PARITY')
    seeds=pd.read_csv(ROOT/'experiments/optimal_latency_margin_20261002/T08_FAST_MODE_PROVENANCE.csv')
    seeds=seeds[np.isclose(seeds.tau_ms,40)]
    seeds=seeds.root_real.to_numpy()+1j*seeds.root_imag.to_numpy()
    roots=[];hopfs=[];scans=[];checks=[];nonlin=[]
    for key in ['N','Z','T']:
        m=Model(key)
        found=[]
        if key in ['N','Z']:
            prior=pd.read_csv(ROOT/'experiments/latency_robust_pll_codesign_20261003/L0_REPRODUCTION.csv')
            p=prior[prior.design_id==key].iloc[0]
            s0=complex(p.critical_root_real,p.critical_root_imag)
            s0,_=m.refine(s0,p.tau_crit_local_ms/1000)
            found.append(m.follow(s0,p.tau_crit_local_ms/1000,.04))
        else:
            prior=pd.read_csv(ROOT/'experiments/latency_robust_pll_closure_20261003/evaluations/full20_step15_medium/ROOTS.csv')
            for p in prior.itertuples():
                s0=complex(p.critical_real,p.critical_imag)
                found.append(m.follow(s0,p.local_crossing_ms/1000,.04))
        # Deterministic numerical root-discovery augmentation, same for all designs.
        # Old seeds alone found only 3 N branches; log in NUMERICAL_NOTES.md.
        discovery=list(seeds)+list(-.05+2j*np.pi*np.arange(2.,12.01,.1))
        for s0 in discovery:
            try: s,_=m.refine(s0,.04)
            except RuntimeError: continue
            if s.imag>1 and all(abs(s-r)>1e-4 for r in found): found.append(s)
        assert len(found)>=5, f'insufficient branch discovery {key}: {len(found)}'
        found.sort(key=lambda x:x.imag)
        for branch,s40 in enumerate(found):
            for tau in [.02,.03,.04,.045]:
                s=m.follow(s40,.04,tau);v,w,res=m.mode(s,tau)
                roots.append(dict(design=key,branch=branch,tau_ms=tau*1000,real=s.real,imag=s.imag,
                    frequency_hz=s.imag/(2*np.pi),residual=res,**decomposition(m,s,tau,v)))
            def growth(t): return m.follow(s40,.04,t).real
            try:
                t=optimize.brentq(growth,.02,.07,xtol=1e-12)
                s=m.follow(s40,.04,t)
                assert abs(s.real)<1e-6, 'root crossing not converged'
                v,w,res=m.mode(s,t)
                hopfs.append(dict(design=key,branch=branch,tau_ms=t*1000,real=s.real,imag=s.imag,
                    frequency_hz=s.imag/(2*np.pi),residual=res,**decomposition(m,s,t,v)))
            except (ValueError,RuntimeError) as e:
                hopfs.append(dict(design=key,branch=branch,error=str(e)))
        # Physical transfer and harmonic-work equality away from poles.
        Bd=np.zeros((len(m.A),10));Bd[m.n,np.arange(10)]=1/m.M
        rng=np.random.default_rng(20261003)
        for tau in [0,.02,.03,.04,.045]:
            for f in np.geomspace(.01,20,120):
                Z=m.impedance(2j*np.pi*f,tau);H=(Z+Z.conj().T)/2
                scans.append(dict(design=key,tau_ms=tau*1000,frequency_hz=f,
                    lambda_min_H=float(la.eigvalsh(H)[0]),
                    offdiag_fraction=la.norm(H-np.diag(np.diag(H)))/la.norm(H)))
            for f in [.1,1,4.5,10]:
                s=2j*np.pi*f;Z=m.impedance(s,tau)
                Y=la.solve(m.delta(s,tau),Bd)[m.n]
                q=rng.normal(size=10)+1j*rng.normal(size=10);q/=la.norm(q)
                d=Z@q;H=(Z+Z.conj().T)/2
                phase=np.linspace(0,2*np.pi,10000,endpoint=False)
                work=np.mean(np.sum(np.real(d[:,None]*np.exp(1j*phase))*np.real(q[:,None]*np.exp(1j*phase)),axis=0))
                formula=np.vdot(q,H@q).real/2
                checks.append(dict(design=key,tau_ms=tau*1000,frequency_hz=f,
                    impedance_transfer_error=la.norm(Y@Z-np.eye(10))/np.sqrt(10),
                    harmonic_work_error=abs(work-formula)/max(1,abs(formula)),
                    hidden_condition=np.linalg.cond(m.delta(s,tau)[np.ix_(m.h,m.h)])))
        if key in ['N','T']:
            for tau in [0,.04,.045]:
                if tau==0:
                    ev=la.eigvals(m.A);ev=ev[(ev.imag>0)&(np.abs(ev)>1e-6)]
                    s=ev[np.argmax(ev.real)]
                else: s=max((m.follow(r,.04,tau) for r in found),key=lambda x:x.real)
                v,w,res=m.mode(s,tau)
                p=OUT/'nonlinear_inputs'/f'{key}_{round(tau*1000)}ms';p.mkdir(parents=True,exist_ok=True)
                pd.DataFrame(dict(v_real=v.real,v_imag=v.imag,w_real=w.real,w_imag=w.imag)).to_csv(p/'mode.csv',index=False)
                (p/'case.json').write_text(json.dumps(dict(design=key,tau=tau,real=s.real,imag=s.imag,residual=res),indent=2))
                (p/'case.toml').write_text(f'design = "{key}"\ntau = {tau}\nreal = {s.real}\nimag = {s.imag}\n')
                nonlin.append(dict(case=p.name,design=key,tau_ms=tau*1000,real=s.real,imag=s.imag,residual=res))
        print('DONE',key,'branches',len(found),flush=True)
        for name,rows in [('TABLE_02_PORT_CHECKS',checks),('TABLE_03_TRACKED_MODES',roots),('TABLE_04_HOPF',hopfs),('TABLE_05_FREQUENCY_SCAN',scans),('TABLE_06_NONLINEAR_INPUTS',nonlin)]:
            pd.DataFrame(rows).to_csv(OUT/(name+'.csv'),index=False)
    plot()

def plot():
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':150})
    r=pd.read_csv(OUT/'TABLE_03_TRACKED_MODES.csv')
    fig,axs=plt.subplots(1,3,figsize=(13,3.6),layout='constrained')
    for key,label in [('N','Nominal gains'),('Z','Zero-delay tuning'),('T','Latency tuning')]:
        a=r[r.design==key].sort_values('real').groupby('tau_ms').tail(1).sort_values('tau_ms')
        axs[0].plot(a.tau_ms,a.real,'o-',label=label)
        axs[1].plot(a.tau_ms,a.mode_damping_total,'o-',label=label)
        axs[2].plot(a.tau_ms,a.cross_absolute_fraction,'o-',label=label)
    for ax in axs: ax.set_xlabel('Uniform PLL error delay (ms)');ax.grid(alpha=.2)
    axs[0].axhline(0,color='k',ls=':',lw=1);axs[0].set_ylabel('Tracked PLL mode growth (1/s)')
    axs[1].axhline(0,color='k',ls=':',lw=1);axs[1].set_ylabel('Physical projected damping (MW)')
    axs[2].set_ylabel('|Cross| / (|Nodal| + |Cross|)');axs[0].legend(fontsize=8)
    fig.suptitle('Fixed replacement: 88.455% — physical torque/speed ports, exact delay characteristic')
    fig.savefig(OUT/'FIG_01_MECHANISM.png');fig.savefig(OUT/'FIG_01_MECHANISM.pdf');plt.close(fig)
    a=r[r.tau_ms==45].sort_values('real').groupby('design').tail(1).sort_values('design')
    fig,ax=plt.subplots(figsize=(7,4),layout='constrained');x=np.arange(len(a));width=.25
    for j,(c,l) in enumerate([('mode_damping_nodal','Nodal diagonal'),('mode_damping_cross','Collective cross term'),('mode_damping_total','Total')]):
        ax.bar(x+(j-1)*width,a[c],width,label=l)
    ax.set_xticks(x,a.design);ax.axhline(0,color='k',lw=.8);ax.set_ylabel('Projected harmonic damping (MW)')
    ax.set_title('45 ms: damping on the tracked critical PLL motion');ax.legend(fontsize=8)
    fig.savefig(OUT/'FIG_02_NODAL_VS_COLLECTIVE.png');fig.savefig(OUT/'FIG_02_NODAL_VS_COLLECTIVE.pdf');plt.close(fig)

if __name__=='__main__': run()
