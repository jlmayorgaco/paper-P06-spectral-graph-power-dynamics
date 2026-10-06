import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1'
import sys,json,time,tomllib,copy
from pathlib import Path
import numpy as np,pandas as pd
import cvxpy as cp
from scipy import linalg as la,optimize
OUT=Path(__file__).resolve().parent;ROOT=OUT.parent.parent
sys.path.insert(0,str(OUT.parent/'physical_collective_damping_20261003'))
from analyze_ports import Model as RootModel

def read(name):return pd.read_csv(OUT/'model'/(name+'.csv')).to_numpy()
def write_design(path,p):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text('\n'.join(f'{k} = [{", ".join(format(float(x),".17g") for x in v)}]' for k,v in [('rho',p[:10]),('Kp',np.exp(p[10:20])),('Ki',np.exp(p[20:]))])+'\ntau = 0.04\n')

class Model(RootModel):
    def __init__(self,p):
        self.key='parametric';self.p=np.array(p);self.rho=p[:10];self.kp=np.exp(p[10:20]);self.ki=np.exp(p[20:])
        for name in ['Adev','Bv','Cs','Cf','Ds','Y','Etheta','Hv','Bp','Bi']:setattr(self,name,read(name))
        self.ports=pd.read_csv(OUT/'model/ports.csv');self.n=self.ports.state_index.to_numpy(int)
        self.h=np.setdiff1d(np.arange(len(self.Adev)),self.n);self.I=np.eye(len(self.Adev))
        self.M=self.ports.M.to_numpy()*(1-self.rho)/(1-self.ports.rho.to_numpy())
        weight=np.repeat(1-self.rho,2)
        self.G=self.Y+weight[:,None]*self.Ds
        self.Vx=-la.solve(self.G,weight[:,None]*self.Cs+(1-weight[:,None])*self.Cf)
        self.A0=self.Adev+self.Bv@self.Vx;self.C=self.Etheta+self.Hv@self.Vx
        self.B=self.Bp*self.kp+self.Bi*self.ki;self.A=self.A0+self.B@self.C
    def gradients(self,s,tau):
        v,w,res=self.mode(s,tau);e=np.exp(-s*tau);den=np.vdot(w,(self.I+tau*e*self.B@self.C)@v)
        ds=[]
        for i in range(10):
            Gr=np.zeros_like(self.G);Cr=np.zeros_like(self.Cs);sl=slice(2*i,2*i+2)
            Gr[sl]=-self.Ds[sl];Cr[sl]=self.Cf[sl]-self.Cs[sl]
            vr=-la.solve(self.G,Gr@self.Vx+Cr)
            Dp=-(self.Bv+e*self.B@self.Hv)@vr
            ds.append(-np.vdot(w,Dp@v)/den)
        for Bg in [self.Bp*self.kp,self.Bi*self.ki]:
            for i in range(10):ds.append(e*np.vdot(w,Bg[:,i])*(self.C[i]@v)/den)
        return np.array(ds),res

def roots(m,extra=()):
    seeds=list(extra)
    old=pd.read_csv(OUT.parent/'physical_collective_damping_20261003/TABLE_03_TRACKED_MODES.csv')
    seeds+=list(old[old.tau_ms==40].real+1j*old[old.tau_ms==40].imag)
    seeds+=list(-.05+2j*np.pi*np.arange(2,12.01,.15))
    ev=la.eigvals(m.A);ev=ev[(ev.imag>=-1e-8)&(abs(ev)>1e-5)]
    seeds+=list(ev[np.argsort(ev.real)[-12:]])
    found=[]
    for s in seeds:
        try:r,_=m.refine(s,.04)
        except (RuntimeError,ValueError,la.LinAlgError):continue
        if r.imag<-1e-6:r=r.conjugate()
        if all(abs(r-x)>1e-4 for x in found):found.append(r)
    return sorted(found,key=lambda z:z.real,reverse=True)

def prepare():
    d=tomllib.loads((OUT/'baseline.toml').read_text());p=np.r_[d['rho'],np.log(d['Kp']),np.log(d['Ki'])]
    m=Model(p);parity=[]
    for j in [1,2]:
        pp=p.copy();pp[:10]=read('rho_check'+str(j)).ravel();a=Model(pp).A
        parity.append(dict(case=j,relative_error=la.norm(a-read('A_check'+str(j)))/la.norm(a)))
    pd.DataFrame(parity).to_csv(OUT/'TABLE_01_MODEL_PARITY.csv',index=False)
    assert max(r['relative_error'] for r in parity)<1e-8
    L=read('L');ev,U=la.eigh(L);one=np.ones(10)/np.sqrt(10)
    assert la.norm(L-L.T)<1e-8 and la.norm(L@one)<1e-8 and ev[0]>-1e-8 and (L-np.diag(np.diag(L))).max()<1e-8
    # Orient the DC graph mode consistently.
    U[:,0]=one
    bases={'common':U[:,:1],**{f'low_q{q}':U[:,:q+1] for q in range(1,6)},'full':np.eye(10),'high_q3':np.column_stack([one,U[:,-3:]])}
    rng=np.random.default_rng(20261003);bases['random_q3']=la.qr(np.column_stack([one,rng.normal(size=(10,3))]),mode='economic')[0]
    Z=m.impedance(10j*np.pi,.04);H=(Z+Z.conj().T)/2
    _,V=la.eigh(H.real);bases['port_q3']=la.qr(np.column_stack([one,V[:,:3]]),mode='economic')[0]
    (OUT/'bases').mkdir(exist_ok=True)
    for name,Phi in bases.items():pd.DataFrame(Phi).to_csv(OUT/'bases'/f'{name}.csv',index=False)
    rr=roots(m);rrows=[];grad=[]
    for j,s in enumerate(rr):
        ds,res=m.gradients(s,.04);grad.append(ds.real)
        rrows.append(dict(root=j,real=s.real,imag=s.imag,frequency_hz=s.imag/(2*np.pi),residual=res))
    pd.DataFrame(rrows).to_csv(OUT/'TABLE_02_BASELINE_ROOTS.csv',index=False)
    pd.DataFrame(grad).to_csv(OUT/'ROOT_GRADIENTS.csv',index=False)
    assert len(rr)>=8 and max(z.real for z in rr)<-.05,'BASELINE_SPECTRAL_GATE_FAILED'
    # Withheld centered differences in one rho, logKp and logKi at three buses.
    fd=[];s=max((z for z in rr if z.imag>10),key=lambda z:z.real);ds,_=m.gradients(s,.04)
    for j in [0,4,8,10,14,18,20,24,28]:
        h=1e-5;pp=p.copy();pm=p.copy();pp[j]+=h;pm[j]-=h
        sp=Model(pp).refine(s,.04)[0];sm=Model(pm).refine(s,.04)[0];f=(sp-sm)/(2*h)
        fd.append(dict(parameter=j,analytic_real=ds[j].real,fd_real=f.real,absolute_error=abs(ds[j]-f),relative_error=abs(ds[j]-f)/max(abs(ds[j]),1e-9)))
    pd.DataFrame(fd).to_csv(OUT/'TABLE_03_SPECTRAL_DERIVATIVES.csv',index=False)
    assert max(r['relative_error'] for r in fd)<.01
    # Graph-modal interaction expansion: retain every registered frequency/order.
    erows=[];grows=[]
    for tau in [0,.02,.04]:
        for f in [.1,.5,1,2,4,5,8,12]:
            Z=m.impedance(2j*np.pi*f,tau);ZH=U.T@Z@U;diag=np.diag(np.diag(ZH));E=ZH-diag;invD=np.diag(1/np.diag(ZH));K=invD@E;r=la.norm(K,2)
            target=la.inv(ZH);approx=invD.copy();power=np.eye(10,dtype=complex)
            H=(Z+Z.conj().T)/2;Hh=U.T@H@U;comm=U.T@(L@H-H@L)@U
            grows.append(dict(tau_ms=1000*tau,frequency_hz=f,modal_cross_fraction=la.norm(Hh-np.diag(np.diag(Hh)))/la.norm(Hh),
                commutator_identity_error=la.norm(comm-(ev[:,None]-ev[None,:])*Hh)/max(1,la.norm(comm)),neumann_norm=r))
            for k in range(9):
                if k>0:power=-power@K;approx+=power@invD
                if k in [0,1,2,3,5,8]:
                    erows.append(dict(tau_ms=1000*tau,frequency_hz=f,order=k,relative_inverse_error=la.norm(approx-target,2)/la.norm(target,2),
                        norm_convergence_condition=r<1,remainder_bound=(r**(k+1)/(1-r)*la.norm(invD,2)) if r<1 else np.nan))
    pd.DataFrame(erows).to_csv(OUT/'TABLE_04_GRAPH_EXPANSIONS.csv',index=False)
    pd.DataFrame(grows).to_csv(OUT/'TABLE_05_GRAPH_MIXING.csv',index=False)
    print('PREPARE_DONE',parity,'roots',len(rr),'alpha',max(z.real for z in rr),flush=True)

def propose():
    d=tomllib.loads((OUT/'baseline.toml').read_text());p=np.r_[d['rho'],np.log(d['Kp']),np.log(d['Ki'])];P0=pd.read_csv(OUT/'model/ports.csv').P0.to_numpy()
    grad=pd.read_csv(OUT/'grad/baseline/gradients.csv').to_numpy()
    fdcheck=pd.read_csv(OUT/'TABLE_06B_EVENT_DERIVATIVE_VALIDATION.csv')
    assert ((fdcheck.relative_error<=.01)|(fdcheck.absolute_error<=1e-7)).all(),'EVENT_DERIVATIVE_GATE_FAILED'
    adaptive=pd.read_csv(OUT/'eval/baseline/events.csv');tg=pd.read_csv(OUT/'grad/baseline/events.csv')
    assert adaptive['pass'].all(),'BASELINE_EVENT_GATE_FAILED'
    vals=adaptive[['F','R','Vmin','Vmax','slack']].to_numpy().ravel()
    tv=tg[['F','R','Vmin','Vmax','slack']].to_numpy().ravel()
    mismatch=abs(vals-tv).reshape(5,5);limits=np.array([.002,.002,.002,.002,2e-5])
    pd.DataFrame(mismatch,columns=['F','R','Vmin','Vmax','slack']).to_csv(OUT/'TABLE_06_TANGENT_MESH_CHECK.csv',index=False)
    assert np.all(mismatch<=limits),'TANGENT_REFINEMENT_REQUIRED'
    signs=np.tile([1,1,-1,1,-1],5);scale=np.tile([.5,.5,.1,.1,.01],5);bounds=np.tile([.5,.5,.9,1.1,.002],5)
    eventg=signs*(vals-bounds)/scale;eventJ=signs[:,None]*grad/scale[:,None]
    rr=pd.read_csv(OUT/'TABLE_02_BASELINE_ROOTS.csv');rj=pd.read_csv(OUT/'ROOT_GRADIENTS.csv').to_numpy()
    g=np.r_[(rr.real.to_numpy()+.05)/.05,eventg];J=np.vstack([rj/.05,eventJ])
    guard=np.r_[np.full(len(rr),.0001/.05),np.tile([.0005/.5,.0005/.5,.0001/.1,.0001/.1,5e-6/.01],5)]
    rows=[]
    for file in sorted((OUT/'bases').glob('*.csv')):
        family=file.stem;Phi=pd.read_csv(file).to_numpy();k=Phi.shape[1]
        T=la.block_diag(np.eye(10),Phi,Phi);Jq=J@T
        # Trust box imposed in nodal coordinates, identical for every basis.
        radius=np.r_[np.full(10,.015),np.full(20,.08)]
        A=np.vstack([Jq,T,-T]);b=np.r_[-g-guard,radius,radius]
        cost=-(P0/P0.sum())@T[:10]
        start=time.perf_counter();sol=optimize.linprog(cost,A_ub=A,b_ub=b,bounds=[(None,None)]*T.shape[1],method='highs')
        if not sol.success:
            rows.append(dict(family=family,gain_variables=2*k,status=sol.message));continue
        # Lexicographic tie break: smallest nodal log-gain change.
        A2=np.vstack([A,cost]);b2=np.r_[b,sol.fun+1e-7]
        H=T[10:].T@T[10:]+1e-12*np.eye(T.shape[1])
        zv=cp.Variable(T.shape[1])
        secondary=cp.Problem(cp.Minimize(.5*cp.quad_form(zv,cp.psd_wrap(H))),[A2@zv<=b2])
        secondary.solve(solver='CLARABEL',tol_gap_abs=1e-11,tol_gap_rel=1e-11,tol_feas=1e-11,max_iter=300)
        assert secondary.status in ['optimal','optimal_inaccurate'] and zv.value is not None,'SECONDARY_QP_FAILED'
        z=zv.value
        assert np.max(A2@z-b2)<1e-7,'SECONDARY_QP_INFEASIBLE'
        dp=T@z
        assert np.max(A@z-b)<1e-7,'LOCAL_SUBPROBLEM_INFEASIBLE'
        for fraction in [1,.5,.25,.125]:
            write_design(OUT/'proposals'/f'{family}_f{fraction:g}.toml',p+fraction*dp)
        capture=la.norm(rj[:,10:]@la.block_diag(Phi@Phi.T,Phi@Phi.T))**2/max(la.norm(rj[:,10:])**2,1e-30)
        eventcapture=la.norm(eventJ[:,10:]@la.block_diag(Phi@Phi.T,Phi@Phi.T))**2/max(la.norm(eventJ[:,10:])**2,1e-30)
        rows.append(dict(family=family,gain_variables=2*k,predicted_GFL_MW=P0@(p[:10]+dp[:10]),predicted_improvement_MW=P0@dp[:10],
            lp_upper_in_trust_region_MW=P0@p[:10]-sol.fun*P0.sum(),gain_step_norm=la.norm(dp[10:]),
            captured_root_gradient_energy=capture,captured_event_gradient_energy=eventcapture,
            linear_constraint_violation=max(0.,np.max(A@z-b)),secondary_success=secondary.status,
            subproblem_seconds=time.perf_counter()-start,status='LOCAL_PREDICTOR_NOT_NONLINEAR_FEASIBLE'))
    pd.DataFrame(rows).to_csv(OUT/'TABLE_07_LOCAL_PREDICTORS.csv',index=False)
    print(pd.DataFrame(rows).to_string(index=False),flush=True)

if __name__=='__main__':prepare() if sys.argv[1]=='prepare' else propose()
