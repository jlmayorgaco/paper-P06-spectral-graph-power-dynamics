"""Independent numerical checks on the explicitly realized 9-bus model."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
from itertools import combinations
import json,time
import numpy as np,pandas as pd
from scipy.linalg import svd,solve
from model import *
OUT=Path(__file__).resolve().parents[1]

def finite_schur(grid,s,kp0,ki0,tau0,kp1,ki1,tau1,site):
    F=grid.char(s,kp0,ki0,tau0)
    p=grid.input_columns(kp0,ki0)[:,site]*np.exp(-s*tau0[site])-grid.input_columns(kp1,ki1)[:,site]*np.exp(-s*tau1[site])
    c=grid.C[site];r=np.arange(3,6);h=np.setdiff1d(np.arange(24),r)
    D=F[np.ix_(h,h)];Hr=F[np.ix_(h,r)];Rh=F[np.ix_(r,h)]
    dhp=solve(D,p[h]);den=1+c[h]@dhp
    a=BASE_MVA*grid.M*(p[r]-Rh@dhp)/den
    vrow=c[r]-c[h]@solve(D,Hr)
    return np.outer(a,vrow),a,vrow.conj(),np.linalg.cond(D)


def contour_count(grid,kp,ki,tau,sigma=.1,density=1):
    """Numerical argument principle with exact exponentials and an analytic tail bound.
    This is not an interval quadrature certificate. Remove one known gauge root.
    """
    B=grid.input_columns(kp,ki)
    Gabs=np.abs(grid.A0)+np.abs(B)@np.abs(grid.C)
    ez,ev=np.linalg.eig(Gabs)
    d=np.maximum(abs(ev[:,np.argmax(ez.real)]),1e-12);d/=max(d)
    A=grid.A0*d[None,:]/d[:,None];bb=B/d[:,None];cc=grid.C*d[None,:]
    Ais=[np.outer(bb[:,i],cc[i]) for i in range(3)]
    tail=np.linalg.norm(A,np.inf)+sum(np.linalg.norm(Ai,np.inf)*np.exp(sigma*tau[i]) for i,Ai in enumerate(Ais))
    R=1.02*tail+1.
    roots,_=root_catalog(grid,kp,ki,tau,order=8,real_floor=-15,max_imag=500)
    # Dense near finite physical roots; logarithmic tails up to the rigorous-radius formula.
    ymax=min(R,120.)
    yy=np.r_[np.linspace(-ymax,ymax,1601*density),np.geomspace(ymax,R,130*density),-np.geomspace(ymax,R,130*density)]
    for lam in roots:
        dist=max(abs(lam.real+sigma),1e-5)
        for sign in [-1,1]:
            yy=np.r_[yy,sign*lam.imag+dist*np.array([-32,-16,-8,-4,-2,-1,-.5,-.25,0,.25,.5,1,2,4,8,16,32])]
    yy=np.unique(np.clip(yy,-R,R))
    xx=np.linspace(-sigma,R,150*density)
    zz=np.r_[xx[:-1]-1j*R,R+1j*np.linspace(-R,R,250*density)[:-1],xx[::-1][:-1]+1j*R,-sigma+1j*yy[::-1]]
    if abs(zz[-1]-zz[0])>1e-10:zz=np.r_[zz,zz[0]]
    def phase(z):
        Fs=z[:,None,None]*np.eye(24)-A[None,:,:]
        for i,Ai in enumerate(Ais):Fs-=np.exp(-z*tau[i])[:,None,None]*Ai[None,:,:]
        sign,log=np.linalg.slogdet(Fs)
        # det(Delta)/s removes the known rotational zero from the winding.
        return np.angle(sign/(z/abs(z)))
    for _ in range(10):
        ph=phase(zz)
        changes=np.angle(np.exp(1j*np.diff(ph)))
        mids=(zz[:-1]+zz[1:])/2
        pm=phase(mids)
        d1=np.angle(np.exp(1j*(pm-ph[:-1])));d2=np.angle(np.exp(1j*(ph[1:]-pm)))
        bad=(abs(d1)>np.pi/8)|(abs(d2)>np.pi/8)|(abs(d1+d2-changes)>.1)
        if not np.any(bad):break
        znew=[]
        for j in range(len(zz)-1):
            znew.append(zz[j])
            if bad[j]:znew.append(mids[j])
        znew.append(zz[-1]);zz=np.array(znew)
    winding=float(np.sum(np.angle(np.exp(1j*np.diff(ph))))/(2*np.pi))
    return dict(count=int(round(winding)),winding=winding,integer_error=abs(winding-round(winding)),
                max_phase_step=float(max(abs(changes))),nodes=len(zz),tail_R=float(R),sigma=sigma)


def run():
    rng=np.random.default_rng(20261005)
    grid=Grid9(rho=.75);cfg=json.load(open(OUT/'data/config.json'));kp=np.array(cfg['kp0']);ki=np.array(cfg['ki0']);tau=np.array(cfg['tau'])
    ds=json.load(open(OUT/'results/designs.json'))
    tests=[]
    V=grid.volts(grid.x0)
    E=grid.E*np.exp(1j*grid.x0[:3]);I=(grid.x0[18:21]+1j*grid.x0[21:])*np.exp(1j*grid.x0[9:12])
    eqerr=np.max(abs(grid.rhs(grid.x0,grid.detector(grid.x0),kp,ki)))
    tests.append(dict(test='equilibrium_rhs',value=eqerr,tolerance=1e-8,passed=eqerr<1e-8))
    verr=max(abs(V-grid.V0));tests.append(dict(test='voltage_matches_AC_powerflow',value=verr,tolerance=1e-10,passed=verr<1e-10))
    KCL=grid.Ytot@V-np.r_[grid.ys*E+I,np.zeros(6)]
    tests.append(dict(test='KCL_complex',value=max(abs(KCL)),tolerance=1e-11,passed=max(abs(KCL))<1e-11))
    for step in [1e-4,1e-5,1e-6,1e-7]:
        J=np.column_stack([(grid.rhs(grid.x0+step*np.eye(24)[j],np.zeros(3),kp,ki)-grid.rhs(grid.x0-step*np.eye(24)[j],np.zeros(3),kp,ki))/(2*step) for j in range(24)])
        C=np.column_stack([(grid.detector(grid.x0+step*np.eye(24)[j])-grid.detector(grid.x0-step*np.eye(24)[j]))/(2*step) for j in range(24)])
        er=np.linalg.norm(J-grid.A0)/np.linalg.norm(grid.A0);ec=np.linalg.norm(C-grid.C)/np.linalg.norm(grid.C)
        tests+=[dict(test=f'Jacobian_rhs_FD_h{step:g}',value=er,tolerance=2e-7,passed=er<2e-7),dict(test=f'detector_FD_h{step:g}',value=ec,tolerance=2e-7,passed=ec<2e-7)]
    gauge=max(np.linalg.norm(grid.A0@grid.gauge),np.linalg.norm(grid.C@grid.gauge))
    tests.append(dict(test='rotational_gauge',value=gauge,tolerance=1e-10,passed=gauge<1e-10))
    pd.DataFrame(tests).to_csv(OUT/'results/model_checks.csv',index=False)
    # Finite physical laws over real model configurations and frequencies.
    rr=[]
    for rho in [.5,.75,.85]:
        g=Grid9(rho=rho)
        for tms in [1.,10.,21.]:
            t0=np.full(3,tms/1000)
            for f in [.5,1.,2.,5.,8.,12.]:
                s=1j*2*np.pi*f
                for site in range(3):
                    for kind in ['delay','kp','ki','both']:
                        p=kp.copy();ii=ki.copy();tt=t0.copy()
                        if kind=='delay':tt[site]+=.001
                        if kind in ['kp','both']:p[site]*=1.05
                        if kind in ['ki','both']:ii[site]*=.95
                        exact=g.impedance(s,p,ii,tt)-g.impedance(s,kp,ki,t0)
                        pred,a,v,cond=finite_schur(g,s,kp,ki,t0,p,ii,tt,site)
                        er=np.linalg.norm(exact-pred)/max(np.linalg.norm(exact),1e-20)
                        D=(exact+exact.conj().T)/2;ev=np.linalg.eigvalsh(D);dsing=svd(exact,compute_uv=False)
                        sine=max(0.,1-abs(np.vdot(a,v))**2/(np.linalg.norm(a)**2*np.linalg.norm(v)**2))
                        product=-.25*(np.linalg.norm(a)**2*np.linalg.norm(v)**2-abs(np.vdot(v,a))**2)
                        rr.append(dict(rho=rho,tau_ms=tms,f_hz=f,site=site+1,kind=kind,relative_error=er,
                          eig_min=ev[0],eig_middle=ev[1],eig_max=ev[2],rank_one_ratio=dsing[1]/dsing[0],
                          sin2_paths=sine,hidden_cond=cond,eig_product_error=abs(ev[0]*ev[2]-product)/max(abs(product),1e-30)))
    ranks=pd.DataFrame(rr);ranks.to_csv(OUT/'results/finite_rank_laws.csv',index=False)
    # Multi-site determinant identity, random admissible gains, support sizes 1..3.
    detrows=[]
    for m in [1,2,3]:
        for trial in range(40):
            support=rng.choice(3,m,replace=False);p=kp.copy();ii=ki.copy();p[support]*=rng.uniform(.7,1.3,m);ii[support]*=rng.uniform(.7,1.3,m)
            s=complex(-.2,2*np.pi*rng.uniform(.3,20));F=grid.char(s,kp,ki,tau);Fnew=grid.char(s,p,ii,tau)
            B0=grid.input_columns(kp,ki);B1=grid.input_columns(p,ii)
            P=(B0[:,support]-B1[:,support])*np.exp(-s*tau[support])[None,:]
            little=np.eye(m)+grid.C[support]@solve(F,P)
            ratio=np.linalg.det(Fnew)/np.linalg.det(F);pred=np.linalg.det(little)
            detrows.append(dict(support_size=m,trial=trial,error=abs(ratio-pred)/max(abs(ratio),1e-15)))
    pd.DataFrame(detrows).to_csv(OUT/'results/determinant_checks.csv',index=False)
    # Reconstruct / place one pole using the exact whole-network return.
    roots,_=root_catalog(grid,kp,ki,tau,8)
    osci=roots[(roots.imag>1)&(roots.real>-6)][:5]
    ar=[];sr=[];tr=[]
    for ri,lam in enumerate(osci):
        gp,gi,gt,_,_=grid.gain_sensitivity(lam,kp,ki,tau)
        for site in range(3):
            pr,ir=grid.assign_one(lam,site,kp,ki,tau)
            ar.append(dict(kind='reconstruction',root=ri,site=site+1,shift=0.,relative_gain_error=max(abs(pr[site]/kp[site]-1),abs(ir[site]/ki[site]-1)),root_residual=svd(grid.char(lam,pr,ir,tau),compute_uv=False)[-1],admissible=True))
            for shift in [-.02,-.2]:
                tar=lam+shift;pa,ia=grid.assign_one(tar,site,kp,ki,tau)
                ratios=np.r_[pa/kp,ia/ki]
                ar.append(dict(kind='new_target',root=ri,site=site+1,shift=shift,relative_gain_error=np.nan,
                   root_residual=svd(grid.char(tar,pa,ia,tau),compute_uv=False)[-1],admissible=bool(np.all((ratios>=.25)&(ratios<=2)))))
            for field in ['kp','ki','tau']:
                p=kp.copy();q=ki.copy();tt=tau.copy()
                base={'kp':kp[site],'ki':ki[site],'tau':tau[site]}[field];eps=1e-4*base
                plus={'kp':p,'ki':q,'tau':tt};plus[field][site]+=eps
                zp,_,op=refine_root(grid,lam,p,q,tt)
                p=kp.copy();q=ki.copy();tt=tau.copy();minus={'kp':p,'ki':q,'tau':tt};minus[field][site]-=eps
                zm,_,om=refine_root(grid,lam,p,q,tt)
                fd=(zp-zm)/(2*eps);an={'kp':gp[site],'ki':gi[site],'tau':gt[site]}[field]
                sr.append(dict(root=ri,site=site+1,parameter=field,relative_error=abs(fd-an)/max(abs(an),1e-10),absolute_error=abs(fd-an),ok=bool(op and om),ratio_identity_error=abs(gp[site]-lam*gi[site])))
            for dh in [.0001,.001,.005]:
                alpha,om=lam.real,lam.imag
                pp=kp.copy();ii=ki.copy();tt=tau.copy();tt[site]+=dh
                pp[site]=np.exp(alpha*dh)*(kp[site]*np.cos(om*dh)+(ki[site]+alpha*kp[site])/om*np.sin(om*dh))
                ii[site]=np.exp(alpha*dh)*(ki[site]*np.cos(om*dh)-(abs(lam)**2*kp[site]+alpha*ki[site])/om*np.sin(om*dh))
                error=svd(grid.char(lam,pp,ii,tt),compute_uv=False)[-1]
                tr.append(dict(root=ri,site=site+1,delay_increment=dh,root_residual=error,kp=pp[site],ki=ii[site]))
    pd.DataFrame(ar).to_csv(OUT/'results/analytic_gain_checks.csv',index=False)
    pd.DataFrame(sr).to_csv(OUT/'results/sensitivity_checks.csv',index=False)
    pd.DataFrame(tr).to_csv(OUT/'results/delay_transport_checks.csv',index=False)
    # Full-region counts with two independently refined contours.
    cr=[];catalog=[]
    for name,d in ds.items():
        p=np.array(d['kp']);ii=np.array(d['ki']);tt=np.array(d.get('tau',tau))
        for order in [4,6,8,10]:
            zz,ee=root_catalog(grid,p,ii,tt,order,real_floor=-30,max_imag=500)
            for l,e in zip(zz,ee):catalog.append(dict(design=name,pade_seed_order=order,real=l.real,imag=l.imag,f_hz=l.imag/2/np.pi,residual=e))
        for sigma in [.0,.1]:
            # A tiny positive shift avoids zero on the undeﬂated contour when sigma=0.
            use_sig=sigma if sigma>0 else 1e-8
            c1=contour_count(grid,p,ii,tt,use_sig,1);c2=contour_count(grid,p,ii,tt,use_sig,2)
            cr.append(dict(design=name,guard=sigma,**c2,coarse_count=c1['count'],refinement_agrees=c1['count']==c2['count']))
        print(name,'counts',cr[-2]['count'],cr[-1]['count'],flush=True)
    pd.DataFrame(cr).to_csv(OUT/'results/exact_DDE_contour_counts.csv',index=False)
    pd.DataFrame(catalog).to_csv(OUT/'results/root_catalog.csv',index=False)
    # Save precisely the exported model for high-precision/interval validation.
    np.savez(OUT/'data/linearization.npz',A0=grid.A0,C=grid.C,M=grid.M,tf=grid.tf,x0=grid.x0,V0=grid.V0,Y=grid.Y,ys=grid.ys,Jv=grid.Jv)
    summary=dict(equilibrium_rhs=eqerr,model_tests=len(tests),model_passed=int(sum(t['passed'] for t in tests)),
                 rank_law_cases=len(ranks),rank_relative_max=float(ranks.relative_error.max()),
                 rank_signed_cases=int(((ranks.eig_min<0)&(ranks.eig_max>0)).sum()),
                 rank_middle_relative_max=float((abs(ranks.eig_middle)/ranks[['eig_min','eig_max']].abs().max(axis=1)).max()),
                 determinant_cases=len(detrows),determinant_relative_max=float(max(x['error'] for x in detrows)),
                 gain_cases=len(ar),gain_max_residual=float(max(x['root_residual'] for x in ar)),
                 sensitivity_cases=len(sr),sensitivity_relative_max=float(max(x['relative_error'] for x in sr)),
                 transport_cases=len(tr),transport_max_residual=float(max(x['root_residual'] for x in tr)),
                 contour_designs=len(ds),contour_refinement_all=bool(all(x['refinement_agrees'] for x in cr)))
    (OUT/'results/validation_summary.json').write_text(json.dumps(summary,indent=2));print(summary,flush=True)
if __name__=='__main__':run()
