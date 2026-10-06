"""Exact moment equality with local nonlinear-eigenvalue guard, at fixed rho/Ki."""
from moments import *
import cvxpy as cp
import time

TARGET=-.0600001
KPLO=.25*2*np.pi*5;KPHI=4*2*np.pi*5

def track(pa,ta,pb,tb,seeds):
    roots=seeds.copy();f=0.;step=.1
    while f<1-1e-12:
        df=min(step,1-f);p=pa+f*(pb-pa);t=ta+f*(tb-ta)
        eig=[root_gradient(p,t,z) for z in roots]
        predicted=np.array([e[0]+df*(e[1]@(pb[10:20]-pa[10:20])+e[2]@(tb-ta)) for e in eig])
        mm=Model(pa+(f+df)*(pb-pa),ta+(f+df)*(tb-ta))
        new=np.array([mm.refine(z)[0] for z in predicted])
        sep=min(abs(new[i]-new[j]) for i in range(len(new)) for j in range(i))
        if sep<1e-5 or max(abs(new-predicted))>.03:
            step=df/2
            if step<1e-6:raise RuntimeError('unresolved catalog continuation')
            continue
        roots=new;f+=df;step=min(.2,step*1.5)
    return roots

def root_gradient(p,tau,seed):
    m=Model(p,tau);s,res=m.refine(seed)
    _,(scale,_)=la.matrix_balance(m.A,permute=False,separate=True)
    D=m.delta(s);Db=D*scale[None,:]/scale[:,None]
    U,sv,Vh=la.svd(Db);v=scale*Vh[-1].conj();w=U[:,-1]/scale
    den=np.vdot(w,m.derivative(s)@v)
    gp=np.exp(-s*tau)*(w.conj()@M['Bp'])*(m.C@v)/den
    gt=-s*np.exp(-s*tau)*(w.conj()@m.B)*(m.C@v)/den
    return s,gp,gt,res

def collective(m,dt,seed):
    kp=P[10:20];ki=P[20:];tau=.04+dt
    c=m.weights/ki**2;cn=c*kp;bb=m.weights@(dt/ki)
    eqscale=la.norm(cn);cn=cn/eqscale;bb=bb/eqscale
    zcur=cn*bb/(cn@cn);rows=[]
    pinit=P.copy();pinit[10:20]=kp*(1+zcur)
    roots=track(P,np.full(10,.04),pinit,tau,seed)
    for k in range(60):
        p=P.copy();p[10:20]=kp*(1+zcur)
        eig=[root_gradient(p,tau,s) for s in roots]
        roots=np.array([e[0] for e in eig]);gz=np.array([e[1].real*kp for e in eig])
        s=roots[np.argmax(roots.real)];res=max(e[3] for e in eig)
        z=cp.Variable(10);slack=cp.Variable(nonneg=True)
        cons=[cn@z==bb,gz@z<=gz@zcur+TARGET-roots.real+slack,
              z>=KPLO/kp-1,z<=KPHI/kp-1,z>=zcur-.05,z<=zcur+.05]
        prob=cp.Problem(cp.Minimize(.5*cp.sum_squares(z)+1000*slack),cons)
        prob.solve(solver='CLARABEL',tol_gap_abs=1e-12,tol_gap_rel=1e-12,tol_feas=1e-12,max_iter=200)
        if z.value is None:raise RuntimeError('QP '+prob.status)
        nxt=z.value;step=float(la.norm(nxt-zcur))
        rows.append(dict(iteration=k,alpha=s.real,imag=s.imag,relative_gain_norm=float(la.norm(zcur)),
                         update_norm=step,root_residual=res,QP_status=prob.status,restoration_slack=float(slack.value)))
        dump('ITERATION_'+CURRENT_PATTERN+'.json',rows)
        if step<1e-8 and s.real<=TARGET+1e-8:break
        pn=P.copy();pn[10:20]=kp*(1+nxt)
        roots=track(p,tau,pn,tau,roots)
        zcur=nxt
    else:
        rows[-1]['QP_status']+=';MAX_ITERATIONS_NO_CONVERGENCE_CLAIM'
        dump('ITERATION_'+CURRENT_PATTERN+'.json',rows)
    p=P.copy();p[10:20]=kp*(1+zcur)
    return p,s,rows

def write_design(name,p,tau):
    dest=OUT/'designs';dest.mkdir(exist_ok=True)
    (dest/f'{name}.toml').write_text('\n'.join(f'{k} = [{", ".join(format(float(x),".17g") for x in v)}]'
        for k,v in [('rho',p[:10]),('Kp',p[10:20]),('Ki',p[20:]),('tau_vector',tau)])+'\n',encoding='utf8')

def main():
    global CURRENT_PATTERN
    lock=OUT/'MULTIMODE_SOURCE_LOCK_V2.json';info={'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'pattern_sha256':hashlib.sha256((OUT/'ADDENDUM_02_MULTIMODE.md').read_bytes()).hexdigest(),
        'target_guard':TARGET,'gain_bounds':[KPLO,KPHI],'trust_region_relative':.05,'max_iterations':60}
    if lock.exists():assert json.loads(lock.read_text())==info
    else:dump(lock.name,info)
    m=Moments();kp=P[10:20];ki=P[20:];tau=np.full(10,.04)
    cat=pd.read_csv(OLD/'TABLE_08_MULTIBRANCH_CONTINUATION.csv');cat=cat[(cat.design=='corrected')&(cat.step==20)]
    seeds=cat.real.to_numpy()+1j*cat.imag.to_numpy();seed=seeds[np.argmax(seeds.real)]
    s,g,gt,res=root_gradient(P,tau,seed);deriv=[]
    for i in [0,3,7]:
        h=1e-3;pp=P.copy();pm=P.copy();pp[10+i]+=h;pm[10+i]-=h
        fd=(Model(pp,tau).refine(s)[0]-Model(pm,tau).refine(s)[0])/(2*h)
        ht=1e-7;tp=tau.copy();tm=tau.copy();tp[i]+=ht;tm[i]-=ht
        fdt=(Model(P,tp).refine(s)[0]-Model(P,tm).refine(s)[0])/(2*ht)
        for var,an,numerical in [('Kp',g[i],fd),('tau',gt[i],fdt)]:
            deriv.append(dict(bus=i+30,variable=var,analytic_real=an.real,analytic_imag=an.imag,
                FD_real=numerical.real,FD_imag=numerical.imag,relative_error=abs(an-numerical)/abs(numerical)))
    pd.DataFrame(deriv).to_csv(OUT/'TABLE_05_POLE_DERIVATIVES.csv',index=False)
    rows=[];iterations=[];rr=[];vectors=[]
    write_design('base',P,tau)
    for name,dt in patterns().items():
        CURRENT_PATTERN=name
        t0=time.perf_counter();pc,sc,it=collective(m,dt,seeds)
        iterations += [dict(pattern=name,**x) for x in it]
        ps=P.copy();ps[10:20]+=ki*dt
        for method,p in [('unchanged',P.copy()),('sitewise',ps),('collective',pc)]:
            roots=track(P,tau,p,tau+dt,seeds)
            mdl=Model(p,tau+dt);residuals=[mdl.refine(z)[1] for z in roots]
            critical=roots[np.argmax(roots.real)]
            kappa=float(m.weights@(dt/ki-(p[10:20]-kp)/ki**2))
            row=dict(pattern=name,method=method,critical_real=critical.real,critical_imag=critical.imag,
                kappa=kappa,relative_Kp_norm=float(la.norm((p[10:20]-kp)/kp)),max_relative_Kp=float(max(abs(p[10:20]/kp-1))),
                gain_bounds_pass=bool(min(p[10:20])>=KPLO and max(p[10:20])<=KPHI),
                catalog_margin_pass=bool(critical.real<=-.05),root_residual=max(residuals),seconds=time.perf_counter()-t0)
            rows.append(row);print(name,method,critical,'kappa',kappa,'effort',row['relative_Kp_norm'],flush=True)
            for j,z in enumerate(roots):rr.append(dict(pattern=name,method=method,mode=j,real=z.real,imag=z.imag,residual=residuals[j]))
            ident=name+'_'+method;write_design(ident,p,tau+dt)
            for i in range(10):vectors.append(dict(pattern=name,method=method,bus=30+i,rho=p[i],Kp=p[10+i],Ki=p[20+i],tau=tau[i]+dt[i]))
        pd.DataFrame(rows).to_csv(OUT/'TABLE_06_CODESIGN.csv',index=False)
        pd.DataFrame(iterations).to_csv(OUT/'TABLE_07_SQP_ITERATIONS.csv',index=False)
        pd.DataFrame(rr).to_csv(OUT/'TABLE_08_CATALOG_ROOTS.csv',index=False)
        pd.DataFrame(vectors).to_csv(OUT/'TABLE_09_DESIGN_VECTORS.csv',index=False)
    dump('CODESIGN_STATUS.json',{'status':'NUMERICALLY_VALIDATED_LOCAL_CATALOG_ONLY',
        'full_spectrum':'PENDING','nonlinear':'PENDING','fixed_GFL_MW':float(PORTS.P0@P[:10]),
        'fixed_GFL_percent':float(100*PORTS.P0@P[:10]/PORTS.P0.sum()),
        'max_derivative_relative_error':max(r['relative_error'] for r in deriv),
        'all_collective_catalog_margin_pass':all(r['catalog_margin_pass'] for r in rows if r['method']=='collective')})

if __name__=='__main__':main()
