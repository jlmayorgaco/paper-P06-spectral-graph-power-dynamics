from codesign_multimode import track,root_gradient,TARGET,KPLO,KPHI
from moments import *
import cvxpy as cp
import time
GOUT=OUT/'globalized';GOUT.mkdir(exist_ok=True)
def saveg(name,obj):(GOUT/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf8')

def merit(z,roots):return .5*float(z@z)+1000*max(0.,max(roots.real)-TARGET)

def solve(m,dt,seeds,method,ident):
    kp=P[10:20];ki=P[20:];tau=.04+dt
    c=m.weights/ki**2;cn=c*kp;bb=m.weights@(dt/ki)
    es=la.norm(cn);cn=cn/es;bb=bb/es
    zcur=np.zeros(10) if method=='modal_only' else cn*bb/(cn@cn)
    p=P.copy();p[10:20]=kp*(1+zcur)
    roots=track(P,np.full(10,.04),p,tau,seeds)
    rows=[];status='MOMENT_ONLY'
    if method=='moment_only':return p,roots,rows,status
    for k in range(60):
        eig=[root_gradient(p,tau,s) for s in roots]
        roots=np.array([e[0] for e in eig]);gz=np.array([e[1].real*kp for e in eig])
        z=cp.Variable(10);slack=cp.Variable(nonneg=True)
        constraints=[gz@z<=gz@zcur+TARGET-roots.real+slack,
            z>=KPLO/kp-1,z<=KPHI/kp-1,z>=zcur-.05,z<=zcur+.05]
        if method=='collective':constraints.append(cn@z==bb)
        prob=cp.Problem(cp.Minimize(.5*cp.sum_squares(z)+1000*slack),constraints)
        prob.solve(solver='CLARABEL',tol_gap_abs=1e-12,tol_gap_rel=1e-12,tol_feas=1e-12,max_iter=200)
        if z.value is None:status='QP_'+prob.status;break
        step=z.value-zcur;pred=merit(zcur,roots)-float(prob.value)
        if la.norm(step)<1e-6 and max(roots.real)<=TARGET+1e-7:
            status='CONVERGED_NUMERICALLY';break
        accepted=False
        for power in range(11):
            beta=2.**(-power);zn=zcur+beta*step;pn=P.copy();pn[10:20]=kp*(1+zn)
            rn=track(p,tau,pn,tau,roots)
            if merit(zn,rn)<=merit(zcur,roots)-1e-4*beta*max(pred,0)+1e-8:
                accepted=True;break
        rows.append(dict(iteration=k,alpha=float(max(roots.real)),objective=float(.5*zcur@zcur),
            merit=merit(zcur,roots),full_step_norm=float(la.norm(step)),beta=beta,
            accepted=accepted,QP_status=prob.status,restoration_slack=float(slack.value)))
        saveg('ITER_'+ident+'.json',rows)
        if not accepted:status='LINE_SEARCH_STALLED';break
        p=pn;zcur=zn;roots=rn
    else:status='MAX_ITERATIONS'
    return p,roots,rows,status

def main():
    info={'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'addendum_sha256':hashlib.sha256((OUT/'ADDENDUM_03_GLOBALIZATION_ABLATION.md').read_bytes()).hexdigest()}
    lock=GOUT/'SOURCE_LOCK.json'
    if lock.exists():assert json.loads(lock.read_text())==info
    else:saveg(lock.name,info)
    m=Moments();cat=pd.read_csv(OLD/'TABLE_08_MULTIBRANCH_CONTINUATION.csv')
    cat=cat[(cat.design=='corrected')&(cat.step==20)];seeds=cat.real.to_numpy()+1j*cat.imag.to_numpy()
    rows=[];designs=[];catalog=[];failures=[]
    for name,dt in patterns().items():
        methods=['collective']+(['modal_only','moment_only'] if name=='uniform_1ms' else [])
        for method in methods:
            ident=name+'_'+method;start=time.perf_counter()
            try:
                p,roots,iterations,status=solve(m,dt,seeds,method,ident)
                kp=p[10:20];tau=.04+dt;critical=roots[np.argmax(roots.real)]
                kap=float(m.weights@(dt/P[20:]-(kp-P[10:20])/P[20:]**2))
                row=dict(pattern=name,method=method,status=status,critical_real=critical.real,critical_imag=critical.imag,
                    kappa=kap,relative_Kp_norm=float(la.norm(kp/P[10:20]-1)),max_relative_Kp=float(max(abs(kp/P[10:20]-1))),
                    catalog_margin_pass=bool(critical.real<=-.05),gain_bounds_pass=bool(min(kp)>=KPLO and max(kp)<=KPHI),
                    iterations=len(iterations),seconds=time.perf_counter()-start)
                rows.append(row)
                for j,z in enumerate(roots):catalog.append(dict(pattern=name,method=method,mode=j,real=z.real,imag=z.imag))
                for i in range(10):designs.append(dict(pattern=name,method=method,bus=30+i,rho=p[i],Kp=kp[i],Ki=p[20+i],tau=tau[i]))
                dest=GOUT/'designs';dest.mkdir(exist_ok=True)
                (dest/(ident+'.toml')).write_text('\n'.join(f'{k} = [{", ".join(format(float(x),".17g") for x in v)}]'
                    for k,v in [('rho',p[:10]),('Kp',kp),('Ki',p[20:]),('tau_vector',tau)])+'\n',encoding='utf8')
                print(ident,status,critical,'effort',row['relative_Kp_norm'],'kappa',kap,flush=True)
            except Exception as exc:
                failures.append(dict(pattern=name,method=method,error=repr(exc)))
                print('FAILURE',ident,repr(exc),flush=True)
            pd.DataFrame(rows).to_csv(GOUT/'TABLE_13_GLOBALIZED.csv',index=False)
            pd.DataFrame(designs).to_csv(GOUT/'TABLE_14_DESIGNS.csv',index=False)
            pd.DataFrame(catalog).to_csv(GOUT/'TABLE_15_CATALOG.csv',index=False)
            saveg('FAILURES.json',failures)

if __name__=='__main__':main()
