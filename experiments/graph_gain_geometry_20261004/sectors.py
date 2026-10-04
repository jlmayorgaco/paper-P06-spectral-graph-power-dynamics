from geometry import *

def sectors(G,s,lo,hi,tau=None):
    tau=np.full(10,.04) if tau is None else np.asarray(tau)
    a=s*s*(1+s*TF)*np.exp(s*tau)
    pc=(lo[:10]+hi[:10])/2;dp=(hi[:10]-lo[:10])/2
    ic=(lo[10:]+hi[10:])/2;di=(hi[10:]-lo[10:])/2
    assert np.all(di>0)
    Q=[];E=[];labels=[]
    def prod(x,y): return np.outer(x.conj(),y)
    def imag(x): return (x-x.conj().T)/(2j)
    for i,g in enumerate(G):
        e=np.eye(10)[i]
        Y=np.r_[g,np.zeros(10)]
        U=np.r_[np.zeros(10),e]
        V=np.r_[(a[i]*e-(s*pc[i]+ic[i])*g)/di[i],-s*dp[i]*e/di[i]]
        # Congruence scaling must be shared across rows, performed below.
        raw=[prod(Y,Y)-prod(U,U),prod(Y,Y)-prod(V,V)]
        raw += [herm(prod(Y+eps*U,Y+zeta*V)) for eps in [-1,1] for zeta in [-1,1]]
        Q+=raw;E += [imag(prod(Y,U)),imag(prod(Y,V)),imag(prod(U,V))]
        labels += [f'{30+i}:p_sector',f'{30+i}:i_sector']+[f'{30+i}:rlt{eps},{zeta}' for eps in [-1,1] for zeta in [-1,1]]
    Q=np.array(Q);E=np.array(E)
    # Normalize coordinates using aggregate diagonal energy, no modal pattern.
    diag=sum(abs(np.diagonal(x)) for x in Q)+sum(abs(np.diagonal(x)) for x in E)
    S=1/np.sqrt(np.maximum(diag,1e-20));S/=max(S)
    Q=Q*S[None,:,None]*S[None,None,:];E=E*S[None,:,None]*S[None,None,:]
    Q=np.array([herm(q)/max(la.norm(q,'fro'),1e-30) for q in Q])
    E=np.array([herm(q)/max(la.norm(q,'fro'),1e-30) for q in E])
    return Q,E,S,labels

def solve_sector(Q,E):
    n=Q.shape[1];weights=cp.Variable(len(Q),nonneg=True);eq=cp.Variable(len(E));t=cp.Variable()
    W=sum(weights[i]*Q[i] for i in range(len(Q)))+sum(eq[i]*E[i] for i in range(len(E)))
    problem=cp.Problem(cp.Maximize(t),[cp.sum(weights)+cp.norm1(eq)<=1,-W-t*np.eye(n)>>0])
    start=time.perf_counter()
    problem.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_feas=1e-9,tol_gap_rel=1e-9,max_iter=200)
    if weights.value is None: raise RuntimeError(problem.status)
    w=np.maximum(weights.value,0);e=eq.value
    mat=np.einsum('i,ijk->jk',w,Q)+np.einsum('i,ijk->jk',e,E)
    slack=-la.eigvalsh(mat)[-1]
    return dict(slack=float(slack),status=problem.status,seconds=time.perf_counter()-start,
                weights=w.tolist(),equalities=e.tolist())

def main():
    lock={n:hashlib.sha256((OUT/n).read_bytes()).hexdigest() for n in ['sectors.py','ADDENDUM_01_REAL_SECTORS.md','geometry.py']}
    target=OUT/'SECTOR_SOURCE_LOCK.json'
    if target.exists():assert json.loads(target.read_text())==lock
    else:save(target.name,lock)
    rows=[];witnesses={}
    for rho in [.875,.90,.95,.975,.99]:
        p=P0.copy();p[:10]=rho;m=Model(p)
        for hz in [.25,.5,1,2,3,4,5,6,8,10]:
            s=-.05+2j*np.pi*hz;G=return_matrix(m,s)
            Q,E,S,lab=sectors(G,s,LO,HI);res=solve_sector(Q,E)
            key=f'rho{rho}_hz{hz}';witnesses[key]={'s':[s.real,s.imag],'rho':rho,'scaling':S.tolist(),'labels':lab,**res}
            rows.append(dict(id=key,rho=rho,hz=hz,slack=res['slack'],seconds=res['seconds'],status=res['status']))
            print(key,'slack',res['slack'],flush=True)
        pd.DataFrame(rows).to_csv(OUT/'TABLE_05_SECTOR_POINT_SCREEN.csv',index=False);save('SECTOR_POINT_WITNESSES.json',witnesses)
    checks=[]
    for name in ['anchor','single30','single37','joint','complex_pair','corrected']:
        d=json.loads((OLD/f'CERT_{name}.json').read_text());p=np.array(d['p']);m=Model(p);s=complex(*d['center'])
        G=return_matrix(m,s);Q,E,S,_=sectors(G,s,LO,HI);res=solve_sector(Q,E)
        F=np.diag(s*s*(1+s*TF))-np.diag((s*p[10:20]+p[20:])*np.exp(-s*.04))@G
        _,_,vh=la.svd(F);q=vh[-1].conj();dp=(HI[:10]-LO[:10])/2;pc=(HI[:10]+LO[:10])/2
        x=np.r_[q,(p[10:20]-pc)/dp*(G@q)]/S;x/=la.norm(x)
        vals=np.einsum('i,kij,j->k',x.conj(),Q,x).real;eqs=np.einsum('i,kij,j->k',x.conj(),E,x).real
        checks.append(dict(design=name,slack=res['slack'],min_inequality=float(min(vals)),max_equality=float(max(abs(eqs))),
                           passed=bool(min(vals)>-1e-7 and max(abs(eqs))<1e-7 and res['slack']<1e-7)))
    pd.DataFrame(checks).to_csv(OUT/'TABLE_06_SECTOR_ROOT_CHECK.csv',index=False)
    assert all(r['passed'] for r in checks)

if __name__=='__main__':main()
