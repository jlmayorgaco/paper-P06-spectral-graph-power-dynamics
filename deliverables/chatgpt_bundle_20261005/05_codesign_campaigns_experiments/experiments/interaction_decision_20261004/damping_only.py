from compensated import *
from functools import lru_cache

@lru_cache(None)
def local(i,t,h=.01):
    ki=ANCHOR[20+i];s=TARGET
    for frac in np.linspace(.1,1,10):
        p=ANCHOR.copy();p[i]+=h*frac;p[10+i]*=1+t*frac
        for it in range(20):
            p[20+i]=ki;m=Model(p)
            lu=la.lu_factor(s*np.eye(204)-m.A0,check_finite=False)
            rb=la.lu_solve(lu,m.B,check_finite=False);e=np.exp(-s*TAU)
            K=np.eye(10)-(m.C@rb)*e[None,:]
            left,sv,vh=la.svd(K,check_finite=False);u=left[:,-1];v=vh[-1].conj()
            Ds=(m.C@la.lu_solve(lu,rb,check_finite=False))*e[None,:]+(m.C@rb)*(TAU*e)[None,:]
            dki=-e[i]*np.vdot(u,m.C@la.lu_solve(lu,M['Bi'][:,i],check_finite=False))*v[i]
            ds=1j*np.vdot(u,Ds@v);r=np.vdot(u,K@v)
            update=la.solve(np.array([[ds.real,dki.real],[ds.imag,dki.imag]]),-np.array([r.real,r.imag]))
            s+=1j*update[0];ki+=update[1]
            if abs(update[0])<1e-10 and abs(update[1])<1e-7:break
        else:raise RuntimeError('local corrector failed')
    p[20+i]=ki
    assert KI_BOUNDS[0]<=ki<=KI_BOUNDS[1],'Ki bound'
    assert abs(s.imag-TARGET.imag)<2,'frequency excursion >2'
    root,res=Model(p).refine(s)
    assert abs(root.real-TARGET.real)<1e-8,'real-part compensation failed'
    return p,root,res

def main():
    import hashlib
    save('ADDENDUM_02_LOCK.json',{'sha256':hashlib.sha256((OUT/'ADDENDUM_02_DAMPING_ONLY.md').read_bytes()).hexdigest()})
    rows=[];single=[];start=time.perf_counter()
    for t in [-.10,-.05,-.02,.02,.05,.10]:
        for i in range(10):
            try:
                p,z,res=local(i,t);single.append(dict(bus=30+i,t=t,valid=True,Kp=p[10+i],Ki=p[20+i],real=z.real,imag=z.imag,residual=res))
            except Exception as ex:single.append(dict(bus=30+i,t=t,valid=False,error=str(ex)))
    pd.DataFrame(single).to_csv(OUT/'TABLE_04_DAMPING_SINGLETONS.csv',index=False)
    for mag in [.02,.05,.10]:
        for t,u in [(mag,mag),(-mag,-mag),(mag,-mag)]:
            for i,j in combinations(range(10),2):
                row=dict(bus_i=30+i,bus_j=30+j,t_i=t,t_j=u,valid=False,error='')
                try:
                    pi,zi,_=local(i,t);pj,zj,_=local(j,u);p=ANCHOR+(pi-ANCHOR)+(pj-ANCHOR)
                    z=TARGET;maxjump=0
                    for f in np.linspace(.1,1,10):
                        # Interpolate actual actions; validate endpoints independently.
                        pp=ANCHOR+f*(p-ANCHOR);new,res=Model(pp).refine(z)
                        maxjump=max(maxjump,abs(new-z));assert abs(new-z)<.5,'joint tracking jump'
                        z=new
                    row.update(valid=True,joint_real=z.real,joint_imag=z.imag,nodal_real=(zi+zj-TARGET).real,
                        nodal_imag=(zi+zj-TARGET).imag,collective_real=(z-zi-zj+TARGET).real,
                        added_GFL_MW=float(.01*(PORTS.P0.iloc[i]+PORTS.P0.iloc[j])),residual=res,max_jump=maxjump,
                        false_accept=bool(z.real>-.0499))
                    if row['false_accept']:save(f'damping_candidate_{len(rows)}.json',{'p':p.tolist(),'root':[z.real,z.imag],**row})
                except Exception as ex:row['error']=str(ex)
                rows.append(row)
            pd.DataFrame(rows).to_csv(OUT/'TABLE_05_DAMPING_COMBINATIONS.csv',index=False)
        print('mag',mag,'violations',sum(r.get('false_accept',False) for r in rows),'seconds',time.perf_counter()-start,flush=True)
    save('DAMPING_SUMMARY.json',{'cases':len(rows),'valid':sum(r['valid'] for r in rows),'violations':sum(r.get('false_accept',False) for r in rows),
        'scope':'adaptive discovery of single-mode decision reversal, not nonlinear feasibility or optimality'})
if __name__=='__main__':main()
