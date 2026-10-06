"""Frozen, inexpensive spectral discovery; no nonlinear simulations."""
from model import *
import time

def main():
    freeze();rows=[];fail=[];start=time.perf_counter()
    for family,count in [('R',18),('P',16),('I',16)]:
        joint=ROOTS.copy();singles=np.tile(ROOTS,(10,1))
        alive=np.ones(11,dtype=bool);singlealive=np.ones((10,11),dtype=bool)
        for k in range(count):
            p=P0.copy()
            if family=='R':p[:10]+=.005*k
            elif family=='P':p[10:20]*=1+.02*k
            else:p[20:]*=1+.02*k
            designs=[p]+[singleton(p,i) for i in range(10)]
            for site,pp in enumerate(designs):
                m=Model(pp);prev=joint if site==0 else singles[site-1]
                valid=alive if site==0 else singlealive[site-1]
                for j in range(11):
                    if not valid[j]:continue
                    try:
                        r,res=m.refine(prev[j]);jump=abs(r-prev[j])
                        if jump>2:raise RuntimeError(f'large branch jump {jump}')
                        if r.imag<0:raise RuntimeError('imaginary branch crossed')
                        prev[j]=r
                    except (RuntimeError,ValueError,la.LinAlgError) as ex:
                        valid[j]=False;fail.append(dict(family=family,k=k,site=site,mode=j,error=str(ex)))
                for j in range(11):
                    for l in range(j):
                        if valid[j] and valid[l] and abs(prev[j]-prev[l])<1e-5:
                            valid[j]=False;valid[l]=False
                            fail.append(dict(family=family,k=k,site=site,mode=j,error=f'collision with {l}'))
            for j in range(11):
                valid=alive[j] and singlealive[:,j].all()
                nodal=ROOTS[j]+sum(singles[:,j]-ROOTS[j])
                rows.append(dict(family=family,k=k,mode=j,rho=p[0],Kp_factor=p[10]/P0[10],Ki_factor=p[20]/P0[20],
                    joint_real=joint[j].real,joint_imag=joint[j].imag,nodal_real=nodal.real,nodal_imag=nodal.imag,
                    collective_real=(joint[j]-nodal).real,valid=bool(valid),
                    false_accept=bool(valid and nodal.real<=-.0501 and joint[j].real>=-.0499),
                    false_reject=bool(valid and nodal.real>=-.0499 and joint[j].real<=-.0501)))
            pd.DataFrame(rows).to_csv(OUT/'TABLE_01_DISCOVERY.csv',index=False)
            save('DISCOVERY_FAILURES.json',fail)
            np.savez(OUT/f'state_{family}_{k:02}.npz',p=p,joint=joint,singles=singles,alive=alive,singlealive=singlealive)
            print(f'{family} {k}/{count-1} valid={sum(alive)} alpha={max(joint.real):.5g} seconds={time.perf_counter()-start:.1f}',flush=True)
    df=pd.DataFrame(rows);save('DISCOVERY_SUMMARY.json',{'rows':len(df),'false_accepts':int(df.false_accept.sum()),
        'false_rejects':int(df.false_reject.sum()),'invalid':int((~df.valid).sum()),'seconds':time.perf_counter()-start,
        'status':'EXPLORATORY_NUMERICAL_ONLY','scope':'11 tracked physical modes, not exhaustive DDE spectrum'})
if __name__=='__main__':main()
