from model import *
from itertools import combinations
import time

ANCHOR=P0.copy();ANCHOR[20:]*=1.14
TARGET=complex(np.load(OUT/'state_I_07.npz')['joint'][2])
KP_BOUNDS=(7.853981633974483,125.66370614359172)
KI_BOUNDS=(61.68502750680849,986.9604401089358)

def compensate(i,h):
    p=ANCHOR.copy();p[i]+=h;m=Model(p)
    s=TARGET;e=np.exp(-s*TAU);B=m.B.copy();B[:,i]=0
    D=s*np.eye(204)-m.A0-(B*e)@m.C
    R=la.solve(D,np.column_stack([M['Bp'][:,i],M['Bi'][:,i]]),check_finite=False)
    a=e[i]*(m.C[i]@R);gain=la.solve(np.vstack([a.real,a.imag]),[1.,0.])
    p[10+i],p[20+i]=gain
    assert KP_BOUNDS[0]<=gain[0]<=KP_BOUNDS[1] and KI_BOUNDS[0]<=gain[1]<=KI_BOUNDS[1], 'gain bound'
    root,res=Model(p).refine(TARGET)
    assert abs(root-TARGET)<1e-7,'compensation failed'
    return p,root,res,float(np.linalg.cond(np.vstack([a.real,a.imag])))

def combined(sites,h):
    p=ANCHOR.copy()
    for i in sites:
        pp,*_=compensate(i,h);ix=[i,10+i,20+i];p[ix]=pp[ix]
    return p

def main():
    import hashlib
    save('ADDENDUM_01_LOCK.json',{'sha256':hashlib.sha256((OUT/'ADDENDUM_01_COMPENSATED_REPLACEMENT.md').read_bytes()).hexdigest(),
        'anchor':ANCHOR.tolist(),'target':[TARGET.real,TARGET.imag]})
    rows=[];sr=[];start=time.perf_counter()
    for h in [.0025,.005,.01,.02,.04]:
        for i in range(10):
            try:
                p,z,res,cond=compensate(i,h)
                sr.append(dict(h=h,bus=30+i,Kp=p[10+i],Ki=p[20+i],root_real=z.real,root_imag=z.imag,residual=res,condition=cond,passed=True,error=''))
            except Exception as ex:sr.append(dict(h=h,bus=30+i,passed=False,error=str(ex)))
        for sites in list(combinations(range(10),2))+[tuple(range(10))]:
            row=dict(h=h,sites=','.join(str(i+30) for i in sites),count=len(sites),valid=False,error='')
            try:
                z=TARGET;maxjump=0
                for t in np.linspace(0,h,9)[1:]:
                    p=combined(sites,t);new,res=Model(p).refine(z)
                    maxjump=max(maxjump,abs(z-new));assert abs(z-new)<.5,'tracking jump'
                    z=new
                row.update(valid=True,joint_real=z.real,joint_imag=z.imag,nodal_real=TARGET.real,
                    interaction_real=z.real-TARGET.real,residual=res,max_jump=maxjump,
                    added_GFL_MW=float(sum(PORTS.P0.iloc[i]*h for i in sites)),
                    false_accept=bool(z.real>-.0499),GFL_MW=float(PORTS.P0@p[:10]))
                if row['false_accept']:save('candidate_'+str(len(rows))+'.json',{'p':p.tolist(),'root':[z.real,z.imag],**row})
            except Exception as ex:row['error']=str(ex)
            rows.append(row)
        pd.DataFrame(rows).to_csv(OUT/'TABLE_02_COMPENSATED_COMBINATIONS.csv',index=False)
        pd.DataFrame(sr).to_csv(OUT/'TABLE_03_SINGLETON_COMPENSATION.csv',index=False)
        print('h',h,'valid',sum(r['valid'] for r in rows),'violations',sum(r.get('false_accept',False) for r in rows),'seconds',time.perf_counter()-start,flush=True)
    save('COMPENSATED_SUMMARY.json',{'cases':len(rows),'valid':sum(r['valid'] for r in rows),'violations':sum(r.get('false_accept',False) for r in rows),
        'scope':'one tracked pole; individually compensated interventions, exploratory, no nonlinear security or optimality'})
if __name__=='__main__':main()
