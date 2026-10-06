"""Interaction compensation, baseline comparisons and selected branch audit."""
from compensated import *
from damping_only import local

SELECTED=json.loads((OUT/'damping_candidate_6.json').read_text())
PAIR=np.array(SELECTED['p']);SITES=[0,7]

def singleton_budget(i,t):
    p,z,_=local(i,.02);ki=p[20+i];alpha=TARGET.real-t/2
    for f in np.linspace(.1,1,10):
        s=complex(TARGET.real-f*t/2,z.imag)
        for it in range(12):
            p[20+i]=ki;m=Model(p);lu=la.lu_factor(s*np.eye(204)-m.A0,check_finite=False)
            rb=la.lu_solve(lu,m.B,check_finite=False);e=np.exp(-s*TAU)
            K=np.eye(10)-(m.C@rb)*e[None,:]
            left,sv,vh=la.svd(K,check_finite=False);u=left[:,-1];v=vh[-1].conj()
            Ds=(m.C@la.lu_solve(lu,rb,check_finite=False))*e[None,:]+(m.C@rb)*(TAU*e)[None,:]
            dki=-e[i]*np.vdot(u,m.C@la.lu_solve(lu,M['Bi'][:,i],check_finite=False))*v[i]
            ds=1j*np.vdot(u,Ds@v);r=np.vdot(u,K@v)
            update=la.solve([[ds.real,dki.real],[ds.imag,dki.imag]],-np.array([r.real,r.imag]))
            s+=1j*update[0];ki+=update[1]
            assert abs(update[0])<.5 and abs(update[1])<20,'budget corrector jump'
            if abs(update[0])<1e-11 and abs(update[1])<1e-8:break
        else:raise RuntimeError('budget corrector failed')
        z=s
    p[20+i]=ki;assert KI_BOUNDS[0]<=ki<=KI_BOUNDS[1]
    z,res=Model(p).refine(z);assert abs(z.real-alpha)<1e-8
    return p,z

def budget(t):
    pi,zi=singleton_budget(0,t);pj,zj=singleton_budget(7,t)
    p=pi+pj-ANCHOR;z,res=Model(p).refine(complex(*SELECTED['root']))
    interaction=(z-zi-zj+TARGET).real
    return p,z,interaction

def write_design(name,p):
    dest=OUT/'designs';dest.mkdir(exist_ok=True)
    (dest/f'{name}.toml').write_text('\n'.join(f'{k} = [{", ".join(format(float(x),".17g") for x in v)}]' for k,v in [('rho',p[:10]),('Kp',p[10:20]),('Ki',p[20:])])+'\ntau = 0.04\n')

def main():
    # Target -0.06: an additional 0.01/s guard relative to the inherited -0.05.
    target=-.06;b=-TARGET.real+target;t=0.;iterations=[]
    for k in range(30):
        p,z,I=budget(t);nxt=max(0,I-b)
        iterations.append(dict(k=k,t=t,interaction=I,joint_real=z.real,next_t=nxt,residual=nxt-t))
        if abs(nxt-t)<1e-10:break
        t=nxt
    else:raise RuntimeError('fixed point did not converge')
    pd.DataFrame(iterations).to_csv(OUT/'TABLE_06_COMPENSATION_ITERATION.csv',index=False)
    samples=[]
    for tt in np.linspace(0,.04,21):
        _,zz,ii=budget(float(tt));samples.append(dict(t=tt,joint_real=zz.real,interaction=ii))
    pd.DataFrame(samples).to_csv(OUT/'TABLE_07_BUDGET_CURVE.csv',index=False)
    slopes=np.diff([x['interaction'] for x in samples])/.002
    save('BUDGET_RESULT.json',{'p':p.tolist(),'root':[z.real,z.imag],'t':t,'target_real':target,
        'iterations':len(iterations),'max_sampled_abs_interaction_slope':float(max(abs(slopes))),
        'contraction_status':'SUPPORTED_EMPIRICAL_ONLY_NO_UNIFORM_DERIVATIVE_CERTIFICATE',
        'optimality':'not gain effort optimum; scalar fixed-point correction under fixed replacement, Kp and equal budget weights'})
    designs={'anchor':ANCHOR,'single30':local(0,.02)[0],'single37':local(7,.02)[0],
             'joint':PAIR,'complex_pair':combined((0,7),.01),'corrected':p}
    rows=[]
    anchor_roots=np.load(OUT/'state_I_07.npz')['joint']
    for name,pp in designs.items():
        write_design(name,pp);roots=anchor_roots.copy()
        for step in range(21):
            xx=ANCHOR+(pp-ANCHOR)*(step/20);m=Model(xx)
            new=np.array([m.refine(q)[0] for q in roots]);
            separation=min(abs(new[j]-new[k]) for j in range(11) for k in range(j))
            assert separation>1e-5,'catalog collision'
            assert max(abs(new-roots))<.5,'catalog branch jump'
            for j,s in enumerate(new):rows.append(dict(design=name,step=step,mode=j,real=s.real,imag=s.imag,min_separation=separation))
            roots=new
    pd.DataFrame(rows).to_csv(OUT/'TABLE_08_MULTIBRANCH_CONTINUATION.csv',index=False)
    save('DESIGN_SUMMARY.json',{'added_GFL_MW':float(PORTS.P0@(PAIR[:10]-ANCHOR[:10])),
        'total_GFL_MW':float(PORTS.P0@PAIR[:10]),'GFL_percent':float(100*(PORTS.P0@PAIR[:10])/PORTS.P0.sum()),
        'retained_SG_MW':float(PORTS.P0@(1-PAIR[:10])),'prescribed_not_maximized':True})
    print(json.dumps(json.loads((OUT/'BUDGET_RESULT.json').read_text()),indent=2))
if __name__=='__main__':main()
