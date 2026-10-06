"""Numerical contour diagnostics for the selected pair; not certification."""
from compensated import *
from damping_only import local

def points(box,n):
    x0,x1,y0,y1=box;c=[x0+1j*y0,x1+1j*y0,x1+1j*y1,x0+1j*y1,x0+1j*y0]
    for a,b in zip(c,c[1:]):
        for k in range(n):yield a+(b-a)*(k+.5)/n,(b-a)/n

def main():
    selected=json.loads((OUT/'damping_candidate_6.json').read_text());p=np.array(selected['p'])
    i,j=selected['bus_i']-30,selected['bus_j']-30
    a=Model(ANCHOR);ii,jj=GROUPS[i],GROUPS[j];ix=np.r_[ii,jj];d=action(p-ANCHOR)[ix]
    pi=ANCHOR.copy();pj=ANCHOR.copy()
    for site,pp in [(i,pi),(j,pj)]:
        idx=[site,10+site,20+site];pp[idx]=p[idx]
    models=[a,Model(pi),Model(pj),Model(p)]
    rows=[]
    for box in [(-.105,.01,30.40,31.05),(-.10,.05,30.35,31.10),(-.09,0,30.45,31.02),(-.11,.1,30.30,31.20)]:
        for n in [64,128]:
            total=0j;beta=0.;pr=0.;kr=0.;counts=np.zeros(4,complex);mom=counts.copy();deterr=0.;bad=0
            for s,ds in points(box,n):
                T=a.T(s)[np.ix_(ix,ix)];H=d[:,None]*T
                K12=la.solve(np.eye(4)+H[:4,:4],H[:4,4:]);K21=la.solve(np.eye(4)+H[4:,4:],H[4:,:4])
                PP=K12@K21;MP=abs(K12)@abs(K21)
                r=max(abs(la.eigvals(MP)));pr=max(pr,float(r));kr=max(kr,float(max(abs(la.eigvals(PP)))))
                total+=np.trace(PP)*ds/(2j*np.pi)
                if r<1:
                    bound=(-np.linalg.slogdet(np.eye(4)-MP)[1]-np.trace(MP))*abs(ds)/(2*np.pi)
                    beta+=float(bound)
                else:bad+=1
                for k,m in enumerate(models):
                    D=m.delta(s);tr=np.trace(la.solve(D,m.derivative(s),check_finite=False))
                    counts[k]+=tr*ds/(2j*np.pi);mom[k]+=s*tr*ds/(2j*np.pi)
            row=dict(box=box,n=n,counts=[[z.real,z.imag] for z in counts],moments=[[z.real,z.imag] for z in mom],
                C2=[total.real,total.imag],remainder_pointwise_trapezoid=beta if bad==0 else None,
                max_pair_majorant_radius=pr,max_pair_complex_radius=kr,bad_points=bad,
                exact_interaction_from_moments=[(mom[3]-mom[1]-mom[2]+mom[0]).real,(mom[3]-mom[1]-mom[2]+mom[0]).imag])
            rows.append(row);save('CONTOUR_DIAGNOSTICS.json',rows)
            print('box',box,'n',n,'counts',counts.real,'C2',total,'beta',beta if not bad else None,'r',pr,flush=True)
if __name__=='__main__':main()
