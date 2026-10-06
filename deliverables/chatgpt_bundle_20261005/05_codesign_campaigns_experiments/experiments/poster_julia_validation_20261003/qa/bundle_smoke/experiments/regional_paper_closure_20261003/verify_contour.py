"""Continuous square-contour enclosures for the frozen exported model.

All tests that determine validity use Arb balls. Float midpoints choose positive
weights only; each weighted row bound is subsequently evaluated with Arb.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
from pathlib import Path
import sys,time,json,hashlib,shutil,tomllib,platform
OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
sys.path.insert(0,str(OUT/'vendor'))
from flint import acb,arb,acb_mat,ctx
import flint
import numpy as np
import pandas as pd
from scipy import linalg as la
ctx.prec=128
from paths import resolve
MODEL=resolve('experiments/graph_gsp_codesign_20261003/model')
NAMES=['Adev','Bv','Cs','Cf','Ds','Y','Etheta','Hv','Bp','Bi']
M={k:pd.read_csv(MODEL/(k+'.csv')).to_numpy() for k in NAMES}
base=tomllib.loads((MODEL.parent/'baseline.toml').read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def save(n,o): (OUT/n).write_text(json.dumps(o,indent=2,allow_nan=False)+'\n')
def freeze():
    p=OUT/'THEORY_before_closure.tex'
    if not p.exists(): shutil.copyfile(ROOT/'experiments/theory_collective_damping_20261003/THEORY.tex',p)
    paths=['experiments/graph_gsp_codesign_20261003/model/'+n+'.csv' for n in NAMES]
    paths+=['experiments/graph_gsp_codesign_20261003/baseline.toml','experiments/regional_paper_closure_20261003/PROTOCOL.txt']
    manifest={'inputs':{str(Path(rel)):sha(resolve(rel)) for rel in paths},
        'theory_snapshot_sha256':sha(p),'python':platform.python_version(),
        'flint':flint.__version__,'arb_precision_bits':ctx.prec,
        'model_semantics':'CSV values parsed as exact binary64, no physical uncertainty enclosure'}
    old=OUT/'INPUT_MANIFEST.json'
    if old.exists():
        before=json.loads(old.read_text())
        before['inputs']={k.replace('\\','/'):v for k,v in before['inputs'].items()}
        candidate=dict(manifest)
        candidate['inputs']={k.replace('\\','/'):v for k,v in manifest['inputs'].items()}
        assert before==candidate
    else: save('INPUT_MANIFEST.json',manifest)

def conv(a): return acb_mat([[acb(float(x)) for x in row] for row in a])
Mb={k:conv(v) for k,v in M.items()}
rho=[arb(float(x)) for x in base['rho']]
kp=[arb(float(x)) for x in base['Kp']]
ki=[arb(float(x)) for x in base['Ki']]
tau=arb('0.04')
radii=[arb('.001')]*20+[arb('.01')*v for v in kp]+[arb('.01')*v for v in ki]
G=acb_mat(20,20); NN=acb_mat(20,204)
for i in range(20):
    w=1-rho[i//2]
    for j in range(20): G[i,j]=Mb['Y'][i,j]+w*Mb['Ds'][i,j]
    for j in range(204): NN[i,j]=w*Mb['Cs'][i,j]+(1-w)*Mb['Cf'][i,j]
B=acb_mat(204,10)
for i in range(204):
    for j in range(10): B[i,j]=Mb['Bp'][i,j]*kp[j]+Mb['Bi'][i,j]*ki[j]
U=acb_mat(234,40)
for i in range(20): U[204+i,i]=1
for i in range(204):
    for j in range(10): U[i,20+j]=Mb['Bp'][i,j];U[i,30+j]=Mb['Bi'][i,j]

def Fmat(s):
    e=(-s*tau).exp()
    F=acb_mat(234,234)
    for i in range(204):
        for j in range(204): F[i,j]=-Mb['Adev'][i,j]+(s if i==j else 0)
        for j in range(20): F[i,204+j]=-Mb['Bv'][i,j]
        for j in range(10): F[i,224+j]=-B[i,j]*e
    for i in range(20):
        for j in range(204): F[204+i,j]=NN[i,j]
        for j in range(20): F[204+i,204+j]=G[i,j]
    for i in range(10):
        for j in range(204): F[224+i,j]=-Mb['Etheta'][i,j]
        for j in range(20): F[224+i,204+j]=-Mb['Hv'][i,j]
        F[224+i,224+i]=1
    return F,e

def Vmat(e):
    V=acb_mat(40,234)
    for i in range(20):
        for j in range(204): V[i,j]=-(Mb['Cs'][i,j]-Mb['Cf'][i,j])
        for j in range(20): V[i,204+j]=-Mb['Ds'][i,j]
    for i in range(10): V[20+i,224+i]=-e;V[30+i,224+i]=-e
    return V

def upper(x): return float(x.upper())
def bounds(x): return {'real':x.real.str(35),'imag':x.imag.str(35)}

def panel(a,b,depth=0):
    # a,b are exact dyadic contour coordinates relative to the decimal center.
    cx=arb('-0.12229688588491999');cy=arb('0.32408378548166134')
    ar,ai=a;br,bi=b
    sc=acb(cx+(arb(ar)+arb(br))/2,cy+(arb(ai)+arb(bi))/2)
    s=sc+acb(arb(0,abs(arb(br)-arb(ar))/2),arb(0,abs(arb(bi)-arb(ai))/2))
    t=time.perf_counter()
    try:
        F,e=Fmat(s)
        Q=F.inv()
        T=Vmat(e)*Q*U
        pos=np.array([[float(abs(T[i,j]).mid())*float(radii[i].mid()) for j in range(40)] for i in range(40)])
        q0=1.1*max(abs(la.eigvals(pos)))+1e-10
        weights=la.solve(q0*np.eye(40)-pos,np.ones(40))
        if weights.min()<=0: raise ArithmeticError('nonpositive weight proposal')
        weights/=weights.max()
        wb=[arb(float(w)) for w in weights]
        qb=[radii[i]*sum((abs(T[i,j])*wb[j] for j in range(40)),arb(0))/wb[i] for i in range(40)]
        if not all(v<1 for v in qb): raise ArithmeticError('box envelope not below one')
        qs=max(upper(v) for v in qb)
        tr=sum((Q[i,i] for i in range(204)),acb(0))
        for j in range(10):
            tr+=tau*e*sum((Q[224+j,i]*B[i,j] for i in range(204)),acb(0))
        ds=acb(arb(br)-arb(ar),arb(bi)-arb(ai))
        count=tr*ds
        moment=s*tr*ds
        # Linear moment coefficients from trace(Theta(d) T), including repeated rho.
        coeff=[-(T[2*i,2*i]+T[2*i+1,2*i+1])*ds for i in range(10)]
        coeff += [-T[i,i]*ds for i in range(20,40)]
        # Integrating a uniform induced-norm bound gives a rigorous remainder.
        q=arb(max(v.upper() for v in qb))
        length=abs(ds)
        beta=length*40*q*q/(2*(1-q))
        return [(dict(a=list(a),b=list(b),depth=depth,q_upper=qs,seconds=time.perf_counter()-t,
                     count=bounds(count),moment=bounds(moment),coefficients=[bounds(c) for c in coeff],
                     beta=beta.str(35)),count,moment,coeff,beta)]
    except (ZeroDivisionError,ArithmeticError,ValueError) as exc:
        if depth>=12: raise RuntimeError(f'Panel failed at depth {depth}: {a} {b}: {exc}') from exc
        mid=((ar+br)/2,(ai+bi)/2)
        return panel(a,mid,depth+1)+panel(mid,b,depth+1)

def main():
    freeze()
    # A separate interval inverse proves algebraic regularity over the rho box.
    GG=acb_mat(20,20)
    for i in range(20):
        w=1-rho[i//2]+arb(0,.0010000000000000002)
        for j in range(20): GG[i,j]=Mb['Y'][i,j]+w*Mb['Ds'][i,j]
    _=GG.inv()
    print('G regularity proved over box',flush=True)
    # 0.04 represented as binary64 is explicitly the geometric contour half-side.
    h=.04; corners=[(-h,-h),(h,-h),(h,h),(-h,h),(-h,-h)]
    total=[];start=time.perf_counter()
    for edge in range(4):
        a,b=corners[edge:edge+2]
        for j in range(32):
            aa=tuple(a[k]+(b[k]-a[k])*j/32 for k in range(2))
            bb=tuple(a[k]+(b[k]-a[k])*(j+1)/32 for k in range(2))
            rows=panel(aa,bb);total.extend(rows)
            print(f'edge={edge} panel={j} accepted={len(rows)} total={len(total)} seconds={time.perf_counter()-start:.1f}',flush=True)
            save('PANELS.json',[x[0] for x in total])
    denom=2*arb.pi()*acb(0,1)
    count=sum((x[1] for x in total),acb(0))/denom
    moment=sum((x[2] for x in total),acb(0))/denom
    coeff=[sum((x[3][j] for x in total),acb(0))/denom for j in range(30)]
    beta=sum((x[4] for x in total),arb(0))/(2*arb.pi())
    # A complex interval with real part strictly in (0.5,1.5) and containing 1+0i
    # identifies the integer count uniquely. No nearest-integer rounding claim.
    count_one=bool(count.real>arb('.5') and count.real<arb('1.5') and count.contains(acb(1)))
    save('CONTOUR_RESULT.json',{'G_regular':True,'all_panels_regular':True,'box_small_gain':True,
         'panel_count':len(total),'q_max_upper':max(x[0]['q_upper'] for x in total),
         'count':bounds(count),'count_one_proved':count_one,'root_sum':bounds(moment),
         'affine_coefficients':[bounds(x) for x in coeff],
         'remainder_upper_for_count_one':str(beta),
         'seconds':time.perf_counter()-start,'scope':'exported binary64 model, one contour only'})
    print((OUT/'CONTOUR_RESULT.json').read_text())

if __name__=='__main__': main()
