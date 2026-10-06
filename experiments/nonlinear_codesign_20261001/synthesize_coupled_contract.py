"""Search full-network invariant ellipsoids using validated Jacobian enclosures.

The nominal modal matrix only initializes the metric. Every accepted region must
pass a derivative enclosure of the full nonlinear DAE and all 39 instruments.
"""
from pathlib import Path
import hashlib
import json
import time
import argparse
import numpy as np
import scipy.linalg as la
import mpmath as mp
from threadpoolctl import threadpool_limits
from coupled_interval_model import Model, OUT, mm, add, neg, up, down, magnitude, inverse_enclosure, lower, upper, positive_product


def symmetric_upper(a):
    """Upper bound via approximate eigencoordinates and interval congruence.

    The approximate eigenvectors need not be exactly orthogonal. Their Gram
    bounds control the metric distortion, and a strict lower Gram bound ensures
    invertibility. Only directed Gershgorin row bounds are used for acceptance.
    """
    lo,hi=a
    mid=(lo+hi)/2;mid=(mid+mid.T)/2
    _,q=la.eigh(mid)
    qt=(q.T,q.T);qi=(q,q)
    transformed=mm(qt,mm(a,qi))
    gram=mm(qt,qi)
    n=len(q)
    def bounds(v):
        off=magnitude(v).copy();np.fill_diagonal(off,0)
        sums=up(off.sum(axis=1)*(1+4*n*np.finfo(float).eps))
        return float(np.min(down(np.diag(v[0])-sums))),float(np.max(up(np.diag(v[1])+sums)))
    gl,gu=bounds(gram);_,hu=bounds(transformed)
    assert gl>0
    return float(up(hu/(gl if hu>=0 else gu)))


def norm_upper(a):
    gram=mm((a[0].T,a[1].T),a)
    return float(up(np.sqrt(max(0.,symmetric_upper(gram)))))


def metric(A,shift=.025,B=None,kind="observability"):
    ab,d=la.matrix_balance(A,permute=False,scale=True)
    if kind=="reachability":
        bb=np.linalg.solve(d,B)
        x=la.solve_continuous_lyapunov(ab+shift*np.eye(len(A)),-(bb@bb.T+1e-12*np.eye(len(A))))
        x=(x+x.T)/2
        h=np.linalg.solve(la.cholesky(x,lower=True),np.eye(len(A)))@np.linalg.inv(d)
    else:
        p=la.solve_continuous_lyapunov((ab+shift*np.eye(len(A))).T,-np.eye(len(A)))
        p=(p+p.T)/2
        h=la.cholesky(p,lower=False)@np.linalg.inv(d)
    t=np.linalg.inv(h)
    h*=np.max(np.linalg.norm(t,axis=1))
    ti,invcheck=inverse_enclosure((h,h))
    return h,ti,invcheck


def search(metric_kind="modal"):
    center=json.loads((OUT/"refined_center.json").read_text(encoding="utf-8"))
    m=Model(center=center["offset"]);A=np.array(m.raw["A"])
    if metric_kind=="modal":
        from probe_coupled_modal_metric import modal_metric
        h,ti,invcheck=modal_metric(A,len(m.keep))
    else:
        h,ti,invcheck=metric(A,B=np.array(m.raw["B"]),kind="observability")
    ht=(h,h)
    # This is only a starting matrix. J(0) is interval recomputed independently.
    nominal=mm(ht,mm((A,A),ti))
    nominal_sym=add(nominal,(nominal[0].T,nominal[1].T))
    nominal_decay=-symmetric_upper(nominal_sym)/2
    trows=up(np.sqrt(up(np.sum(magnitude(ti)**2,axis=1)*(1+8*m.nx*np.finfo(float).eps))))
    print(json.dumps(dict(stage="metric",dimension=m.nx,condition_H=float(np.linalg.cond(h)),
                          metric_inverse_check=invcheck,nominal_decay=nominal_decay)),flush=True)
    start=time.perf_counter()
    baseline=m.evaluate([mp.iv.mpf(0)]*m.nx,[mp.iv.mpf(0)]*m.nw)[0]
    f0=(np.array([lower(v) for v in baseline])[:,None],np.array([upper(v) for v in baseline])[:,None])
    hf0=mm(ht,f0)
    drift=float(up(np.sqrt(up(np.sum(magnitude(hf0)**2)*(1+8*m.nx*np.finfo(float).eps)))))
    print(json.dumps(dict(stage="baseline",seconds=time.perf_counter()-start,drift=drift,
                          max_actual_rhs=float(magnitude(f0).max()))),flush=True)
    attempts=[];best=None
    radii=np.logspace(-8,-5,10) if metric_kind=="modal" else np.logspace(-10,-7,10)
    for r in radii:
        epsilon=0.
        for iteration in range(3):
            out,v,ff,rr=m.intervals(up(r*trows),np.full(m.nw,epsilon))
            jl=np.array([x.d[0][:m.nx] for x in out]);jh=np.array([x.d[1][:m.nx] for x in out])
            dl=down(jl-A);dh=up(jh-A)
            remainder=mm(ht,mm((dl,dh),ti))
            sy=add(nominal_sym,add(remainder,(remainder[0].T,remainder[1].T)))
            decay=-symmetric_upper(sy)/2
            bj=(np.array([x.d[0][m.nx:] for x in out]),np.array([x.d[1][m.nx:] for x in out]))
            inputgain=norm_upper(mm(ht,bj))
            allowed=max(0.,(decay*r-drift)/inputgain)
            if iteration<2:
                epsilon=.8*allowed
            else:
                inward=decay*r-inputgain*epsilon-drift
                fmax=max(magnitude(x.v) for x in ff);rmax=max(magnitude(x.v) for x in rr)
                vmlo=[];vmhi=[]
                for b in range(39):
                    vm=v[2*b]**2+v[2*b+1]**2
                    vmlo.append(np.sqrt(max(0.,vm.v[0])));vmhi.append(np.sqrt(max(0.,vm.v[1])))
                passed=epsilon>0 and inward>0 and fmax<.5 and rmax<.5 and min(vmlo)>.9 and max(vmhi)<1.1
                row=dict(radius=float(r),epsilon=float(epsilon),decay=float(decay),input_gain=float(inputgain),
                         trim_drift=drift,inward_margin=float(inward),f_bound=float(fmax),rocof_bound=float(rmax),
                         voltage_min=float(min(vmlo)),voltage_max=float(max(vmhi)),passed=bool(passed),
                         inverse_neumann_check=m.last_inverse_residual)
                attempts.append(row);print(json.dumps(row),flush=True)
                if passed and (best is None or epsilon>best["epsilon"]):
                    best=row
                    np.savez_compressed(OUT/f"certificate_{metric_kind}.npz",H=h,Tlo=ti[0],Thi=ti[1],
                        Jlo=jl,Jhi=jh,Blo=bj[0],Bhi=bj[1],f0lo=f0[0],f0hi=f0[1],
                        x_radius=r*trows,radius=r,epsilon=epsilon)
    report=dict(status="COUPLED_ELLIPSOID_FOUND_PENDING_INDEPENDENT_VERIFICATION" if best else "NO_COUPLED_ELLIPSOID_FOUND",metric=metric_kind,
                best=best,attempts=attempts,state_count=m.nx,input_count=m.nw,
                nominal_decay=nominal_decay,metric_condition=float(np.linalg.cond(h)),
                replacement_fraction=float(np.dot(m.raw["rho"],m.raw["power_MW"])/sum(m.raw["power_MW"])),
                model_sha256=hashlib.sha256((OUT/"model.toml").read_bytes()).hexdigest(),
                center_sha256=hashlib.sha256((OUT/"refined_center.json").read_bytes()).hexdigest(),
                implementation_sha256=hashlib.sha256(Path(__file__).with_name("coupled_interval_model.py").read_bytes()).hexdigest(),
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                scope="Fixed design; all measurable additive currents at all 39 buses and DC supply inputs, in a joint Euclidean envelope; no optimized capacity claim")
    (OUT/f"synthesis_{metric_kind}.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":report["status"],"best":best}),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--metric",choices=("modal","observability"),default="modal")
    args=parser.parse_args()
    with threadpool_limits(limits=1):search(args.metric)
