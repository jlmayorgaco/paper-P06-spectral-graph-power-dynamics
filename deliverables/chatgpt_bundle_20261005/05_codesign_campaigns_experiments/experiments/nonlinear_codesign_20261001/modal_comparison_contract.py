"""Nonlinear block-modal invariant regions with separate time scales.

V(z)=max_k ||(H z)_k||/r_k. Interval Jacobians give a Metzler comparison system.
This preserves mode-specific decay/input gains instead of multiplying the slowest
decay by the largest (typically fast-control) disturbance gain.
"""
from pathlib import Path
import argparse
import hashlib
import importlib
import json
import time
import numpy as np
import scipy.linalg as la
from threadpoolctl import threadpool_limits
from coupled_interval_model import mm,add,up,down,magnitude,inverse_enclosure,positive_product
from synthesize_coupled_contract import symmetric_upper


def block_metric(A,nplant):
    ap=A[:nplant,:nplant];am=A[nplant:,nplant:];c=A[nplant:,:nplant]
    ev,vec=la.eig(ap);columns=[];groups=[];modes=[]
    for i,e in enumerate(ev):
        if e.imag>1e-7:
            pair=np.column_stack([vec[:,i].real,vec[:,i].imag]);pair/=la.norm(pair)
            groups.append(list(range(len(columns),len(columns)+2)))
            columns.extend([pair[:,0],pair[:,1]]);modes.append([float(e.real),float(e.imag)])
        elif abs(e.imag)<=1e-7:
            groups.append([len(columns)]);columns.append(vec[:,i].real/la.norm(vec[:,i].real))
            modes.append([float(e.real),0.])
    tp=np.column_stack(columns);assert tp.shape==ap.shape
    cross=la.solve_sylvester(am,-ap,-c)
    tm=np.diag(np.r_[np.ones(39),5*np.ones(39)])
    T=np.block([[tp,np.zeros((nplant,78))],[cross@tp,tm]])
    for bus in range(39):
        groups.append([nplant+bus,nplant+39+bus]);modes.append([-10.,0.])
    H=np.linalg.inv(T);Ti,check=inverse_enclosure((H,H))
    return H,Ti,groups,modes,check


def vector_norm_bound(values):
    q=np.asarray(values)
    return float(up(np.sqrt(up(np.sum(q*q)*(1+8*q.size*np.finfo(float).eps)))))


class Comparison:
    def __init__(self,coordinate="dq"):
        module=importlib.import_module("coupled_dq_model" if coordinate=="dq" else "coupled_interval_model")
        self.out=module.OUT;self.coordinate=coordinate;self.module=module
        self.center=json.loads((self.out/"refined_center.json").read_text(encoding="utf-8"))
        self.model=module.Model(center=self.center["offset"])
        m=self.model
        self.H,self.Ti,self.groups,self.modes,self.inverse_check=block_metric(np.array(m.raw["A"]),len(m.keep))
        self.ng=len(self.groups)
        self.U=np.column_stack([up(np.sqrt(up(np.sum(magnitude(self.Ti)[:,g]**2,axis=1)*(1+8*len(g)*np.finfo(float).eps)))) for g in self.groups])
        fi=(np.array(self.center["f0lo"])[:,None],np.array(self.center["f0hi"])[:,None])
        hi=mm((self.H,self.H),fi)
        self.drift=np.array([vector_norm_bound(magnitude(hi)[g]) for g in self.groups])
        zi=(-m.center[:,None],-m.center[:,None]);initial=mm((self.H,self.H),zi)
        self.initial=np.array([vector_norm_bound(magnitude(initial)[g]) for g in self.groups])

    def evaluate(self,r,epsilon):
        m=self.model
        xr=positive_product(self.U,r[:,None]).ravel()
        out,v,ff,rr=m.intervals(xr,np.full(m.nw,epsilon))
        J=(np.array([x.d[0][:m.nx] for x in out]),np.array([x.d[1][:m.nx] for x in out]))
        B=(np.array([x.d[0][m.nx:] for x in out]),np.array([x.d[1][m.nx:] for x in out]))
        AJ=mm((self.H,self.H),mm(J,self.Ti));HB=mm((self.H,self.H),B)
        M=np.zeros((self.ng,self.ng));a=np.zeros(self.ng);b=np.zeros(self.ng)
        aa=magnitude(AJ);bb=magnitude(HB)
        for i,g in enumerate(self.groups):
            diag=(AJ[0][np.ix_(g,g)],AJ[1][np.ix_(g,g)])
            sym=add(diag,(diag[0].T,diag[1].T))
            a[i]=down(-symmetric_upper(sym)/2)
            b[i]=vector_norm_bound(bb[g])
            for j,h in enumerate(self.groups):
                if i!=j:M[i,j]=vector_norm_bound(aa[np.ix_(g,h)])
        rhs=up(positive_product(M,r[:,None]).ravel()+up(b*epsilon))
        rhs=up(rhs+self.drift)
        margins=down(down(a*r)-rhs)
        fmax=max(magnitude(q.v) for q in ff);rmax=max(magnitude(q.v) for q in rr)
        vv=[v[2*i]**2+v[2*i+1]**2 for i in range(39)]
        vmin=min(down(np.sqrt(max(0.,q.v[0]))) for q in vv)
        vmax=max(up(np.sqrt(max(0.,q.v[1]))) for q in vv)
        return dict(a=a,M=M,b=b,rhs=rhs,margins=margins,x_radius=xr,
                    frequency=float(fmax),rocof=float(rmax),voltage_min=float(vmin),voltage_max=float(vmax),
                    output_pass=bool(fmax<.5 and rmax<.5 and vmin>.9 and vmax<1.1),
                    inverse_network_check=m.last_inverse_residual)

    def trial(self,epsilon,base,max_iterations=12):
        r=1.25*np.maximum.reduce([base["b"]*epsilon/base["a"],self.initial,self.drift/base["a"],np.full(self.ng,1e-26)])
        trace=[];result=None
        for k in range(max_iterations):
            try:result=self.evaluate(r,epsilon)
            except (AssertionError,ValueError) as e:
                return dict(status="DOMAIN_ENCLOSURE_FAILED",epsilon=epsilon,trace=trace,error=str(e)),None
            a,M,b=result["a"],result["M"],result["b"]
            if a.min()<=0:
                trace.append(dict(iteration=k,min_decay=float(a.min())))
                return dict(status="NONPOSITIVE_BLOCK_DECAY",epsilon=epsilon,trace=trace),None
            ratio=float(np.max(result["rhs"]/(a*r)))
            # Perron root is a diagnostic/iteration guide; acceptance uses directed
            # inequalities for this explicit positive radius vector.
            spectral=float(max(abs(la.eigvals(M/a[:,None]))))
            trace.append(dict(iteration=k,min_decay=float(a.min()),max_comparison_usage=ratio,
                              perron_root_diagnostic=spectral,max_coordinate_radius=float(result["x_radius"].max()),
                              frequency=result["frequency"],rocof=result["rocof"],output_pass=result["output_pass"]))
            if result["margins"].min()>0 and np.all(r>=self.initial) and result["output_pass"]:
                return dict(status="BLOCK_MODAL_REGION_FOUND",epsilon=epsilon,trace=trace,
                            min_margin=float(result["margins"].min()),max_comparison_usage=ratio,
                            frequency=result["frequency"],rocof=result["rocof"],
                            voltage_min=result["voltage_min"],voltage_max=result["voltage_max"]),dict(r=r,**result)
            if spectral>=.999:
                return dict(status="COMPARISON_SMALL_GAIN_FAILED",epsilon=epsilon,trace=trace),None
            if not result["output_pass"]:
                return dict(status="OUTPUT_ENCLOSURE_FAILED",epsilon=epsilon,trace=trace),None
            proposed=np.linalg.solve(np.diag(a)-M,b*epsilon+self.drift)
            assert proposed.min()>=-1e-16
            r=1.15*np.maximum.reduce([proposed,self.initial,np.full(self.ng,1e-26)])
        return dict(status="RADIUS_FIXED_POINT_NOT_CONVERGED",epsilon=epsilon,trace=trace),None


def run(coordinate):
    cc=Comparison(coordinate);m=cc.model
    base=cc.evaluate(np.maximum(cc.initial*1.01,1e-27),0.)
    print(json.dumps(dict(stage="block_metric",coordinate=coordinate,blocks=cc.ng,states=m.nx,
                          condition=float(np.linalg.cond(cc.H)),minimum_block_decay=float(base["a"].min()),
                          maximum_block_decay=float(base["a"].max()),maximum_input_gain=float(base["b"].max()))),flush=True)
    assert base["a"].min()>0
    attempts=[];best=None;best_certificate=None
    for eps in np.logspace(-12,-1,12):
        row,cert=cc.trial(float(eps),base);attempts.append(row)
        print(json.dumps(row),flush=True)
        if cert is None:
            if best is not None:
                lo=best["epsilon"];hi=float(eps)
                for _ in range(5):
                    mid=np.sqrt(lo*hi);rr,qq=cc.trial(mid,base);attempts.append(rr)
                    print(json.dumps(rr),flush=True)
                    if qq is not None:lo=mid;best=rr;best_certificate=qq
                    else:hi=mid
            break
        best=row;best_certificate=cert
    if best is not None:
        cert=best_certificate
        np.savez_compressed(cc.out/"block_modal_certificate.npz",H=cc.H,Tlo=cc.Ti[0],Thi=cc.Ti[1],
            U=cc.U,radii=cert["r"],a=cert["a"],M=cert["M"],b=cert["b"],drift=cc.drift,
            initial=cc.initial,x_radius=cert["x_radius"],epsilon=best["epsilon"])
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    report=dict(status="BLOCK_MODAL_CANDIDATE_PENDING_RECHECK" if best else "NO_BLOCK_MODAL_REGION_FOUND",
        coordinate=coordinate,blocks=cc.ng,groups=cc.groups,nominal_modes=cc.modes,best=best,attempts=attempts,
        state_count=m.nx,input_count=m.nw,metric_inverse_check=cc.inverse_check,
        model_sha256=sha(cc.out/"model.toml"),center_sha256=sha(cc.out/"refined_center.json"),
        implementation_sha256=sha(Path(cc.module.__file__)),script_sha256=sha(Path(__file__)),
        scope="Fixed physical candidate; nonlinear block-modal invariant region; all declared measurable current/DC input signals, not an optimized replacement capacity")
    (cc.out/"block_modal_synthesis.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(dict(status=report["status"],best=best)),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--coordinates",choices=("dq","ri"),default="dq")
    args=parser.parse_args()
    with threadpool_limits(limits=1):run(args.coordinates)
