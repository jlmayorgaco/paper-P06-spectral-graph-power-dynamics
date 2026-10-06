"""Optimize modal radii/input size in a frozen nonlinear-domain comparison LP.

This is an inner certificate problem at fixed rho/Kp/Ki, not physical replacement
optimality. The LP table and its primal/dual bounds are independently rationally
checked, and the resulting region is checked against unnormalized coefficients.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.optimize import linprog
from threadpoolctl import threadpool_limits
from modal_comparison_contract import Comparison
from coupled_interval_model import positive_product,up


def dot(a,b):return sum((F(float(x))*F(float(y)) for x,y in zip(a,b)),F(0))


with threadpool_limits(limits=1):
    cc=Comparison("dq");out=cc.out
    old=np.load(out/"block_modal_certificate.npz")
    assert np.max(np.abs(cc.H-old["H"]))<1e-8
    rref=old["radii"];epsref=float(old["epsilon"])
    epsbox=1e-6
    envelope=cc.evaluate(rref,epsbox)
    assert envelope["output_pass"]
    a,M,b=envelope["a"],envelope["M"],envelope["b"]
    assert a.min()>0
    xbar=envelope["x_radius"]
    guard=1e-6;n=cc.ng
    Ac=np.column_stack([(M-np.diag(a))*rref[None,:]/(a*rref)[:,None],b*epsref/(a*rref)])
    bc=-cc.drift/(a*rref)-guard
    Ad=np.column_stack([cc.U*rref[None,:]/xbar[:,None],np.zeros(cc.model.nx)])
    bd=np.full(cc.model.nx,1-guard)
    Acap=np.zeros((1,n+1));Acap[0,-1]=1
    A=np.vstack([Ac,Ad,Acap]);bound=np.r_[bc,bd,epsbox/epsref]
    lower=np.r_[cc.initial/rref,np.float64(0)]
    # Nonnegative shifted variables expose a simple dual residual certificate.
    rhs=bound-A@lower
    cost=np.zeros(n+1);cost[-1]=-1
    result=linprog(cost,A_ub=A,b_ub=rhs,bounds=(0,None),method="highs",
                   options={"dual_feasibility_tolerance":1e-9,"primal_feasibility_tolerance":1e-9})
    assert result.success,result.message
    dual=np.maximum(-result.ineqlin.marginals,0.)
    # Repair primal tolerances by solving a strictly interior auxiliary problem.
    # Retain the original LP table and dual for the optimality upper bound.
    repair_margin=1e-7
    repaired=linprog(cost,A_ub=A,b_ub=rhs-repair_margin,bounds=(0,None),method="highs",
                    options={"dual_feasibility_tolerance":1e-9,"primal_feasibility_tolerance":1e-9})
    assert repaired.success,repaired.message
    candidate=np.maximum(repaired.x,0.)
    slacks=[F(float(bb))-dot(row,candidate) for row,bb in zip(A,rhs)]
    assert min(slacks)>=0,float(min(slacks))
    x=candidate+lower;r=rref*x[:-1];epsilon=epsref*x[-1]
    margins=[]
    for i in range(n):
        m=F(float(a[i]))*F(float(r[i]))-dot(M[i],r)-F(float(b[i]))*F(epsilon)-F(float(cc.drift[i]))
        assert m>0,(i,float(m));margins.append(m)
    assert np.all(r>=cc.initial)
    # Exact unnormalized domain containment; entries of U already upper bounds.
    for row,hi in zip(cc.U,xbar):assert dot(row,r)<=F(float(hi))
    assert epsilon<=epsbox
    # Explicit finite variable bounds: each nonnegative radius is limited by
    # at least one positive domain coefficient. The cap bounds epsilon.
    upper=[]
    domain_start=n;domain_end=n+cc.model.nx
    for j in range(n):
        ratios=[F(float(rhs[i]))/F(float(A[i,j])) for i in range(domain_start,domain_end) if A[i,j]>0 and rhs[i]>=0]
        upper.append(min(ratios))
    upper.append(F(float(rhs[-1])))
    residual=[]
    for j in range(n+1):
        val=F(1 if j==n else 0)-dot(A[:,j],dual)
        residual.append(max(F(0),val))
    dual_upper=dot(rhs,dual)+sum(e*u for e,u in zip(residual,upper))
    primal=F(float(candidate[-1]));assert dual_upper>=primal
    gap=(dual_upper-primal)/dual_upper
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    report=dict(status="FROZEN_COMPARISON_LP_PRIMAL_DUAL_AND_NONLINEAR_CONTAINMENT_VERIFIED",
        input_radius=epsilon,input_radius_before=epsref,input_envelope_used_for_derivative_bounds=epsbox,
        improvement=epsilon/epsref,lp_objective_verified=float(primal),lp_objective_upper_verified=float(dual_upper),
        lp_relative_gap_verified=float(gap),exact_dual_upper=str(dual_upper),
        min_boundary_margin=float(min(margins)),
        frequency_bound=envelope["frequency"],rocof_bound=envelope["rocof"],
        voltage_min=envelope["voltage_min"],voltage_max=envelope["voltage_max"],
        domain_guard=guard,primal_repair_margin=repair_margin,
        exact_LP_scope="Binary64 normalized LP table is treated as the fixed problem; rational dual residual bound. Physical region also checked against unnormalized coefficients.",
        scope="Optimality bound for the frozen comparison LP only; no optimum over nonlinear domains, metrics, gains, replacement or physical network.",
        source_hashes={name:sha(out/name) for name in ("model.toml","refined_center.json","block_modal_certificate.npz")},
        model_implementation_sha256=sha(Path(cc.module.__file__)),
        comparison_implementation_sha256=sha(Path(__file__).with_name("modal_comparison_contract.py")),script_sha256=sha(Path(__file__)))
    np.savez_compressed(out/"optimized_comparison_LP.npz",A=A,rhs=rhs,primal=candidate,dual=dual,
                        lower=lower,H=cc.H,radii=r,epsilon=epsilon,a=a,M=M,b=b,drift=cc.drift,
                        U=cc.U,domain_x=xbar,input_box=epsbox)
    (out/"optimized_comparison_LP.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in report.items() if k not in ("exact_dual_upper","source_hashes")},indent=2))
