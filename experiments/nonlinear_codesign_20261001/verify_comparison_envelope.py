"""Recompute nonlinear bounds and check the saved fixed-domain LP exactly."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import numpy as np
from threadpoolctl import threadpool_limits
from modal_comparison_contract import Comparison


def dot(a,b):
    return sum((F(float(x))*F(float(y)) for x,y in zip(a,b)),F(0))


with threadpool_limits(limits=1):
    cc=Comparison("dq"); out=cc.out
    saved=np.load(out/"optimized_comparison_LP.npz")
    base=np.load(out/"block_modal_certificate.npz")
    assert np.array_equal(cc.H,saved["H"])
    env=cc.evaluate(base["radii"],float(saved["input_box"]))
    assert env["output_pass"]
    for key in ("a","M","b"):
        assert np.array_equal(env[key],saved[key]),key
    for key,value in (("U",cc.U),("drift",cc.drift),("domain_x",env["x_radius"])):
        assert np.array_equal(saved[key],value),key
    a,M,b,r=saved["a"],saved["M"],saved["b"],saved["radii"]
    epsilon=float(saved["epsilon"])
    margins=[]; rates=[]; perron=[]
    for i in range(cc.ng):
        coupling=dot(M[i],r)
        ar=F(float(a[i]))*F(float(r[i]))
        margin=ar-coupling-F(float(b[i]))*F(epsilon)-F(float(cc.drift[i]))
        assert margin>0 and F(float(r[i]))>=F(float(cc.initial[i]))
        margins.append(margin);rates.append((ar-coupling)/F(float(r[i])))
        perron.append(coupling/ar)
    for row,hi in zip(cc.U,saved["domain_x"]):
        assert dot(row,r)<=F(float(hi))
    assert epsilon<=float(saved["input_box"])
    A,rhs,p,dual=saved["A"],saved["rhs"],saved["primal"],saved["dual"]
    # Reconstruct the exact binary64 LP definition and its physical variable map.
    definition=json.loads((out/"optimized_comparison_LP.json").read_text())
    rref=base["radii"];epsref=float(base["epsilon"]);guard=definition["domain_guard"]
    n=cc.ng;lower=np.r_[cc.initial/rref,0.]
    Ac=np.column_stack([(M-np.diag(a))*rref[None,:]/(a*rref)[:,None],b*epsref/(a*rref)])
    Ad=np.column_stack([cc.U*rref[None,:]/env["x_radius"][:,None],np.zeros(cc.model.nx)])
    cap=np.zeros((1,n+1));cap[0,-1]=1
    expected_A=np.vstack([Ac,Ad,cap])
    bound=np.r_[-cc.drift/(a*rref)-guard,np.full(cc.model.nx,1-guard),float(saved["input_box"])/epsref]
    assert np.array_equal(A,expected_A) and np.array_equal(rhs,bound-A@lower)
    assert np.array_equal(lower,saved["lower"])
    assert np.array_equal(r,rref*(p+lower)[:-1]) and epsilon==epsref*(p+lower)[-1]
    assert p.min()>=0 and dual.min()>=0
    slacks=[F(float(h))-dot(row,p) for row,h in zip(A,rhs)]
    assert min(slacks)>=0
    upper=[]
    for j in range(n):
        upper.append(min(F(float(rhs[i]))/F(float(A[i,j]))
                         for i in range(n,n+cc.model.nx) if A[i,j]>0 and rhs[i]>=0))
    upper.append(F(float(rhs[-1])))
    residual=[max(F(0),F(int(j==n))-dot(A[:,j],dual)) for j in range(n+1)]
    bound=dot(rhs,dual)+sum(u*e for u,e in zip(upper,residual))
    primal=F(float(p[-1])); assert bound>=primal
    sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
    report=dict(status="FROZEN_LP_REENCLOSURE_AND_RATIONAL_RECHECK_PASSED",
        epsilon=epsilon,min_boundary_margin=float(min(margins)),
        contraction_rate_lower=float(np.nextafter(float(min(rates)),-np.inf)),
        comparison_perron_upper=float(np.nextafter(float(max(perron)),np.inf)),
        exact_contraction_rate_lower=str(min(rates)),min_LP_slack=float(min(slacks)),
        LP_primal=float(primal),LP_dual_upper=float(np.nextafter(float(bound),np.inf)),
        LP_relative_gap_upper=float(np.nextafter(float((bound-primal)/bound),np.inf)),
        frequency_bound=env["frequency"],rocof_bound=env["rocof"],
        voltage_min=env["voltage_min"],voltage_max=env["voltage_max"],
        scope="Fixed-domain comparison LP at fixed gains/replacement; shared interval backend, independently recomputed bounds and rational scalar inequalities. Not global physical optimality.",
        source_hashes={name:sha(out/name) for name in ("model.toml","refined_center.json","block_modal_certificate.npz","optimized_comparison_LP.npz","optimized_comparison_LP.json")},
        script_sha256=sha(Path(__file__)))
    (out/"optimized_comparison_LP_verification.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in report.items() if k not in ("source_hashes","exact_contraction_rate_lower")},indent=2))
