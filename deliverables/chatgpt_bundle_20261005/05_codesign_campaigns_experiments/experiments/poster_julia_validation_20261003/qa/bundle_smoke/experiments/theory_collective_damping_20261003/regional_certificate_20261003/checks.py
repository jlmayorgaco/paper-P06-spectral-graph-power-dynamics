"""Small theory checks only. No grid optimization, trajectories or DDE sweep."""
import os
os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["MKL_NUM_THREADS"]="1"
from pathlib import Path
import hashlib
import json
import shutil
import platform
import tomllib
import numpy as np
import pandas as pd
import scipy
from scipy import linalg as la
from scipy.optimize import linprog
import sympy as sp
import mpmath as mp

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]
THEORY=OUT.parent/"THEORY.tex"
MODEL=ROOT/"experiments/graph_gsp_codesign_20261003/model"

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def save(name,obj):
    (OUT/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")

def rel(a,b):
    return float(la.norm(a-b)/max(la.norm(a),la.norm(b),1e-30))

def phi_bound(q):
    """Conservative analytic upper expression; float output is not a certificate."""
    if not 0<=q<1:
        raise ValueError("Norm condition unresolved: require a verified 0<=q<1.")
    return q*q/(2*(1-q))

def freeze():
    snapshot=OUT/"source_snapshot"
    snapshot.mkdir(exist_ok=True)
    target=snapshot/"THEORY_before_regional.tex"
    manifest_path=OUT/"INPUT_MANIFEST.json"
    if manifest_path.exists():
        old=json.loads(manifest_path.read_text())
        assert sha(target)==old["theory_before_sha256"]
        for name,expected in old["inputs"].items():
            assert sha(ROOT/name)==expected,name
        assert sha(OUT/"PROTOCOL.txt")==old["protocol_sha256"]
        return
    shutil.copyfile(THEORY,target)
    inputs=list(MODEL.glob("*.csv"))+[
        ROOT/"experiments/graph_gsp_codesign_20261003/baseline.toml",
        ROOT/"experiments/graph_gsp_codesign_20261003/graph_design.py",
        ROOT/"experiments/all_pll_gain_map_validation_20261003/SUMMARY.json",
        ROOT/"experiments/all_pll_gain_map_validation_20261003/POSTER_CLAIMS.txt"]
    save("INPUT_MANIFEST.json",{
        "protocol_sha256":sha(OUT/"PROTOCOL.txt"),
        "theory_before_sha256":sha(target),
        "inputs":{str(p.relative_to(ROOT)):sha(p) for p in inputs},
        "python":platform.python_version(),"numpy":np.__version__,
        "scipy":scipy.__version__,"sympy":sp.__version__,"mpmath":mp.__version__})

def symbolic():
    checks=[]
    def add(name,value):
        passed=bool(value)
        checks.append({"check":name,"pass":passed,"status":"EXACT_IDENTITY"})
        assert passed,name
    t=sp.symbols("t")
    G0=sp.Matrix([[3,1],[1,2]])
    dG=t*sp.Matrix([[1,2],[-1,1]])
    N0=sp.Matrix([[1,0,2],[0,2,1]])
    dN=t*sp.Matrix([[1,-1,0],[2,0,1]])
    V0=-G0.inv()*N0
    J=G0.inv()*dG
    V1=-G0.inv()*(dN+dG*V0)
    dV=-(G0+dG).inv()*(N0+dN)-V0
    simpl=lambda M:M.applyfunc(sp.factor)
    add("exact algebraic voltage increment",simpl(dV-(sp.eye(2)+J).inv()*V1)==sp.zeros(2,3))
    add("exact algebraic inverse remainder",simpl(dV-V1+(sp.eye(2)+J).inv()*J*V1)==sp.zeros(2,3))
    a,b,c,d,z=sp.symbols("a b c d z")
    H=sp.Matrix([[a,b],[c,d]])
    add("closed-walk decomposition order two",sp.expand(sp.trace(H**2)-(a*a+d*d+2*b*c))==0)
    D=sp.diag(a,d)
    T=(sp.eye(2)+D).inv()*(H-D)
    add("off-diagonal correction has zero trace",sp.trace(T)==0)
    add("block-diagonal determinant separation",
        sp.factor((sp.eye(2)+H).det()-(sp.eye(2)+D).det()*(sp.eye(2)+T).det())==0)
    K=sp.Matrix([[1,2,0],[0,1,1],[2,0,3]])
    U=sp.Matrix([[1,0],[1,1],[0,2]])
    W=z*sp.Matrix([[1,2,0],[0,-1,1]])
    add("rectangular determinant lemma",
        sp.factor((K+U*W).det()-K.det()*(sp.eye(2)+W*K.inv()*U).det())==0)
    # Boundary integration sign: Delta=s-a-p, radius R, contour around a.
    # The first-order integrand is -p/(s-a); minus its contour integral is +p.
    s,p,aa=sp.symbols("s p aa")
    add("root-sum contour integration sign",-sp.residue(-p/(s-aa),s,aa)==p)
    eta=sp.Matrix([sp.Rational(1,2),sp.Rational(1,2)])
    A=sp.Matrix([1,1]);C=sp.Matrix([1,-1]);rhs=sp.Matrix([sp.Rational(11,10)]*2)
    add("dual witness cancels controller",(eta.T*C)[0]==0)
    add("dual witness matches replacement objective",(eta.T*A)[0]==1)
    add("dual objective exact rational",(eta.T*rhs)[0]==sp.Rational(11,10))
    add("proposed step excluded independently of control",
        (eta.T*(A*sp.Rational(6,5)-rhs))[0]==sp.Rational(1,10))
    save("SYMBOLIC_CHECKS.json",checks)
    return checks

def dde_fixture():
    # These constants bound the entire contour and parameter interval:
    # ||Delta0^-1||<=4, exp(-Re(s)/10)<=exp(1/8)<=8/7,
    # ||[[1,1],[1,-1]]||=sqrt(2)<=3/2, |p|<=1/100.
    q=sp.Rational(4)*sp.Rational(1,100)*sp.Rational(8,7)*sp.Rational(3,2)
    beta=sp.Rational(1,4)*2*q*q/(2*(1-q))
    assert q==sp.Rational(12,175) and beta==sp.Rational(36,28525)
    rows=[]
    mp.mp.dps=70
    B=np.array([[1.,1.],[1.,-1.]])
    for p in [-.01,-.005,0.,.005,.01]:
        pm=mp.mpf(str(p))
        determinant=lambda s:(s+1)*(s+2)-pm*mp.exp(-s/10)-2*pm**2*mp.exp(-s/5)
        root=mp.findroot(determinant,(-1.01,-.99)) if p else mp.mpf(-1)
        exact_move=float(root+1)
        tangent=float(pm*mp.exp(mp.mpf("0.1")))
        for count in [256,512]:
            angles=2*np.pi*np.arange(count)/count
            circle=np.exp(1j*angles)
            vals=[];tracevals=[];maxnorm=0.
            for point in -1+.25*circle:
                D0=np.diag([point+1,point+2])
                H=-p*np.exp(-point/10)*B@la.inv(D0)
                maxnorm=max(maxnorm,float(la.norm(H,2)))
                vals.append(np.log1p(la.eigvals(H)).sum())
                tracevals.append(np.trace(H))
            moment=-.25*np.mean(np.array(vals)*circle)
            first=-.25*np.mean(np.array(tracevals)*circle)
            err=abs(moment-exact_move)
            remainder=abs(exact_move-tangent)
            passed=err<1e-11 and abs(first-tangent)<1e-11 and remainder<=float(beta) and maxnorm<=float(q)
            rows.append(dict(parameter=p,quadrature_nodes=count,root=float(root),
                contour_moment_real=moment.real,contour_moment_imag=moment.imag,
                root_contour_discrepancy=err,tangent=tangent,actual_remainder=remainder,
                analytic_uniform_remainder=float(beta),sampled_H_norm=maxnorm,
                analytic_uniform_H_norm=float(q),pass_check=passed))
            assert passed,rows[-1]
    pd.DataFrame(rows).to_csv(OUT/"TABLE_01_TOY_DDE.csv",index=False)
    save("TOY_ANALYTIC_BOUND.json",{
        "status":"PROVED_REDUCED_MODEL",
        "scope":"two-state scalar-parameter DDE fixture only",
        "parameter_interval":["-1/100","1/100"],"contour_center":-1,"contour_radius":"1/4",
        "reference_root_count":1,"small_matrix_dimension":2,
        "uniform_kappa_exact":str(q),"uniform_remainder_exact":str(beta),
        "proof":"On the entire circle min|s+1|=1/4 and min|s+2|>=3/4. exp(-Re(s)/10)<=exp(1/8)<=sum_(k>=0)(1/8)^k=8/7. Matrix norm sqrt(2)<=3/2. Thus kappa<=12/175<1. Baseline has one enclosed root; homotopy preserves its count. The matrix perturbation is affine in p, so H_rem=0. The log-series remainder is bounded by radius*dimension*kappa^2/[2(1-kappa)]=36/28525.",
        "IEEE39_certificate":False})
    return rows

def physical_identity():
    M={name:pd.read_csv(MODEL/(name+".csv")).to_numpy() for name in
       ["Adev","Bv","Cs","Cf","Ds","Y","Etheta","Hv","Bp","Bi"]}
    baseline=tomllib.loads((ROOT/"experiments/graph_gsp_codesign_20261003/baseline.toml").read_text())
    rho=np.array(baseline["rho"]);kp=np.array(baseline["Kp"]);ki=np.array(baseline["Ki"])
    rng=np.random.default_rng(20261003)
    dr=rng.uniform(-.001,.001,10)
    dp=rng.uniform(-1,1,10)*.01*kp
    di=rng.uniform(-1,1,10)*.01*ki
    tau=np.linspace(0,.04,10)
    def closure(rr):
        w=np.repeat(1-rr,2)
        G=M["Y"]+w[:,None]*M["Ds"]
        N=w[:,None]*M["Cs"]+(1-w[:,None])*M["Cf"]
        return G,N,-la.solve(G,N)
    G0,N0,V0=closure(rho)
    G1,N1,Vp=closure(rho+dr)
    C0=M["Etheta"]+M["Hv"]@V0
    B0=M["Bp"]*kp+M["Bi"]*ki
    B1=M["Bp"]*(kp+dp)+M["Bi"]*(ki+di)
    dV=Vp-V0
    J=la.solve(G0,G1-G0)
    V1=-la.solve(G0,(N1-N0)+(G1-G0)@V0)
    VR=-la.solve(np.eye(20)+J,J@V1)
    # Uniform algebraic-block bounds over the preregistered rho box.
    a=v=0.
    for i in range(10):
        Gr=np.zeros_like(G0);Nr=np.zeros_like(N0)
        sl=slice(2*i,2*i+2)
        Gr[sl]=-M["Ds"][sl];Nr[sl]=M["Cf"][sl]-M["Cs"][sl]
        a+=.001*la.norm(la.solve(G0,Gr),2)
        v+=.001*la.norm(la.solve(G0,Nr+Gr@V0),2)
    rows=[]
    for freq in [.5,5.,10.]:
        s=-.05+2j*np.pi*freq
        E=np.diag(np.exp(-s*tau))
        D0=s*np.eye(204)-M["Adev"]-M["Bv"]@V0-B0@E@C0
        C1=M["Etheta"]+M["Hv"]@Vp
        D1=s*np.eye(204)-M["Adev"]-M["Bv"]@Vp-B1@E@C1
        U=np.column_stack([M["Bv"]+B0@E@M["Hv"],M["Bp"]@E,M["Bi"]@E])
        W=-np.row_stack([dV,dp[:,None]*(C0+M["Hv"]@dV),di[:,None]*(C0+M["Hv"]@dV)])
        Wlin=-np.row_stack([V1,dp[:,None]*C0,di[:,None]*C0])
        Wrem=-np.row_stack([VR,dp[:,None]*(M["Hv"]@dV),di[:,None]*(M["Hv"]@dV)])
        resolvent_U=la.solve(D0,U)
        H=W@resolvent_U
        Hrem=Wrem@resolvent_U
        order=np.array([j for i in range(10) for j in [2*i,2*i+1,20+i,30+i]])
        Hb=H[np.ix_(order,order)]
        local=0j;cross=0j
        for i in range(10):
            for j in range(10):
                product=np.trace(Hb[4*i:4*i+4,4*j:4*j+4]@Hb[4*j:4*j+4,4*i:4*i+4])
                if i==j:local+=product
                else:cross+=product
        walk_error=abs(np.trace(H@H)-local-cross)/max(abs(np.trace(H@H)),1e-30)
        errors=[rel(D1-D0,U@W),rel(dV,V1+VR),rel(W-Wlin,Wrem),float(walk_error)]
        kappa_bound=eps_bound=None
        if a<1:
            vb=v/(1-a);vr=a*v/(1-a);hv=la.norm(M["Hv"],2)
            cb=la.norm(C0,2)+hv*vb
            pbar=max(.01*kp);ibar=max(.01*ki)
            wb=np.sqrt(vb**2+(pbar*cb)**2+(ibar*cb)**2)
            rb=np.sqrt(vr**2+(pbar*hv*vb)**2+(ibar*hv*vb)**2)
            kappa_bound=float(wb*la.norm(resolvent_U,2))
            eps_bound=float(40*rb*la.norm(resolvent_U,2))
            assert la.norm(H,2)<=kappa_bound*(1+1e-12)
            assert la.norm(Hrem,"nuc")<=eps_bound*(1+1e-12)
        row=dict(frequency_Hz=freq,characteristic_factor_error=errors[0],
            voltage_remainder_error=errors[1],port_remainder_error=errors[2],
            closed_walk_identity_error=errors[3],point_H_norm=float(la.norm(H,2)),
            algebraic_box_a=float(a),algebraic_box_v=float(v),
            box_H_bound_at_this_s=kappa_bound,box_Hrem_nuclear_bound_at_this_s=eps_bound,
            q_condition_at_this_s=bool(kappa_bound is not None and kappa_bound<1),
            uniform_contour_verified=False,pass_check=max(errors)<1e-9)
        rows.append(row)
        assert row["pass_check"],row
    pd.DataFrame(rows).to_csv(OUT/"TABLE_02_IEEE39_ALGEBRA.csv",index=False)
    save("FROZEN_ALGEBRA_POINT.json",{"seed":20261003,
        "rho":(rho+dr).tolist(),"Kp":(kp+dp).tolist(),"Ki":(ki+di).tolist(),
        "fixed_tau_seconds":tau.tolist(),"scope":"matrix identity point, not a validated design"})
    return rows

def lp_fixture():
    result=linprog([-1.,0.],A_ub=[[1.,1.],[1.,-1.]],b_ub=[1.1,1.1],
                   bounds=[(0.,3.),(-2.,2.)],method="highs")
    assert result.success and abs(result.x[0]-1.1)<1e-10 and abs(result.x[1])<1e-10
    save("LP_FIXTURE.json",{"status":"PROVED_REDUCED_MODEL",
        "scope":"abstract two-constraint relaxation, not MW on IEEE-39",
        "exact_primal_optimum":"11/10","exact_dual_weights":["1/2","1/2"],
        "exact_controller_aggregate":"0","excluded_trial":"6/5",
        "numerical_primal":result.x.tolist(),"pass":True})

def main():
    freeze()
    sy=symbolic();toy=dde_fixture();grid=physical_identity();lp_fixture()
    summary={"status":"THEORY_AND_SMALL_CHECKS_COMPLETE",
        "symbolic_checks":len(sy),"symbolic_passed":sum(r["pass"] for r in sy),
        "toy_DDE_checks":len(toy),"toy_DDE_passed":sum(r["pass_check"] for r in toy),
        "toy_max_moment_error":max(r["root_contour_discrepancy"] for r in toy),
        "toy_max_actual_remainder":max(r["actual_remainder"] for r in toy),
        "toy_uniform_remainder":toy[0]["analytic_uniform_remainder"],
        "IEEE39_algebra_points":len(grid),"IEEE39_algebra_passed":sum(r["pass_check"] for r in grid),
        "IEEE39_max_factor_error":max(r["characteristic_factor_error"] for r in grid),
        "IEEE39_scalar_envelope_max":max(r["box_H_bound_at_this_s"] for r in grid),
        "IEEE39_scalar_envelope_min":min(r["box_H_bound_at_this_s"] for r in grid),
        "IEEE39_scalar_envelope_under_one_points":sum(r["q_condition_at_this_s"] for r in grid),
        "IEEE39_region_certificate":"BLOCKED_UNIFORM_ENCLOSURES_NOT_EVALUATED",
        "IEEE39_replacement_upper_MW":None,"IEEE39_new_feasible_design":False,
        "grid_simulations_run":0,"grid_optimization_runs":0,
        "novelty_status":"NOT_ESTABLISHED"}
    save("SUMMARY.json",summary)
    print(json.dumps(summary,indent=2))

if __name__=="__main__":
    main()
