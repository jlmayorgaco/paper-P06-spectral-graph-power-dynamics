"""One preregistered all-PLL gain-map validation. Exact delays; no optimization."""
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
import sys
sys.dont_write_bytecode = True
import hashlib, json, platform, shutil, time, tomllib
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import linalg as la

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent.parent
OLD = ROOT/"experiments/graph_gsp_codesign_20261003"
sys.path.insert(0, str(OLD))
from graph_design import Model

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save_json(name, obj):
    (OUT/name).write_text(json.dumps(obj, indent=2, allow_nan=False)+"\n", encoding="utf-8")
def frame(name, rows): pd.DataFrame(rows).to_csv(OUT/name, index=False)
def pack(d): return np.r_[d["rho"], np.log(d["Kp"]), np.log(d["Ki"])]
def design(name, rho, kp, ki):
    p=OUT/"designs"/(name+".toml");p.parent.mkdir(exist_ok=True)
    p.write_text("\n".join(k+" = ["+", ".join(format(float(z),".17g") for z in v)+"]"
                         for k,v in [("rho",rho),("Kp",kp),("Ki",ki)])+"\ntau = 0.04\n",encoding="utf-8")

def pair(m,seed):
    """Balanced full-characteristic Newton; avoid ill-scaled small-return residual."""
    _,(scale,_)=la.matrix_balance(m.A,permute=False,separate=True)
    s=complex(seed)
    for iteration in range(20):
        D=m.delta(s,.04);Db=(D*scale[None,:])/scale[:,None]
        U,S,Vh=la.svd(Db,check_finite=False);v=Vh[-1].conj();w=U[:,-1]
        Ds=m.I+.04*np.exp(-s*.04)*(m.B@m.C)
        Dsb=(Ds*scale[None,:])/scale[:,None]
        step=-np.vdot(w,Db@v)/np.vdot(w,Dsb@v)
        vo=scale*v;wo=w/scale
        residual=la.norm(D@vo)/(la.norm(D)*la.norm(vo))
        if abs(step)<=1e-10*max(1,abs(s)) and residual<=1e-9:
            vo/=max(abs(vo[m.ports.pll_angle_index.to_numpy(int)]))
            return s,vo,wo,float(residual)
        if abs(step)>1:step/=abs(step)
        s+=step
    raise RuntimeError(f"Balanced full Newton failed: {s}, residual={residual}, step={step}")

class Return:
    def __init__(self,m,s):
        self.m=m;self.s=s
        self.theta=m.ports.pll_angle_index.to_numpy(int)
        omega=m.ports.pll_frequency_index.to_numpy(int)
        xi=np.argmax(abs(m.Bi),axis=0)
        removed=np.r_[self.theta,omega,xi]
        self.h=np.setdiff1d(np.arange(len(m.A)),removed)
        assert la.norm(m.A0[np.ix_(self.h,np.r_[omega,xi])])<1e-12
        self.tf=1/m.Bp[omega,np.arange(10)]
        self.T=s*np.eye(len(self.h))-m.A0[np.ix_(self.h,self.h)]
        self.lu=la.lu_factor(self.T)
        self.X=la.lu_solve(self.lu,m.A0[np.ix_(self.h,self.theta)])
        self.G=m.C[:,self.theta]+m.C[:,self.h]@self.X
    def gain(self,q,tau=.04):
        y=self.G@q
        W=self.s**2*(1+self.tf*self.s)*np.exp(self.s*tau)*q/y
        kp=W.imag/self.s.imag
        ki=W.real-self.s.real*kp
        return kp,ki,W,y
    def gain_derivative(self,q,bus_index):
        m=self.m;i=bus_index;sl=slice(2*i,2*i+2)
        Gr=np.zeros_like(m.G);Cr=np.zeros_like(m.Cs)
        Gr[sl]=-m.Ds[sl];Cr[sl]=m.Cf[sl]-m.Cs[sl]
        vr=-la.solve(m.G,Gr@m.Vx+Cr)
        Ap=m.Bv@vr;Cp=m.Hv@vr
        Xp=la.lu_solve(self.lu,Ap[np.ix_(self.h,self.theta)]+Ap[np.ix_(self.h,self.h)]@self.X)
        Gp=Cp[:,self.theta]+Cp[:,self.h]@self.X+m.C[:,self.h]@Xp
        kp,ki,W,y=self.gain(q)
        Wp=-W*(Gp@q)/y
        return np.r_[Wp.imag/self.s.imag,Wp.real-self.s.real*Wp.imag/self.s.imag]

def main():
    started=time.perf_counter()
    frozen=[OLD/"baseline.toml",OLD/"TABLE_02_BASELINE_ROOTS.csv",OLD/"graph_design.py",
            OLD/"DelayedEvents.jl",OLD/"eval/baseline/events.csv",OLD/"spectral/baseline_contour.csv",
            ROOT/"experiments/physical_collective_damping_20261003/analyze_ports.py",
            ROOT/"experiments/nonlinear_codesign_20261001/ReducedDAE.jl",
            ROOT/"src/bnd_model_expN/PDExactDesignN.jl",
            ROOT/"experiments/theory_collective_damping_20261003/THEORY.tex",
            ROOT/"experiments/analytical_delay_codesign_mega_20261002/m3_a_trace_integral.jl"]
    frozen+=sorted((OLD/"model").glob("*.csv"))
    manifest={"protocol_sha256":sha(OUT/"PROTOCOL.txt"),"python":platform.python_version(),
              "numpy":np.__version__,"inputs":{str(p.relative_to(ROOT)):sha(p) for p in frozen},
              "scope":"one preregistered finite replacement, no optimization"}
    save_json("INPUT_MANIFEST.json",manifest)
    d=tomllib.loads((OLD/"baseline.toml").read_text())
    p=pack(d);m=Model(p);catalog=pd.read_csv(OLD/"TABLE_02_BASELINE_ROOTS.csv")
    modes={};parity=[]
    for row in catalog.itertuples():
        if row.imag<=0:continue
        s,v,_,res=pair(m,complex(row.real,row.imag))
        th=m.ports.pll_angle_index.to_numpy(int);q=v[th];q=q/q[np.argmax(abs(q))]
        ret=Return(m,s);kp,ki,W,y=ret.gain(q)
        err=np.max(abs(np.r_[kp/m.kp-1,ki/m.ki-1]))
        modes[row.root]=(s,q,ret)
        parity.append(dict(root=row.root,real=s.real,imag=s.imag,gain_max_relative_error=err,
                           full_residual=res,hidden_condition=np.linalg.cond(ret.T),
                           min_relative_pattern=min(abs(q))/max(abs(q)),
                           min_relative_detector=min(abs(y))/max(abs(y)),
                           gain_pass=bool(err<=1e-5 and res<=1e-9)))
    frame("TABLE_01_GAIN_RECONSTRUCTION.csv",parity)
    if not all(r["gain_pass"] for r in parity):
        save_json("STATUS.json",{"status":"BLOCKED_GAIN_RECONSTRUCTION"});return
    fast=catalog[(catalog.frequency_hz>=2)&(catalog.frequency_hz<=12)].sort_values("real",ascending=False)
    target=int(fast.iloc[0].root);second=int(fast.iloc[1].root)
    slow=int(catalog[catalog.frequency_hz<2].sort_values("real",ascending=False).iloc[0].root)
    checks=[]
    for j in [slow,target,second]:
        s,q,r=modes[j]
        for i in [0,4,8]:
            h=1e-5;pp=p.copy();pm=p.copy();pp[i]+=h;pm[i]-=h
            fp=np.r_[Return(Model(pp),s).gain(q)[:2]].reshape(-1)
            fm=np.r_[Return(Model(pm),s).gain(q)[:2]].reshape(-1)
            fd=(fp-fm)/(2*h);an=r.gain_derivative(q,i)
            absolute=la.norm(an-fd,np.inf);relative=absolute/max(la.norm(an,np.inf),1e-12)
            checks.append(dict(root=j,bus=i+30,absolute_inf_error=absolute,relative_inf_error=relative,
                               derivative_pass=bool(relative<=.01 or absolute<=1e-6)))
    frame("TABLE_02_GAIN_DERIVATIVES.csv",checks)
    if not all(r["derivative_pass"] for r in checks):
        save_json("STATUS.json",{"status":"BLOCKED_GAIN_DERIVATIVES"});return
    s,q,ret=modes[target];rhonew=p[:10]+.001
    pfixed=p.copy();pfixed[:10]=rhonew
    rf=Return(Model(pfixed),s);kp,ki,W,y=rf.gain(q)
    gain_j=np.column_stack([ret.gain_derivative(q,i) for i in range(10)])
    dg=gain_j@np.full(10,.001)
    pnew=np.r_[rhonew,np.log(kp),np.log(ki)]
    kpnom=2*np.pi*5;kinom=kpnom**2/4
    gain_pass=bool(np.all((kp>=.25*kpnom)&(kp<=4*kpnom)) and np.all((ki>=.25*kinom)&(ki<=4*kinom)))
    predictions=[]
    dp_tangent=np.r_[np.full(10,.001),dg[:10]/m.kp,dg[10:]/m.ki]
    for j,(sj,qj,rj) in modes.items():
        grad,_=m.gradients(sj,.04)
        for name,delta in [("analytic",dp_tangent),("fixed",np.r_[np.full(10,.001),np.zeros(20)])]:
            pred=sj+grad@delta
            predictions.append(dict(design=name,root=j,baseline_real=sj.real,baseline_imag=sj.imag,
                                    predicted_real=pred.real,predicted_imag=pred.imag))
    frame("PREDICTED_ROOTS.csv",predictions)
    frame("TABLE_03_GAIN_PREDICTION.csv",[dict(bus=i+30,rho_before=p[i],rho_after=rhonew[i],
        Kp_before=m.kp[i],Ki_before=m.ki[i],Kp_tangent=m.kp[i]+dg[i],Ki_tangent=m.ki[i]+dg[i+10],
        Kp_exact=kp[i],Ki_exact=ki[i]) for i in range(10)])
    frame("TARGET_PATTERN.csv",[dict(bus=i+30,real=z.real,imag=z.imag) for i,z in enumerate(q)])
    design("baseline",p[:10],m.kp,m.ki);design("analytic",rhonew,kp,ki);design("fixed",rhonew,m.kp,m.ki)
    # Predictions are sealed before full trial roots, Julia AD or trajectories.
    save_json("PREDICTION_LOCK.json",{"target_root":target,"lambda_real":s.real,"lambda_imag":s.imag,
        "P_total_MW":float(m.ports.P0.sum()),"baseline_GFL_MW":float(m.ports.P0@p[:10]),
        "trial_GFL_MW":float(m.ports.P0@rhonew),"gain_bounds_pass":gain_pass,
        "trial_hashes":{str(x.relative_to(OUT)):sha(x) for x in [OUT/"designs/analytic.toml",
        OUT/"designs/fixed.toml",OUT/"PREDICTED_ROOTS.csv",OUT/"TARGET_PATTERN.csv"]}})
    if not gain_pass:
        save_json("STATUS.json",{"status":"BLOCKED_GAIN_BOUNDS"});return
    rows=[]
    for pred in predictions:
        mdl=Model(pnew if pred["design"]=="analytic" else pfixed)
        oldroot=modes[pred["root"]][0]
        root,v,w,res=pair(mdl,complex(pred["predicted_real"],pred["predicted_imag"]))
        qn=v[ret.theta]
        mac=abs(np.vdot(q,qn))**2/(np.vdot(q,q).real*np.vdot(qn,qn).real)
        predicted=complex(pred["predicted_real"],pred["predicted_imag"])
        rows.append(dict(**pred,actual_real=root.real,actual_imag=root.imag,
            first_order_error=abs(root-predicted),actual_movement=abs(root-oldroot),
            full_residual=res,PLL_pattern_MAC=mac,
            target_preservation_pass=bool(abs(root-s)<=1e-5 and mac>=.999999) if pred["root"]==target else None))
    frame("TABLE_04_ROOT_PREDICTION.csv",rows)
    # Two prescribed patterns; this is not an impossibility test for stabilization.
    C=np.array([[z for z in row] for j in [target,second]
                for row in [[modes[j][0].real,1.],[modes[j][0].imag,0.]]])
    projection=np.eye(4)-C@la.pinv(C)
    compat=[]
    for name,pm in [("baseline",p),("replacement",pfixed)]:
        Wrs=[Return(Model(pm),modes[j][0]).gain(modes[j][1])[2] for j in [target,second]]
        for i in range(10):
            ww=np.array([v for wr in Wrs for v in [wr[i].real,wr[i].imag]])
            conflict=projection@ww
            compat.append(dict(design=name,bus=i+30,absolute_residual=la.norm(conflict),
                               relative_residual=la.norm(conflict)/max(la.norm(ww),1e-30)))
    frame("TABLE_05_FIXED_PATTERN_COMPATIBILITY.csv",compat)
    targetrow=next(x for x in rows if x["design"]=="analytic" and x["root"]==target)
    save_json("STATUS.json",{"status":"PREDICTOR_VALIDATED_PENDING_INDEPENDENT_FULL_MODEL",
        "gain_reconstruction_pass":True,"derivatives_pass":True,
        "target_pass":targetrow["target_preservation_pass"],"seconds":time.perf_counter()-started,
        "maximum_gain_reconstruction_error":max(x["gain_max_relative_error"] for x in parity),
        "maximum_derivative_error":max(x["relative_inf_error"] for x in checks),
        "nonlinear_validation":"NOT_RUN","optimality":"NOT_OPTIMIZED"})
    print(json.dumps(json.loads((OUT/"STATUS.json").read_text())),flush=True)

if __name__=="__main__":main()
