"""Explore finite extra replacement with fixed rho and all twenty PLL gains.

No failed search is labeled an impossibility certificate. All rejected proposals
are retained by the Julia oracle. Only full nonlinear evaluations accept steps.
"""
import csv,json,os,socket,subprocess,threading,tomllib,time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize,linprog

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"reports/retuning_boundary_20261002"
OUT.mkdir(exist_ok=True)
seed=tomllib.loads((ROOT/"reports/analytic_iteration_20261001/refined/joint_final.toml").read_text())
p0=np.r_[seed["rho"],np.log(seed["Kp"]),np.log(seed["Ki"])]
power=np.array([float(r["P_gen_MW"]) for r in csv.DictReader(
    (ROOT/"reports/experiment_N/TABLE_N01_original_operating_point.csv").open())])
kp=2*np.pi*5;ki=kp*kp/4
lo=np.r_[np.full(10,np.log(kp/4)),np.full(10,np.log(ki/4))]
hi=np.r_[np.full(10,np.log(kp*4)),np.full(10,np.log(ki*4))]
protocol={"seed":str(ROOT/"reports/analytic_iteration_20261001/refined/joint_final.toml"),
          "extra_MW_grid":[2,5,10,20],"rho_direction":"uniform increment at all ten generator buses",
          "all_20_PLL_gains_free":True,"iterations_per_target":5,"primary_event":[8,100],
          "limits":{"F_Hz":.5,"R_Hz_s":.5,"Vmin":.9,"Vmax":1.1,"alpha":-.05,"actuator_slack":.002},
          "claim":"exploratory finite-event feasibility; no global or continuous-disturbance guarantee"}
(OUT/"search_protocol.json").write_text(json.dumps(protocol,indent=2))
proc=subprocess.Popen(["julia","--startup-file=no","--project=.",
    "experiments/retuning_boundary_20261002/oracle.jl"],cwd=ROOT,stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,text=True,bufsize=1,
    creationflags=getattr(subprocess,"CREATE_NO_WINDOW",0))
log=(OUT/"oracle.log").open("w",encoding="utf8")
port=None
while True:
    line=proc.stdout.readline()
    if not line: raise RuntimeError("Julia exited before READY")
    log.write(line);log.flush();print(line,end="",flush=True)
    if line.startswith("READY "):port=int(line.split()[1]);break
def drain():
    for line in proc.stdout:
        log.write(line);log.flush();print(line,end="",flush=True)
thread=threading.Thread(target=drain,daemon=True);thread.start()
sock=socket.create_connection(("127.0.0.1",port));stream=sock.makefile("rw")
calls=0;history=[]
def ask(p,grad=False,horizon=60,bus=8,delta=100):
    global calls
    calls+=1
    stream.write(f'{"GRAD" if grad else "EVAL"};{horizon};{bus};{delta};'+",".join(map(repr,map(float,p)))+"\n")
    stream.flush();line=stream.readline().strip()
    if line.startswith("ERR"):raise RuntimeError(line)
    out=np.fromstring(line,sep=",")
    if grad:return out[:6],out[6:].reshape((6,30),order="F")
    return out
def save_candidate(label,p,g,extra):
    lines=[f'rho = {p[:10].tolist()}',f'Kp = {np.exp(p[10:20]).tolist()}',
           f'Ki = {np.exp(p[20:]).tolist()}','dc_convention = "physical_supply"',
           f'constraints = {g.tolist()}',f'extra_MW = {extra}',
           f'replacement_percent = {100*np.dot(power,p[:10])/power.sum()}',
           'status = "FINITE_EVENT_EVALUATION_NOT_ROBUST_CERTIFICATE"']
    (OUT/(label+".toml")).write_text("\n".join(lines)+"\n")
try:
    g0=ask(p0);save_candidate("anchor",p0,g0,0)
    for extra in protocol["extra_MW_grid"]:
        p=p0.copy();p[:10]+=extra/power.sum()
        if np.any(p[:10]>=.99):continue
        g=ask(p);fixed=g.copy();save_candidate(f"extra_{extra}_fixed",p,g,extra)
        start_g=g.copy()
        print("TARGET",extra,"fixed_maxg",g.max(),flush=True)
        if g.max()<=-1e-4:
            history.append({"extra_MW":extra,"fixed_maxg":float(g.max()),"tuned_maxg":float(g.max()),
                            "iterations":0,"recovered":False,"already_feasible":True})
            continue
        trust=.2;iterations=0
        for iteration in range(1,protocol["iterations_per_target"]+1):
            iterations=iteration
            tg,J=ask(p,grad=True,horizon=20)
            if np.any((g>-.08)&(np.abs(g-tg)>.002)):
                tg,J=ask(p,grad=True,horizon=60)
            if np.any((g>-.08)&(np.abs(g-tg)>.003)):
                print("STOP_GRADIENT_HORIZON_MISMATCH",extra,flush=True);break
            A=J[:,10:]
            lower=np.maximum(lo-p[10:],-trust);upper=np.minimum(hi-p[10:],trust)
            # Minimax linear proposal; the slack is not interpreted as a finite certificate.
            lp=linprog(np.r_[np.zeros(20),1.],A_ub=np.c_[A,-np.ones(6)],
                       b_ub=-g-.0002,bounds=list(zip(lower,upper))+[(-.01,None)],method="highs")
            if not lp.success:break
            proposal=lp.x[:20]
            predicted=g+A@proposal
            # Minimum-norm correction when the linearized system is feasible.
            if predicted.max()<=-.0001:
                opt=minimize(lambda z:.5*np.dot(z,z),proposal,jac=lambda z:z,
                    bounds=list(zip(lower,upper)),
                    constraints={"type":"ineq","fun":lambda z:-g-.0001-A@z,"jac":lambda z:-A},
                    method="SLSQP",options={"maxiter":150,"ftol":1e-12})
                if opt.success:proposal=opt.x
            accepted=False
            for bt in range(4):
                trial=p.copy();trial[10:]+=proposal*2**(-bt)
                gt=ask(trial)
                if gt.max()<g.max()-1e-7:
                    p=trial;g=gt;accepted=True;break
            if not accepted:
                trust*=.5
                if trust<.025:break
            if g.max()<=-1e-5:break
        save_candidate(f"extra_{extra}_tuned",p,g,extra)
        history.append({"extra_MW":extra,"fixed_maxg":float(fixed.max()),"tuned_maxg":float(g.max()),
                        "iterations":iterations,"recovered":bool(fixed.max()>0 and g.max()<=0),
                        "already_feasible":False,"initial_constraints":start_g.tolist(),
                        "final_constraints":g.tolist()})
        (OUT/"search_result.json").write_text(json.dumps(history,indent=2))
    (OUT/"search_result.json").write_text(json.dumps(history,indent=2))
finally:
    try:stream.write("QUIT\n");stream.flush()
    finally:stream.close();sock.close()
    proc.wait(timeout=30);thread.join(timeout=5);log.close()
print(json.dumps(history,indent=2),flush=True)
