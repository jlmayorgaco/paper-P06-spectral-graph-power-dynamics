"""Trust-region SQP proposals, with independent full nonlinear acceptance.
The QP uses nonlinear discrete trajectory tangents. Modal gradients define a
positive regularization metric and an equilibrium-margin constraint. Local search
does not certify a global maximum; the exploratory operation envelope is explicit.
"""
import csv,json,socket,subprocess,threading,time,argparse,os,tomllib
from pathlib import Path
import numpy as np
from scipy.optimize import minimize,nnls
ROOT=Path(__file__).resolve().parents[2]
parser=argparse.ArgumentParser();parser.add_argument('--iterations',type=int,default=8)
parser.add_argument('--seed');parser.add_argument('--label',default='search')
parser.add_argument('--full-horizon',action='store_true')
parser.add_argument('--radius',type=float)
parser.add_argument('--corrections',type=int,default=0)
parser.add_argument('--modal-bundle',type=int,default=1)
args=parser.parse_args()
cfg=tomllib.loads(Path(__file__).with_name('problem.toml').read_text())
nc=args.modal_bundle+5*len(cfg['scenarios'])
cfg['modal_bundle_size']=args.modal_bundle
OUT=ROOT/'reports/nonlinear_codesign_20261001'/args.label;OUT.mkdir(exist_ok=True)
weights=np.array([float(r['P_gen_MW']) for r in csv.DictReader((ROOT/'reports/experiment_N/TABLE_N01_original_operating_point.csv').open())])
kp0=10*np.pi;ki0=kp0**2/4
x=np.r_[np.full(10,.8),np.full(10,np.log(kp0)),np.full(10,np.log(ki0))]
if args.seed:x=np.array(json.loads(Path(args.seed).read_text())['x'])
lo=np.r_[np.full(10,.001),np.full(10,np.log(kp0/4)),np.full(10,np.log(ki0/4))]
hi=np.r_[np.full(10,.999),np.full(10,np.log(kp0*4)),np.full(10,np.log(ki0*4))]
cost=np.r_[-weights/weights.sum(),np.zeros(20)]
scales=np.r_[np.ones(10),np.full(20,6.)]
proc=subprocess.Popen(['julia','--startup-file=no','--project=.',str(Path(__file__).with_name('search_oracle.jl'))],
    cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,encoding='utf-8',creationflags=subprocess.CREATE_NO_WINDOW,
    env={**os.environ,'BND_SEARCH_SUBDIRECTORY':args.label,'BND_MODAL_BUNDLE_SIZE':str(args.modal_bundle)})
start=time.monotonic();history=[]
try:
    while True:
        line=proc.stdout.readline()
        if not line:raise RuntimeError('Julia stopped before READY')
        print(line.rstrip(),flush=True)
        if line.startswith('READY '):port=int(line.split()[1]);break
    # Drain progress independently so long oracle calls cannot fill the pipe.
    def drain():
        with (OUT/'oracle.log').open('w',encoding='utf-8') as log:
            for line in proc.stdout:
                log.write(line);log.flush();print(line.rstrip(),flush=True)
    threading.Thread(target=drain,daemon=True).start()
    with socket.create_connection(('127.0.0.1',port),timeout=1200) as sock:
        f=sock.makefile('rw')
        def oracle(p,gradient=False,full_horizon=False):
            command=('GRADFULL' if full_horizon else 'GRAD') if gradient else 'EVAL'
            f.write(command+';'+','.join(format(v,'.17g') for v in p)+'\n');f.flush()
            answer=f.readline().strip()
            if answer.startswith('ERR'):raise RuntimeError(answer)
            v=np.fromstring(answer,sep=',')
            if not np.all(np.isfinite(v)) or len(v)!=(nc*31+1 if gradient else nc):
                raise RuntimeError(f'Invalid oracle length {len(v)}')
            return (v[:nc],v[nc:-1].reshape(nc,30,order='F'),v[-1]) if gradient else v
        g=oracle(x)
        if max(g)>0:raise RuntimeError('Anchor does not satisfy declared search constraints')
        def save(label,**kw):
            data={'x':x.tolist(),'rho':x[:10].tolist(),'Kp':np.exp(x[10:20]).tolist(),'Ki':np.exp(x[20:]).tolist(),
                'retained_MW':float(weights@(1-x[:10])),'replacement_percent':float(100*weights@x[:10]/weights.sum()),
                'constraints':g.tolist(),'status':'BEST_FEASIBLE_FOUND_NOT_CERTIFIED_OPTIMUM',
                'problem_instance':cfg,**kw}
            (OUT/(label+'.json')).write_text(json.dumps(data,indent=2))
        save('anchor');radius=args.radius if args.radius is not None else (.001 if args.full_horizon else .025)
        full_horizon=args.full_horizon
        for iteration in range(1,args.iterations+1):
            cg,A,condition=oracle(x,True,full_horizon)
            # If a relevant 60 s constraint is absent from the 20 s tangent window,
            # do not mislabel that derivative as a full-horizon derivative.
            mismatch=g-cg
            relevant=(g>-.10)&(mismatch>.005)
            if np.any(relevant) and not full_horizon:
                print('ENRICHING_GRADIENT_HORIZON',np.flatnonzero(relevant).tolist(),flush=True)
                full_horizon=True
                cg,A,condition=oracle(x,True,True)
            if full_horizon and np.any((g>-.10)&(abs(g-cg)>.005)):
                save('best',termination='GRADIENT_MESH_REFINEMENT_REQUIRED')
                print('GRADIENT_MESH_REFINEMENT_REQUIRED',flush=True);break
            modal=A[0]*scales
            H=.002*np.eye(30)+.02*np.outer(modal,modal)/(1+modal@modal)
            accepted=False
            for trial in range(6):
                lower=np.maximum((lo-x)/scales,-radius);upper=np.minimum((hi-x)/scales,radius)
                B=A*scales
                # Small guard covers measured SDIRK/adaptive peak discrepancy.
                safety=.002
                q=minimize(lambda z:float((cost*scales)@z+.5*z@H@z),np.zeros(30),
                    jac=lambda z:cost*scales+H@z,method='SLSQP',bounds=list(zip(lower,upper)),
                    constraints=[{'type':'ineq','fun':lambda z:-g-safety-B@z,'jac':lambda z:-B}],
                    options={'maxiter':150,'ftol':1e-11})
                if not q.success:
                    radius*=.5;continue
                dx=scales*q.x
                if cost@dx>=-1e-7:break
                candidate=np.clip(x+dx,lo,hi)
                candidate_g=oracle(candidate)
                corrections=0
                # A bounded feasibility correction can follow curvature of active
                # constraints without discarding an otherwise useful SQP direction.
                # Every correction is evaluated in the complete nonlinear model.
                for correction in range(args.corrections):
                    if max(candidate_g)<=0 or max(candidate_g)>=9.:break
                    clower=np.maximum((lo-candidate)/scales,-radius-(candidate-x)/scales)
                    cupper=np.minimum((hi-candidate)/scales,radius-(candidate-x)/scales)
                    minimum_decrease=.1*min(float(cost@dx),-1e-7)
                    repair=minimize(lambda z:.5*float(z@z),np.zeros(30),jac=lambda z:z,
                        method='SLSQP',bounds=list(zip(clower,cupper)),constraints=[
                            {'type':'ineq','fun':lambda z:-candidate_g-safety-B@z,'jac':lambda z:-B},
                            {'type':'ineq','fun':lambda z:minimum_decrease-cost@(candidate-x+scales*z),
                             'jac':lambda z:-cost*scales}],options={'maxiter':150,'ftol':1e-12})
                    if not repair.success:break
                    candidate=np.clip(candidate+scales*repair.x,lo,hi)
                    candidate_g=oracle(candidate);corrections+=1
                row={'iteration':iteration,'trial':trial,'radius':radius,'retained_MW':float(weights@(1-candidate[:10])),
                    'max_constraint':float(max(candidate_g)),'accepted':bool(max(candidate_g)<=0),
                    'modal_condition':float(condition),'corrections':corrections,'elapsed_s':time.monotonic()-start}
                history.append(row);(OUT/'history.json').write_text(json.dumps(history,indent=2));print('TRIAL',row,flush=True)
                if row['accepted']:
                    x=candidate;g=candidate_g;accepted=True
                    save('best',iteration=iteration);radius=min(.04,radius*1.3);break
                radius*=.5
            if not accepted:
                save('best',termination='NO_ACCEPTED_TRUST_STEP');break
        else:save('best',termination='ITERATION_LIMIT')
        # A fresh gradient is required for a KKT claim. Never reuse a previous-point
        # Jacobian to assert optimality of the newly accepted point.
        f.write('QUIT\n');f.flush()
    proc.wait(timeout=30)
finally:
    if proc.poll() is None:proc.terminate()
