"""Numerical optimizer driving only the frozen Julia physical evaluator.
SLSQP is a local falsification/search tool, never a global certificate.
"""
import csv,json,os,socket,subprocess,time,tomllib,argparse
from pathlib import Path
import numpy as np
from scipy.optimize import minimize,nnls
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'reports/codesign_validation_20261001'
parser=argparse.ArgumentParser();parser.add_argument('--seed');parser.add_argument('--label');parser.add_argument('--prune',action='store_true');parser.add_argument('--iterations',type=int,default=65);parser.add_argument('--dynamic-jacobian',action='store_true')
args=parser.parse_args()
if args.label:OUT=OUT/args.label;OUT.mkdir(exist_ok=True)
d=tomllib.loads((Path(__file__).parent/'frozen/original_candidate.toml').read_text())
weights=np.array([float(r['P_gen_MW']) for r in csv.DictReader((ROOT/'reports/experiment_N/TABLE_N01_original_operating_point.csv').open())])
kp0=10*np.pi;ki0=kp0**2/4
x0=np.r_[d['epsilon'],(np.array(d['Kp'])-.25*kp0)/(3.75*kp0),(np.array(d['Ki'])-.25*ki0)/(3.75*ki0)]
c=np.r_[weights/1000,np.zeros(20)]
lo=np.r_[np.full(10,1e-5),np.zeros(20)];hi=np.r_[np.full(10,1-1e-5),np.ones(20)]
seed=x0.copy();seed[:10]*=1.15
if args.seed:
  seed=np.array(json.loads(Path(args.seed).read_text())['x'])
if args.prune:
  removed=np.flatnonzero(seed[:10]<5e-5);seed[removed]=0;lo[removed]=0;hi[removed]=0
  print('PRUNED_SG_BUSES',(removed+30).tolist(),flush=True)
proc=subprocess.Popen(['julia','--startup-file=no','--project=.',str(Path(__file__).with_name('oracle.jl'))],cwd=ROOT,
 stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,encoding='utf-8',creationflags=subprocess.CREATE_NO_WINDOW)
try:
  while True:
    line=proc.stdout.readline()
    if not line:raise RuntimeError('Julia exited before READY')
    print(line.rstrip(),flush=True)
    if line.startswith('READY '):port=int(line.split()[1]);break
  with socket.create_connection(('127.0.0.1',port),timeout=600) as sock:
    f=sock.makefile('rw');cache={};calls=0;start=time.monotonic();best=None;hist=[]
    def model(x,fine=False):
      global calls,best
      key=(np.asarray(x).tobytes(),fine)
      if key in cache:return cache[key]
      f.write(('FINE' if fine else 'EVAL')+';'+','.join(format(v,'.17g') for v in x)+'\n');f.flush()
      ans=f.readline().strip()
      if ans.startswith('ERR'):raise RuntimeError(ans)
      a=np.fromstring(ans,sep=',')
      if len(a)!=21 or not np.all(np.isfinite(a)):raise RuntimeError('Bad response: '+ans)
      cache[key]=a;calls+=1
      if np.max(a[13:])<=1e-7 and (best is None or a[0]<best['J']):
        best={'J':float(a[0]),'x':x.tolist(),'metrics':a.tolist(),'evaluations':calls}
        (OUT/'best_search.json').write_text(json.dumps(best,indent=2))
      return a
    def gjac(x,h=2e-5):
      out=np.zeros((8,30))
      for j in range(30):
        if hi[j]==lo[j]:continue
        xp=x.copy();xm=x.copy();xp[j]=min(hi[j],x[j]+h);xm[j]=max(lo[j],x[j]-h)
        out[:,j]=(model(xp)[13:]-model(xm)[13:])/(xp[j]-xm[j])
      if args.dynamic_jacobian:
        f.write('DYNJAC;'+','.join(format(v,'.17g') for v in x)+'\n');f.flush()
        ans=f.readline().strip()
        if ans.startswith('ERR'):raise RuntimeError(ans)
        out[:2]=np.fromstring(ans,sep=',').reshape((2,30),order='F')
      return out
    def callback(x):
      a=model(x);row={'iteration':len(hist)+1,'J':float(a[0]),'maxg':float(max(a[13:])),
        'alpha':float(a[1]),'beta_sampled':float(a[2]),'F':a[4:7].tolist(),'evaluations':calls,'elapsed_s':time.monotonic()-start}
      hist.append(row);(OUT/'optimization_history.json').write_text(json.dumps(hist,indent=2))
      print('ITER',row,flush=True)
    result=minimize(lambda x:float(c@x),seed,jac=lambda x:c,bounds=list(zip(lo,hi)),
      constraints=[{'type':'ineq','fun':lambda x:-model(x)[13:],'jac':lambda x:-gjac(x)}],
      method='SLSQP',callback=callback,options={'maxiter':args.iterations,'ftol':2e-9,'disp':True})
    candidates=[result.x]
    if best is not None:candidates.append(np.array(best['x']))
    audited=[]
    for ix,x in enumerate(candidates):
      a=model(x,True);J=gjac(x,h=5e-6);active=np.where(a[13:]>=-2e-4)[0]
      rows=[J[k] for k in active];names=['g'+str(k) for k in active]
      for j in range(30):
        if x[j]-lo[j]<1e-6:rows.append(-np.eye(30)[j]);names.append('lower'+str(j))
        if hi[j]-x[j]<1e-6:rows.append(np.eye(30)[j]);names.append('upper'+str(j))
      M=np.array(rows).T
      mu=nnls(M,-c,maxiter=10000)[0] if len(rows) else np.zeros(0)
      residual=float(np.max(np.abs(c+M@mu))) if len(rows) else float(max(abs(c)))
      item={'x':x.tolist(),'metrics_fine':a.tolist(),'J':float(a[0]),'stationarity_inf':residual,
        'active':names,'multipliers':mu.tolist(),'max_constraint_fine':float(max(a[13:])),
        'solver_success':bool(result.success),'solver_message':str(result.message),
        'global_optimum_certified':False,'SOSC_tested':False}
      audited.append(item)
      (OUT/f'candidate_search_{ix}.json').write_text(json.dumps(item,indent=2))
      print('AUDIT',item,flush=True)
    (OUT/'optimization_results.json').write_text(json.dumps({'audits':audited,'elapsed_s':time.monotonic()-start,'evaluations':calls},indent=2))
    f.write('QUIT\n');f.flush()
finally:
  if proc.poll() is None:
    try:proc.wait(timeout=20)
    except subprocess.TimeoutExpired:proc.terminate();proc.wait(timeout=20)
