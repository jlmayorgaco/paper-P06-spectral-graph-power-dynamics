"""Read-only frozen model; full exact exponential characteristic, no Pade."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
from pathlib import Path
import json,hashlib,tomllib,sys,platform,subprocess
import numpy as np
import pandas as pd
from scipy import linalg as la

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
SOURCE=ROOT/'experiments/graph_gsp_codesign_20261003'
NAMES=['Adev','Bv','Cs','Cf','Ds','Y','Etheta','Hv','Bp','Bi']
M={k:pd.read_csv(SOURCE/'model'/f'{k}.csv').to_numpy() for k in NAMES}
BASE=tomllib.loads((SOURCE/'baseline.toml').read_text())
P0=np.r_[BASE['rho'],BASE['Kp'],BASE['Ki']]
PORTS=pd.read_csv(SOURCE/'model/ports.csv')
CAT=pd.read_csv(SOURCE/'TABLE_02_BASELINE_ROOTS.csv')
ROOTS=CAT.real.to_numpy()+1j*CAT.imag.to_numpy()
TAU=np.full(10,.04)
U=np.block([[np.zeros((204,20)),M['Bp'],M['Bi']],
            [np.eye(20),np.zeros((20,20))],[np.zeros((10,40))]])
GROUPS=[np.array([2*i,2*i+1,20+i,30+i]) for i in range(10)]

def save(name,obj): (OUT/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def freeze():
    paths=[SOURCE/'baseline.toml',SOURCE/'TABLE_02_BASELINE_ROOTS.csv',OUT/'PROTOCOL.md']
    paths+=sorted((SOURCE/'model').glob('*.csv'))
    head=subprocess.run(['git','rev-parse','HEAD'],cwd=ROOT,text=True,capture_output=True)
    dirty=subprocess.run(['git','status','--porcelain'],cwd=ROOT,text=True,capture_output=True)
    manifest={'inputs':{str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
              'python':platform.python_version(),'numpy':np.__version__,
              'head':head.stdout.strip() if head.returncode==0 else None,
              'dirty':dirty.stdout if dirty.returncode==0 else None}
    target=OUT/'INPUT_MANIFEST.json'
    if target.exists(): assert json.loads(target.read_text())['inputs']==manifest['inputs']
    else: save('INPUT_MANIFEST.json',manifest)

class Model:
    def __init__(self,p=P0,tau=TAU):
        self.p=np.array(p);self.tau=np.broadcast_to(tau,(10,)).copy()
        w=np.repeat(1-self.p[:10],2)
        self.G=M['Y']+w[:,None]*M['Ds']
        self.N=w[:,None]*M['Cs']+(1-w[:,None])*M['Cf']
        self.Vx=-la.solve(self.G,self.N,check_finite=False)
        self.A0=M['Adev']+M['Bv']@self.Vx
        self.C=M['Etheta']+M['Hv']@self.Vx
        self.B=M['Bp']*p[10:20]+M['Bi']*p[20:30]
        self.A=self.A0+self.B@self.C;self.eye=np.eye(204)
    def delta(self,s):return s*self.eye-self.A0-(self.B*np.exp(-s*self.tau))@self.C
    def derivative(self,s):return self.eye+(self.B*(self.tau*np.exp(-s*self.tau)))@self.C
    def refine(self,s):
        s=complex(s)
        for it in range(35):
            D=s*self.eye-self.A0;lu=la.lu_factor(D,check_finite=False)
            rb=la.lu_solve(lu,self.B,check_finite=False)
            e=np.exp(-s*self.tau)
            K=np.eye(10)-(self.C@rb)*e[None,:]
            left,sv,vh=la.svd(K,check_finite=False);v=vh[-1].conj();u=left[:,-1]
            Ks=(self.C@la.lu_solve(lu,rb,check_finite=False))*e[None,:]+(self.C@rb)*(self.tau*e)[None,:]
            step=-np.vdot(u,K@v)/np.vdot(u,Ks@v)
            if abs(step)<1e-10*max(1,abs(s)):
                fullv=rb@(e*v);res=la.norm(self.delta(s)@fullv)/(la.norm(self.delta(s))*la.norm(fullv))
                if res<1e-10:return s,float(res)
            if abs(step)>1:step/=abs(step)
            s+=step
        raise RuntimeError(f'root failed seed={s} step={step}')
    def descriptor(self,s):
        E=np.diag(np.exp(-s*self.tau))
        return np.block([[s*np.eye(204)-M['Adev'],-M['Bv'],-self.B@E],
                 [self.N,self.G,np.zeros((20,10))],[-M['Etheta'],-M['Hv'],np.eye(10)]])
    def T(self,s):
        E=np.diag(np.exp(-s*self.tau))
        V=np.block([[-(M['Cs']-M['Cf']),-M['Ds'],np.zeros((20,10))],
                    [np.zeros((10,224)),-E],[np.zeros((10,224)),-E]])
        return V@la.solve(self.descriptor(s),U,check_finite=False)

def action(d):return np.r_[np.repeat(d[:10],2),d[10:]]
def singleton(p,i):
    out=P0.copy();ix=[i,10+i,20+i];out[ix]=p[ix];return out
