"""IEEE/WSCC 9-bus electrical data with an explicitly reduced SG/GFL DDAE.

Electrical numbers: MATPOWER case9 (accessed 2026-10-05).
Dynamic parameters are declared experimental assumptions, NOT standard case9 data.
State blocks (3 per block): delta, nu, pm, theta, wpll, xi, id, iq.
nu is per-unit SG speed; PLL angles are in the nominal synchronous frame.
All powers/currents use a 100-MVA system base.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from functools import lru_cache
from math import factorial
import numpy as np
from scipy.linalg import eig, eigvals, matrix_balance, svd, solve
from scipy.optimize import root
from scipy.signal import tf2ss

BASE_MVA = 100.0
F0 = 60.0
WB = 2*np.pi*F0
# f, t, r, x, line charging, MVA rating. Bus indices here are one-based.
BRANCH = np.array([
 [1,4,0,.0576,0,250], [4,5,.017,.092,.158,250],
 [5,6,.039,.17,.358,150], [3,6,0,.0586,0,300],
 [6,7,.0119,.1008,.209,150], [7,8,.0085,.072,.149,250],
 [8,2,0,.0625,0,250], [8,9,.032,.161,.306,250],
 [9,4,.01,.085,.176,250]], dtype=float)
LOAD = np.zeros(9,complex); LOAD[[4,6,8]] = np.array([90+30j,100+35j,125+50j])/BASE_MVA
VG = np.array([1.04,1.025,1.025])
PG = np.array([np.nan,1.63,.85])


def ybus(branch_scale: float = 1.0) -> np.ndarray:
    """Series R/X and shunt charging; no ideal slack retained dynamically."""
    Y=np.zeros((9,9),complex)
    for f,t,r,x,b,_ in BRANCH:
        a,bus=int(f)-1,int(t)-1
        y=1/(branch_scale*(r+1j*x))
        Y[a,a]+=y+1j*b/2;Y[bus,bus]+=y+1j*b/2
        Y[a,bus]-=y;Y[bus,a]-=y
    return Y


def powerflow(load_scale: float=1.0, branch_scale: float=1.0):
    """PV buses 2,3; reference angle at 1; voltage setpoints from case9 gen."""
    Y=ybus(branch_scale); L=LOAD*load_scale
    def unpack(q):
        phi=np.r_[0.0,q[:8]];v=np.r_[VG,q[8:]]
        return v*np.exp(1j*phi)
    pspec=-L.real.copy();pspec[1:3]=PG[1:3]
    qspec=-L.imag
    def eq(q):
        v=unpack(q);s=v*np.conj(Y@v)
        return np.r_[s.real[1:]-pspec[1:],s.imag[3:]-qspec[3:]]
    sol=root(eq,np.r_[np.zeros(8),np.ones(6)],tol=1e-11)
    # Tight residual cleanup, independent of scipy.root's step-based exit flag.
    for _ in range(4):
        residual=eq(sol.x)
        if np.max(np.abs(residual))<2e-14:break
        hh=1e-5;J=np.column_stack([(eq(sol.x+hh*np.eye(14)[j])-eq(sol.x-hh*np.eye(14)[j]))/(2*hh) for j in range(14)])
        sol.x-=np.linalg.solve(J,residual)
    res=np.max(np.abs(eq(sol.x)))
    if res>1e-9: raise RuntimeError(f'powerflow failed: {sol.message}, {res:g}')
    V=unpack(sol.x);S=V*np.conj(Y@V)+L
    return V,S[:3],Y,L,res

@dataclass
class Grid9:
    rho: np.ndarray = field(default_factory=lambda:np.full(3,.75))
    # These ratings/H/xd are assumptions on explicit bases, not case9 dynamics.
    H: np.ndarray = field(default_factory=lambda:np.array([4.,3.,2.5]))
    rating: np.ndarray = field(default_factory=lambda:np.array([250.,200.,150.]))
    xd: np.ndarray = field(default_factory=lambda:np.array([.12,.15,.18]))
    d0: np.ndarray = field(default_factory=lambda:np.array([2.,2.,2.]))
    tg: np.ndarray = field(default_factory=lambda:np.array([.45,.55,.65]))
    droop: np.ndarray = field(default_factory=lambda:np.full(3,.05))
    tc: np.ndarray = field(default_factory=lambda:np.full(3,.01))
    tf: np.ndarray = field(default_factory=lambda:np.full(3,.003))
    load_scale: float = 1.0
    branch_scale: float = 1.0

    def __post_init__(self):
        self.rho=np.broadcast_to(self.rho,(3,)).astype(float).copy()
        if np.any(self.rho<=0) or np.any(self.rho>=1):
            raise ValueError('Fixed support requires 0 < rho_i < 1.')
        self.V0,self.S0,self.Y,self.load,self.pf_res=powerflow(self.load_scale,self.branch_scale)
        I0=np.conj(self.S0/self.V0[:3])
        E0=self.V0[:3]+1j*self.xd*I0
        self.E=np.abs(E0)
        self.ys=(1-self.rho)/(1j*self.xd)
        self.M=2*self.H*(self.rating/BASE_MVA)*(1-self.rho)
        self.D=self.d0*(1-self.rho)
        self.Kg=(self.rating/BASE_MVA)*(1-self.rho)/self.droop
        self.pc=self.rho*self.S0.real;self.qc=self.rho*self.S0.imag
        self.pm0=(1-self.rho)*self.S0.real
        self.Yload=np.conj(self.load)/np.abs(self.V0)**2
        Ytot=self.Y+np.diag(self.Yload)
        Ytot[np.arange(3),np.arange(3)]+=self.ys
        self.Ytot=Ytot;self.Yinv=np.linalg.inv(Ytot)
        self.Kv=self.Yinv[:,:3]
        self.x0=np.zeros(24)
        self.x0[:3]=np.angle(E0);self.x0[6:9]=self.pm0
        self.x0[9:12]=np.angle(self.V0[:3])
        self.x0[18:21]=self.pc/np.abs(self.V0[:3])
        self.x0[21:24]=-self.qc/np.abs(self.V0[:3])
        self.A0,self.C,self.Jv=self.linearize_analytic()
        self.gauge=np.zeros(24);self.gauge[:3]=1;self.gauge[9:12]=1

    def volts(self,x,event_bus=-1,event_y=0j):
        E=self.E*np.exp(1j*x[:3])
        I=(x[18:21]+1j*x[21:24])*np.exp(1j*x[9:12])
        V=self.Kv@(self.ys*E+I)
        if event_bus>=0 and event_y!=0:
            V-=self.Yinv[:,event_bus]*event_y*V[event_bus]/(1+event_y*self.Yinv[event_bus,event_bus])
        return V

    def detector(self,x,event_bus=-1,event_y=0j):
        return np.imag(self.volts(x,event_bus,event_y)[:3]*np.exp(-1j*x[9:12]))

    def rhs(self,x,delayed_error,kp,ki,event_bus=-1,event_y=0j):
        V=self.volts(x,event_bus,event_y)
        Em=self.E*np.exp(1j*x[:3]);Is=self.ys*(Em-V[:3])
        pe=(Em*np.conj(Is)).real
        vd=(V[:3]*np.exp(-1j*x[9:12])).real
        if np.min(vd)<.15: raise FloatingPointError('PLL-frame voltage outside reduced model domain')
        f=np.zeros(24)
        f[:3]=WB*x[3:6]
        f[3:6]=(x[6:9]-pe-self.D*x[3:6])/self.M
        f[6:9]=(self.pm0-x[6:9]-self.Kg*x[3:6])/self.tg
        f[9:12]=x[12:15]
        f[12:15]=(-x[12:15]+x[15:18]+kp*delayed_error)/self.tf
        f[15:18]=ki*delayed_error
        f[18:21]=(-x[18:21]+self.pc/vd)/self.tc
        f[21:24]=(-x[21:24]-self.qc/vd)/self.tc
        return f

    def linearize_analytic(self):
        x=self.x0;Em=self.E*np.exp(1j*x[:3]);rot=np.exp(1j*x[9:12])
        I=(x[18:21]+1j*x[21:24])*rot
        dV=np.zeros((9,24),complex)
        dV[:,:3]=self.Kv*np.expand_dims(self.ys*1j*Em,0)
        dV[:,9:12]=self.Kv*np.expand_dims(1j*I,0)
        dV[:,18:21]=self.Kv*np.expand_dims(rot,0)
        dV[:,21:24]=self.Kv*np.expand_dims(1j*rot,0)
        dEm=np.zeros((3,24),complex);dEm[:,:3]=np.diag(1j*Em)
        V=self.volts(x);Is=self.ys*(Em-V[:3])
        dIs=self.ys[:,None]*(dEm-dV[:3])
        dpe=(dEm*np.conj(Is[:,None])+Em[:,None]*np.conj(dIs)).real
        dvdq=np.conj(rot[:,None])*dV[:3]
        dvdq[:,9:12]+=np.diag(-1j*np.conj(rot)*V[:3])
        dvd=dvdq.real;C=dvdq.imag
        vd=(V[:3]*np.conj(rot)).real
        A=np.zeros((24,24));idx=np.arange(3)
        A[idx,3+idx]=WB
        A[3:6]=-dpe/self.M[:,None]
        A[3+idx,3+idx]-=self.D/self.M;A[3+idx,6+idx]+=1/self.M
        A[6+idx,6+idx]=-1/self.tg;A[6+idx,3+idx]=-self.Kg/self.tg
        A[9+idx,12+idx]=1
        A[12+idx,12+idx]=-1/self.tf;A[12+idx,15+idx]=1/self.tf
        A[18:21]=-self.pc[:,None]/vd[:,None]**2*dvd/self.tc[:,None]
        A[21:24]=self.qc[:,None]/vd[:,None]**2*dvd/self.tc[:,None]
        A[18+idx,18+idx]-=1/self.tc;A[21+idx,21+idx]-=1/self.tc
        return A,C,dV

    def input_columns(self,kp,ki):
        B=np.zeros((24,3));ii=np.arange(3)
        B[12+ii,ii]=np.asarray(kp)/self.tf
        B[15+ii,ii]=ki
        return B

    def char(self,s,kp,ki,tau):
        B=self.input_columns(kp,ki)
        return s*np.eye(24)-self.A0-(B*np.exp(-s*np.asarray(tau))[None,:])@self.C

    def char_s(self,s,kp,ki,tau):
        B=self.input_columns(kp,ki);tau=np.asarray(tau)
        return np.eye(24)+(B*(tau*np.exp(-s*tau))[None,:])@self.C

    def impedance(self,s,kp,ki,tau):
        F=self.char(s,kp,ki,tau);r=np.arange(3,6);h=np.setdiff1d(np.arange(24),r)
        Z=F[np.ix_(r,r)]-F[np.ix_(r,h)]@solve(F[np.ix_(h,h)],F[np.ix_(h,r)])
        return BASE_MVA*self.M[:,None]*Z

    def gain_sensitivity(self,lam,kp,ki,tau):
        F=self.char(lam,kp,ki,tau);U,ss,Vh=svd(F);v=Vh[-1].conj();w=U[:,-1]
        den=w.conj()@self.char_s(lam,kp,ki,tau)@v
        ce=self.C@v;expo=np.exp(-lam*np.asarray(tau))
        gp=(w[12:15].conj()/self.tf)*ce*expo/den
        gi=w[15:18].conj()*ce*expo/den
        gt=-(w.conj()@self.input_columns(kp,ki))*(lam*expo)*ce/den
        return gp,gi,gt,ss[-1],abs(den)

    def single_pll_return(self,s,site,kp,ki,tau):
        # Other PLLs stay closed, including their exact delays.
        F=self.char(s,kp,ki,tau)
        triple=np.array([9+site,12+site,15+site]);h=np.setdiff1d(np.arange(24),triple)
        theta=9+site
        z=solve(F[np.ix_(h,h)],-F[h,theta])
        return self.C[site,theta]+self.C[site,h]@z

    def assign_one(self,target,site,kp,ki,tau):
        if abs(target.imag)<1e-12:raise ValueError('nonreal target required')
        G=self.single_pll_return(target,site,kp,ki,tau)
        if abs(G)<1e-12:raise ValueError('zero detector return')
        W=target**2*(1+self.tf[site]*target)*np.exp(target*tau[site])/G
        p=W.imag/target.imag;q=W.real-target.real*p
        pp=np.array(kp,copy=True);ii=np.array(ki,copy=True);pp[site]=p;ii[site]=q
        return pp,ii

@lru_cache(maxsize=None)
def pade_unit(n):
    # [n/n] Padé of exp(-z), used ONLY to generate initial eigenvalue guesses.
    c=np.array([factorial(2*n-k)*factorial(n)/(factorial(2*n)*factorial(k)*factorial(n-k)) for k in range(n+1)])
    den=c[::-1];num=(c*(-1.)**np.arange(n+1))[::-1]
    A,B,C,D=tf2ss(num,den)
    Ab,T=matrix_balance(A,permute=False)
    return Ab,solve(T,B),C@T,float(D[0,0])

def pade_seeds(grid,kp,ki,tau,order=6):
    positive=[i for i,t in enumerate(tau) if t>1e-12];N=24+order*len(positive)
    A=np.zeros((N,N));A[:24,:24]=grid.A0
    B=grid.input_columns(kp,ki)
    for i,t in enumerate(tau):
        if t<=1e-12:A[:24,:24]+=np.outer(B[:,i],grid.C[i])
    ap,bp,cp,dp=pade_unit(order)
    for k,i in enumerate(positive):
        ss=slice(24+k*order,24+(k+1)*order)
        A[:24,:24]+=dp*np.outer(B[:,i],grid.C[i])
        A[:24,ss]=np.outer(B[:,i],cp.ravel())
        A[ss,:24]=bp@grid.C[i:i+1]/tau[i]
        A[ss,ss]=ap/tau[i]
    return eigvals(A)

def refine_root(grid,s,kp,ki,tau,tol=2e-11,maxiter=20):
    s=complex(s)
    for it in range(maxiter):
        F=grid.char(s,kp,ki,tau)
        U,ss,Vh=svd(F);v=Vh[-1].conj();w=U[:,-1]
        den=w.conj()@grid.char_s(s,kp,ki,tau)@v
        if abs(den)<1e-13:return s,ss[-1],False
        step=(w.conj()@F@v)/den
        if abs(step)>50:step*=50/abs(step)
        s-=step
        if abs(step)<tol*(1+abs(s)):
            err=svd(grid.char(s,kp,ki,tau),compute_uv=False)[-1]
            return s,float(err),bool(err<2e-7)
    return s,float(svd(grid.char(s,kp,ki,tau),compute_uv=False)[-1]),False

def root_catalog(grid,kp,ki,tau,order=6,real_floor=-30,max_imag=300):
    seeds=pade_seeds(grid,kp,ki,tau,order)
    seeds=seeds[(seeds.real>real_floor)&(seeds.imag>=-1e-8)&(np.abs(seeds.imag)<max_imag)]
    roots=[];errs=[]
    for s in sorted(seeds,key=lambda x:-x.real):
        z,e,ok=refine_root(grid,s,kp,ki,tau)
        if not ok or abs(z)<1e-6 or z.imag<-1e-5:continue
        if not any(abs(z-r)<1e-5 for r in roots):roots.append(z);errs.append(e)
    ix=np.argsort([-z.real for z in roots])
    return np.array(roots)[ix],np.array(errs)[ix]

def alpha_fast(grid,kp,ki,tau,order=5):
    z=pade_seeds(grid,kp,ki,tau,order)
    z=z[np.abs(z)>1e-6]
    return float(np.max(z.real))
