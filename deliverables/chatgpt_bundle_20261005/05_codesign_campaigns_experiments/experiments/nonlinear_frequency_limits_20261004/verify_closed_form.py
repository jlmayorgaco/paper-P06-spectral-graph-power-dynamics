"""Evaluate exact constructions, not ODE simulations or IEEE-39 predictions."""
from pathlib import Path
import hashlib
import json
import platform
import csv
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import root
import sympy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent
P=1.0
zstar=np.arcsin(P/2)
floor=2*zstar/P
m_bad=41/40
m_good=1.5
k=2/m_good
ell_min=zstar*k*k/(P-zstar*k)
ell=2*ell_min
a=zstar*k*(k+ell)/ell
Omega=P/2

def trajectory(t):
    e=np.exp(-k*t)
    e2=np.exp(-(k+ell)*t)
    w=Omega*(1-e)
    wd=Omega*k*e
    z=a*((1-e)/k-(1-e2)/(k+ell))
    zd=a*(e-e2)
    zdd=a*(-k*e+(k+ell)*e2)
    nu1=w+zd/2
    nu2=w-zd/2
    nud1=wd+zdd/2
    nud2=wd-zdd/2
    u1=m_good*zdd/4+zd/2+np.sin(z)-P/2
    return w,z,nu1,nu2,nud1,nud2,u1

t=np.linspace(0,18,3601)
w,z,nu1,nu2,nud1,nud2,u1=trajectory(t)
r1=m_good*nud1/2+nu1+np.sin(z)-P-u1
r2=m_good*nud2/2+nu2-np.sin(z)+u1
assert np.max(np.abs(r1))<2e-14 and np.max(np.abs(r2))<2e-14
assert 0<a<P
assert np.max(nu1)<=Omega+1e-12 and np.min(nu2)>=-Omega-1e-12
assert np.min(nu1)>=-1e-12 and np.max(nu2)<=Omega+1e-12
q1_area=quad(lambda s: Omega-trajectory(s)[2],0,40,epsabs=1e-12)[0]
q2_area=quad(lambda s: Omega-trajectory(s)[3],0,40,epsabs=1e-12)[0]
assert abs(q1_area-(P*m_good/4-zstar/2))<1e-11
assert abs(q2_area-(P*m_good/4+zstar/2))<1e-11

m,zs,pp=sp.symbols('m z P', positive=True)
crit=2*zs/pp
assert sp.simplify(pp*crit/4-zs/2)==0
bad_area=sp.Rational(41,160)-sp.pi/12
assert sp.simplify(bad_area-(123-40*sp.pi)/480)==0
assert float(bad_area)<0
barrier=1-sp.pi/4-sp.Rational(41,320)
assert float(barrier)>0

# Frozen static checks: six-node graph, every site, both signs, two amplitudes.
# These solve endpoints only; they are not a dynamic validation campaign.
edges=[(0,1),(1,2),(2,3),(3,4),(4,5),(5,0),(0,3),(1,4)]
B=np.zeros((6,len(edges)))
for e,(i,j) in enumerate(edges):
    B[i,e],B[j,e]=1,-1
gamma=np.array([2.,3.,1.5,2.5,4.,2.2,1.7,2.8])
D=np.array([1.,2.,1.5,.8,1.2,.7])
c=D/D.sum()
q0=np.array([.08,-.12,.09,-.07,.03,-.01])
F=lambda q: B@(gamma*np.sin(B.T@q))
rows=[]
for bus in range(6):
    for step in [-.5,-.2,.2,.5]:
        rhs=step*(np.eye(6)[bus]-c)
        sol=root(lambda x: (F(q0+np.r_[x,0])-F(q0)-rhs)[:5],
                 np.zeros(5),tol=1e-11)
        delta=np.r_[sol.x,0]
        residual=float(np.linalg.norm(F(q0+delta)-F(q0)-rhs))
        assert residual<1e-10
        eta0=B.T@q0
        eta1=B.T@(q0+delta)
        assert max(abs(eta0).max(),abs(eta1).max())<np.pi/2
        wbar=gamma*np.cos((eta0+eta1)/2)*np.sinc((eta1-eta0)/(2*np.pi))
        Lbar=(B*wbar)@B.T
        g=np.eye(6)[bus]-c
        R=float(g@np.linalg.pinv(Lbar,hermitian=True)@g)
        phi=delta/step-c@(delta/step)
        assert int(np.argmax(phi))==bus
        assert abs(phi[bus]-R)<1e-10
        vals,vecs=np.linalg.eigh(Lbar)
        Rg=float(np.sum((vecs[:,1:].T@g)**2/vals[1:]))
        cutmax=0.
        for mask in range(1,2**6-1):
            if not (mask>>bus)&1: continue
            ix=np.array([(mask>>j)&1 for j in range(6)],dtype=bool)
            crossing=np.array([ix[i]!=ix[j] for i,j in edges])
            cutmax=max(cutmax,float((1-c[ix].sum())**2/wbar[crossing].sum()))
        assert cutmax<=R+1e-11 and abs(Rg-R)<1e-10
        rows.append({'bus':bus+1,'P':step,'endpoint_residual':residual,
          'max_potential_bus':int(np.argmax(phi))+1,'R_secant':R,
          'max_potential':float(phi.max()),'R_GSP':Rg,'max_cut_lower_bound':cutmax,
          'unfiltered_inertia_floor':D.sum()**2*R,
          'area_window_W0p5_floor':D.sum()**2*R-D.sum()*.5/2})
with (OUT/'STATIC_GRAPH_CHECKS.csv').open('w',newline='',encoding='utf-8') as f:
    writer=csv.DictWriter(f,fieldnames=rows[0].keys())
    writer.writeheader();writer.writerows(rows)

results={
 'status':'PROVED_REDUCED_MODEL; exact formulas with numerical algebra checks',
 'dynamic_ode_dde_simulations':0,
 'model':'Dimensionless two-node nonlinear lossless sine network; D=I, M=(m/2)I, d=(P,0)',
 'P':P,'Omega':Omega,'z_star':float(zstar),
 'sharp_total_inertia_infimum':float(floor),
 'sharp_infimum_exact':'2 asin(P/2)/P; at P=1: pi/3',
 'linearized_floor':1.0,
 'below_threshold':{'m':m_bad,'COI_exact':'0.5*(1-exp(-2*t/m))',
   'bus1_area_exact':'(123-40*pi)/480','bus1_area':float(bad_area),
   'bus2_area':m_bad/4+zstar/2,'COI_area':m_bad/4,
   'initial_rotating_energy':41/320,'cohesive_boundary_energy':float(1-sp.pi/4),
   'energy_barrier_margin':float(barrier),
   'claim':'Convergence proved for u=0; some local endpoint crossing is unavoidable for every convergent energy-neutral controller.'},
 'constructive_case':{'m':m_good,'k':k,'ell_min':ell_min,'ell':ell,'a':a,
   'bus1_area_formula':m_good/4-zstar/2,'bus1_area_quadrature':q1_area,
   'bus2_area_formula':m_good/4+zstar/2,'bus2_area_quadrature':q2_area,
   'max_equation_residual':float(max(abs(r1).max(),abs(r2).max())),
   'u1_initial':float(u1[0]),
   'scope':'Exact bounded feedforward with u2=-u1, no actuator limit. Absolute node frequency <= final-frequency magnitude is proved algebraically; grid values are diagnostics only.'},
 'window_note':'For a causal rectangular frequency window W, each area adds Omega*W/2. The necessary inertia floor subtracts d0*W/2. Sharpness theorem above is W=0 only.',
 'static_graph_checks':{'cases':len(rows),'all_maxima_at_event_bus':True,
   'max_endpoint_residual':max(x['endpoint_residual'] for x in rows),
   'max_secant_potential_error':max(abs(x['R_secant']-x['max_potential']) for x in rows),
   'scope':'Six-node illustrative graph, endpoint calculations only; no trajectory convergence asserted.'},
 'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'sympy':sp.__version__},
 'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
}
(OUT/'CLOSED_FORM_RESULTS.json').write_text(json.dumps(results,indent=2),encoding='utf-8')

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(1,2,figsize=(11.7,4.5),constrained_layout=True)
pg=np.linspace(.001,1.93,350)
fg=2*np.arcsin(pg/2)/pg
ax[0].plot(pg,fg,color='#146b70',lw=2.5,label='Exact nonlinear threshold')
ax[0].axhline(1,color='#b36b1d',ls='--',lw=1.8,label='Initial-Jacobian threshold')
ax[0].fill_between(pg,.92,fg,color='#ce5c4a',alpha=.12,label='No convergent neutral control avoids crossing')
ax[0].scatter([1],[m_bad],color='#a6382e',s=48,zorder=5)
ax[0].annotate('Monotone COI,\nlocal crossing forced',xy=(1,m_bad),xytext=(.10,1.22),
 arrowprops={'arrowstyle':'->','color':'.35'},fontsize=10)
ax[0].set(xlabel='Localized step P',ylabel='Total inertia m',title='A. A sharp limit for a nonlinear network',ylim=(.92,1.42))
ax[0].legend(loc='upper left',fontsize=8,frameon=False)
ax[1].plot(t,nu1,lw=2,color='#146b70',label='Bus 1')
ax[1].plot(t,nu2,lw=2,color='#94633e',label='Bus 2')
ax[1].plot(t,w,lw=1.8,ls='--',color='#34475a',label='COI')
ax[1].axhline(Omega,color='.5',ls=':',label='Final frequency')
ax[1].axhline(-Omega,color='.65',ls=':',lw=1)
ax[1].set(xlim=(0,5),xlabel='Time',ylabel='Frequency deviation',title='B. An exact feasible construction above it')
ax[1].legend(loc='lower right',frameon=False,fontsize=9)
for axis in ax: axis.grid(alpha=.18)
fig.suptitle('Beyond aggregate damping: the network imposes a local frequency limit',fontsize=14)
fig.supxlabel('Dimensionless analytical example. No ODE simulation, no IEEE-39 claim, no actuator constraints.',fontsize=9,color='.3')
fig.savefig(OUT/'FIG_SHARP_NETWORK_LIMIT.png',dpi=240)
fig.savefig(OUT/'FIG_SHARP_NETWORK_LIMIT.svg')
plt.close(fig)

print(json.dumps({'floor':float(floor),'bad_local_area':float(bad_area),
 'constructed_equation_residual':results['constructive_case']['max_equation_residual'],
 'dynamic_simulations':0},indent=2))
