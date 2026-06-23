"""
RESCUE TEST — why is the substitutable fraction flat, and is there a metric that discriminates?

Hypotheses:
  H1. The Frobenius-norm fraction averages out node structure. A LOCAL defect metric
      (defect AT the injected node, relative to the inertia effect there) may discriminate.
  H2. The right design question is not "fraction of dL" but: to cancel the SHUNT defect
      with real lines, how much extra reinforcement is needed? That cost varies by node
      degree/connectivity.
  H3. The discriminating quantity is the defect projected onto the CRITICAL mode, since
      only the critical mode sets the damping margin. A node whose shunt defect is
      orthogonal to q_c is effectively substitutable FOR DAMPING even if not globally.
"""

import numpy as np
np.set_printoptions(precision=4, suppress=True)

edges = [(0,1,1.2),(0,2,0.8),(1,2,1.5),(1,3,0.6),(2,4,1.0),(3,4,0.9),(3,5,1.3),(4,5,0.7)]
N=6
def build_L(edges,N):
    L=np.zeros((N,N))
    for(i,j,w) in edges:
        L[i,j]-=w;L[j,i]-=w;L[i,i]+=w;L[j,j]+=w
    return L
L=build_L(edges,N)
M=np.diag([2.0,0.5,3.0,0.8,1.5,0.4])
def Ltilde(L,M):
    ms=np.diag(1/np.sqrt(np.diag(M)));return ms@L@ms
Lt0=Ltilde(L,M)
nu,Q=np.linalg.eigh(Lt0)
# critical damping mode: smallest nonzero nu (Fiedler) AND highest nu both candidates;
# for damping ratio the critical is problem-dependent. Test against Fiedler (mode 1) and top.
q_fiedler=Q[:,1]
q_top=Q[:,-1]

eps=0.05
deg=np.array([sum(w for (i,j,w) in edges if i==n or j==n) for n in range(N)])
print("node degrees (weighted):",deg)
print("inertia:",np.diag(M))

print("\n%-4s %-10s %-12s %-12s %-14s %-14s"%("node","localdef","shunt/effect","defproj_Fdlr","defproj_top","reinf_to_cancel"))
H1=[];H3F=[];H3T=[];H2=[]
for node in range(N):
    e=np.zeros(N);e[node]=eps;E=np.diag(e)
    dL=0.5*(E@L+L@E)
    defect=0.5*(L@e)                       # shunt defect vector
    # H1: LOCAL defect = defect at injected node / inertia effect magnitude at that node
    inertia_effect=abs(0.5*(L@e)[node])+np.sum(np.abs(dL[node,:]))  # total change touching node
    local_def=abs(defect[node])/(inertia_effect+1e-12)
    # H2: to cancel shunt defect with real lines you must add grounding-equivalent
    #     reinforcement ~ |defect|/degree (cheaper if node well connected)
    reinf_to_cancel=np.sum(np.abs(defect))/(deg[node]+1e-12)
    # H3: projection of shunt defect onto critical mode (only this matters for that mode)
    projF=abs(defect@q_fiedler)/(np.linalg.norm(defect)+1e-12)
    projT=abs(defect@q_top)/(np.linalg.norm(defect)+1e-12)
    H1.append(local_def);H3F.append(projF);H3T.append(projT);H2.append(reinf_to_cancel)
    print("%-4d %-10.4f %-12.4f %-12.4f %-14.4f %-14.4f"%(node,local_def,np.linalg.norm(defect)/np.linalg.norm(dL),projF,projT,reinf_to_cancel))

def spread(x):
    x=np.array(x);return x.max()-x.min(), x.max()/(x.min()+1e-12)
for name,arr in [("H1 local defect",H1),("H2 reinf-to-cancel",H2),("H3 proj Fiedler",H3F),("H3 proj top",H3T)]:
    s,r=spread(arr)
    print(f"\n{name}: spread={s:.4f}  ratio max/min={r:.2f}  -> {'DISCRIMINATES' if r>1.5 else 'flat'}")
