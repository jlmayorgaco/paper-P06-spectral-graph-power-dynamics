"""
FULL RESULT PACKAGE — Modal inertia-topology substitutability
Generates the numbers and the classifier figure data for the poster/paper.
"""
import numpy as np
import json
np.set_printoptions(precision=4, suppress=True)

edges=[(0,1,1.2),(0,2,0.8),(1,2,1.5),(1,3,0.6),(2,4,1.0),(3,4,0.9),(3,5,1.3),(4,5,0.7)]
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

# Pick the critical mode. For damping-ratio with homogeneous d/m, zeta_k ~ 1/sqrt(nu_k),
# so the HIGHEST nu is least damped -> critical. (Documented choice; problem-dependent.)
crit=np.argmax(nu)
q_c=Q[:,crit]
print("Critical mode index:",crit,"nu_c=",nu[crit])

# ---- RESULT 1: modal substitutability classifier ----
eps=0.05
subst=[]
for node in range(N):
    e=np.zeros(N);e[node]=eps
    defect=0.5*(L@e)
    proj=abs(defect@q_c)/(np.linalg.norm(defect)+1e-12)
    s_i=1.0-proj   # high = substitutable by topology for the critical mode
    subst.append(s_i)
subst=np.array(subst)
print("\nRESULT 1 - Modal substitutability per node (high=topology-cheap):")
for n in range(N):
    tag="TOPOLOGY (cheap)" if subst[n]>subst.mean() else "INERTIA (needed)"
    print(f"  node {n}: s={subst[n]:.3f}  -> {tag}")
print(f"  spread={subst.max()-subst.min():.3f}, ratio={ (1-subst.min())/(1-subst.max()+1e-12):.2f}x")

# ---- RESULT 2: finite-difference validation of the equivalence ----
# Show that injecting inertia eps at node, vs applying the equivalent dL, give the
# SAME critical-mode shift to first order.
print("\nRESULT 2 - Equivalence validation (inertia vs equivalent topology), nu_c shift:")
print(f"  {'node':>4} {'d_nu_c (inertia)':>18} {'d_nu_c (equiv topo)':>20} {'rel.err':>10}")
fd_data=[]
for node in range(N):
    e=np.zeros(N);e[node]=eps;E=np.diag(e)
    # Path A: real inertia injection
    M2=M.copy();M2[node,node]*=(1+eps)
    nuA,_=np.linalg.eigh(Ltilde(L,M2))
    dnuA=nuA[crit]-nu[crit]
    # Path B: equivalent topology change dL (lines part only, the realizable part)
    dL=0.5*(E@L+L@E)
    dL_lines=dL.copy();np.fill_diagonal(dL_lines,np.diag(dL))  # keep full dL for fair eq test
    nuB,_=np.linalg.eigh(Ltilde(L+dL,M))   # apply topology, keep M
    dnuB=nuB[crit]-nu[crit]
    relerr=abs(dnuA-dnuB)/(abs(dnuA)+1e-12)
    fd_data.append((node,dnuA,dnuB,relerr))
    print(f"  {node:>4} {dnuA:>18.6f} {dnuB:>20.6f} {relerr:>10.2e}")

# ---- RESULT 3: the design trade-off witness ----
# Most-substitutable node: same damping benefit from lines as from inertia.
# Least-substitutable node: inertia is modally irreplaceable.
best=subst.argmax();worst=subst.argmin()
print(f"\nRESULT 3 - Design witness:")
print(f"  Node {best}: MOST substitutable (s={subst[best]:.3f}) -> reinforce incident lines instead of inertia")
print(f"  Node {worst}: LEAST substitutable (s={subst[worst]:.3f}) -> inertia (or shunt/SynCon) is needed")

# ---- export for figures ----
out={
 "nodes":list(range(N)),
 "inertia":list(np.diag(M)),
 "substitutability":list(np.round(subst,4)),
 "critical_mode":int(crit),
 "nu":list(np.round(nu,4)),
 "edges":edges,
 "best_node":int(best),"worst_node":int(worst),
 "fd_validation":[{"node":n,"dnu_inertia":round(a,6),"dnu_topo":round(b,6),"relerr":float(f'{e:.2e}')} for (n,a,b,e) in fd_data],
}
with open("poster_data.json","w") as f:
    json.dump(out,f,indent=2)
print("\n[exported poster_data.json]")
