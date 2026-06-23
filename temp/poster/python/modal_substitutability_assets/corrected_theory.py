"""
CORRECTED THEORY — find the exact statement that survives finite-difference.

The naive 'add inertia = reinforce lines' is WRONG in sign (verified).
Correct direction: dLtilde from inertia is NEGATIVE-definite-ish; from line reinforcement POSITIVE.
So the real, defensible objects are:

(A) FIRST-ORDER MODAL SHIFT decomposition: how much does nu_c move from
    (i) inertia at node i, vs (ii) reinforcing line (i,j), in the SAME currency (d nu_c).
    Both are exact Rayleigh first-order formulas. Compare them -> the cross/lever question.

(B) The SHUNT DEFECT and its projection on q_c = the modal substitutability classifier (survives).

Goal: verify (A) to machine precision (first-order Rayleigh), and state the honest theorem.
"""
import numpy as np
np.set_printoptions(precision=6,suppress=True)

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
crit=np.argmax(nu); q=Q[:,crit]
print("critical mode",crit,"nu_c=",nu[crit])

# ---- (A) First-order Rayleigh shift of nu_c from each lever ----
# Inertia at node i:  d Ltilde / d(ln M_i) = -1/2 (Ei Lt + Lt Ei),  Ei=e_i e_i^T
# d nu_c = q^T (dLtilde) q
print("\n(A) ANALYTIC vs FD first-order shift of nu_c")
print(f"{'lever':>16} {'analytic':>14} {'finite-diff':>14} {'rel.err':>10}")

eps=1e-4
# inertia levers
inertia_shift={}
for i in range(N):
    Ei=np.zeros((N,N));Ei[i,i]=1.0
    dLt=-0.5*(Ei@Lt0+Lt0@Ei)        # d Ltilde / d ln M_i  (analytic, congruence deriv)
    dnu_an=q@dLt@q
    # FD
    M2=M.copy();M2[i,i]*=np.exp(eps)
    nu2,_=np.linalg.eigh(Ltilde(L,M2))
    dnu_fd=(nu2[crit]-nu[crit])/eps
    inertia_shift[i]=dnu_an
    print(f"{'inertia n'+str(i):>16} {dnu_an:>14.6f} {dnu_fd:>14.6f} {abs(dnu_an-dnu_fd)/(abs(dnu_fd)+1e-12):>10.2e}")

# line levers: d Ltilde/d w_ij = M^-1/2 (ei-ej)(ei-ej)^T M^-1/2
print()
line_shift={}
ms=np.diag(1/np.sqrt(np.diag(M)))
for (i,j,w) in edges:
    b=np.zeros(N);b[i]=1;b[j]=-1
    dLt=ms@np.outer(b,b)@ms          # per unit w
    dnu_an=q@dLt@q
    # FD: increase this line's weight
    edges2=[(a,c,ww*(1+eps) if (a,c)==(i,j) else ww) for (a,c,ww) in edges]
    nu2,_=np.linalg.eigh(Ltilde(build_L(edges2,N),M))
    dnu_fd=(nu2[crit]-nu[crit])/(eps*w)
    line_shift[(i,j)]=dnu_an
    print(f"{'line '+str(i)+'-'+str(j):>16} {dnu_an:>14.6f} {dnu_fd:>14.6f} {abs(dnu_an-dnu_fd)/(abs(dnu_fd)+1e-12):>10.2e}")

print("\n>>> If rel.err ~ 1e-4 or better, the first-order shift formulas are EXACT.")
print(">>> These are the honest, sign-correct objects. Inertia LOWERS nu_c, lines RAISE it.")

# ---- (B) classifier (sign-independent, survives) ----
print("\n(B) Modal substitutability classifier (defect projection on q_c):")
subst=[]
for i in range(N):
    e=np.zeros(N);e[i]=1.0
    defect=0.5*(L@e)
    s=1-abs(defect@q)/(np.linalg.norm(defect)+1e-12)
    subst.append(s)
subst=np.array(subst)
for i in range(N):
    print(f"  node {i}: s={subst[i]:.3f}")
print(f"  ratio max/min of (1-s): {(1-subst.min())/(1-subst.max()):.2f}x  -> classifier strength")
