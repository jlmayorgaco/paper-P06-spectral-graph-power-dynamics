"""
GATE TEST — Inertia-Topology Equivalence Theorem (first order)
================================================================
Claim:  delta_L = 1/2 (E L + L E),  E = diag(dM_i / M_i)
        off-diag (realizable as lines):  dw_ij = 1/2 (e_i + e_j) w_ij
        row-sum defect (NOT a line, = shunt):  1/2 L e
        substitutable fraction per node = 1 - ||defect_i|| / ||delta_L row i||

Tests:
  T1. Analytic delta_L vs finite-difference: does adding inertia at node i,
      then computing the equivalent delta_L, actually leave L_tilde invariant?
  T2. Is the off-diagonal part a physical (positive-weight) line reinforcement?
  T3. Does the substitutable fraction VARY across nodes? (need variation -> classifier)
  T4. Spectral check: does delta_L keep the modal geometry (nu_k, q_k) fixed to O(eps)?
"""

import numpy as np
np.set_printoptions(precision=4, suppress=True)

# ----------------------------------------------------------------------
# 6-bus test network. Weighted graph Laplacian from line susceptances.
# Heterogeneous inertia (the IBR-relevant case: virtual inertia uneven).
# ----------------------------------------------------------------------
# Edges (i,j,weight) -- a small meshed 6-bus system
edges = [
    (0,1,1.2), (0,2,0.8), (1,2,1.5), (1,3,0.6),
    (2,4,1.0), (3,4,0.9), (3,5,1.3), (4,5,0.7),
]
N = 6

def build_L(edges, N, scale=None):
    L = np.zeros((N,N))
    for (i,j,w) in edges:
        ww = w if scale is None else w*scale.get((i,j), 1.0)
        L[i,j] -= ww; L[j,i] -= ww
        L[i,i] += ww; L[j,j] += ww
    return L

L = build_L(edges, N)

# Heterogeneous inertia (deliberately uneven, as virtual inertia would be)
M = np.diag([2.0, 0.5, 3.0, 0.8, 1.5, 0.4])
Minv_sqrt = np.diag(1.0/np.sqrt(np.diag(M)))

def Ltilde(L, M):
    ms = np.diag(1.0/np.sqrt(np.diag(M)))
    return ms @ L @ ms

Lt0 = Ltilde(L, M)
nu0, Q0 = np.linalg.eigh(Lt0)   # modal frequencies^2 and modes of the inertial Laplacian

print("="*64)
print("6-BUS GATE TEST: Inertia-Topology Equivalence Theorem")
print("="*64)
print("\nInertia M (diag):", np.diag(M))
print("Modal stiffness nu_k (eig of L_tilde):", nu0)

# ----------------------------------------------------------------------
# For each node, inject a small relative inertia eps at that node,
# compute the analytic equivalent delta_L = 1/2 (E L + L E),
# its off-diagonal (line) part, its row-sum defect (shunt part),
# and the substitutable fraction.
# ----------------------------------------------------------------------
eps_mag = 0.05   # 5% relative inertia injection (first-order regime)

print("\n" + "-"*64)
print("PER-NODE ANALYSIS (inject {:.0%} relative inertia at each node)".format(eps_mag))
print("-"*64)
print(f"{'node':>4} {'subst.frac':>11} {'||defect||':>11} {'||dL||':>9} {'lines>0?':>9} {'Ltilde inv?':>12}")

results = []
for node in range(N):
    e = np.zeros(N)
    e[node] = eps_mag                     # relative inertia change at this node
    E = np.diag(e)

    # Analytic equivalent topology change (Jordan product)
    dL = 0.5 * (E @ L + L @ E)

    # Off-diagonal = line-realizable part; build the "lines-only" delta
    dL_lines = dL.copy()
    np.fill_diagonal(dL_lines, 0.0)
    # Row-sum defect = shunt part (what lines cannot realize)
    defect = dL @ np.ones(N)              # = 1/2 L e  (theorem prediction)
    defect_pred = 0.5 * (L @ e)

    # Substitutable fraction = 1 - ||defect|| / ||dL||  (Frobenius-ish, use row norms)
    norm_dL = np.linalg.norm(dL)
    norm_defect = np.linalg.norm(defect)
    subst_frac = 1.0 - norm_defect / (norm_dL + 1e-12)

    # Are the equivalent line reinforcements physical (positive)? dw_ij = 1/2(e_i+e_j)w_ij
    lines_positive = True
    for (i,j,w) in edges:
        dw = 0.5*(e[i]+e[j])*w
        if dw < -1e-12:
            lines_positive = False

    # SPECTRAL CHECK: apply the *exact* equivalent and see if L_tilde is preserved
    # Exact: M2 = M(1+e) at node, L2 = L + dL ; compare L_tilde(L2,M2) to Lt0
    M2 = M.copy()
    M2[node,node] = M[node,node]*(1+eps_mag)
    L2 = L + dL
    Lt2 = Ltilde(L2, M2)
    inv_err = np.linalg.norm(Lt2 - Lt0) / np.linalg.norm(Lt0)

    # verify defect matches theorem prediction
    defect_match = np.linalg.norm(defect - defect_pred) < 1e-10

    results.append((node, subst_frac, norm_defect, norm_dL, lines_positive, inv_err))
    print(f"{node:>4} {subst_frac:>11.4f} {norm_defect:>11.4f} {norm_dL:>9.4f} {str(lines_positive):>9} {inv_err:>12.2e}")

# ----------------------------------------------------------------------
# VERDICT
# ----------------------------------------------------------------------
fracs = np.array([r[1] for r in results])
inv_errs = np.array([r[5] for r in results])

print("\n" + "="*64)
print("VERDICT")
print("="*64)
print(f"Substitutable fraction range: [{fracs.min():.3f}, {fracs.max():.3f}]  spread = {fracs.max()-fracs.min():.3f}")
print(f"Most substitutable node (topology cheap):   node {fracs.argmax()}  (frac={fracs.max():.3f})")
print(f"Least substitutable node (inertia needed):  node {fracs.argmin()}  (frac={fracs.min():.3f})")
print(f"First-order invariance error (should be O(eps^2)~{eps_mag**2:.0e}): max={inv_errs.max():.2e}")

if fracs.max()-fracs.min() > 0.05:
    print("\n>>> CLASSIFIER EXISTS: substitutable fraction varies across nodes.")
    print(">>> The theorem distinguishes 'topology-substitutable' from 'inertia-inevitable' nodes.")
else:
    print("\n>>> WARNING: fraction nearly constant. Classifier weak; theorem may be vacuous here.")
