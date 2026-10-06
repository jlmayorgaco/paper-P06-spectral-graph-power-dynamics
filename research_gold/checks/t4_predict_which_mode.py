"""T4: does the toll formula predict, from ONE baseline eigen-analysis, which protected mode is safe?
For each candidate protected mode lam: first-order prediction of every other catalogued root after +4 ms at all sites
versus the exact finite all-site transport (tracked roots)."""
import os, sys, pathlib
os.environ.setdefault("OPENBLAS_NUM_THREADS","1"); os.environ.setdefault("MKL_NUM_THREADS","1")
import numpy as np, pandas as pd
from scipy import linalg as la
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'experiments' / 'interaction_decision_20261004'))
import model as MM
from model import Model, P0, TAU
from t1_spillover_ieee39 import transport, roots
base = Model(P0, TAU); H = 0.004
def muKi(mu):                                   # analytic d mu / d Ki_i for all sites (ordinary gain)
    D = base.delta(mu); U, S, Vh = la.svd(D); v = Vh[-1].conj(); l = U[:, -1]
    den = np.vdot(l, base.derivative(mu) @ v)
    return np.array([np.vdot(l, MM.M['Bi'][:, i]) * np.exp(-mu * TAU[i]) * (base.C[i] @ v) / den for i in range(10)])
sens = {mu: muKi(mu) for mu in roots}
rows = []
for lam in roots:
    if lam.imag < 0.2: continue
    pred = {mu: mu.real + H * np.sum((-P0[10:20] * (mu - lam) * (mu - lam.conjugate()) * sens[mu]).real) for mu in roots if abs(mu - lam) > 1e-9}
    cur = list(roots); ok = True
    for h in np.linspace(0, H, 17)[1:]:
        p = P0.copy()
        for i in range(10): p[10 + i], p[20 + i] = transport(P0[10 + i], P0[20 + i], lam, h)
        m = Model(p, TAU + h)
        cur = [m.refine(z)[0] for z in cur]
    exact = {mu0: z.real for mu0, z in zip(roots, cur) if abs(mu0 - lam) > 1e-9}
    rows.append(dict(protected_Hz=lam.imag / 2 / np.pi, protected_re=lam.real, predicted_worst=max(pred.values()), exact_worst=max(exact.values()),
                     Ki_at_4ms=p[20], pred_vs_exact_corr=np.corrcoef(list(pred.values()), list(exact.values()))[0, 1]))
df = pd.DataFrame(rows); df.to_csv('T4_predict_which_mode.csv', index=False); print(df.round(4).to_string(index=False))
from scipy.stats import spearmanr
print('Spearman(predicted worst, exact worst) =', spearmanr(df.predicted_worst, df.exact_worst))
