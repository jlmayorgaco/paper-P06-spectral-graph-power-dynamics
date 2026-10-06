"""T1: validate the delay-transport law and the modal spillover identity on the exported IEEE-39 DDE
(204 states, exact exponentials, 40 ms at every PLL). Read-only use of experiments/interaction_decision_20261004/model.py."""
import os, sys, json, pathlib
os.environ.setdefault("OPENBLAS_NUM_THREADS","1"); os.environ.setdefault("MKL_NUM_THREADS","1")
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'experiments' / 'interaction_decision_20261004'))
from model import Model, P0, TAU, ROOTS

def transport(kp, ki, lam, h):
    a, w = lam.real, lam.imag; e = np.exp(a * h)
    return (e * (kp * np.cos(w * h) + (ki + a * kp) / w * np.sin(w * h)),
            e * (ki * np.cos(w * h) - (abs(lam) ** 2 * kp + a * ki) / w * np.sin(w * h)))

base = Model(P0, TAU)
roots = []
for s in ROOTS:
    try:
        r, _ = base.refine(s)
        if r.imag > 1e-6 and all(abs(r - q) > 1e-6 for q in roots): roots.append(r)
    except Exception: pass
roots = sorted(roots, key=lambda z: z.imag)
if __name__ == '__main__': print('catalog roots (rad/s):', [f'{z.real:+.4f}{z.imag:+.3f}j' for z in roots])
def main():
    out = []
    for site in [0, 7, 3]:                       # buses 30, 37, 33
        for lam in roots:
            if lam.imag < 1: continue          # protected pole must be oscillatory and well scaled
            kp, ki = P0[10 + site], P0[20 + site]
            def model_h(h):
                p = P0.copy(); tau = TAU.copy(); tau[site] += h
                p[10 + site], p[20 + site] = transport(kp, ki, lam, h)
                return Model(p, tau)
            def model_ki(d):
                p = P0.copy(); p[20 + site] += d; return Model(p, TAU)
            h = 2e-6; d = 1e-3
            # protected pole stays?
            drift = abs(model_h(1e-3).refine(lam)[0] - lam)
            for mu in roots:
                if abs(mu - lam) < 1e-6: continue
                lhs = (model_h(h).refine(mu)[0] - model_h(-h).refine(mu)[0]) / (2 * h)
                muki = (model_ki(d).refine(mu)[0] - model_ki(-d).refine(mu)[0]) / (2 * d)
                rhs = -kp * (mu - lam) * (mu - lam.conjugate()) * muki
                out.append(dict(site=30 + site, lam_Hz=lam.imag / 2 / np.pi, mu_Hz=mu.imag / 2 / np.pi, lhs_re=lhs.real, rhs_re=rhs.real,
                                abs_err=abs(lhs - rhs), rel_err=abs(lhs - rhs) / max(abs(lhs), 1e-12), protected_drift_1ms=drift,
                                dist2=abs((mu - lam) * (mu - lam.conjugate()))))
    import pandas as pd
    df = pd.DataFrame(out); df.to_csv('T1_spillover_ieee39.csv', index=False)
    big = df[df.lhs_re.abs() > 1e-6]
    print('pairs', len(df), '| well-scaled', len(big), '| median rel err', big.rel_err.median(), '| max rel err', big.rel_err.max())
    print('max protected-pole drift after +1 ms transport:', df.protected_drift_1ms.max())
    print('sign agreement Re:', (np.sign(big.lhs_re) == np.sign(big.rhs_re)).mean())
    print(big.sort_values('lhs_re').tail(6)[['site', 'lam_Hz', 'mu_Hz', 'lhs_re', 'rhs_re', 'rel_err']].to_string())
    print(big.sort_values('lhs_re').head(4)[['site', 'lam_Hz', 'mu_Hz', 'lhs_re', 'rhs_re', 'rel_err']].to_string())

if __name__ == "__main__":
    main()
