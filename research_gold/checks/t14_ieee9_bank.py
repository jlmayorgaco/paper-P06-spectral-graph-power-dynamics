"""T14: nine-bus bank. Delay transport of the three PLL gains 20->21 ms toward a chosen oscillatory root (assumed dynamic data)."""
import os, sys, json, pathlib
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np
B = pathlib.Path(__file__).resolve().parents[1] / 'external_gpt_kp_only' / 'BND_IEEE9_Python'
sys.path.insert(0, str(B / 'src'))
from model import *
cfg = json.load(open(B / 'data/config.json'))
g = Grid9(rho=.75); kp = np.array(cfg['kp0']); ki = np.array(cfg['ki0'])
t20 = np.full(3, .020); t21 = np.full(3, .021)
def tr(kp, ki, lam, h):
    a, w = lam.real, lam.imag; e = np.exp(a*h)
    return (e*(kp*np.cos(w*h)+(ki+a*kp)/w*np.sin(w*h)), e*(ki*np.cos(w*h)-(abs(lam)**2*kp+a*ki)/w*np.sin(w*h)))
for name, tt in (('20ms', t20), ('21ms', t21)):
    r, _ = root_catalog(g, kp, ki, tt, 8, real_floor=-15, max_imag=500)
    print(name, 'cfg gains, rightmost', r[:4], flush=True)
r20, _ = root_catalog(g, kp, ki, t20, 8, real_floor=-15, max_imag=500)
for lam in [z for z in r20 if z.imag > 1][:8]:
    p, i = tr(kp, ki, lam, .001)
    r, _ = root_catalog(g, p, i, t21, 8, real_floor=-15, max_imag=500)
    print('protect', lam, '-> rightmost after 21ms', r[0].real, flush=True)
