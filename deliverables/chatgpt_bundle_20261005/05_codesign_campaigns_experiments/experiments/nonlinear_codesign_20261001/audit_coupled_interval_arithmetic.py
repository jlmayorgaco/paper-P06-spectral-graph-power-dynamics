"""Independent rational/MP checks of custom outward-rounded interval primitives.

Finite implementation tests, not a substitute for the arithmetic error analysis.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import numpy as np
import mpmath as mp
from threadpoolctl import threadpool_limits
from coupled_interval_model import add,mul,inv,mm,inverse_enclosure,OUT

rng=np.random.default_rng(415891)
checks=0
def contains(interval,value):
    global checks
    assert F(float(interval[0]))<=value<=F(float(interval[1]))
    checks+=1

with threadpool_limits(limits=1):
    for k in range(250):
        a=np.sort(rng.normal(size=2)*10**rng.uniform(-10,10))
        b=np.sort(rng.normal(size=2)*10**rng.uniform(-10,10))
        for x in a:
            for y in b:
                contains(add(a,b),F(float(x))+F(float(y)))
                contains(mul(a,b),F(float(x))*F(float(y)))
        if a[0]*a[1]>0:
            for x in a:contains(inv(a),1/F(float(x)))
    for k in range(30):
        am=rng.normal(size=(5,7));bm=rng.normal(size=(7,4))
        ar=rng.uniform(0,.01,size=am.shape);br=rng.uniform(0,.01,size=bm.shape)
        ai=(am-ar,am+ar);bi=(bm-br,bm+br);ci=mm(ai,bi)
        for p in range(3):
            aa=np.where(rng.random(am.shape)>.5,ai[0],ai[1]);bb=np.where(rng.random(bm.shape)>.5,bi[0],bi[1])
            for i in range(5):
                for j in range(4):
                    value=sum(F(float(aa[i,l]))*F(float(bb[l,j])) for l in range(7))
                    contains((ci[0][i,j],ci[1][i,j]),value)
    mp.mp.dps=80
    for k in range(20):
        g=rng.normal(size=(5,5))+.5*np.eye(5)
        radius=np.full((5,5),1e-7)
        ii,check=inverse_enclosure((g-radius,g+radius))
        for p in range(3):
            gg=g+radius*rng.choice([-1,1],size=g.shape)
            exact=mp.inverse(mp.matrix(gg.tolist()))
            for i in range(5):
                for j in range(5):
                    assert mp.mpf(float(ii[0][i,j]))<=exact[i,j]<=mp.mpf(float(ii[1][i,j]))
                    checks+=1
    report=dict(status="ARITHMETIC_IMPLEMENTATION_CHECKS_PASSED",checks=checks,
                seed=415891,scope="Finite exact-rational and 80-digit checks of primitive containment; no independent proof of the entire interval implementation",
                implementation_sha256=hashlib.sha256(Path(__file__).with_name("coupled_interval_model.py").read_bytes()).hexdigest(),
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/"arithmetic_audit.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))
