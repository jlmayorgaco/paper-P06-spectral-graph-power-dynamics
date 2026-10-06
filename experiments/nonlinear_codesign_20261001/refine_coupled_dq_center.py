"""Refine the equilibrium in the exact angular quotient, retaining common speed.

The center is stored as original frozen coordinates plus a separate binary64
offset. It is never rounded back to one Float64 sum for the enclosure evaluation.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import mpmath as mp
from threadpoolctl import threadpool_limits
from coupled_dq_model import Model, OUT, lower, upper

with threadpool_limits(limits=1):
    m=Model();A=np.array(m.raw["A"]);offset=np.zeros(m.nx);history=[]
    for k in range(6):
        m=Model(center=offset)
        f=m.evaluate([mp.iv.mpf(0)]*m.nx,[mp.iv.mpf(0)]*m.nw)[0]
        flo=np.array([lower(q) for q in f]);fhi=np.array([upper(q) for q in f])
        norm=float(np.maximum(np.abs(flo),np.abs(fhi)).max())
        history.append(norm);print(k,norm,flush=True)
        if norm<1e-21:break
        offset-=np.linalg.solve(A,(flo+fhi)/2)
    report=dict(status="REFINED_QUOTIENT_EQUILIBRIUM_WITH_RESIDUAL_BOUND",
                offset=offset.tolist(),residual_history=history,
                maximum_coordinate_offset=float(np.max(np.abs(offset))),
                f0lo=flo.tolist(),f0hi=fhi.tolist(),
                model_sha256=hashlib.sha256((OUT/"model.toml").read_bytes()).hexdigest(),
                implementation_sha256=hashlib.sha256(Path(__file__).with_name("coupled_dq_model.py").read_bytes()).hexdigest(),
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/"refined_center.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    assert norm<1e-19
