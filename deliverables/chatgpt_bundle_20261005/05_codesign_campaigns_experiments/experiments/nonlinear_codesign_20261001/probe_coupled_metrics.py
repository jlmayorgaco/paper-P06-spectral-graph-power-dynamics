import json
import numpy as np
from threadpoolctl import threadpool_limits
from coupled_interval_model import Model,OUT,mm,add,up,down,magnitude
from synthesize_coupled_contract import metric,symmetric_upper,norm_upper

with threadpool_limits(limits=1):
    cc=json.loads((OUT/"refined_center.json").read_text(encoding="utf-8"))
    m=Model(center=cc["offset"]);A=np.array(m.raw["A"]);B=np.array(m.raw["B"])
    rows=[]
    for kind in ("plant_observability",):
        n=len(m.keep)
        A=A[:n,:n];B=B[:n,:]
        h,ti,iv=metric(A,B=B,kind=kind)
        nominal=mm((h,h),mm((A,A),ti));ns=add(nominal,(nominal[0].T,nominal[1].T))
        tr=up(np.sqrt(up(np.sum(magnitude(ti)**2,axis=1)*(1+8*m.nx*np.finfo(float).eps))))
        for r in (1e-9,1e-7,1e-5):
            fo,*_=m.intervals(np.r_[r*tr,np.zeros(78)],np.zeros(m.nw))
            jl=np.array([q.d[0][:n] for q in fo[:n]]);jh=np.array([q.d[1][:n] for q in fo[:n]])
            re=mm((h,h),mm((down(jl-A),up(jh-A)),ti))
            decay=-symmetric_upper(add(ns,add(re,(re[0].T,re[1].T))))/2
            gain=norm_upper(mm((h,h),(B,B)))
            row=dict(kind=kind,radius=r,condition=float(np.linalg.cond(h)),decay=decay,
                     input_gain=gain,linear_radius_per_input=1/gain)
            rows.append(row);print(json.dumps(row),flush=True)
    (OUT/"plant_metric_probes.json").write_text(json.dumps(rows,indent=2)+"\n",encoding="utf-8")
