"""Point checks to diagnose enclosure conservatism; not safety certification."""
import json
import numpy as np
from threadpoolctl import threadpool_limits
from coupled_interval_model import Model,OUT,mm,add,magnitude,inverse_enclosure
from synthesize_coupled_contract import symmetric_upper

with threadpool_limits(limits=1):
    center=np.array(json.loads((OUT/"refined_center.json").read_text(encoding="utf-8"))["offset"])
    saved=np.load(OUT/"certificate_modal.npz");h=saved["H"];t=np.linalg.inv(h)
    ti,_=inverse_enclosure((h,h));rng=np.random.default_rng(605693)
    rows=[]
    for r in (float(saved["radius"]),1e-4):
        for k in range(12):
            y=rng.normal(size=len(h));y*=r/np.linalg.norm(y)
            delta=t@y;m=Model(center=center+delta)
            out,*_=m.intervals(np.zeros(m.nx),np.zeros(m.nw))
            ji=(np.array([q.d[0][:m.nx] for q in out]),np.array([q.d[1][:m.nx] for q in out]))
            a=mm((h,h),mm(ji,ti));s=add(a,(a[0].T,a[1].T))
            decay=-symmetric_upper(s)/2
            rows.append(dict(radius=r,trial=k,point_decay_lower_bound=decay))
        print(json.dumps(dict(radius=r,min_point_decay=min(q["point_decay_lower_bound"] for q in rows if q["radius"]==r))),flush=True)
    report=dict(status="POINT_DIAGNOSTICS_ONLY_NO_REGIONAL_GUARANTEE",seed=605693,rows=rows,
                interpretation="Sampled point decay can remain positive where the regional box enclosure fails; this does not certify unsampled states or inputs")
    (OUT/"wrapping_diagnostic.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
