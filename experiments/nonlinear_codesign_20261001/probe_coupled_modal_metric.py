import json
import numpy as np
import scipy.linalg as la
from threadpoolctl import threadpool_limits
from coupled_interval_model import Model,OUT,mm,add,up,down,magnitude,inverse_enclosure


def modal_metric(A,np_):
    ap=A[:np_,:np_]; am=A[np_:,np_:]; c=A[np_:,:np_]
    ev,v=la.eig(ap)
    cols=[]
    for i,e in enumerate(ev):
        if e.imag>1e-7:
            pair=np.column_stack([v[:,i].real,v[:,i].imag])
            pair/=np.linalg.norm(pair)
            cols.extend([pair[:,0],pair[:,1]])
        elif abs(e.imag)<=1e-7:
            cols.append(v[:,i].real/np.linalg.norm(v[:,i].real))
    tp=np.column_stack(cols)
    assert tp.shape==ap.shape
    cross=la.solve_sylvester(am,-ap,-c)
    tm=np.diag(np.r_[np.ones(39),5*np.ones(39)])
    t=np.block([[tp,np.zeros((np_,78))],[cross@tp,tm]])
    h=np.linalg.inv(t)
    h*=np.max(np.linalg.norm(t,axis=1))
    ti,iv=inverse_enclosure((h,h))
    return h,ti,iv


if __name__=="__main__":
    from synthesize_coupled_contract import symmetric_upper,norm_upper
    with threadpool_limits(limits=1):
        cc=json.loads((OUT/"refined_center.json").read_text(encoding="utf-8"))
        m=Model(center=cc["offset"]);A=np.array(m.raw["A"]);B=np.array(m.raw["B"])
        h,ti,iv=modal_metric(A,len(m.keep))
        nominal=mm((h,h),mm((A,A),ti));ns=add(nominal,(nominal[0].T,nominal[1].T))
        tr=up(np.sqrt(up(np.sum(magnitude(ti)**2,axis=1)*(1+8*m.nx*np.finfo(float).eps))))
        rows=[]
        for r in (0.,1e-9,1e-7,1e-5):
            fo,*_=m.intervals(r*tr,np.zeros(m.nw))
            jl=np.array([q.d[0][:m.nx] for q in fo]);jh=np.array([q.d[1][:m.nx] for q in fo])
            re=mm((h,h),mm((down(jl-A),up(jh-A)),ti))
            decay=-symmetric_upper(add(ns,add(re,(re[0].T,re[1].T))))/2
            gain=norm_upper(mm((h,h),(B,B)))
            row=dict(kind="modal_plus_exact_sensor_cascade",radius=r,condition=float(np.linalg.cond(h)),
                     inverse_check=iv,decay=decay,input_gain=gain)
            rows.append(row);print(json.dumps(row),flush=True)
        (OUT/"modal_metric_probes.json").write_text(json.dumps(rows,indent=2)+"\n",encoding="utf-8")
        np.savez_compressed(OUT/"modal_metric.npz",H=h,Tlo=ti[0],Thi=ti[1])
