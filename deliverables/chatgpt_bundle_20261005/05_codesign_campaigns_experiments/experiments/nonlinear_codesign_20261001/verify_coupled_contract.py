"""Recompute the nonlinear enclosures and check the proposed invariant ellipsoid.

This checker does not optimize. It shares the audited interval-model evaluator;
the finite tests do not replace its universal derivative enclosure assumptions.
"""
from pathlib import Path
from fractions import Fraction
import hashlib
import json
import time
import numpy as np
from threadpoolctl import threadpool_limits
from coupled_interval_model import Model,OUT,mm,add,up,down,magnitude,inverse_enclosure
from synthesize_coupled_contract import symmetric_upper,norm_upper

sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
with threadpool_limits(limits=1):
    start=time.perf_counter()
    syn=json.loads((OUT/"synthesis_modal.json").read_text(encoding="utf-8"))
    center=json.loads((OUT/"refined_center.json").read_text(encoding="utf-8"))
    assert syn["model_sha256"]==center["model_sha256"]==sha(OUT/"model.toml")
    assert syn["center_sha256"]==sha(OUT/"refined_center.json")
    assert syn["implementation_sha256"]==center["implementation_sha256"]==sha(Path(__file__).with_name("coupled_interval_model.py"))
    assert syn["script_sha256"]==sha(Path(__file__).with_name("synthesize_coupled_contract.py"))
    m=Model(center=center["offset"]);saved=np.load(OUT/"certificate_modal.npz")
    h=saved["H"];r=float(saved["radius"]);eps=float(saved["epsilon"])
    ti,invcheck=inverse_enclosure((h,h))
    rowbounds=up(np.sqrt(up(np.sum(magnitude(ti)**2,axis=1)*(1+8*m.nx*np.finfo(float).eps))))
    xr=up(r*rowbounds)
    fo,v,ff,rr=m.intervals(xr,np.full(m.nw,eps))
    ji=(np.array([x.d[0][:m.nx] for x in fo]),np.array([x.d[1][:m.nx] for x in fo]))
    bi=(np.array([x.d[0][m.nx:] for x in fo]),np.array([x.d[1][m.nx:] for x in fo]))
    # Direct congruence of the freshly recomputed full Jacobian interval, without
    # the nominal-A decomposition used during synthesis.
    aj=mm((h,h),mm(ji,ti))
    sym=add(aj,(aj[0].T,aj[1].T))
    decay=float(down(-symmetric_upper(sym)/2))
    b=norm_upper(mm((h,h),bi))
    f0=(np.array(center["f0lo"])[:,None],np.array(center["f0hi"])[:,None])
    hf0=mm((h,h),f0)
    drift=float(up(np.sqrt(up(np.sum(magnitude(hf0)**2)*(1+8*m.nx*np.finfo(float).eps)))))
    margin=Fraction(decay)*Fraction(r)-Fraction(b)*Fraction(eps)-Fraction(drift)
    assert decay>0 and margin>0
    initial=mm((h,h),(-m.center[:,None],-m.center[:,None]))
    initial_norm=float(up(np.sqrt(up(np.sum(magnitude(initial)**2)*(1+8*m.nx*np.finfo(float).eps)))))
    assert initial_norm<r
    fmax=max(magnitude(q.v) for q in ff);rmax=max(magnitude(q.v) for q in rr)
    assert fmax<.5 and rmax<.5
    vlo=[];vhi=[]
    for k in range(39):
        v2=v[2*k]**2+v[2*k+1]**2
        vlo.append(float(down(np.sqrt(max(0.,v2.v[0])))));vhi.append(float(up(np.sqrt(max(0.,v2.v[1])))))
    assert min(vlo)>.9 and max(vhi)<1.1
    domains=[]
    for i,ix in enumerate(m.gf):
        vals=[]
        for index in ix:
            k=int(np.where(m.keep==index)[0][0]);base=m.x0[index]+m.center[k]
            vals.append((float(down(base-xr[k])),float(up(base+xr[k]))))
        imax=float(up(np.sqrt(sum(max(abs(a),abs(b))**2 for a,b in vals[5:7]))))
        domains.append(dict(bus=i+30,current_magnitude_model_base_required=imax,
                            vdc_min=vals[8][0],vdc_max=vals[8][1]))
    report=dict(status="NONLINEAR_COUPLED_ENCLOSURE_RECHECK_PASSED_MICROSCOPIC_REGION",
        verification_type="Directed binary64 matrix enclosures plus mpmath interval elementary functions; not a proof-assistant certificate",
        state_count=m.nx,input_count=m.nw,radius=r,input_radius=eps,
        certified_decay_lower_bound=decay,input_gain_upper_bound=b,trim_drift_upper_bound=drift,
        inward_margin_exact=str(margin),inward_margin_float=float(margin),
        initial_original_trim_metric_norm_bound=initial_norm,
        original_trim_inside=bool(initial_norm<r),frequency_bound=float(fmax),rocof_bound=float(rmax),
        voltage_min=min(vlo),voltage_max=max(vhi),hardware_requirements_not_validated=domains,
        metric_inverse_neumann_check=invcheck,network_inverse_neumann_check=m.last_inverse_residual,
        seconds=time.perf_counter()-start,
        source_hashes={name:sha(OUT/name) for name in ("model.toml","refined_center.json","synthesis_modal.json","certificate_modal.npz")},
        implementation_sha256=sha(Path(__file__).with_name("coupled_interval_model.py")),
        script_sha256=sha(Path(__file__)),
        scope="All 39-bus current channels and 10 DC-source channels in the declared joint instantaneous ball; fixed design; practically negligible envelope, no capacity optimum or hardware feasibility claim")
    (OUT/"verification.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in report.items() if k not in ("hardware_requirements_not_validated","source_hashes","inward_margin_exact")},indent=2))
