from pathlib import Path
import hashlib
import json
import time
import tomllib
import numpy as np
from threadpoolctl import threadpool_limits
from coupled_interval_model import Model, OUT, Jet, magnitude

with threadpool_limits(limits=1):
    m = Model()
    data = tomllib.loads((OUT/"julia_probes.toml").read_text(encoding="utf-8"))
    assert data["model_sha256"] == hashlib.sha256((OUT/"model.toml").read_bytes()).hexdigest()
    error = np.zeros(4)
    for probe in data["probes"]:
        result = m.evaluate(probe["z"], probe["w"])
        ref = [probe[k] for k in ("f","v","frequency","rocof")]
        error = np.maximum(error,[np.max(np.abs(np.array(a)-b)) for a,b in zip(result,ref)])
    assert error[0]<2e-7 and error[1]<1e-10 and error[2]<1e-9 and error[3]<1e-8, error
    t = time.perf_counter()
    out,v,ff,rr = m.intervals(np.zeros(m.nx),np.zeros(m.nw))
    lo = np.array([x.d[0] for x in out]); hi = np.array([x.d[1] for x in out])
    ref = np.hstack([np.array(m.raw["A"]),np.array(m.raw["B"])])
    excess = np.maximum(lo-ref,ref-hi).max()
    midpoint_error = np.abs((lo+hi)/2-ref).max()
    # Independent Julia Float64 evaluation need not lie inside the exact enclosure;
    # its own rounding is not interval validated. Record rather than erase excess.
    assert midpoint_error<1e-5
    report = dict(status="FULL_NETWORK_IMPLEMENTATION_AUDITED_NOT_YET_CERTIFIED",
        state_count=m.nx,input_count=m.nw,physical_frequency_retained=True,
        julia_probes=len(data["probes"]),max_absolute_mapping_errors=error.tolist(),
        interval_evaluation_seconds=time.perf_counter()-t,
        reference_jacobian_midpoint_error=float(midpoint_error),
        reference_Float64_jacobian_enclosure_excess=float(excess),
        maximum_point_jacobian_interval_width=float(np.max(hi-lo)),
        network_inverse_neumann_residual=m.last_inverse_residual,
        point_rhs_bound=float(max(magnitude(x.v) for x in out)),
        model_sha256=hashlib.sha256((OUT/"model.toml").read_bytes()).hexdigest(),
        implementation_sha256=hashlib.sha256(Path(__file__).with_name("coupled_interval_model.py").read_bytes()).hexdigest())
    (OUT/"mapping_audit.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    np.savez_compressed(OUT/"point_interval.npz",lo=lo,hi=hi,
                        flo=np.array([x.v[0] for x in out]),fhi=np.array([x.v[1] for x in out]))
    print(json.dumps(report,indent=2))
