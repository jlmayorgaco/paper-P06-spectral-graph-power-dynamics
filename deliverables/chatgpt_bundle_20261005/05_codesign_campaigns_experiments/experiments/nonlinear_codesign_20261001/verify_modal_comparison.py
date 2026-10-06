"""Re-enclose a saved block-modal region and check every boundary row exactly."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import numpy as np
from threadpoolctl import threadpool_limits
from modal_comparison_contract import Comparison,vector_norm_bound
from coupled_interval_model import mm,magnitude,inverse_enclosure,up


def run(coordinate):
    cc=Comparison(coordinate);out=cc.out
    syn=json.loads((out/"block_modal_synthesis.json").read_text(encoding="utf-8"))
    saved=np.load(out/"block_modal_certificate.npz")
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    assert syn["model_sha256"]==sha(out/"model.toml")
    assert syn["center_sha256"]==sha(out/"refined_center.json")
    assert syn["implementation_sha256"]==sha(Path(cc.module.__file__))
    assert syn["script_sha256"]==sha(Path(__file__).with_name("modal_comparison_contract.py"))
    cc.H=saved["H"];cc.groups=syn["groups"];cc.ng=len(cc.groups)
    cc.Ti,iv=inverse_enclosure((cc.H,cc.H))
    cc.U=np.column_stack([up(np.sqrt(up(np.sum(magnitude(cc.Ti)[:,g]**2,axis=1)*(1+8*len(g)*np.finfo(float).eps)))) for g in cc.groups])
    fi=(np.array(cc.center["f0lo"])[:,None],np.array(cc.center["f0hi"])[:,None])
    drift=mm((cc.H,cc.H),fi)
    cc.drift=np.array([vector_norm_bound(magnitude(drift)[g]) for g in cc.groups])
    initial=mm((cc.H,cc.H),(-cc.model.center[:,None],-cc.model.center[:,None]))
    cc.initial=np.array([vector_norm_bound(magnitude(initial)[g]) for g in cc.groups])
    r=saved["radii"];epsilon=float(saved["epsilon"])
    result=cc.evaluate(r,epsilon)
    a,M,b=result["a"],result["M"],result["b"]
    margins=[];perron=[]
    for i in range(cc.ng):
        interaction=sum(F(float(M[i,j]))*F(float(r[j])) for j in range(cc.ng))
        ar=F(float(a[i]))*F(float(r[i]))
        margin=ar-interaction-F(float(b[i]))*F(epsilon)-F(float(cc.drift[i]))
        assert margin>0,(i,float(margin))
        assert F(float(r[i]))>=F(float(cc.initial[i]))
        margins.append(margin);perron.append(interaction/ar)
    assert result["output_pass"]
    perron_bound=max(perron);assert perron_bound<1
    report=dict(status="BLOCK_MODAL_NONLINEAR_REGION_RECHECK_PASSED",coordinate=coordinate,
        state_count=cc.model.nx,input_count=cc.model.nw,blocks=cc.ng,epsilon=epsilon,
        min_boundary_margin=float(min(margins)),exact_min_boundary_margin=str(min(margins)),
        certified_comparison_Perron_upper=float(up(float(perron_bound))),
        max_initial_region_usage=float(np.max(cc.initial/r)),
        frequency_bound=result["frequency"],rocof_bound=result["rocof"],
        voltage_min=result["voltage_min"],voltage_max=result["voltage_max"],
        min_block_decay=float(a.min()),max_block_decay=float(a.max()),
        metric_inverse_check=iv,network_inverse_check=result["inverse_network_check"],
        implementation_type="Nonlinear interval model enclosure; rational check of all boundary inequalities; shares interval model with synthesis",
        source_hashes={n:sha(out/n) for n in ("model.toml","refined_center.json","block_modal_synthesis.json","block_modal_certificate.npz")},
        implementation_sha256=sha(Path(cc.module.__file__)),script_sha256=sha(Path(__file__)),
        scope="Fixed design and declared current/DC ball. Not global optimality, useful operating-capacity proof, or independently validated hardware limits.")
    (out/"block_modal_verification.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    # Plain CSV bridge for an independent Julia evaluation of the nonlinear model.
    np.savetxt(out/"validation_H.csv",cc.H,delimiter=",")
    np.savetxt(out/"validation_T.csv",np.linalg.inv(cc.H),delimiter=",")
    np.savetxt(out/"validation_center.csv",cc.model.center,delimiter=",")
    np.savetxt(out/"validation_radii.csv",r,delimiter=",")
    group_of=np.zeros(cc.model.nx,dtype=int)
    for k,g in enumerate(cc.groups):group_of[g]=k+1
    np.savetxt(out/"validation_groups.csv",group_of,delimiter=",",fmt="%d")
    np.savetxt(out/"validation_epsilon.csv",[epsilon],delimiter=",")
    print(json.dumps({k:v for k,v in report.items() if k not in ("source_hashes","exact_min_boundary_margin")},indent=2))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--coordinates",choices=("dq","ri"),default="dq")
    args=p.parse_args()
    with threadpool_limits(limits=1):run(args.coordinates)
