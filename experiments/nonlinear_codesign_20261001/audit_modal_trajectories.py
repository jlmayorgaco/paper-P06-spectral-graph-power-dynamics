"""Audit finite Julia runs and tolerance sensitivity separately from the proof."""
from pathlib import Path
import hashlib
import json
import tomllib
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"reports/nonlinear_codesign_20261001/coupled_dq_contract"
E=ROOT/"experiments/nonlinear_codesign_20261001"
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
model=tomllib.loads((OUT/"model.toml").read_text())
for name,value in model["source_hashes"].items():assert sha(ROOT/name)==value,name
records=[]
for tag in ("standard","refined"):
    record=tomllib.loads((OUT/f"nonlinear_validation_julia_{tag}.toml").read_text())
    assert record["status"]=="THREE_TRAJECTORY_CHECKS_PASSED"
    assert len(record["cases"])==3 and all(r["complete"] and r["passed"] for r in record["cases"])
    for key,name in (("script_sha256","validate_modal_comparison_julia.jl"),("reference_sha256","CoupledDQReference.jl"),("jacobian_sha256","CoupledDQJacobian.jl")):
        assert record[key]==sha(E/name)
    assert record["certificate_sha256"]==sha(OUT/"block_modal_certificate.npz")
    records.append(record)
cert=json.loads((OUT/"block_modal_verification.json").read_text())
scale=np.array([1.,cert["frequency_bound"],cert["rocof_bound"],.1,.1])
rows=[]
for case in records[0]["cases"]:
    name=case["name"]
    first=np.loadtxt(OUT/f"trajectory_standard_{name}.csv",delimiter=",")
    second=np.loadtxt(OUT/f"trajectory_refined_{name}.csv",delimiter=",")
    assert first.shape==second.shape and np.array_equal(first[:,0],second[:,0])
    assert first[-1,0]==second[-1,0]==20.
    assert np.max(np.abs(first[:,6]-.8))<1e-12 and np.max(np.abs(second[:,6]-.8))<1e-12
    diff=np.max(np.abs(first[:,1:6]-second[:,1:6]),axis=0)
    normalized=diff/scale
    assert normalized.max()<1e-3,(name,normalized)
    rows.append(dict(case=name,samples=len(first),maximum_absolute_difference=diff.tolist(),
                     maximum_scaled_difference=normalized.tolist()))
report=dict(status="SIX_FINITE_JULIA_RUNS_AND_SAMPLED_OUTPUT_REFINEMENT_CHECK_PASSED",
    tolerance_pair=[r["normalized_coordinate_abstol_reltol"] for r in records],
    horizon_seconds=20,epsilon=cert["epsilon"],columns=["region_usage","max_abs_frequency","max_abs_rocof","voltage_min","voltage_max"],
    criterion="Each sampled-output difference divided by [1, certified f bound, certified RoCoF bound, 0.1 pu, 0.1 pu] is below 1e-3.",
    cases=rows,model_source_hashes_verified=model["source_hashes"],
    scope="Finite sample checks on the pre-LP modal region, not exhaustive time/input validation or a trajectory discretization error proof. The all-input claim rests on the separate nonlinear certificate.",
    source_hashes={p.name:sha(p) for p in OUT.glob("trajectory_*.csv")},script_sha256=sha(Path(__file__)))
for tag in ("standard","refined"):
    p=OUT/f"nonlinear_validation_julia_{tag}.toml";report["source_hashes"][p.name]=sha(p)
(OUT/"trajectory_refinement_audit.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
print(json.dumps({k:v for k,v in report.items() if k not in ("source_hashes","model_source_hashes_verified")},indent=2))
