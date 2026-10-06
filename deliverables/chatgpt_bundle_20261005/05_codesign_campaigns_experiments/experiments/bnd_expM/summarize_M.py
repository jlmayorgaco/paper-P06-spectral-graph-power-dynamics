"""Assemble Experiment M evidence without imputing unmeasured gates."""
from __future__ import annotations

import csv
from pathlib import Path


ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"reports"/"experiment_M"
TABLES=OUT/"tables"
CASES=("all_SG","ExpG_candidate","ExpK_nominal")


def read(path:Path)->list[dict]:
    if not path.is_file():return []
    with path.open(newline="",encoding="utf-8") as fh:return list(csv.DictReader(fh))


def write(name:str,rows:list[dict])->None:
    if not rows:return
    TABLES.mkdir(parents=True,exist_ok=True)
    with (TABLES/name).open("w",newline="",encoding="utf-8") as fh:
        writer=csv.DictWriter(fh,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def main()->None:
    pq=[];params=[];dims=[];closure=[];alignment=[]
    for case in CASES:
        p=OUT/"matrices"/case
        pq+=read(p/"component_PQ_sharing.csv")
        dims+=read(p/"component_dimensions.csv")
        alignment+=read(TABLES/f"TABLE_M10_state_alignment_{case}.csv")
        for row in read(p/"PD_parameter_values.csv"):
            params.append(dict(case=case,parameter=row["parameter"],initialized_value=row["value"],
                source="PowerDynamics_NWState_pflat"))
        poles=read(p/"PD_poles.csv")
        cpoles=read(p/"PD_closure_positive_poles.csv")
        if poles and cpoles:
            import numpy as np
            from scipy.optimize import linear_sum_assignment
            direct=np.array([float(r["real"])+1j*float(r["imag"]) for r in poles])
            rebuilt=np.array([float(r["real"])+1j*float(r["imag"]) for r in cpoles])
            if len(direct)==len(rebuilt):
                row,col=linear_sum_assignment(np.abs(direct[:,None]-rebuilt[None,:]))
                errors=np.abs(direct[row]-rebuilt[col])
                closure.append(dict(case=case,status="EVALUATED_POSITIVE_FEEDBACK",PD_poles=len(direct),
                    closure_poles=len(rebuilt),max_pole_error=float(errors.max()),
                    median_pole_error=float(np.median(errors)),pass_1e_8=bool(errors.max()<1e-8)))
            else:
                closure.append(dict(case=case,status="DIMENSION_MISMATCH",PD_poles=len(direct),
                    closure_poles=len(rebuilt),max_pole_error="",median_pole_error="",pass_1e_8=False))
    write("TABLE_M04_component_PQ_sharing.csv",pq)
    write("TABLE_M06_component_parameter_dump.csv",params)
    write("TABLE_M07_component_dimensions.csv",dims)
    write("TABLE_M10_state_alignment.csv",alignment)
    write("TABLE_M08_PD_open_loop_closure_identity.csv",closure)
    write("TABLE_M13_base_frequency_ablation.csv",[dict(status="NOT_APPLICABLE",
        reason="Measured original and replacement components all store 60 Hz; no 50 Hz mismatch was found.",
        alpha_50="",alpha_60="",explains_mismatch="NO")])
    write("TABLE_M18_parameter_identification.csv",[dict(status="NOT_USED",
        reason="Installed parameters are exposed; grey-box fitting is not justified.",
        parameters_estimated="",identifiability="NOT_APPLICABLE")])
    print("assembled",len(pq),"P/Q rows,",len(params),"actual parameters,",len(dims),"component dimensions")


if __name__=="__main__":main()
