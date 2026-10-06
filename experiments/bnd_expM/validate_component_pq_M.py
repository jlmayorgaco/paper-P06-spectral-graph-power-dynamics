"""Check SG and GFL output currents against initialized bus and ZIP power."""
from __future__ import annotations

import csv
import re
from pathlib import Path


ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"reports"/"experiment_M"


def read(path:Path)->list[dict]:
    with path.open(newline="",encoding="utf-8") as fh:return list(csv.DictReader(fh))


def main()->None:
    rows=[]
    for case in ("all_SG","ExpG_candidate","ExpK_nominal"):
        p=OUT/"matrices"/case
        states={r["state_name"]:float(r["equilibrium_value"]) for r in read(p/"PD_state_map.csv")}
        params={r["parameter"]:float(r["value"]) for r in read(p/"PD_parameter_values.csv")}
        for share in read(p/"PD_share_audit.csv"):
            bus=int(share["bus"])
            def v(name:str)->float:return states[f"VIndex({bus}, :{name})"]
            ur,ui=v("busbar₊u_r"),v("busbar₊u_i")
            keyr=f"VIndex({bus}, :gfl₊filter₊i_f_r)"
            keyi=f"VIndex({bus}, :gfl₊filter₊i_f_i)"
            if keyr in states and keyi in states:
                scale=params.get(f"VIndex({bus}, :gfl₊filter₊PortScale)",1.0)
                ir,ii=scale*states[keyr],scale*states[keyi]
                gflp=100*(ur*ir+ui*ii)
                gflq=100*(ui*ir-ur*ii)
            else:
                scale=0.0;gflp=0.0;gflq=0.0
            sgp=float(share["sg_p_MW"]);sgq=float(share["sg_q_Mvar"])
            zip_p=float(share["zip_p_MW"]);zip_q=float(share["zip_q_Mvar"])
            netp=float(share["net_p_MW"]);netq=float(share["net_q_Mvar"])
            rows.append(dict(case=case,bus=bus,sg_direct_p_MW=sgp,sg_direct_q_Mvar=sgq,
                gfl_direct_p_MW=gflp,gfl_direct_q_Mvar=gflq,gfl_port_scale=scale,
                zip_observed_p_MW=zip_p,zip_observed_q_Mvar=zip_q,
                net_observed_p_MW=netp,net_observed_q_Mvar=netq,
                p_balance_MW=sgp+gflp+zip_p-netp,
                q_balance_Mvar=sgq+gflq+zip_q-netq,
                expected_frozen_sg_p_MW=share["expected_analytic_sg_p_MW"],
                sg_dispatch_difference_MW=share["sg_p_difference_MW"]))
    path=OUT/"tables"/"TABLE_M04_direct_device_PQ_balance.csv"
    with path.open("w",newline="",encoding="utf-8") as fh:
        writer=csv.DictWriter(fh,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print("rows",len(rows),"max P residual",max(abs(r["p_balance_MW"]) for r in rows),
          "max Q residual",max(abs(r["q_balance_Mvar"]) for r in rows))
    for row in rows:
        if row["case"]=="ExpK_nominal" and row["bus"] in (38,39):print(row)


if __name__=="__main__":main()
