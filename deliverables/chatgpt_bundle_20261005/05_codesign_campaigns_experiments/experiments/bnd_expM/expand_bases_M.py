"""Audit every initialized component base, including bound subdevice bases."""
from __future__ import annotations

import csv
import math
import re
from pathlib import Path


ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"reports"/"experiment_M"
TABLES=OUT/"tables"
SOURCE=Path.home()/".julia"/"packages"/"PowerDynamics"/"VzOiZ"/"docs"/"examples"/"ieee39data"
CASES=("all_SG","ExpG_candidate","ExpK_nominal")


def read(path:Path)->list[dict]:
    with path.open(newline="",encoding="utf-8") as fh:return list(csv.DictReader(fh))


def main()->None:
    buses={int(r["bus"]):r for r in read(SOURCE/"bus.csv")}
    lines=read(SOURCE/"branch.csv")
    definitions={(r["case"],int(r["bus"])):r for r in read(TABLES/"TABLE_M05_case_definitions.csv")}
    rows=[]
    for case in CASES:
        params=read(OUT/"matrices"/case/"PD_parameter_values.csv")
        bycomponent={}
        for row in params:
            match=re.fullmatch(r"([VE])Index\((\d+), :(.*)\)",row["parameter"])
            if match:
                bycomponent.setdefault((match[1],int(match[2])),{})[match[3]]=float(row["value"])
        for bus in range(1,40):
            d=bycomponent[("V",bus)]
            s=d["systembase₊Sbase"]; w=d["systembase₊ωbase"]
            frame=d["systembase₊ωframe"];v=d["busbar₊Vbase"]
            expected_v=float(buses[bus]["base_kv"])
            good=math.isclose(s,100,abs_tol=1e-12) and math.isclose(w,2*math.pi*60,abs_tol=1e-10) and math.isclose(frame,1,abs_tol=1e-12) and math.isclose(v,expected_v,abs_tol=1e-12)
            kinds=[("bus",buses[bus]["category"])]
            if buses[bus]["has_load"]=="true":kinds.append(("ZIPLoad","bound_to_bus_systembase"))
            if bus>=30:
                architecture=definitions[(case,bus)]["architecture"]
                if "SG" in architecture:kinds.append(("SG","bound_to_bus_systembase"))
                if "GFL" in architecture:kinds.append(("GFL","bound_to_bus_systembase"))
            for name,kind in kinds:
                rows.append(dict(case=case,bus=bus,component=f"bus{bus}_{name}",component_type=kind,
                    Sbase=s,fbase_equivalent=w/(2*math.pi),omega_base=w,omega_frame=frame,
                    Vbase=v,src_Vbase="",dst_Vbase="",expected_Sbase=100,
                    expected_omega_base=2*math.pi*60,expected_Vbase=expected_v,
                    base_scope="compiled_bus" if name=="bus" else "device_bound_to_compiled_bus",
                    pass_=good))
        for edge,row in enumerate(lines,1):
            d=bycomponent[("E",edge)]
            s=d["systembase₊Sbase"];w=d["systembase₊ωbase"]
            frame=d["systembase₊ωframe"]
            sv=d["src₊Vbase"];dv=d["dst₊Vbase"]
            expected_s=float(buses[int(row["src_bus"])]["base_kv"])
            expected_d=float(buses[int(row["dst_bus"])]["base_kv"])
            good=all((math.isclose(s,100,abs_tol=1e-12),math.isclose(w,2*math.pi*60,abs_tol=1e-10),
                math.isclose(frame,1,abs_tol=1e-12),math.isclose(sv,expected_s,abs_tol=1e-12),
                math.isclose(dv,expected_d,abs_tol=1e-12)))
            rows.append(dict(case=case,bus="",component=f"line{edge}",component_type="PiLine_fault",
                Sbase=s,fbase_equivalent=w/(2*math.pi),omega_base=w,omega_frame=frame,
                Vbase="",src_Vbase=sv,dst_Vbase=dv,expected_Sbase=100,
                expected_omega_base=2*math.pi*60,expected_Vbase="",
                base_scope="compiled_line",pass_=good))
    old=read(TABLES/"TABLE_M01_component_bases.csv")
    for row in old:
        if row["case"] in ("GFL_template","mixed_bus38"):
            rows.append(dict(case=row["case"],bus=row["bus"],component=row["component"],
                component_type=row["component_type"],Sbase=row["Sbase"],
                fbase_equivalent=row["fbase_equivalent"],omega_base=row["omega_base"],
                omega_frame=row["omega_frame"],Vbase=row["Vbase"],src_Vbase="",dst_Vbase="",
                expected_Sbase=row["expected_Sbase"],expected_omega_base=row["expected_omega_base"],
                expected_Vbase=row["expected_Vbase"],base_scope="compiled_template",pass_=row["pass"].lower()=="true"))
    with (TABLES/"TABLE_M01_component_bases.csv").open("w",newline="",encoding="utf-8") as fh:
        writer=csv.DictWriter(fh,fieldnames=["pass" if k=="pass_" else k for k in rows[0]])
        writer.writeheader()
        for row in rows:
            writer.writerow({("pass" if k=="pass_" else k):v for k,v in row.items()})
    print("component base rows",len(rows),"pass",sum(bool(r["pass_"]) for r in rows))


if __name__=="__main__":main()
