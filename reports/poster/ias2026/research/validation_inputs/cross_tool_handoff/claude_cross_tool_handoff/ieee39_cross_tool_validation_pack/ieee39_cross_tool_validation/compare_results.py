
from pathlib import Path
import json, sys

root=Path(__file__).resolve().parent
r=root/"results"
def load(name):
    p=r/name
    return json.loads(p.read_text()) if p.exists() else None

py=load("python_reference.json")
pp=load("pandapower_parity.json")
an=load("andes_dynamic_summary.json")

rows=[]
def add(check, value, gate, ok):
    rows.append({"check":check,"measured":value,"gate":gate,"verdict":"PASS" if ok else "FAIL"})

if py:
    add("Python PF mismatch",py["pf_max_mismatch"],1e-10,py["pf_max_mismatch"]<1e-10)
    add("Python vs ANDES Ybus",py["ybus_max_abs_error"],1e-9,py["ybus_max_abs_error"]<1e-9)
    add("Python vs ANDES Vm",py["vm_max_abs_error_pu"],1e-5,py["vm_max_abs_error_pu"]<1e-5)
    add("Python vs ANDES Va",py["va_max_gauge_aligned_error_rad"],1e-5,py["va_max_gauge_aligned_error_rad"]<1e-5)
if pp:
    add("pandapower vs Python Vm",pp["vm_max_abs_error_pu"],1e-5,pp["vm_max_abs_error_pu"]<1e-5)
    add("pandapower vs Python Va",pp["va_max_gauge_aligned_error_rad"],1e-5,pp["va_max_gauge_aligned_error_rad"]<1e-5)
    add("pandapower branch terminal flows",pp["canonical_branch_terminal_flow_max_error_MVA"],1e-3,pp["canonical_branch_terminal_flow_max_error_MVA"]<1e-3)
if an:
    add("ANDES dynamic alpha",an["R3_max_alpha_error"],1e-4,an["R3_max_alpha_error"]<1e-4)
    add("ANDES dynamic frequency",an["R3_max_frequency_error_hz"],1e-3,an["R3_max_frequency_error_hz"]<1e-3)
    add("ANDES dynamic RHP count",an["R3_all_RHP_counts_agree"],True,bool(an["R3_all_RHP_counts_agree"]))

import pandas as pd
df=pd.DataFrame(rows)
df.to_csv(r/"CROSS_TOOL_EVIDENCE_TABLE.csv",index=False)
print(df.to_string(index=False))
if pp is None or an is None:
    print("\nINCOMPLETE: native pandapower/ANDES result files are missing.")
    sys.exit(2)
sys.exit(0 if (df.verdict=="PASS").all() else 1)
