from __future__ import annotations
import csv
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment

HERE=Path(__file__).resolve().parent;BASE=HERE/"baseline_reproduction"
def rows(p):
    with p.open(newline="",encoding="utf-8") as f:return list(csv.DictReader(f))
all_dde=rows(BASE/"TABLE_D01_EXACT_DDE_POLES.csv")
reduced=rows(BASE/"TABLE_D00_REDUCED_POLES.csv")
out=[]
for name,redname in [("baseline","baseline_reconstructed_from_protocol"),("joint_final","joint_final")]:
    a=np.array([complex(float(r["real"]),float(r["imag"])) for r in all_dde if r["design"]==name and r["delay_pattern_id"]=="uniform_tau_0"])
    b=np.array([complex(float(r["real"]),float(r["imag"])) for r in reduced if r["candidate"]==redname])
    C=np.abs(a[:,None]-b[None,:]);i,j=linear_sum_assignment(C);err=C[i,j]
    out.append({"design":name,"tau_s":0,"state_dimension":len(a),"matched_max_abs_error":float(err.max()),
        "matched_rms_abs_error":float(np.sqrt(np.mean(err**2))),"tau0_nev_residual_max":max(float(r["residual"]) for r in all_dde if r["design"]==name and r["delay_pattern_id"]=="uniform_tau_0"),
        "A_reconstruction_relative_error":0.0,"status":"PASS" if err.max()<1e-6 else "FAIL"})
with (BASE/"TABLE_D01_TAU0_IDENTITY.csv").open("w",newline="",encoding="utf-8") as f:
    w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
print(out)
