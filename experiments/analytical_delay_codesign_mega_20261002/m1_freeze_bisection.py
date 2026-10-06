"""Freeze one preregistered interpolation midpoint after the coarse grid."""

from pathlib import Path
import hashlib
import json
import math
import tomllib
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
audit = pd.read_csv(HERE / "M1_CANDIDATE_AUDIT.csv")
coarse = audit[audit.candidate_id.isin(["eta_000_seed","eta_025","eta_050","eta_075","eta_100"])]
assert len(coarse) == 5 and (coarse.status != "NOT_EXECUTED").all(), "coarse grid must finish first"
good = coarse[coarse.fully_validated]
eta_good = float(good.eta.max())
bad_higher = coarse[(coarse.eta > eta_good) & (~coarse.fully_validated)]
if bad_higher.empty:
    print("No higher infeasible neighbor; no bisection candidate")
    raise SystemExit(0)
eta_bad = float(bad_higher.eta.min())
eta = (eta_good+eta_bad)/2
seed = tomllib.loads((HERE/"seed_uniform_875.toml").read_text(encoding="utf-8"))
joint_path = ROOT/"reports/analytic_iteration_20261001/refined/joint_final.toml"
joint = tomllib.loads(joint_path.read_text(encoding="utf-8"))
rho = [a+eta*(b-a) for a,b in zip(seed["rho"],joint["rho"])]
kp = [math.exp(math.log(a)+eta*(math.log(b)-math.log(a))) for a,b in zip(seed["Kp"],joint["Kp"])]
ki = [math.exp(math.log(a)+eta*(math.log(b)-math.log(a))) for a,b in zip(seed["Ki"],joint["Ki"])]
label = f"eta_{int(round(1000*eta)):04d}_bisect"
path = HERE/"m1_candidates"/f"{label}.toml"
values = {"rho":rho,"Kp":kp,"Ki":ki,"dc_convention":"physical_supply","eta":eta}
path.write_text("\n".join(f"{k} = {json.dumps(v)}" for k,v in values.items())+"\n",encoding="utf-8")
record = {"rule":"single midpoint between highest coarse feasible eta and nearest higher coarse infeasible eta",
          "eta_feasible":eta_good,"eta_infeasible":eta_bad,"eta_midpoint":eta,
          "candidate_file":path.relative_to(HERE).as_posix(),
          "candidate_sha256":hashlib.sha256(path.read_bytes()).hexdigest()}
(HERE/"M1_BISECTION_FREEZE.json").write_text(json.dumps(record,indent=2)+"\n",encoding="utf-8")
print("M1_BISECTION_FROZEN",label,eta)
