"""Predefine nodal rho path with PLL gains fixed at stored nominal seed values."""
from pathlib import Path
import hashlib, json, tomllib

HERE=Path(__file__).resolve().parent
MEGA=HERE.parent/"analytical_delay_codesign_mega_20261002"
seed=tomllib.loads((MEGA/"seed_uniform_875.toml").read_text(encoding="utf-8"))
best=tomllib.loads((MEGA/"M1_ZERO_DELAY_DESIGN.toml").read_text(encoding="utf-8"))
low,high=87.5,float(best["GFL_percent"])
levels=(87.5,87.75,88.0,88.25,high)
out=HERE/"fixed_gain_designs";out.mkdir(exist_ok=True)
manifest=[]
for j,target in enumerate(levels):
    eta=(target-low)/(high-low)
    rho=[a+eta*(b-a) for a,b in zip(seed["rho"],best["rho"])]
    label=f"rho_{j:02d}"
    path=out/f"{label}.toml"
    data={"rho":rho,"Kp":seed["Kp"],"Ki":seed["Ki"],
          "target_GFL_percent":target,"path_eta":eta,"label":label}
    path.write_text("\n".join(f"{k} = {json.dumps(v)}" for k,v in data.items())+"\n",encoding="utf-8")
    manifest.append({"label":label,"target_GFL_percent":target,"path_eta":eta,
                     "file":path.relative_to(HERE).as_posix(),
                     "sha256":hashlib.sha256(path.read_bytes()).hexdigest()})
(HERE/"FIXED_GAIN_DESIGNS_FROZEN.json").write_text(json.dumps({
    "description":"rho line from stored 87.5% seed to stored 88.455% zero-delay design; all Kp/Ki frozen at seed values",
    "designs":manifest},indent=2)+"\n",encoding="utf-8")
print("FROZEN_FIXED_GAIN",len(manifest),"designs")
