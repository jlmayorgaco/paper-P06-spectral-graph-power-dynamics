"""Freeze the zero-delay search direction from seed toward the prior joint design."""

from pathlib import Path
import hashlib
import json
import math
import tomllib

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SEED = tomllib.loads((HERE / "seed_uniform_875.toml").read_text(encoding="utf-8"))
JOINT_PATH = ROOT / "reports/analytic_iteration_20261001/refined/joint_final.toml"
JOINT = tomllib.loads(JOINT_PATH.read_text(encoding="utf-8"))
GRID = (0.25, 0.50, 0.75, 1.00)


def candidate(eta):
    rho = [a + eta * (b-a) for a,b in zip(SEED["rho"], JOINT["rho"])]
    kp = [math.exp(math.log(a) + eta*(math.log(b)-math.log(a))) for a,b in zip(SEED["Kp"], JOINT["Kp"])]
    ki = [math.exp(math.log(a) + eta*(math.log(b)-math.log(a))) for a,b in zip(SEED["Ki"], JOINT["Ki"])]
    return {"rho": rho, "Kp": kp, "Ki": ki, "dc_convention": "physical_supply", "eta": eta}


out = HERE / "m1_candidates"
out.mkdir(exist_ok=True)
manifest = {
    "search_contract": "M1_SEARCH_CONTRACT.md",
    "seed_sha256": hashlib.sha256((HERE / "seed_uniform_875.toml").read_bytes()).hexdigest(),
    "joint_source": JOINT_PATH.relative_to(ROOT).as_posix(),
    "joint_sha256": hashlib.sha256(JOINT_PATH.read_bytes()).hexdigest(),
    "coarse_eta_grid": GRID,
    "candidates": {},
}
for eta in GRID:
    d = candidate(eta)
    label = f"eta_{int(round(100*eta)):03d}"
    path = out / f"{label}.toml"
    content = "\n".join(f"{key} = {json.dumps(val)}" for key,val in d.items()) + "\n"
    path.write_text(content, encoding="utf-8")
    manifest["candidates"][label] = {
        "eta": eta,
        "file": path.relative_to(HERE).as_posix(),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }
(HERE / "M1_FROZEN_CANDIDATES.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
print("frozen", len(GRID), "candidate designs before evaluation")
