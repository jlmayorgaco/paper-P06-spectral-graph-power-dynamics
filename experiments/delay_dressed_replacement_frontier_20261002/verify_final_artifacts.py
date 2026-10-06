"""Integrity and gate-consistency checks for the final handoff directory."""
from pathlib import Path
import csv
import hashlib
import json
import pandas as pd
from PIL import Image

HERE = Path(__file__).resolve().parent
required = [
    "README.md","STATUS.md","MODEL_PROVENANCE.md","DELAY_MODEL_CONTRACT.md",
    "THEORY_DELAY_REPLACEMENT.md","POSTER_CLAIMS.md","EXPERIMENT_MANIFEST.json","REPRODUCE.md",
    *[f"TABLE_D{i:02d}_{name}.csv" for i,name in [
        (0,"BASELINE_REPRODUCTION"),(1,"EXACT_DDE_POLES"),(2,"DERIVATIVE_VALIDATION"),
        (3,"TAYLOR_NETWORK_VALIDATION"),(4,"GRAPH_DAMPING_OPERATORS"),(5,"GSP_DELAY_METRICS"),
        (6,"GSP_GAIN_REDUCTION"),(7,"REPLACEMENT_FRONTIER"),(8,"DELAY_SHADOW_PRICES"),
        (9,"NONLINEAR_DDE_VALIDATION"),(10,"OPTIMALITY_BOUNDS")]],
    *[f"FIG_D{i:02d}_{name}.png" for i,name in enumerate([
        "THEORY_PIPELINE","ANALYTIC_VS_FD","OPERATOR_APPROXIMATION_ERROR",
        "REPLACEMENT_VS_GRAPH_MODES","REPLACEMENT_FRONTIER_VS_DELAY",
        "SAME_DELAY_DIFFERENT_PLACEMENT","DELAY_SHADOW_PRICE_MAP",
        "FULL_MODEL_VALIDATION","OPTIMALITY_GAP"],1)],
]
missing = [name for name in required if not (HERE/name).is_file()]
assert not missing, f"Missing required artifacts: {missing}"
for path in HERE.glob("*.json"):
    json.loads(path.read_text(encoding="utf-8"))
for path in HERE.glob("TABLE_D*.csv"):
    with path.open(newline="",encoding="utf-8") as f:
        list(csv.DictReader(f))
for path in HERE.glob("FIG_D*.png"):
    with Image.open(path) as im:
        im.verify()

deriv = pd.read_csv(HERE/"TABLE_D02_DERIVATIVE_VALIDATION.csv")
assert len(deriv)==27 and deriv.relative_error.max()<=0.01
taylor = pd.read_csv(HERE/"TABLE_D03_TAYLOR_NETWORK_VALIDATION.csv")
assert len(taylor)==60 and taylor.third_order_relative_error.max()<1e-5
schur = pd.read_csv(HERE/"TABLE_D04_GRAPH_DAMPING_OPERATORS.csv")
assert len(schur)==112 and set(schur.operator_status)=={"BLOCKED_SCHUR_SINGULAR_AFTER_GAUGE_DEFLATION"}
gsp = pd.read_csv(HERE/"TABLE_D05_GSP_DELAY_METRICS.csv")
assert len(gsp)==112 and gsp.commutator_identity_relative_error.max()<1e-12
frontier = pd.read_csv(HERE/"TABLE_D07_REPLACEMENT_FRONTIER.csv")
assert len(frontier)==1 and frontier.status.iloc[0]=="BLOCKED_EXACT_DDE_SPECTRUM"

hash_manifest = json.loads((HERE/"EXPERIMENT_ARTIFACT_HASHES.json").read_text(encoding="utf-8"))
for item in hash_manifest["artifacts"]:
    path = HERE/item["path"]
    assert path.stat().st_size==item["size_bytes"]
    assert hashlib.sha256(path.read_bytes()).hexdigest()==item["sha256"], f"Hash mismatch: {item['path']}"
print(f"PASS: {len(required)} required outputs; 27 sensitivities; 60 Taylor checks; "
      f"112 blocked Schur rows; graph-commutator identity; {len(hash_manifest['artifacts'])} SHA-256 entries")
