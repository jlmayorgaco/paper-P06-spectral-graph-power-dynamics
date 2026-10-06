import hashlib
import json
import platform
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCES = [
    "Project.toml", "Manifest.toml",
    "reports/analytic_iteration_20261001/refined/protocol.toml",
    "reports/analytic_iteration_20261001/refined/joint_final.toml",
    "reports/experiment_N/TABLE_N01_original_operating_point.csv",
    "reports/experiment_N/TABLE_N01_original_operating_point.csv.sha256",
    "src/bnd_model_expN/PDExactDesignN.jl",
    "src/bnd_design_e/CollectiveModel.jl",
    "experiments/nonlinear_codesign_20261001/PDPhysicalReference.jl",
    "experiments/nonlinear_codesign_20261001/ReducedDAE.jl",
    "experiments/delay_dressed_replacement_frontier_20261002/DelayCharacteristic.jl",
    "experiments/delay_dressed_replacement_frontier_20261002/DELAY_MODEL_CONTRACT.md",
    "experiments/delay_dressed_replacement_frontier_20261002/physical_holdout_reproduction.jl",
    "experiments/codesign_validation_20261001/validate_pd.jl",
    "experiments/bnd_expQ2B/certified_search/ConsistentEventObservation.jl",
]

def run(*args):
    return subprocess.check_output(args, cwd=ROOT, text=True, stderr=subprocess.STDOUT).strip()

hashes = {}
for rel in SOURCES:
    p = ROOT / rel
    hashes[rel] = hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None

try:
    head = run("git", "rev-parse", "HEAD")
    branch = run("git", "branch", "--show-current")
    porcelain = run("git", "status", "--porcelain", "--untracked-files=normal").splitlines()
    preexisting = [line for line in porcelain if HERE.relative_to(ROOT).as_posix() not in line.replace("\\", "/")]
except Exception as e:
    head, branch, preexisting = None, None, [f"status_read_error: {e}"]

try:
    julia_version = run("julia", "--startup-file=no", "-e", "print(VERSION)")
except Exception as e:
    julia_version = f"unavailable: {e}"

manifest = {
    "experiment": "IEEE-39 exact action-space SG-to-GFL co-design",
    "experiment_date": "2026-10-02",
    "git_head": head,
    "git_branch": branch,
    "repository_dirty_before_this_experiment": bool(preexisting),
    "preexisting_status_entry_count_excluding_this_directory": len(preexisting),
    "all_new_artifacts_directory": HERE.relative_to(ROOT).as_posix(),
    "commits_or_pushes": False,
    "julia_version": julia_version,
    "python_version": platform.python_version(),
    "packages": {
        "PowerDynamics": "5.0.0 (repository environment, verified in prior frozen run)",
        "NetworkDynamics": "1.3.0 (repository environment, verified in prior frozen run)",
        "SciMLBase": "3.53.2 (repository environment, verified in prior frozen run)",
    },
    "frozen_events": [
        {"event": "bus8_minus100", "bus": 8, "load_step_MW": -100.0},
        {"event": "bus16_plus100", "bus": 16, "load_step_MW": 100.0},
        {"event": "bus16_minus100", "bus": 16, "load_step_MW": -100.0},
        {"event": "bus29_plus100", "bus": 29, "load_step_MW": 100.0},
        {"event": "bus29_minus100", "bus": 29, "load_step_MW": -100.0},
    ],
    "source_sha256": hashes,
    "experiment_artifact_sha256": {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(HERE.iterdir())
        if p.is_file() and p.name != "EXPERIMENT_MANIFEST.json"
    },
    "gate_decisions": {
        "Q0": "PASS: fully validated zero-delay feasibility witness; not optimized",
        "Q1_Q2": "PASS at stated source-identity and numerical-reconstruction scope",
        "Q3": "INDETERMINATE at 20 and 40 ms; positive-delay winding indices are unverified diagnostics",
        "Q4": "conditional target-root equation validated; nearest-pole boundary not established",
        "Q5": "rank-two conditional law validated; no physical roots in tested grid",
        "Q6": "fixed-coordinate action-space factorization validated; rank bound <=30",
        "Q7_Q8": "NOT RUN because the Q3 gate did not pass",
    },
    "note": "Repository has substantial pre-existing dirty/untracked content. New experiment outputs are confined to this directory; hashes identify the read-only source snapshot used here.",
}
(HERE / "EXPERIMENT_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
print(f"wrote {HERE / 'EXPERIMENT_MANIFEST.json'}; source hashes={sum(v is not None for v in hashes.values())}/{len(hashes)}")
