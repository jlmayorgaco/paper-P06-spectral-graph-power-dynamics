"""Record source snapshot and the preregistered Mega experiment contract."""

import hashlib
import json
import platform
import subprocess
import tomllib
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCES = [
    "Project.toml",
    "Manifest.toml",
    "src/bnd_model_expN/PDExactDesignN.jl",
    "experiments/nonlinear_codesign_20261001/PDPhysicalReference.jl",
    "experiments/bnd_expQ2B/certified_search/ConsistentEventObservation.jl",
    "reports/analytic_iteration_20261001/refined/protocol.toml",
    "reports/analytic_iteration_20261001/refined/joint_final.toml",
    "experiments/analytic_iteration_20261001/run_experiment.jl",
    "experiments/delay_dressed_replacement_frontier_20261002/DelayCharacteristic.jl",
    "experiments/delay_dressed_replacement_frontier_20261002/DDE_NEV.jl",
    "experiments/delay_dressed_replacement_frontier_20261002/FROZEN_DELAY_PATTERNS.json",
    "experiments/ieee39_exact_action_space_codesign_20261002/q0_validate_five_events.jl",
    "experiments/ieee39_exact_action_space_codesign_20261002/q12_validate_low_rank.jl",
    "experiments/ieee39_exact_action_space_codesign_20261002/q4_validate_pi_boundary.jl",
    "experiments/ieee39_exact_action_space_codesign_20261002/q56_validate_action_space.jl",
]


def run(*args):
    return subprocess.check_output(args, cwd=ROOT, text=True, stderr=subprocess.STDOUT).strip()


manifest_data = tomllib.loads((ROOT / "Manifest.toml").read_text(encoding="utf-8"))
deps = manifest_data.get("deps", {})
package_versions = {}
for name in ("PowerDynamics", "NetworkDynamics", "SciMLBase", "OrdinaryDiffEqRosenbrock"):
    entries = deps.get(name, [])
    if isinstance(entries, dict):
        entries = [entries]
    package_versions[name] = entries[0].get("version") if entries else None

status = run("git", "status", "--porcelain", "--untracked-files=normal").splitlines()
relative_here = HERE.relative_to(ROOT).as_posix()
preexisting = [line for line in status if relative_here not in line.replace("\\", "/")]
frozen = json.loads(
    (ROOT / "experiments/delay_dressed_replacement_frontier_20261002/FROZEN_DELAY_PATTERNS.json").read_text(encoding="utf-8")
)

result = {
    "experiment": "analytical delay-aware SG-to-GFL co-design mega experiment",
    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    "git_head": run("git", "rev-parse", "HEAD"),
    "git_branch": run("git", "branch", "--show-current"),
    "repository_dirty_before_new_outputs": bool(preexisting),
    "preexisting_status_entry_count_excluding_this_directory": len(preexisting),
    "new_artifact_directory": relative_here,
    "commit_or_push_performed": False,
    "julia_version": run("julia", "--startup-file=no", "-e", "print(VERSION)"),
    "python_version": platform.python_version(),
    "packages": package_versions,
    "source_sha256": {
        name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        for name in SOURCES
    },
    "frozen_contract": {
        "candidate_buses": list(range(30, 40)),
        "gain_bounds_multiplier_of_nominal": [0.25, 4.0],
        "nominal_Kp": 10 * 3.141592653589793,
        "nominal_Ki": (10 * 3.141592653589793) ** 2 / 4,
        "equilibrium_residual_inf_limit": 1e-8,
        "spectral_abscissa_limit_s_inv": -0.05,
        "frequency_limit_Hz": 0.5,
        "RoCoF_limit_Hz_s": 0.5,
        "voltage_range_pu": [0.9, 1.1],
        "min_SG_actuator_fraction_slack": 0.002,
        "measurement_window_s": 0.5,
        "event_time_s": 1.0,
        "simulation_end_s": 61.0,
        "event_solver": "Rodas5P",
        "event_abstol": 1e-9,
        "event_reltol": 1e-9,
        "event_saveat_s": 0.01,
        "hard_converter_current_and_dc_energy_limits_certified": False,
        "events": [
            {"id": "bus8_minus100", "bus": 8, "load_step_MW": -100.0},
            {"id": "bus16_plus100", "bus": 16, "load_step_MW": 100.0},
            {"id": "bus16_minus100", "bus": 16, "load_step_MW": -100.0},
            {"id": "bus29_plus100", "bus": 29, "load_step_MW": 100.0},
            {"id": "bus29_minus100", "bus": 29, "load_step_MW": -100.0},
        ],
        "uniform_delay_grid_ms": [0, 5, 10, 20, 30, 40, 50],
        "heterogeneous_frozen_patterns_path": "experiments/delay_dressed_replacement_frontier_20261002/FROZEN_DELAY_PATTERNS.json",
        "heterogeneous_frozen_multiset_sha256": frozen["multiset_sha256"],
        "heterogeneous_seed": frozen["heuristic_search"]["seed"],
        "experiment_random_seed": 20261002,
    },
    "model_variants": {
        "nonlinear_events": "PowerDynamics IEEE-39 PhysicalGFLDC, fixed interior architecture",
        "finite_spectrum": "zero-delay PowerDynamics network linearization, gauge mode removed",
        "DDE": "exact delayed PLL phase-error channel in DelayCharacteristic.jl; no Pade primary truth",
        "action_space": "fixed-coordinate augmented descriptor in PDExactDesignN.jl",
    },
}
(HERE / "MEGA_EXPERIMENT_MANIFEST.json").write_text(
    json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
)
print("wrote MEGA_EXPERIMENT_MANIFEST.json; source hashes", len(result["source_sha256"]))
