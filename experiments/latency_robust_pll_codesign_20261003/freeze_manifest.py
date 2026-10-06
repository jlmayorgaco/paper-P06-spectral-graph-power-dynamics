"""Freeze the L-stage contract and hash all read-only source dependencies."""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
import sys
import tomllib
from datetime import datetime, timezone
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
SOURCES = [
    "src/bnd_model_expN/PDExactDesignN.jl",
    "experiments/nonlinear_codesign_20261001/ReducedDAE.jl",
    "experiments/nonlinear_codesign_20261001/PDPhysicalReference.jl",
    "experiments/delay_dressed_replacement_frontier_20261002/DelayCharacteristic.jl",
    "experiments/analytical_delay_codesign_mega_20261002/m3_a_trace_integral.jl",
    "experiments/analytical_delay_codesign_mega_20261002/q56_validate_action_space.jl",
    "experiments/optimal_latency_margin_20261002/run_fast_root_locus.jl",
    "experiments/analytical_delay_codesign_mega_20261002/m1_validate_candidate.jl",
    "experiments/latency_robust_pll_codesign_20261003/run_event_local.jl",
    "experiments/latency_robust_pll_codesign_20261003/designs/Z_zero_delay_tuned.toml",
    "experiments/latency_robust_pll_codesign_20261003/designs/N_nominal.toml",
]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def output(command: list[str]) -> str:
    return subprocess.check_output(command, cwd=ROOT, text=True).strip()


def main() -> None:
    target = HERE / "EXPERIMENT_MANIFEST.json"
    if target.exists():
        raise SystemExit("Manifest already frozen; will not overwrite it")
    z = tomllib.loads((HERE / "designs/Z_zero_delay_tuned.toml").read_text())
    n = tomllib.loads((HERE / "designs/N_nominal.toml").read_text())
    assert z["rho"] == n["rho"] and len(z["rho"]) == 10
    kp0 = 10 * math.pi
    ki0 = (10 * math.pi) ** 2 / 4
    state = {
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_head": output(["git", "rev-parse", "HEAD"]),
        "git_dirty": bool(output(["git", "status", "--porcelain"])),
        "julia_version": output(["julia", "--version"]),
        "python_version": sys.version.split()[0],
        "model_sources_sha256": {name: sha(ROOT / name) for name in SOURCES},
        "physical_replacement_percent": z["GFL_percent"],
        "rho_vector": z["rho"],
        "gain_bounds": {
            "Kp_min": 0.25 * kp0,
            "Kp_max": 4 * kp0,
            "Ki_min": 0.25 * ki0,
            "Ki_max": 4 * ki0,
            "source": "m1_validate_candidate.jl frozen bounds",
        },
        "spectral_margin_alpha_max_s_inv": -0.05,
        "zero_delay_alpha_limit_s_inv": -0.05,
        "uniform_delay_definition": "same exogenous PLL error/measurement delay at all ten active GFL PLLs",
        "heterogeneous_delay_definition": "fixed exogenous vector on the same ten PLL error paths; not optimized",
        "numerical_precision": "Float64/ComplexF64; not interval certified",
        "frozen_events": [
            {"id": "bus8_minus100", "bus": 8, "delta_load_MW": -100},
            {"id": "bus16_plus100", "bus": 16, "delta_load_MW": 100},
            {"id": "bus16_minus100", "bus": 16, "delta_load_MW": -100},
            {"id": "bus29_plus100", "bus": 29, "delta_load_MW": 100},
            {"id": "bus29_minus100", "bus": 29, "delta_load_MW": -100},
        ],
        "event_limits": {
            "max_frequency_deviation_Hz": 0.5,
            "max_RoCoF_Hz_s": 0.5,
            "min_bus_voltage_pu": 0.9,
            "max_bus_voltage_pu": 1.1,
            "min_SG_actuator_fraction_slack": 0.002,
            "event_duration_s": 61.0,
            "sampling_dt_s": 0.01,
            "frequency_RoCoF_window_s": 0.5,
        },
        "optimizer_condition": "Only designs passing equilibrium, zero-delay physical spectrum and all five frozen nonlinear events are feasible; no global-optimum claim without proof.",
        "critical_delay_accuracy_target_ms": 0.05,
        "do_not_commit_or_push": True,
    }
    target.write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")
    print(f"frozen {target}; git={state['git_head']}; GFL={state['physical_replacement_percent']}%")


if __name__ == "__main__":
    main()
