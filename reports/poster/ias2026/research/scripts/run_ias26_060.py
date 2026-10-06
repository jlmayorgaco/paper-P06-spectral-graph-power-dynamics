#!/usr/bin/env python
"""Validate frozen IAS26-060 inputs, run its Julia smoke, then the full MC."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[5]
CAMPAIGN = ROOT / "reports/poster/ias2026/research/bnd_h4_mechanism"
CONFIG_REPO = CAMPAIGN / "configs/IAS26-060_MC_OPERATING_V1.json"
MANIFEST_REPO = CAMPAIGN / "inputs/IAS26-060_SCENARIOS_V1.csv"
SIDECAR_REPO = CAMPAIGN / "inputs/IAS26-060_PRE_RUN_SHA256.json"
CONFIG_050 = CAMPAIGN / "configs/IAS26-050_INTERVENTION_V1.json"
JULIA_RUNNER = ROOT / "reports/poster/ias2026/research/src/ias26_060_julia_runner.jl"
TX4_RUNNER = ROOT / "code/tx4/run_tx4_exact_p4_ieee39.jl"
TX4_DATA = ROOT / "code/tx4/frozen_ieee39_data.jl"
JULIA_PROJECT = ROOT / "research/ias2026_last_validation/env/julia/Project.toml"
JULIA_MANIFEST = ROOT / "research/ias2026_last_validation/env/julia/Manifest.toml"
EXACT_REPORT = ROOT / "reports/TX4_EXACT_P4_JULIA_FINAL_REPORT.md"
PD_ALT_REPORT = ROOT / "research/ias2026_last_validation/reports/P2_POWERDYNAMICS_V4_STATUS.md"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_frozen_inputs(run_root: Path) -> tuple[dict, list[dict[str, str]]]:
    config_path = run_root / "config/IAS26-060_MC_OPERATING_V1.json"
    manifest_path = run_root / "inputs/IAS26-060_SCENARIOS_V1.csv"
    sidecar_path = run_root / "inputs/IAS26-060_PRE_RUN_SHA256.json"
    sidecar = json.loads(sidecar_path.read_text(encoding="utf-8"))
    declared = sidecar["sha256"]
    rel_config = "reports/poster/ias2026/research/bnd_h4_mechanism/configs/IAS26-060_MC_OPERATING_V1.json"
    rel_manifest = "reports/poster/ias2026/research/bnd_h4_mechanism/inputs/IAS26-060_SCENARIOS_V1.csv"
    if sha256(config_path) != declared[rel_config]:
        raise RuntimeError("copied IAS26-060 config does not match frozen SHA256")
    if sha256(manifest_path) != declared[rel_manifest]:
        raise RuntimeError("copied scenario manifest does not match frozen SHA256")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    with manifest_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    expected = [f"IAS26-060-S{i:04d}" for i in range(1, 1001)]
    if [row["scenario_id"] for row in rows] != expected:
        raise RuntimeError("scenario manifest IDs or order changed")
    if config["scenario_generation"]["n_scenarios"] != 1000 or config["scenario_generation"]["seed"] != 20260926:
        raise RuntimeError("frozen N/seed mismatch")
    if any(row["power_flow_DAE_eigenvalue_outcomes"] != "NOT_RUN" for row in rows):
        raise RuntimeError("manifest unexpectedly contains outcomes")
    if any(row["outcome_dependent_filtering_applied"].lower() != "false" for row in rows):
        raise RuntimeError("manifest contains outcome-dependent filtering")
    if sha256(CONFIG_050) != config["source_provenance"]["sha256"][
        "reports/poster/ias2026/research/bnd_h4_mechanism/configs/IAS26-050_INTERVENTION_V1.json"
    ]:
        raise RuntimeError("IAS26-050 freeze hash mismatch")
    for rel, expected_hash in config["source_provenance"]["sha256"].items():
        path = ROOT / rel
        if not path.is_file() or sha256(path) != expected_hash:
            raise RuntimeError(f"frozen model/evidence hash mismatch: {rel}")
    return config, rows


def verify_backend() -> dict:
    exact_text = EXACT_REPORT.read_text(encoding="utf-8")
    pd_text = PD_ALT_REPORT.read_text(encoding="utf-8")
    if "exact 11-state GFL" not in exact_text or "86=4*11+6*7" not in exact_text:
        raise RuntimeError("the custom Julia exact-GFL11 regression report is not affirmative")
    if "FRESH_ALTERNATIVE_DYNAMIC_MODEL_V4_NEGATIVE_HOLDOUT" not in pd_text:
        raise RuntimeError("PowerDynamics alternative-model status could not be verified")
    return {
        "julia_backend": "custom exact TX4/P4/GFL11 Julia implementation",
        "powerdynamics_exact_model_available": False,
        "powerdynamics_reason": "Existing PowerDynamics IEEE-39 portfolio census is explicitly an alternative dynamic model; it does not reproduce frozen TX4/no-governor machine physics. The frozen GFL11 injector is exact at device level, but its available network census is not the canonical TX4 DAE.",
        "custom_julia_canonical_regression": "PASS per TX4_EXACT_P4_JULIA_FINAL_REPORT.md: 86 H4 states, 16/16 verdict agreement, alpha/frequency and H4 mode MAC pass.",
        "custom_julia_runner_sha256": sha256(TX4_RUNNER),
        "custom_julia_data_sha256": sha256(TX4_DATA),
        "campaign_runner_sha256": sha256(JULIA_RUNNER),
        "project_sha256": sha256(JULIA_PROJECT),
        "manifest_sha256": sha256(JULIA_MANIFEST),
        "exact_regression_report_sha256": sha256(EXACT_REPORT),
        "powerdynamics_alternative_report_sha256": sha256(PD_ALT_REPORT),
    }


def record_execution_overlay(run_root: Path, config: dict, backend: dict, attempt: str) -> dict:
    julia = shutil.which("julia")
    if not julia:
        raise RuntimeError("Julia executable not found")
    numerical = config["numerical_policy"]
    overlay = {
        "ticket": "IAS26-060B",
        "attempt": attempt,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "authorization": "The current user instruction explicitly requires Julia primary when an exact GFL11 runner exists; frozen scenario/config files are used unchanged.",
        **backend,
        "julia_executable": str(Path(julia).resolve()),
        "julia_project": str(JULIA_PROJECT.relative_to(ROOT)).replace("\\", "/"),
        "julia_version_required": "1.11.9 per existing exact-run record",
        "frozen_input_hashes": {
            "config": sha256(run_root / "config/IAS26-060_MC_OPERATING_V1.json"),
            "scenario_manifest": sha256(run_root / "inputs/IAS26-060_SCENARIOS_V1.csv"),
            "intervention_config": sha256(run_root / "config/IAS26-050_INTERVENTION_V1.json"),
        },
        "numerical_backend_note": {
            "frozen_config_backend": {
                "power_flow": numerical["power_flow"]["solver"],
                "eigenvalues": numerical["eigenvalues"]["solver"],
            },
            "executed_backend": "Julia custom TX4 analytic-equilibrium path; finite-difference DAE Jacobian and LinearAlgebra.eigen; frozen P4_G binding is made mutable only in-memory for the two preregistered g treatments.",
            "thresholds_preserved": {
                "power_flow_residual_acceptance": 1e-9,
                "dae_f": numerical["dae_equilibrium"]["max_abs_differential_residual"],
                "dae_g": numerical["dae_equilibrium"]["max_abs_algebraic_residual"],
                "eigenpair_relative_residual": numerical["eigenvalues"]["eigenpair_relative_residual_max"],
                "tau_dec": numerical["eigenvalues"]["tau_dec_s-1"],
            },
            "physics_changed": False,
            "parameters_or_scenarios_changed": False,
            "note": "This execution-backend choice is user-directed and is recorded separately; the frozen V1 config and manifest are not edited. Numerical tolerances and model equations are unchanged. Results are tied to the Julia source hashes above.",
        },
        "scenario_order": "frozen ascending manifest order for full MC; five technical smoke selectors in the order min L, max L, nearest L=1, S0001, S1000.",
        "no_outcome_dependent_filtering": True,
        "no_retries_or_replacements": True,
    }
    path = run_root / f"environment/EXECUTION_BACKEND_OVERLAY_ATTEMPT_{attempt}.json"
    if path.exists():
        old = json.loads(path.read_text(encoding="utf-8"))
        if old.get("frozen_input_hashes") != overlay["frozen_input_hashes"] or old.get("custom_julia_runner_sha256") != overlay["custom_julia_runner_sha256"] or old.get("campaign_runner_sha256") != overlay["campaign_runner_sha256"]:
            raise RuntimeError(f"existing overlay conflicts with the current immutable inputs/source: {path}")
        return old
    path.write_text(json.dumps(overlay, indent=2), encoding="utf-8")
    return overlay


def smoke_ids(rows: list[dict[str, str]]) -> list[str]:
    scales = [float(row["load_scale_L_s"]) for row in rows]
    indices = [min(range(len(scales)), key=scales.__getitem__),
               max(range(len(scales)), key=scales.__getitem__),
               min(range(len(scales)), key=lambda i: abs(scales[i] - 1.0)), 0, len(rows) - 1]
    ids = list(dict.fromkeys(rows[i]["scenario_id"] for i in indices))
    if len(ids) != 5:
        raise RuntimeError(f"smoke selector expected five distinct IDs, got {ids}")
    return ids


def run_julia(run_root: Path, phase: str, attempt: str, config: dict) -> None:
    julia = shutil.which("julia")
    assert julia is not None
    (run_root / "reports").mkdir(exist_ok=True)
    (run_root / "src").mkdir(exist_ok=True)
    np = config["numerical_policy"]
    physical = config["physical_feasibility"]
    command = [
        julia, f"--project={JULIA_PROJECT}", str(JULIA_RUNNER),
        "--campaign-root", str(run_root), "--run-root", str(run_root),
        "--manifest", str(run_root / "inputs/IAS26-060_SCENARIOS_V1.csv"),
        "--repo-root", str(ROOT), "--phase", phase, "--attempt", attempt,
        "--v-min", str(physical["bus_voltage_pu"]["minimum"]),
        "--v-max", str(physical["bus_voltage_pu"]["maximum"]),
        "--limit-tolerance", str(physical["limit_comparison_tolerance_pu"]),
        "--slack-p-min", str(config["slack_policy"]["pmin_pu"]),
        "--slack-p-max", str(config["slack_policy"]["pmax_pu"]),
        "--pf-residual-tolerance", "1e-9",
        "--pf-timeout", str(np["power_flow"]["timeout_seconds_per_case"]),
        "--dae-f-tolerance", str(np["dae_equilibrium"]["max_abs_differential_residual"]),
        "--dae-g-tolerance", str(np["dae_equilibrium"]["max_abs_algebraic_residual"]),
        "--eigenpair-tolerance", str(np["eigenvalues"]["eigenpair_relative_residual_max"]),
    ]
    logname = f"julia_{phase}_attempt_{attempt}.log"
    logfile = run_root / "logs" / logname
    if logfile.exists():
        raise RuntimeError(f"refusing to overwrite immutable execution log: {logfile}")
    print("COMMAND:", subprocess.list2cmdline(command), flush=True)
    with logfile.open("w", encoding="utf-8", newline="") as stream:
        process = subprocess.Popen(command, cwd=ROOT, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True, bufsize=1)
        assert process.stdout is not None
        for line in process.stdout:
            stream.write(line)
            stream.flush()
            print(line, end="", flush=True)
        code = process.wait()
    if code != 0:
        raise RuntimeError(f"Julia {phase} failed with exit code {code}; see {logfile}")


def validate_smoke(run_root: Path, ids: list[str], attempt: str) -> None:
    result_path = run_root / f"derived/MC_CASES_SMOKE_ATTEMPT_{attempt}.csv"
    with result_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 85:
        raise RuntimeError(f"smoke case table has {len(rows)} rows, expected 85")
    if set(row["scenario_id"] for row in rows) != set(ids):
        raise RuntimeError("smoke scenario IDs differ from the frozen selector")
    for scenario_id in ids:
        cases = [row for row in rows if row["scenario_id"] == scenario_id]
        if len(cases) != 17 or len({(row["portfolio_id"], row["g"]) for row in cases}) != 17:
            raise RuntimeError(f"smoke lattice incomplete for {scenario_id}")
        spec = run_root / "raw/julia" / f"smoke_attempt_{attempt}" / f"{scenario_id}.jld2"
        if not spec.is_file() or spec.stat().st_size < 128:
            raise RuntimeError(f"smoke spectrum serialization missing for {scenario_id}")
    completed_spectra = [row for row in rows if row["alpha_perp_complete_spectrum"].lower() == "true"]
    if not completed_spectra:
        raise RuntimeError("smoke did not complete any full transverse spectrum; do not start full MC")
    bad_dims = [row for row in completed_spectra if row["portfolio_id"] == "30+33+35+37" and int(row["state_count"]) != 86]
    if bad_dims:
        raise RuntimeError("smoke H4 state count deviates from exact GFL11 canonical 86-state gate")
    print(f"SMOKE_GATE_PASS scenarios={len(ids)} cases=85 spectra={len(completed_spectra)}", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--phase", choices=("smoke", "full", "all"), default="all")
    parser.add_argument("--attempt", default="01")
    args = parser.parse_args()
    run_root = args.run_root.resolve()
    if not run_root.is_relative_to(CAMPAIGN.resolve() / "results"):
        raise RuntimeError("RUN_ID must be located under the IAS26-060 results directory")
    config, rows = verify_frozen_inputs(run_root)
    backend = verify_backend()
    record_execution_overlay(run_root, config, backend, args.attempt)
    ids = smoke_ids(rows)
    if args.phase in ("smoke", "all"):
        run_julia(run_root, "smoke", args.attempt, config)
        validate_smoke(run_root, ids, args.attempt)
    if args.phase in ("full", "all"):
        if args.phase == "full":
            smoke_csvs = list((run_root / "derived").glob("MC_CASES_SMOKE_ATTEMPT_*.csv"))
            if not smoke_csvs:
                raise RuntimeError("full MC is disallowed before a verified smoke on this RUN_ID")
        run_julia(run_root, "full", args.attempt, config)
    print("IAS26_060_EXECUTION_PHASES_COMPLETE", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"IAS26_060_EXECUTION_BLOCKED: {exc}", file=sys.stderr, flush=True)
        raise SystemExit(2)
