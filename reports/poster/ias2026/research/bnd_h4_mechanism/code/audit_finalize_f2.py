"""Finalize IAS26-030 from frozen TX4 evidence; do not rerun the model.

This script reconciles the frozen alpha sweep, pre-normalized physical local
factors, and exact-boundary collective factor. It performs no equilibrium
solve, eigenvalue solve, parameter continuation, or root solve.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


BASELINE_RUN_ID = "20260925T204017_549c3c07_h4_crossmode_v3"
PORT_ORDER = (30, 33, 35, 37)
FROZEN_SHA256 = {
    "derived/TX4_CONTEXTUAL_RETURN_SWEEP.csv": "2dd4559acdd847349f68faa26818687eb3a6194f5a2387bab17205955525429b",
    "derived/TX4_PHYSICAL_LOCAL_FACTORS.csv": "4a7b1b64a4bfe39467856d36d0da798ae9106822cef69ab80073bb9e255f49d3",
    "derived/TX4_CONTEXTUAL_RETURN_AT_BOUNDARY.csv": "6bac61f319efb4216e46ec78156c1c4446c1e2cfe15c922cdb0b461ae6138e2e",
    "derived/TX4_CONTEXTUAL_RETURN_CORE_SUMMARY.json": "8308409f4edce2e657a6f9bd67517c3422acfce961578540896eb7bb3aa70b24",
    "claims/gates.json": "a06b65d2141c16a6f972ef2e7fd7ebebde9fd1a61bd8ae92db948e54294de1e4",
    "environment/runtime.json": "e99613e9b5848e21442f91fa33abc341e6bfc11c8db6c6c384507997c2f597b6",
}
FROZEN_CODE_SHA256 = {
    "reports/poster/ias2026/research/bnd_h4_mechanism/code/run_campaign.py": "8d2f54267bc5852ca21735002e34b9e6bf9864068be0255f431d3a00d9f55c09",
    "reports/poster/ias2026/research/bnd_h4_mechanism/code/render_f1_f2.py": "02a41dd9c3c4fc0c0ba9ba33b71639e93dbdfb98cae929c91e3bf2ca83ac6193",
    "reports/poster/ias2026/research/experiments/tx4_contextual_return.py": "77d4ed0f76f8c4a6ec369bc30b306361cff9b789c6f432e287d5f180d9210b74",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise RuntimeError(f"refusing to write an empty table: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, default=float) + "\n", encoding="utf-8")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def one(rows: list[dict[str, str]], key: str, value: str, source: str) -> dict[str, str]:
    matches = [row for row in rows if row[key] == value]
    require(len(matches) == 1, f"expected one {key}={value} in {source}, found {len(matches)}")
    return matches[0]


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("usage: audit_finalize_f2.py <repo-root> <new-run-id>")
    repo = Path(sys.argv[1]).resolve()
    run_id = sys.argv[2]
    require("/" not in run_id and "\\" not in run_id and run_id not in ("", ".", ".."), "RUN_ID must be a single path component")
    result_root = repo / "reports/poster/ias2026/research/bnd_h4_mechanism/results"
    run_root = result_root / run_id
    require(not run_root.exists(), f"immutable RUN_ID already exists: {run_root}")
    baseline = result_root / BASELINE_RUN_ID
    require(baseline.is_dir(), f"frozen baseline is missing: {baseline}")

    # Validate every frozen data input and adapter/source revision before any output is created.
    frozen_input_paths = {name: baseline / name for name in FROZEN_SHA256}
    input_hashes = {name: sha256(path) for name, path in frozen_input_paths.items()}
    require(input_hashes == FROZEN_SHA256, "one or more frozen IAS26-030 input hashes changed")
    code_hashes = {name: sha256(repo / name) for name in FROZEN_CODE_SHA256}
    require(code_hashes == FROZEN_CODE_SHA256, "one or more frozen F2/model source hashes changed")

    sweep = read_csv(frozen_input_paths["derived/TX4_CONTEXTUAL_RETURN_SWEEP.csv"])
    local = read_csv(frozen_input_paths["derived/TX4_PHYSICAL_LOCAL_FACTORS.csv"])
    boundary = read_csv(frozen_input_paths["derived/TX4_CONTEXTUAL_RETURN_AT_BOUNDARY.csv"])
    summary = json.loads(frozen_input_paths["derived/TX4_CONTEXTUAL_RETURN_CORE_SUMMARY.json"].read_text(encoding="utf-8"))
    baseline_gates = json.loads(frozen_input_paths["claims/gates.json"].read_text(encoding="utf-8"))
    runtime = json.loads(frozen_input_paths["environment/runtime.json"].read_text(encoding="utf-8"))

    g_star = float(summary["g_eigen_boundary"])
    sweep_by_g = {float(row["g"]): row for row in sweep}
    require(len(sweep_by_g) == 17, f"frozen alpha/collective sweep must have 17 unique g values, got {len(sweep_by_g)}")
    local_by_g: dict[float, dict[int, dict[str, str]]] = {}
    for row in local:
        g = float(row["g"])
        bus = int(row["device_bus"])
        require(bus in PORT_ORDER, f"unexpected physical-local port {bus}")
        require(bus not in local_by_g.setdefault(g, {}), f"duplicate local factor at g={g}, port={bus}")
        local_by_g[g][bus] = row
    local_sweep_g = set(local_by_g).intersection(sweep_by_g)
    require(local_sweep_g == set(sweep_by_g), "physical-local factors do not exactly cover all frozen sweep points")
    root_local_g = min(local_by_g, key=lambda g: abs(g - g_star))
    root_g_encoding_delta = abs(root_local_g - g_star)
    require(root_local_g not in sweep_by_g and root_g_encoding_delta <= math.ulp(g_star), "no unique physical-local row at the stored boundary to floating-point encoding precision")
    require(set(local_by_g) == set(sweep_by_g) | {root_local_g}, "physical-local factors contain unexpected g points")
    require(all(tuple(sorted(rows)) == PORT_ORDER for rows in local_by_g.values()), "local-factor port ordering/completeness changed")
    require(len(boundary) == len(PORT_ORDER), "boundary table must contain one row per H4 port")
    require(tuple(int(row["device_bus"]) for row in boundary) == PORT_ORDER, "boundary port ordering changed")
    require(all(float(row["g_root"]) == g_star for row in boundary), "boundary rows do not share the frozen g*")
    boundary_frequency = {float(row["frequency_hz"]) for row in boundary}
    boundary_collective = {float(row["collective_sigma_min"]) for row in boundary}
    require(len(boundary_frequency) == len(boundary_collective) == 1, "boundary rows disagree on common operating frequency or collective factor")
    boundary_hz = next(iter(boundary_frequency))
    boundary_collective_smin = next(iter(boundary_collective))

    aligned: list[dict] = []
    frequency_deltas: list[float] = []
    for g in sorted(sweep_by_g):
        op = sweep_by_g[g]
        frequency = float(op["frequency_hz"])
        for bus in PORT_ORDER:
            physical = local_by_g[g][bus]
            delta_f = abs(float(physical["frequency_hz"]) - frequency)
            frequency_deltas.append(delta_f)
            aligned.append(
                {
                    "point_kind": "frozen_sweep",
                    "g": g,
                    "physical_factor_g_recorded": g,
                    "g_abs_delta_local_vs_boundary": 0.0,
                    "device_bus": bus,
                    "alpha_perp_s-1": float(op["alpha"]),
                    "lambda_real_s-1": float(op["critical_real"]),
                    "lambda_imag_rad_s-1": float(op["critical_imag"]),
                    "frequency_hz": frequency,
                    "physical_local_sigma_min_I_plus_Mii": float(physical["physical_local_sigma_min"]),
                    "collective_sigma_min_I_plus_QH": float(op["collective_sigma_min"]),
                    "frequency_abs_delta_local_vs_sweep_hz": delta_f,
                    "equilibrium_residual_max_g": float(op["equilibrium_residual"]),
                    "port_order": "+".join(map(str, PORT_ORDER)),
                    "local_convention": "physical pre-normalized I+M_ii",
                    "collective_convention": "I+Q_H",
                    "alpha_provenance": "frozen TX4_CONTEXTUAL_RETURN_SWEEP.csv",
                }
            )

    # Add the exact stored boundary point; alpha=0 is the defining constraint,
    # not a value interpolated from the nearest sweep samples.
    for bus in PORT_ORDER:
        physical = local_by_g[root_local_g][bus]
        delta_f = abs(float(physical["frequency_hz"]) - boundary_hz)
        frequency_deltas.append(delta_f)
        aligned.append(
            {
                "point_kind": "stored_eigen_boundary",
                    "g": g_star,
                    "physical_factor_g_recorded": root_local_g,
                    "g_abs_delta_local_vs_boundary": root_g_encoding_delta,
                "device_bus": bus,
                "alpha_perp_s-1": 0.0,
                "lambda_real_s-1": 0.0,
                "lambda_imag_rad_s-1": float(boundary[0]["s_imag"]),
                "frequency_hz": boundary_hz,
                "physical_local_sigma_min_I_plus_Mii": float(physical["physical_local_sigma_min"]),
                "collective_sigma_min_I_plus_QH": boundary_collective_smin,
                "frequency_abs_delta_local_vs_sweep_hz": delta_f,
                "equilibrium_residual_max_g": "",
                "port_order": "+".join(map(str, PORT_ORDER)),
                "local_convention": "physical pre-normalized I+M_ii",
                "collective_convention": "I+Q_H",
                "alpha_provenance": "root constraint s_real=0; frozen bisection stop |alpha|<=1e-9",
            }
        )
    require(len(aligned) == 4 * (len(sweep) + 1), "aligned F2 table has an unexpected row count")

    # Check mode-branch continuity only from stored critical eigenvalues. Full
    # eigenvectors and competing-mode gaps were not retained by the baseline.
    ordered_sweep = [sweep_by_g[g] for g in sorted(sweep_by_g)]
    lambdas = [complex(float(row["critical_real"]), float(row["critical_imag"])) for row in ordered_sweep]
    alpha_values = [float(row["alpha"]) for row in ordered_sweep]
    frequencies = [float(row["frequency_hz"]) for row in ordered_sweep]
    adjacent_lambda_steps = [abs(right - left) for left, right in zip(lambdas, lambdas[1:])]
    mode_continuity = {
        "status": "SCALAR_CONTINUITY_SUPPORTED_NOT_EIGENVECTOR_CERTIFIED",
        "selection_rule": "largest real part among positive-frequency transverse modes in 0.3-1.5 Hz band (frozen tx4_contextual_return.py)",
        "sweep_points": len(ordered_sweep),
        "alpha_monotone_decreasing_on_frozen_grid": all(b < a for a, b in zip(alpha_values, alpha_values[1:])),
        "frequency_monotone_increasing_on_frozen_grid": all(b > a for a, b in zip(frequencies, frequencies[1:])),
        "maximum_adjacent_complex_lambda_step": max(adjacent_lambda_steps),
        "minimum_adjacent_complex_lambda_step": min(adjacent_lambda_steps),
        "adjacent_complex_lambda_steps": adjacent_lambda_steps,
        "eigenvectors_or_competing_mode_gaps_available": False,
        "limitation": "The frozen sweep stores only the selected eigenvalue at each g; mode-family/MAC certification cannot be inferred from these artifacts.",
    }

    gate_by_name = {item["name"]: item for item in baseline_gates}
    for name in ("G8 eigenvalue boundary", "G9 port/Schur identities", "G10 collective not local"):
        require(gate_by_name[name]["status"] == "PASS", f"baseline gate is not PASS: {name}")
    g8 = gate_by_name["G8 eigenvalue boundary"]
    g9 = gate_by_name["G9 port/Schur identities"]
    g10 = gate_by_name["G10 collective not local"]
    g8_observed = g8["observed"]
    g9_observed = g9["observed"]
    g10_observed = g10["observed"]
    boundary_difference = abs(float(summary["g_eigen_boundary"]) - float(summary["g_return_minimum"]))
    min_p4_local = min(float(local_by_g[0.03625][bus]["physical_local_sigma_min"]) for bus in PORT_ORDER)
    min_root_local = min(float(local_by_g[root_local_g][bus]["physical_local_sigma_min"]) for bus in PORT_ORDER)
    max_boundary_schur = max(float(row["schur_identity_residual"]) for row in boundary)
    min_collective_on_sweep = min(float(row["collective_sigma_min"]) for row in sweep)
    expected_g10 = g10["expected"]
    gate_results = {
        "boundary_alignment": {
            "status": "PASS" if boundary_difference <= float(g8["tolerance"]) else "BLOCKED",
            "g_eigen_boundary": g_star,
            "g_return_minimum": float(summary["g_return_minimum"]),
            "absolute_g_difference": boundary_difference,
            "baseline_tolerance": float(g8["tolerance"]),
            "nearest_sweep_alpha_abs": abs(float(g8_observed["nearest_grid_alpha"])),
            "baseline_gate": g8["status"],
        },
        "local_physical_regularity": {
            "status": "PASS" if min_p4_local >= float(expected_g10["p4_physical_local_sigma"][2:]) and min_root_local >= float(expected_g10["root_physical_local_sigma"][2:]) else "BLOCKED",
            "p4_min_physical_local_sigma": min_p4_local,
            "p4_frozen_minimum": float(expected_g10["p4_physical_local_sigma"][2:]),
            "boundary_min_physical_local_sigma": min_root_local,
            "boundary_frozen_minimum": float(expected_g10["root_physical_local_sigma"][2:]),
            "baseline_gate": g10["status"],
        },
        "collective_closure": {
            "status": "PASS" if boundary_collective_smin <= float(expected_g10["root_collective_sigma"][2:]) else "BLOCKED",
            "boundary_collective_sigma_min": boundary_collective_smin,
            "sweep_min_collective_sigma_min": min_collective_on_sweep,
            "frozen_boundary_maximum": float(expected_g10["root_collective_sigma"][2:]),
            "baseline_gate": g10["status"],
        },
        "port_schur_identity": {
            "status": "PASS" if max_boundary_schur <= float(g9["tolerance"]) and float(g9_observed["max_port"]) <= float(g9["tolerance"]) else "BLOCKED",
            "max_boundary_schur_residual": max_boundary_schur,
            "max_sweep_port_identity_residual": float(g9_observed["max_port"]),
            "frozen_maximum": float(g9["tolerance"]),
            "baseline_gate": g9["status"],
        },
        "frequency_alignment": {
            "status": "PASS_SAME_G_AND_PORT_ORDER; FLOAT_DIFFERENCES_RECORDED",
            "max_abs_local_vs_collective_frequency_difference_hz": max(frequency_deltas),
            "threshold": "none added; exact frequency differences are reported, not gated",
        },
        "mode_continuity": mode_continuity,
        "m1_reduced_root_comparison": "NOT_PERFORMED; M1 remains BLOCKED_M1_STRICT",
    }
    gate_pass = all(gate_results[name]["status"] == "PASS" for name in (
        "boundary_alignment", "local_physical_regularity", "collective_closure", "port_schur_identity"
    ))
    require(gate_pass, "an IAS26-030 frozen gate failed; see diagnostics before publishing F2")

    # Output creation starts only after all frozen inputs and gates validate.
    for folder in ("config", "raw/inputs", "tables", "figures", "claims", "report", "environment"):
        (run_root / folder).mkdir(parents=True, exist_ok=True)
    for name, path in frozen_input_paths.items():
        shutil.copy2(path, run_root / "raw" / "inputs" / Path(name).name)

    write_csv(run_root / "tables" / "F2_ALIGNED_LOCAL_COLLECTIVE_ALPHA.csv", aligned)
    write_csv(
        run_root / "tables" / "F2_BOUNDARY_PORT_AUDIT.csv",
        [
            {
                **row,
                "physical_factor_g_recorded": root_local_g,
                "g_abs_delta_local_vs_boundary": root_g_encoding_delta,
                "physical_local_sigma_min_I_plus_Mii": float(local_by_g[root_local_g][int(row["device_bus"])]["physical_local_sigma_min"]),
                "local_frequency_hz": float(local_by_g[root_local_g][int(row["device_bus"])]["frequency_hz"]),
                "frequency_abs_delta_local_vs_boundary_hz": abs(float(local_by_g[root_local_g][int(row["device_bus"])]["frequency_hz"]) - float(row["frequency_hz"])),
                "port_order": "+".join(map(str, PORT_ORDER)),
                "local_operator_convention": "physical pre-normalized I+M_ii",
                "collective_operator_convention": "I+Q_H",
                "m1_reduced_comparison": "BLOCKED_NOT_PERFORMED",
            }
            for row in boundary
        ],
    )

    # Two-panel F2: frozen alpha_perp(g), then physical local versus collective
    # closure. The exact root measurements are marked separately from the grid.
    colors = {30: "#386cb0", 33: "#f0027f", 35: "#1b9e77", 37: "#e6ab02"}
    fig, (ax_alpha, ax_sigma) = plt.subplots(2, 1, figsize=(8.2, 6.6), sharex=True, gridspec_kw={"height_ratios": [0.85, 1.25]})
    grid_g = sorted(sweep_by_g)
    ax_alpha.axhline(0.0, color="#20252b", linewidth=0.8)
    ax_alpha.plot(grid_g, [float(sweep_by_g[g]["alpha"]) for g in grid_g], "o-", color="#386cb0", linewidth=1.4, markersize=4, label=r"target-band $\alpha_\perp(g)$")
    ax_alpha.scatter([g_star], [0.0], marker="*", s=105, facecolor="#8b1e3f", edgecolor="#20252b", linewidth=0.6, zorder=4, label=r"stored boundary $\alpha_\perp=0$")
    ax_alpha.axvline(g_star, color="#8b1e3f", linewidth=0.9, linestyle=":")
    ax_alpha.set_ylabel(r"$\alpha_\perp$ [s$^{-1}$]")
    ax_alpha.set_title(r"F2 — Stability boundary and local/collective closure")
    ax_alpha.grid(axis="y", color="#d8dde3", linewidth=0.5, alpha=0.8)
    ax_alpha.legend(frameon=False, fontsize=8, loc="best")

    for bus in PORT_ORDER:
        gvals = sorted(sweep_by_g)
        localvals = [float(local_by_g[g][bus]["physical_local_sigma_min"]) for g in gvals]
        ax_sigma.semilogy(gvals, localvals, "o-", color=colors[bus], linewidth=1.1, markersize=3.5, label=fr"physical local $I+M_{{{bus},{bus}}}$")
        ax_sigma.scatter([g_star], [float(local_by_g[root_local_g][bus]["physical_local_sigma_min"])], marker="o", s=35, facecolor="white", edgecolor=colors[bus], linewidth=1.1, zorder=4)
    collective_grid = [float(sweep_by_g[g]["collective_sigma_min"]) for g in grid_g]
    ax_sigma.semilogy(grid_g, collective_grid, "s--", color="#20252b", linewidth=1.6, markersize=3.6, label=r"collective $I+Q_H$ (sweep)")
    ax_sigma.scatter([g_star], [boundary_collective_smin], marker="D", s=43, facecolor="#8b1e3f", edgecolor="#20252b", linewidth=0.6, zorder=5, label=r"collective at stored $g^*$")
    ax_sigma.axvline(g_star, color="#8b1e3f", linewidth=0.9, linestyle=":")
    ax_sigma.set_xlabel(r"feedback parameter $g$")
    ax_sigma.set_ylabel(r"minimum singular value $\sigma_{\min}$")
    ax_sigma.grid(which="both", color="#d8dde3", linewidth=0.5, alpha=0.8)
    ax_sigma.legend(frameon=False, fontsize=7.5, ncol=2, loc="best")
    fig.text(0.99, 0.01, "local: physical pre-normalized I+M_ii · collective: I+Q_H · fixed port order 30, 33, 35, 37", ha="right", va="bottom", fontsize=7, color="#4b5563")
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    for ext in ("png", "pdf", "svg"):
        fig.savefig(run_root / "figures" / f"F2_local_collective_alpha.{ext}", dpi=320 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)

    commit = __import__("subprocess").check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip()
    dirty = bool(__import__("subprocess").check_output(["git", "status", "--porcelain"], cwd=repo, text=True).strip())
    model_hash = hashlib.sha256(json.dumps(runtime["source_hashes"], sort_keys=True).encode("utf-8")).hexdigest()
    config = {
        "ticket": "IAS26-030",
        "run_id": run_id,
        "figure": "F2 only",
        "scientific_scope": "reconciliation and rendering from the frozen 17-point TX4 sweep plus exact stored boundary; no model/eigenvalue/root recomputation",
        "baseline_run_id": BASELINE_RUN_ID,
        "frozen_inputs_sha256": input_hashes,
        "frozen_adapter_and_model_code_sha256": code_hashes,
        "g_grid": grid_g,
        "g_star": g_star,
        "physical_local_g_recorded": root_local_g,
        "physical_local_g_float_encoding_delta": root_g_encoding_delta,
        "port_order": list(PORT_ORDER),
        "local_operator": "physical pre-normalized I+M_ii",
        "collective_operator": "I+Q_H",
        "frozen_gate_source": "baseline claims/gates.json; no thresholds changed or added",
        "root_solver": "not rerun; stored 26-iteration bisection result reused",
        "mode_identity_limit": mode_continuity["limitation"],
        "equilibrium_fingerprint": None,
        "equilibrium_provenance": "same frozen tx4_case construction, H4=(30,33,35,37), k=1.425, t=1.5, h=1.0, per-g voltage_gain=g; equilibrium state vectors/hashes were not retained, so bitwise state equality is not claimed",
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    write_json(run_root / "config" / "FROZEN_CONFIG.json", config)
    status = {
        "ticket": "IAS26-030",
        "status": "PASS" if gate_pass else "BLOCKED",
        "figure": "F2",
        "gates": gate_results,
        "claim_safe": ["Physical local factors remain regular while collective closure approaches singularity."],
        "claims_unsafe": [
            "Do not claim positive local damping.",
            "Do not compare reduced feedback poles with full DAE poles; IAS26-010 remains BLOCKED_M1_STRICT.",
            "Do not claim eigenvector-certified mode-family identity across g; the frozen sweep lacks eigenvectors and competing-mode gaps.",
        ],
        "boundary_metrics": {
            "g_star": g_star,
            "g_boundary_difference": boundary_difference,
            "root_alpha_residual_bound_from_frozen_solver": 1e-9,
            "physical_local_sigma_min": min_root_local,
            "collective_sigma_min": boundary_collective_smin,
            "max_boundary_schur_residual": max_boundary_schur,
            "max_frequency_alignment_delta_hz": max(frequency_deltas),
        },
    }
    write_json(run_root / "claims" / "IAS26-030_STATUS.json", status)
    (run_root / "report" / "IAS26-030_SUMMARY.md").write_text(
        "# IAS26-030 — F2 closure and boundary audit\n\n"
        f"Status: **{status['status']}**\n\n"
        f"Run ID: `{run_id}`\n\n"
        "F2 was assembled from the frozen baseline only. No equilibrium, eigenvalue, continuation, or root solve was run. The 17 alpha/collective sweep points were exact-joined to their physical pre-normalized local factors in fixed port order `(30, 33, 35, 37)`. The exact stored boundary row was added from the boundary artifact, not interpolated.\n\n"
        f"- `g*`: `{g_star:.17g}`; eigen/return boundary difference: `{boundary_difference:.6g}` (baseline tolerance `{float(g8['tolerance']):g}`).\n"
        f"- At the stored boundary, minimum physical local `sigma_min(I+M_ii)`: `{min_root_local:.9g}`; collective `sigma_min(I+Q_H)`: `{boundary_collective_smin:.9g}`.\n"
        f"- Maximum local-versus-collective frequency difference in matched records: `{max(frequency_deltas):.6g} Hz` (reported, no new threshold).\n"
        f"- Maximum boundary Schur identity residual: `{max_boundary_schur:.6g}`; frozen port identity maximum: `{float(g9_observed['max_port']):.6g}`.\n"
        f"- Mode audit: alpha monotone={mode_continuity['alpha_monotone_decreasing_on_frozen_grid']}; frequency monotone={mode_continuity['frequency_monotone_increasing_on_frozen_grid']}; no eigenvectors/gaps were retained, so mode-family identity is not certified.\n\n"
        "The figure supports only the narrow claim that physical local factors remain regular while collective closure approaches singularity. It does not establish positive local damping or validate reduced poles against the DAE.\n",
        encoding="utf-8",
    )
    write_json(run_root / "environment" / "provenance.json", {
        "run_id": run_id,
        "baseline_run_id": BASELINE_RUN_ID,
        "baseline_git_commit": runtime["git_commit"],
        "current_git_commit": commit,
        "current_worktree_dirty": dirty,
        "baseline_source_hashes": runtime["source_hashes"],
        "frozen_data_input_sha256": input_hashes,
        "frozen_analysis_code_sha256": {**code_hashes, str(Path(__file__).resolve().relative_to(repo)): sha256(Path(__file__).resolve())},
        "python": sys.version,
        "numpy": np.__version__,
        "matplotlib": matplotlib.__version__,
        "equilibrium_state_hash": None,
        "equilibrium_hash_limitation": "not present in the frozen source artifacts; no new model solve performed",
        "execution_policy": "post-processing only; no temporary outputs outside run root; no archives",
    })

    output_hashes = {}
    for path in sorted(p for p in run_root.rglob("*") if p.is_file()):
        output_hashes[str(path.relative_to(run_root)).replace("\\", "/")] = sha256(path)
    write_json(run_root / "environment" / "OUTPUT_SHA256.json", output_hashes)
    print(json.dumps({"run_id": run_id, "status": status["status"], "g_star": g_star, "aligned_rows": len(aligned), "max_frequency_alignment_delta_hz": max(frequency_deltas), "output": str(run_root)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
