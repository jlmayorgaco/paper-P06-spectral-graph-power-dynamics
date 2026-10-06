"""Freeze IAS26-060 operating-point scenarios only; never solves the model.

This input generator reads the frozen IEEE-39 JSON data and the preregistered
IAS26-060 config. It performs no power flow, DAE construction, Jacobian, or
eigenvalue calculation. All 1,000 scenario IDs are retained in the manifest.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.special import ndtr, ndtri


CONFIG_REL = Path("reports/poster/ias2026/research/bnd_h4_mechanism/configs/IAS26-060_MC_OPERATING_V1.json")
NETWORK_REL = Path("reports/poster/ias2026/research/configs/ias2026/ieee39_network.json")
MANIFEST_REL = Path("reports/poster/ias2026/research/bnd_h4_mechanism/inputs/IAS26-060_SCENARIOS_V1.csv")
SUMMARY_REL = Path("reports/poster/ias2026/research/bnd_h4_mechanism/inputs/IAS26-060_SCENARIOS_V1_SUMMARY.json")
HASHES_REL = Path("reports/poster/ias2026/research/bnd_h4_mechanism/inputs/IAS26-060_PRE_RUN_SHA256.json")
GENERATOR_REL = Path("reports/poster/ias2026/research/bnd_h4_mechanism/code/freeze_ias26_060_prerun.py")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_cell(value: object) -> str:
    return json.dumps(value, separators=(",", ":"), ensure_ascii=True, allow_nan=False)


def dispatch_with_base_participation(
    delta_p: float, generators: list[dict], tolerance: float
) -> tuple[dict[int, float], list[int], float]:
    """Apply base-P participation, saturating bounds and redistributing remainder."""
    scheduled = {int(row["bus"]): float(row["p0_pu"]) for row in generators}
    lower = {int(row["bus"]): float(row["pmin_pu"]) for row in generators}
    upper = {int(row["bus"]): float(row["pmax_pu"]) for row in generators}
    base_p = {int(row["bus"]): float(row["p0_pu"]) for row in generators}
    remaining = float(delta_p)
    active = set(scheduled)
    saturated: list[int] = []

    for _ in range(len(scheduled) + 1):
        if abs(remaining) <= tolerance:
            remaining = 0.0
            break
        if not active:
            break
        denom = sum(base_p[bus] for bus in active)
        if denom <= 0.0:
            break
        proposed = {bus: remaining * base_p[bus] / denom for bus in active}
        if remaining > 0.0:
            hit = sorted(bus for bus in active if scheduled[bus] + proposed[bus] > upper[bus])
            bound = upper
        else:
            hit = sorted(bus for bus in active if scheduled[bus] + proposed[bus] < lower[bus])
            bound = lower
        if not hit:
            for bus in sorted(active):
                scheduled[bus] += proposed[bus]
            remaining = 0.0
            break
        for bus in hit:
            delivered = bound[bus] - scheduled[bus]
            scheduled[bus] = bound[bus]
            remaining -= delivered
            active.remove(bus)
            saturated.append(bus)

    residual = float(remaining)
    return scheduled, saturated, residual


def summarize(values: np.ndarray) -> dict[str, float]:
    return {
        "count": int(values.size),
        "mean": float(np.mean(values)),
        "std_population": float(np.std(values, ddof=0)),
        "min": float(np.min(values)),
        "q05": float(np.quantile(values, 0.05)),
        "median": float(np.median(values)),
        "q95": float(np.quantile(values, 0.95)),
        "max": float(np.max(values)),
    }


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: freeze_ias26_060_prerun.py <repo-root>")
    repo = Path(sys.argv[1]).resolve()
    config_path = repo / CONFIG_REL
    network_path = repo / NETWORK_REL
    manifest_path = repo / MANIFEST_REL
    summary_path = repo / SUMMARY_REL
    hashes_path = repo / HASHES_REL

    config = json.loads(config_path.read_text(encoding="utf-8"))
    if config.get("status") != "FROZEN_PRE_RUN_AWAITING_USER_REVIEW":
        raise RuntimeError("config is not frozen and awaiting review")
    if config.get("execution_authorization") != "PRE_RUN_INPUT_GENERATION_ONLY":
        raise RuntimeError("refusing to run outside the pre-run phase")
    if any(path.exists() for path in (manifest_path, summary_path, hashes_path)):
        raise FileExistsError("refusing to overwrite a frozen IAS26-060 pre-run artifact")

    source_hashes = config["source_provenance"]["sha256"]
    for relative, expected in source_hashes.items():
        actual = sha256_file(repo / relative)
        if actual != expected:
            raise RuntimeError(f"frozen source hash changed: {relative}")
    if sha256_file(repo / GENERATOR_REL) != config["source_provenance"]["prerun_generator_sha256"]:
        raise RuntimeError("pre-run generator hash differs from the frozen config")

    network = json.loads(network_path.read_text(encoding="utf-8"))
    if network["source"]["system_base_mva"] != config["system"]["base_mva"]:
        raise RuntimeError("system base differs from the frozen config")
    base_loads = sorted(network["loads"], key=lambda row: int(row["bus"]))
    base_pv = sorted(network["pv"], key=lambda row: int(row["bus"]))
    config_loads = config["base_case_snapshot"]["loads_pu"]
    config_generators = config["base_case_snapshot"]["pv_generators_pu"]
    if [(int(row["bus"]), float(row["p0"]), float(row["q0"])) for row in base_loads] != [
        (int(row["bus"]), float(row["p0_pu"]), float(row["q0_pu"])) for row in config_loads
    ]:
        raise RuntimeError("base load snapshot differs from the frozen IEEE-39 JSON")
    if [(int(row["bus"]), float(row["p0"]), float(row["q0"]), float(row["pmin"]), float(row["pmax"]), float(row["qmin"]), float(row["qmax"]), float(row["Sn"]), float(row["v0"])) for row in base_pv] != [
        (int(row["bus"]), float(row["p0_pu"]), float(row["q0_pu"]), float(row["pmin_pu"]), float(row["pmax_pu"]), float(row["qmin_pu"]), float(row["qmax_pu"]), float(row["sn_mva"]), float(row["v0_pu"]))
        for row in config_generators
    ]:
        raise RuntimeError("base generation snapshot differs from the frozen IEEE-39 JSON")
    slack_source = network["slack"][0]
    slack_frozen = config["base_case_snapshot"]["slack_pu"]
    if (int(slack_source["bus"]), float(slack_source["p0"]), float(slack_source["q0"]), float(slack_source["pmin"]), float(slack_source["pmax"]), float(slack_source["Sn"]), float(slack_source["v0"]), float(slack_source["a0"])) != (
        int(slack_frozen["bus"]), float(slack_frozen["p0_pu"]), float(slack_frozen["q0_pu"]), float(slack_frozen["pmin_pu"]), float(slack_frozen["pmax_pu"]), float(slack_frozen["sn_mva"]), float(slack_frozen["v0_pu"]), float(slack_frozen["a0_rad"])
    ):
        raise RuntimeError("base slack snapshot differs from the frozen IEEE-39 JSON")

    gen_cfg = config["scenario_generation"]
    n = int(gen_cfg["n_scenarios"])
    seed = int(gen_cfg["seed"])
    if n != 1000 or seed != 20260926:
        raise RuntimeError("frozen scenario count or seed changed")
    loads = config_loads
    p0_load = np.asarray([float(row["p0_pu"]) for row in loads], dtype=np.float64)
    q0_load = np.asarray([float(row["q0_pu"]) for row in loads], dtype=np.float64)
    load_buses = [int(row["bus"]) for row in loads]
    load_weights = p0_load / p0_load.sum()
    generators = [
        row for row in config_generators
        if int(row["bus"]) in set(config["redispatch"]["participation_buses"])
    ]
    expected_participants = [int(row["bus"]) for row in base_pv if int(row["bus"]) != int(network["slack"][0]["bus"])]
    if config["redispatch"]["participation_buses"] != expected_participants:
        raise RuntimeError("redispatch participant set is not all scheduled non-slack PV generators")
    gen_buses = [int(row["bus"]) for row in generators]
    if gen_buses != config["redispatch"]["participation_buses"]:
        raise RuntimeError("redispatch participant order differs from the frozen config")
    participation = {int(row["bus"]): float(row["p0_pu"]) / sum(float(x["p0_pu"]) for x in generators) for row in generators}

    epsilon = gen_cfg["spatial_epsilon"]
    sigma = float(epsilon["truncated_normal_scale_pu"])
    lower = float(epsilon["support_pu"][0])
    upper = float(epsilon["support_pu"][1])
    trunc_lo = float(ndtr(lower / sigma))
    trunc_hi = float(ndtr(upper / sigma))
    a = upper / sigma
    phi_a = math.exp(-0.5 * a * a) / math.sqrt(2.0 * math.pi)
    actual_truncated_sd = sigma * math.sqrt(1.0 - 2.0 * a * phi_a / (2.0 * ndtr(a) - 1.0))
    if abs(actual_truncated_sd - float(epsilon["truncated_population_sd_pu"])) > 1e-14:
        raise RuntimeError("truncated Gaussian scale does not produce the frozen population SD")
    rng = np.random.Generator(np.random.PCG64(seed))
    l_min, l_max = (float(v) for v in gen_cfg["load_scale_support"])
    max_attempts = int(epsilon["max_proposals_per_scenario"])
    support_roundoff = float(epsilon["support_comparison_roundoff_pu"])
    dispatch_tol = float(config["redispatch"]["allocation_residual_tolerance_pu"])

    rows: list[dict[str, object]] = []
    load_p_totals, load_q_totals, load_scales, raw_eps_all, final_eps_all = [], [], [], [], []
    gen_delta_by_bus = {bus: [] for bus in gen_buses}
    gen_setpoint_by_bus = {bus: [] for bus in gen_buses}
    saturation_counts = {bus: 0 for bus in gen_buses}
    rejected_proposals = 0
    max_abs_weighted_mean = 0.0
    redispatch_shortfalls = 0

    for scenario_index in range(n):
        scenario_id = f"IAS26-060-S{scenario_index + 1:04d}"
        load_u = float(rng.random())
        load_scale = l_min + (l_max - l_min) * load_u
        proposals: list[dict[str, object]] = []
        accepted = None
        for attempt_index in range(max_attempts):
            uniforms = rng.random(len(loads), dtype=np.float64)
            raw_epsilon = sigma * ndtri(trunc_lo + uniforms * (trunc_hi - trunc_lo))
            weighted_mean = float(np.dot(load_weights, raw_epsilon))
            centered = raw_epsilon - weighted_mean
            centered_max = float(np.max(np.abs(centered)))
            proposal_ok = centered_max <= upper + support_roundoff and float(np.min(centered)) >= lower - support_roundoff
            proposals.append({
                "attempt_index_zero_based": attempt_index,
                "raw_uniform_by_load_bus": {str(bus): float(u) for bus, u in zip(load_buses, uniforms, strict=True)},
                "raw_truncated_epsilon_pu_by_load_bus": {str(bus): float(v) for bus, v in zip(load_buses, raw_epsilon, strict=True)},
                "removed_base_load_weighted_mean_pu": weighted_mean,
                "projected_epsilon_min_pu": float(np.min(centered)),
                "projected_epsilon_max_pu": float(np.max(centered)),
                "accepted_for_scenario": bool(proposal_ok),
            })
            if proposal_ok:
                accepted = (uniforms, raw_epsilon, centered, weighted_mean, centered_max, attempt_index)
                break
        if accepted is None:
            raise RuntimeError(f"bounded epsilon sampler exhausted for {scenario_id}; no outputs written")

        uniforms, raw_epsilon, eps_final, removed_mean, centered_max, accepted_attempt = accepted
        if np.min(eps_final) < lower - support_roundoff or np.max(eps_final) > upper + support_roundoff:
            raise RuntimeError(f"final epsilon support violation in {scenario_id}")
        weighted_residual = float(np.dot(load_weights, eps_final))
        if abs(weighted_residual) > 1e-12:
            raise RuntimeError(f"base-load-weighted epsilon mean is not zero in {scenario_id}")
        max_abs_weighted_mean = max(max_abs_weighted_mean, abs(weighted_residual))
        rejected_proposals += accepted_attempt

        load_multiplier = load_scale * (1.0 + eps_final)
        p_load_pu = p0_load * load_multiplier
        q_load_pu = q0_load * load_multiplier  # exact preservation of each signed base Q/P ratio
        p_load_total = float(np.sum(p_load_pu))
        q_load_total = float(np.sum(q_load_pu))
        requested_delta_p = p_load_total - float(np.sum(p0_load))
        pg_scheduled, saturated, residual = dispatch_with_base_participation(
            requested_delta_p, generators, dispatch_tol
        )
        pg_delta = {bus: pg_scheduled[bus] - float(next(x["p0_pu"] for x in generators if int(x["bus"]) == bus)) for bus in gen_buses}
        capacity_feasible = abs(residual) <= dispatch_tol
        redispatch_shortfalls += int(not capacity_feasible)
        for bus in saturated:
            saturation_counts[bus] += 1
        for bus in gen_buses:
            gen_delta_by_bus[bus].append(pg_delta[bus])
            gen_setpoint_by_bus[bus].append(pg_scheduled[bus])
        load_scales.append(load_scale)
        load_p_totals.append(p_load_total)
        load_q_totals.append(q_load_total)
        raw_eps_all.extend(raw_epsilon.tolist())
        final_eps_all.extend(eps_final.tolist())

        rows.append({
            "scenario_id": scenario_id,
            "scenario_index_zero_based": scenario_index,
            "seed": seed,
            "bit_generator": "NumPy.PCG64",
            "rng_stream_policy": "per-ID sequential: one load-scale uniform, then 19 epsilon uniforms per proposal until bounded weighted-centered vector accepted",
            "load_scale_uniform_u": load_u,
            "load_scale_L_s": load_scale,
            "epsilon_sampler_status": "ACCEPTED_INPUT_DRAW",
            "epsilon_proposal_count": len(proposals),
            "accepted_proposal_index_zero_based": accepted_attempt,
            "epsilon_proposal_history_json": json_cell(proposals),
            "accepted_raw_epsilon_uniforms_by_bus_json": json_cell({str(bus): float(u) for bus, u in zip(load_buses, uniforms, strict=True)}),
            "accepted_raw_truncated_epsilon_pu_by_bus_json": json_cell({str(bus): float(v) for bus, v in zip(load_buses, raw_epsilon, strict=True)}),
            "epsilon_removed_base_load_weighted_mean_pu": removed_mean,
            "epsilon_final_by_bus_pu_json": json_cell({str(bus): float(v) for bus, v in zip(load_buses, eps_final, strict=True)}),
            "epsilon_final_weighted_mean_pu": weighted_residual,
            "epsilon_final_min_pu": float(np.min(eps_final)),
            "epsilon_final_max_pu": float(np.max(eps_final)),
            "epsilon_projected_max_abs_before_acceptance_pu": centered_max,
            "P_L_base_by_bus_pu_json": json_cell({str(bus): float(v) for bus, v in zip(load_buses, p0_load, strict=True)}),
            "Q_L_base_by_bus_pu_json": json_cell({str(bus): float(v) for bus, v in zip(load_buses, q0_load, strict=True)}),
            "P_L_s_by_bus_pu_json": json_cell({str(bus): float(v) for bus, v in zip(load_buses, p_load_pu, strict=True)}),
            "Q_L_s_by_bus_pu_json": json_cell({str(bus): float(v) for bus, v in zip(load_buses, q_load_pu, strict=True)}),
            "P_L_s_by_bus_MW_json": json_cell({str(bus): float(v * config["system"]["base_mva"]) for bus, v in zip(load_buses, p_load_pu, strict=True)}),
            "Q_L_s_by_bus_Mvar_json": json_cell({str(bus): float(v * config["system"]["base_mva"]) for bus, v in zip(load_buses, q_load_pu, strict=True)}),
            "P_L_total_pu": p_load_total,
            "Q_L_total_pu": q_load_total,
            "P_L_total_MW": p_load_total * float(config["system"]["base_mva"]),
            "Q_L_total_Mvar": q_load_total * float(config["system"]["base_mva"]),
            "P_G_base_by_bus_pu_json": json_cell({str(int(x["bus"])): float(x["p0_pu"]) for x in generators}),
            "P_G_participation_factor_by_bus_json": json_cell({str(bus): participation[bus] for bus in gen_buses}),
            "delta_P_G_requested_total_pu": requested_delta_p,
            "delta_P_G_requested_total_MW": requested_delta_p * float(config["system"]["base_mva"]),
            "delta_P_G_allocated_by_bus_pu_json": json_cell({str(bus): pg_delta[bus] for bus in gen_buses}),
            "P_G_scheduled_by_bus_pu_json": json_cell({str(bus): pg_scheduled[bus] for bus in gen_buses}),
            "P_G_scheduled_by_bus_MW_json": json_cell({str(bus): pg_scheduled[bus] * float(config["system"]["base_mva"]) for bus in gen_buses}),
            "Q_G_schedule_status_by_bus_json": json_cell({str(bus): "NOT_PRESET_PV_VOLTAGE_CONTROL_SOLVES_Q" for bus in gen_buses}),
            "delta_P_G_allocated_total_pu": float(sum(pg_delta.values())),
            "delta_P_G_allocated_total_MW": float(sum(pg_delta.values())) * float(config["system"]["base_mva"]),
            "delta_P_G_allocation_residual_pu": residual,
            "redispatch_saturated_buses_json": json_cell(saturated),
            "pre_run_redispatch_capacity_status": "WITHIN_P_LIMITS" if capacity_feasible else "CAPACITY_SHORTFALL_RETAIN_SCENARIO",
            "slack_bus": int(config["slack_policy"]["bus"]),
            "slack_PQ_setpoint": "NOT_PRESET; swing/slack absorbs only solved PF residual",
            "slack_actual_PQ": "NOT_EVALUATED_IN_PRE_RUN",
            "outcome_dependent_filtering_applied": False,
            "power_flow_DAE_eigenvalue_outcomes": "NOT_RUN",
        })

    if len(rows) != n or len({row["scenario_id"] for row in rows}) != n:
        raise RuntimeError("scenario manifest did not produce exactly 1,000 unique IDs")
    if [row["scenario_id"] for row in rows] != [f"IAS26-060-S{i + 1:04d}" for i in range(n)]:
        raise RuntimeError("scenario ID sequence is incomplete or out of order")

    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    eps_by_bus = {
        str(bus): {
            "mean": float(np.mean([json.loads(str(row["epsilon_final_by_bus_pu_json"]))[str(bus)] for row in rows])),
            "std_population": float(np.std([json.loads(str(row["epsilon_final_by_bus_pu_json"]))[str(bus)] for row in rows], ddof=0)),
        }
        for bus in load_buses
    }
    pg_bus_summary = {
        str(bus): {
            "base_pu": float(next(x["p0_pu"] for x in generators if int(x["bus"]) == bus)),
            "scheduled_min_pu": float(np.min(gen_setpoint_by_bus[bus])),
            "scheduled_max_pu": float(np.max(gen_setpoint_by_bus[bus])),
            "mean_redispatch_pu": float(np.mean(gen_delta_by_bus[bus])),
            "at_limit_count": int(saturation_counts[bus]),
        }
        for bus in gen_buses
    }
    summary = {
        "ticket": "IAS26-060",
        "phase": "PRE_RUN_INPUT_GENERATION_ONLY",
        "status": "FROZEN_PRE_RUN_AWAITING_USER_REVIEW",
        "scenario_count": len(rows),
        "first_scenario_id": rows[0]["scenario_id"],
        "last_scenario_id": rows[-1]["scenario_id"],
        "seed": seed,
        "bit_generator": "NumPy.PCG64",
        "no_outcome_dependent_filtering": True,
        "scenario_rows_dropped_or_replaced": 0,
        "outcome_solves": {"power_flow": 0, "DAE": 0, "eigenvalue": 0},
        "epsilon_sampler_rejections_for_declared_final_support_only": rejected_proposals,
        "epsilon_sampler_rejections_used_any_scientific_outcome": False,
        "epsilon_final_max_abs_weighted_mean_pu": max_abs_weighted_mean,
        "epsilon_final_pooled_distribution_pu": summarize(np.asarray(final_eps_all, dtype=np.float64)),
        "epsilon_raw_truncated_pooled_distribution_pu": summarize(np.asarray(raw_eps_all, dtype=np.float64)),
        "epsilon_final_by_load_bus_pu": eps_by_bus,
        "load_scale_distribution": summarize(np.asarray(load_scales, dtype=np.float64)),
        "P_L_total_distribution_pu": summarize(np.asarray(load_p_totals, dtype=np.float64)),
        "Q_L_total_distribution_pu": summarize(np.asarray(load_q_totals, dtype=np.float64)),
        "P_L_total_distribution_MW": summarize(np.asarray(load_p_totals, dtype=np.float64) * float(config["system"]["base_mva"])),
        "Q_L_total_distribution_Mvar": summarize(np.asarray(load_q_totals, dtype=np.float64) * float(config["system"]["base_mva"])),
        "P_G_redispatch_total_distribution_pu": summarize(np.asarray([float(row["delta_P_G_allocated_total_pu"]) for row in rows], dtype=np.float64)),
        "P_G_redispatch_total_distribution_MW": summarize(np.asarray([float(row["delta_P_G_allocated_total_pu"]) for row in rows], dtype=np.float64) * float(config["system"]["base_mva"])),
        "generator_redispatch_by_bus_pu": pg_bus_summary,
        "redispatch_capacity_shortfall_scenario_count": redispatch_shortfalls,
        "slack": {
            "bus": int(config["slack_policy"]["bus"]),
            "base_scheduled_p0_pu_from_case": float(config["base_case_snapshot"]["slack_pu"]["p0_pu"]),
            "actual_scenario_slack_outputs": "not computed; require power flow",
            "role": "absorb only real/reactive residual required by frozen slack-bus PF convention, subject to frozen limits",
        },
        "branch_thermal_feasibility": "NOT_SCREENED: frozen IEEE-39 JSON contains no branch thermal ratings",
        "pre_run_only_note": "This file summarizes generated inputs and arithmetic redispatch only; it contains no PF, DAE, Jacobian, eigenvalue, or stability outcomes.",
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")

    checksums = {
        "ticket": "IAS26-060",
        "phase": "PRE_RUN_INPUT_GENERATION_ONLY",
        "sha256": {
            CONFIG_REL.as_posix(): sha256_file(config_path),
            MANIFEST_REL.as_posix(): sha256_file(manifest_path),
        },
        "summary_path": SUMMARY_REL.as_posix(),
        "generator_sha256": sha256_file(repo / GENERATOR_REL),
    }
    hashes_path.write_text(json.dumps(checksums, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "ticket": "IAS26-060",
        "phase": "PRE_RUN_INPUT_GENERATION_ONLY",
        "status": "FROZEN_PRE_RUN_AWAITING_USER_REVIEW",
        "scenario_count": len(rows),
        "manifest": str(manifest_path),
        "summary": str(summary_path),
        "sha256": checksums["sha256"],
        "power_flow_DAE_eigenvalue_solves": 0,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
