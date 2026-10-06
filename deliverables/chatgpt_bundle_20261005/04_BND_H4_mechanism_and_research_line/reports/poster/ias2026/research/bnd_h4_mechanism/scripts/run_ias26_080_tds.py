"""Frozen IAS26-080 nonlinear phasor DAE trajectories, using audited G2 solver."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
REPO = HERE.parents[6]
RESEARCH = REPO / "reports/poster/ias2026/research"
EXP = RESEARCH / "experiments"
sys.path.insert(0, str(RESEARCH / "src"))
sys.path.insert(0, str(EXP))

import G2_tds as g2  # noqa: E402
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case  # noqa: E402
from ibr_cycles.models.ieee39_devices import ConverterParameters  # noqa: E402
from ibr_cycles.models.ieee39_network import load_network  # noqa: E402


BASE_G = 0.03625
RETUNED_G = 0.25
K = 1.425
T = 1.5
LEAK = 0.05
H4 = (30, 33, 35, 37)
H4_NOMINAL_ALPHA = 0.12700646782830968
H4_RETUNED_NOMINAL_ALPHA = -0.017384676061556643
TRIPLE_NOMINAL_ALPHA = -0.20468769022203584


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def make_scenario_network(base, row):
    active = {int(k): float(v) for k, v in json.loads(row.P_L_s_by_bus_pu_json).items()}
    reactive = {int(k): float(v) for k, v in json.loads(row.Q_L_s_by_bus_pu_json).items()}
    loads = {bus: complex(active[bus], reactive[bus]) for bus in active}
    schedule = {int(k): float(v) for k, v in json.loads(row.P_G_scheduled_by_bus_pu_json).items()}
    pv = {bus: dict(spec) for bus, spec in base.pv.items()}
    for bus, p in schedule.items():
        pv[bus]["p"] = p
    return replace(base, loads=loads, pv=pv)


def solve_python_case(network, members: tuple[int, ...], g: float):
    return solve_case(
        ReplacementPlan.of({bus: 1.0 for bus in members}),
        network=network,
        converter=ConverterParameters(voltage_control=True, voltage_gain=g, voltage_leak=LEAK),
        machine_scaling={"ka": K, "ta": T},
        tol=1e-9,
    )


def critical_lambda(case) -> complex:
    values = np.linalg.eigvals(case.system.A)
    keep = np.flatnonzero(np.abs(values) >= 1e-3)
    if not len(keep):
        raise RuntimeError("empty non-gauge spectrum")
    return complex(values[keep[np.argmax(values[keep].real)]])


def common_observable_signal(sol, case, t_end: float, trace: pd.DataFrame) -> np.ndarray:
    labels = case.system.labels
    i = labels.index("omega_sg34")
    j = labels.index("omega_sg39")
    n = int(min(6000, max(1500, t_end * 50)))
    times = np.linspace(0.0, t_end, n)
    stride = max(1, n // 1500)
    states = sol(times)
    signal = states[i] - states[j]
    if len(trace) != len(signal[::stride]) or not np.allclose(trace.t.to_numpy(), times[::stride], rtol=0, atol=1e-10):
        raise RuntimeError("G2 trace sampling did not match the stored common-observable samples")
    return signal[::stride]


def fit_r2(trace: pd.DataFrame, alpha: float, frequency: float) -> float:
    if not np.isfinite(alpha) or not np.isfinite(frequency) or frequency <= 0:
        return math.nan
    t = trace.t.to_numpy(dtype=float)
    y = trace.signal.to_numpy(dtype=float)
    spread = trace.spread.to_numpy(dtype=float)
    in_linear = np.maximum.accumulate(spread) <= g2.LINEAR
    t_lin = float(t[in_linear][-1]) if np.any(in_linear) else 0.0
    start = min(g2.FIT_START, t_lin / 3.0)
    mask = in_linear & (t >= start)
    if int(mask.sum()) < 20:
        return math.nan
    tau = t[mask] - t[mask][0]
    yy = y[mask]
    phase = 2 * np.pi * frequency * tau
    env = np.exp(alpha * tau)
    design = np.column_stack([env * np.cos(phase), env * np.sin(phase), np.ones_like(tau)])
    coeff, *_ = np.linalg.lstsq(design, yy, rcond=None)
    residual = yy - design @ coeff
    den = float(np.sum((yy - yy.mean()) ** 2))
    return float(1.0 - np.sum(residual**2) / den) if den > 0 else math.nan


def selected_tasks(run: Path, tds_config: dict, manifest: pd.DataFrame) -> list[dict]:
    tasks: list[dict] = []
    base = load_network()
    # Canonical six-point disturbance grid for each of the three frozen cases.
    canonical = tds_config["canonical_cases"]
    for spec in canonical:
        for disturbance in ("D1", "D2"):
            for multiplier in tds_config["canonical_amplitude_multipliers"]:
                tasks.append({
                    "trajectory_id": f"CANONICAL_{spec['case_id']}_{disturbance}_A{multiplier:g}",
                    "phase": "CANONICAL", "scenario_id": "CANONICAL_NOMINAL",
                    "stratum": "NOMINAL", "case_id": spec["case_id"],
                    "members": tuple(spec["members"]), "g": float(spec["g"]),
                    "disturbance": disturbance, "amplitude_multiplier": float(multiplier),
                    "network": base,
                })
    holdout = pd.read_csv(run / "inputs/IAS26-080_TDS_HOLDOUT_V1.csv")
    mc_cases = pd.read_csv(run / "derived/MC_CASES.csv")
    for entry in holdout.itertuples(index=False):
        row = manifest.loc[entry.scenario_id]
        network = make_scenario_network(base, row)
        case_specs = [
            ("H4_ORIGINAL", H4, BASE_G),
            ("H4_RETUNED", H4, RETUNED_G),
            ("MOST_CRITICAL_PROPER_SUBSET", tuple(int(x) for x in entry.most_critical_proper_subset_id.split("+") if x != "BASE"), BASE_G),
        ]
        # An empty subset is encoded as BASE.
        if entry.most_critical_proper_subset_id == "BASE":
            case_specs[2] = ("MOST_CRITICAL_PROPER_SUBSET", (), BASE_G)
        for case_id, members, g in case_specs:
            tasks.append({
                "trajectory_id": f"HOLDOUT_{entry.scenario_id}_{case_id}_D2_A1",
                "phase": "HOLDOUT", "scenario_id": entry.scenario_id,
                "stratum": entry.tds_stratum, "case_id": case_id,
                "members": members, "g": g, "disturbance": "D2",
                "amplitude_multiplier": 1.0, "network": network,
                "most_critical_proper_subset_alpha": float(entry.most_critical_proper_subset_alpha),
            })
    if len(tasks) != 90 or len({t["trajectory_id"] for t in tasks}) != 90:
        raise RuntimeError(f"frozen TDS task set invalid: {len(tasks)} trajectories")
    return tasks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    args = parser.parse_args()
    run = args.run_root.resolve()
    tds_config_path = run / "config/IAS26-080_TDS_OPERATING_V1.json"
    tds_selection_path = run / "inputs/IAS26-080_TDS_HOLDOUT_V1.csv"
    config = json.loads(tds_config_path.read_text(encoding="utf-8"))
    mc_config = json.loads((run / "config/IAS26-060_MC_OPERATING_V1.json").read_text(encoding="utf-8"))
    if sha256(tds_selection_path) != config["selection_sha256"]:
        raise RuntimeError("frozen IAS26-080 selection hash mismatch")
    if sha256(run / "derived/MC_CASES.csv") != config["mc_case_table_sha256"]:
        raise RuntimeError("IAS26-060 case table changed after TDS selection")
    for rel, expected in config["historical_and_model_source_sha256"].items():
        if sha256(REPO / rel) != expected:
            raise RuntimeError(f"audited TDS/model source hash changed: {rel}")

    raw_dir = run / "raw/tds/traces"
    metrics_dir = run / "derived/tds/metrics"
    logs_dir = run / "logs/tds"
    for folder in (raw_dir, metrics_dir, logs_dir):
        folder.mkdir(parents=True, exist_ok=True)
    aggregate = run / "tables/IAS26-080_TDS_RESULTS.csv"
    if aggregate.exists():
        raise RuntimeError("final TDS aggregate exists; refusing overwrite")

    manifest = pd.read_csv(run / "inputs/IAS26-060_SCENARIOS_V1.csv").set_index("scenario_id")
    tasks = selected_tasks(run, config, manifest)
    metric_files = {p.stem for p in metrics_dir.glob("*.json")}
    planned = {t["trajectory_id"] for t in tasks}
    unexpected = metric_files - planned
    if unexpected:
        raise RuntimeError(f"unplanned TDS metric artifacts found: {sorted(unexpected)[:3]}")

    # Canonical direct-DAE regression is a hard software/model gate, not a retuning step.
    nominal = load_network()
    regression = []
    for label, members, g, expected in (
        ("proper", (30, 33, 35), BASE_G, TRIPLE_NOMINAL_ALPHA),
        ("H4_original", H4, BASE_G, H4_NOMINAL_ALPHA),
        ("H4_retuned", H4, RETUNED_G, H4_RETUNED_NOMINAL_ALPHA),
    ):
        case = solve_python_case(nominal, members, g)
        measured = float(np.linalg.eigvals(case.system.A)[np.abs(np.linalg.eigvals(case.system.A)) >= 1e-3].real.max())
        regression.append({"case": label, "alpha_python": measured, "alpha_frozen": expected, "abs_error": abs(measured - expected)})
    if max(x["abs_error"] for x in regression) > 2e-6:
        raise RuntimeError(f"canonical Python model regression failed: {regression}")
    pd.DataFrame(regression).to_csv(run / "tables/IAS26-080_CANONICAL_REGRESSION.csv", index=False)

    cache: dict[tuple[str, tuple[int, ...], float], object] = {}
    current_case = {"case": None}
    original_build = g2.build
    original_simulate = g2.simulate
    captured: dict[str, object] = {}

    def build_from_cache(bench, members, point, condenser):
        if bench != "IEEE-39" or condenser is not None:
            raise RuntimeError("unexpected TDS model request")
        key = (current_case["scenario_id"], tuple(members), float(point["g"]))
        if key not in cache:
            task = current_case["task"]
            cache[key] = solve_python_case(task["network"], tuple(members), float(point["g"]))
        return cache[key]

    def capture_simulate(case, disturbance, horizon, speeds):
        result = original_simulate(case, disturbance, horizon, speeds)
        captured["simulation"] = (result[0], result[1], result[2], case)
        return result

    g2.build = build_from_cache
    g2.simulate = capture_simulate
    started_all = time.perf_counter()
    try:
        for idx, task in enumerate(tasks, 1):
            tid = task["trajectory_id"]
            trace_path = raw_dir / f"{tid}.csv"
            metric_path = metrics_dir / f"{tid}.json"
            if metric_path.exists() and trace_path.exists():
                print(f"IAS26_080_TDS_SKIP_EXISTING {idx}/90 {tid}", flush=True)
                continue
            if metric_path.exists() != trace_path.exists():
                raise RuntimeError(f"partial trajectory artifact for {tid}; preserve and inspect before resuming")
            current_case["scenario_id"] = task["scenario_id"]
            current_case["task"] = task
            if task["phase"] == "CANONICAL":
                scenario_key = "CANONICAL_NOMINAL"
            else:
                scenario_key = task["scenario_id"]
            current_case["scenario_id"] = scenario_key
            # Cache keys for canonical cases are unambiguous and holdout keys are exact IDs.
            multiplier = float(task["amplitude_multiplier"])
            g2.KICK = 1e-4 * multiplier
            g2.STEP = 0.02 * multiplier
            members = task["members"]
            point = {"g": task["g"], "k": K}
            case = build_from_cache("IEEE-39", members, point, None)
            lam = critical_lambda(case)
            task_tuple = ("IEEE-39", task["case_id"], members, point, None, task["disturbance"])
            captured.pop("simulation", None)
            started = time.perf_counter()
            status = "COMPLETED"
            failure = ""
            try:
                result, trace = g2.run(task_tuple)
            except Exception as exc:
                result = {"benchmark": "IEEE-39", "case": task["case_id"], "lambda_re": lam.real, "lambda_f_hz": abs(lam.imag)/(2*np.pi), "run_status": "UNCLASSIFIED_EXCEPTION"}
                trace = pd.DataFrame(columns=["t", "signal", "spread"])
                status = "UNCLASSIFIED_EXCEPTION"
                failure = f"{type(exc).__name__}: {exc}"
            if "simulation" in captured and len(trace):
                sol, t_end, sim_status, case_captured = captured["simulation"]
                status = sim_status
                common_signal = common_observable_signal(sol, case_captured, float(t_end), trace)
                output_trace = pd.DataFrame({
                    "time_s": trace.t.to_numpy(dtype=float),
                    "mode_specific_relative_speed_signal_pu": trace.signal.to_numpy(dtype=float),
                    "common_relative_speed_omega_sg34_minus_omega_sg39_pu": common_signal,
                    "machine_speed_spread_pu": trace.spread.to_numpy(dtype=float),
                })
            else:
                output_trace = pd.DataFrame(columns=["time_s", "mode_specific_relative_speed_signal_pu", "common_relative_speed_omega_sg34_minus_omega_sg39_pu", "machine_speed_spread_pu"])
            alpha_tds = float(result.get("primary_re", math.nan))
            frequency_tds = float(result.get("primary_f_hz", math.nan))
            alpha_eig = float(result.get("lambda_re", lam.real))
            frequency_eig = float(result.get("lambda_f_hz", abs(lam.imag)/(2*np.pi)))
            fit_r2_value = fit_r2(trace, float(result.get("pencil_re", math.nan)), float(result.get("pencil_f_hz", math.nan))) if len(trace) else math.nan
            if status != "COMPLETED":
                tds_class = status
            elif not np.isfinite(alpha_tds):
                tds_class = "NO_MODAL_FIT"
            elif alpha_tds < -1e-8:
                tds_class = "STABLE"
            elif alpha_tds > 1e-8:
                tds_class = "UNSTABLE"
            else:
                tds_class = "INDETERMINATE"
            sign_agreement = bool(np.sign(alpha_tds) == np.sign(alpha_eig)) if np.isfinite(alpha_tds) else None
            trace_meta = {
                "trajectory_id": tid, "run_id": run.name, "scenario_id": task["scenario_id"],
                "phase": task["phase"], "case_id": task["case_id"], "members": list(members),
                "g": task["g"], "disturbance": task["disturbance"],
                "amplitude_multiplier": multiplier,
                "base_pulse_fraction": g2.STEP / multiplier if task["disturbance"] == "D2" else None,
                "base_speed_kick_pu": g2.KICK / multiplier if task["disturbance"] == "D1" else None,
                "pulse_duration_s": g2.PULSE_S if task["disturbance"] == "D2" else 0.0,
                "selected_observable": result.get("observable", "NOT_AVAILABLE"),
                "common_plot_observable": "omega_sg34-omega_sg39",
                "status": status,
                "failure": failure,
                "samples": int(len(output_trace)),
                "alpha_eig_s_inv": alpha_eig,
                "frequency_eig_hz": frequency_eig,
                "alpha_tds_s_inv": alpha_tds,
                "frequency_tds_hz": frequency_tds,
                "alpha_abs_error": abs(alpha_tds-alpha_eig) if np.isfinite(alpha_tds) else math.nan,
                "frequency_abs_error_hz": abs(frequency_tds-frequency_eig) if np.isfinite(frequency_tds) else math.nan,
                "sign_agreement": sign_agreement,
                "fit_window_s": result.get("pencil_window_s", math.nan),
                "fit_r2_diagnostic": fit_r2_value,
                "max_spread_pu": result.get("max_spread_pu", math.nan),
                "envelope_growth_second_half_s_inv": result.get("env_growth_2nd_half", math.nan),
                "fft_frequency_hz": result.get("fft_f_hz", math.nan),
                "t_end_s": result.get("t_end_s", math.nan),
                "predicted_verdict": "UNSTABLE" if alpha_eig > 1e-8 else ("STABLE" if alpha_eig < -1e-8 else "INDETERMINATE"),
                "trajectory_verdict": tds_class,
                "verdict_agreement": ((tds_class == "UNSTABLE") == (alpha_eig > 1e-8)) if status == "COMPLETED" and tds_class in ("STABLE", "UNSTABLE") else None,
                "small_signal_fit_status": "FIT_AVAILABLE" if np.isfinite(alpha_tds) else "FIT_UNAVAILABLE_RETAINED",
                "nonlinear_or_invalid_flag": status if status != "COMPLETED" else ("NO_MODAL_FIT" if not np.isfinite(alpha_tds) else "NONE"),
                "python_pf_residual": case.dae.power_flow.max_mismatch,
                "python_dae_f_residual": float(np.linalg.norm(case.dae.f(case.equilibrium.x, case.equilibrium.z, {}), np.inf)),
                "python_dae_g_residual": float(np.linalg.norm(case.dae.g(case.equilibrium.x, case.equilibrium.z, {}), np.inf)),
                "elapsed_s": time.perf_counter() - started,
                "git_commit": mc_config["source_provenance"]["git_head"],
            }
            output_trace.to_csv(trace_path, index=False)
            metric_path.write_text(json.dumps(trace_meta, indent=2, allow_nan=True) + "\n", encoding="utf-8")
            print(f"IAS26_080_TDS_COMPLETE {idx}/90 id={tid} status={status} alpha_eig={alpha_eig:.6g} alpha_tds={alpha_tds:.6g} freq={frequency_tds:.5g}Hz elapsed={trace_meta['elapsed_s']:.1f}s", flush=True)
    finally:
        g2.build = original_build
        g2.simulate = original_simulate

    metric_paths = sorted(metrics_dir.glob("*.json"))
    metrics = pd.DataFrame([json.loads(p.read_text(encoding="utf-8")) for p in metric_paths])
    metrics = metrics.sort_values("trajectory_id", kind="mergesort")
    if len(metrics) != 90 or metrics.trajectory_id.nunique() != 90:
        raise RuntimeError(f"TDS output integrity failure: {len(metrics)} complete metric records")
    metrics.to_csv(aggregate, index=False)
    fitted = metrics[np.isfinite(metrics.alpha_tds_s_inv)]
    errors_alpha = fitted.alpha_abs_error.to_numpy(dtype=float)
    errors_f = fitted.frequency_abs_error_hz.to_numpy(dtype=float)
    summary = {
        "ticket": "IAS26-080", "run_id": run.name, "N_planned": 90,
        "N_completed_records": int(len(metrics)), "N_canonical": int((metrics.phase == "CANONICAL").sum()),
        "N_holdout": int((metrics.phase == "HOLDOUT").sum()),
        "run_status_counts": metrics.status.value_counts().to_dict(),
        "fit_available_N": int(len(fitted)),
        "sign_agreement_N": int(fitted.sign_agreement.notna().sum()),
        "sign_agreement_rate": float(fitted.sign_agreement.mean()) if len(fitted) else math.nan,
        "alpha_mae_s_inv": float(errors_alpha.mean()) if len(errors_alpha) else math.nan,
        "alpha_median_absolute_error_s_inv": float(np.median(errors_alpha)) if len(errors_alpha) else math.nan,
        "alpha_q95_absolute_error_s_inv": float(np.quantile(errors_alpha, .95)) if len(errors_alpha) else math.nan,
        "frequency_mae_hz": float(errors_f.mean()) if len(errors_f) else math.nan,
        "frequency_median_absolute_error_hz": float(np.median(errors_f)) if len(errors_f) else math.nan,
        "frequency_q95_absolute_error_hz": float(np.quantile(errors_f, .95)) if len(errors_f) else math.nan,
        "mean_fit_r2_diagnostic": float(fitted.fit_r2_diagnostic.mean()) if len(fitted) else math.nan,
        "same_model_tds_not_emt": True,
        "elapsed_wall_s": time.perf_counter() - started_all,
        "no_scenario_replacements": True,
    }
    (run / "derived/tds/IAS26-080_TDS_SUMMARY.json").write_text(json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, allow_nan=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
