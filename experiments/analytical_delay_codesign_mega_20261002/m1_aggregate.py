"""Select the largest fully validated zero-delay candidate actually executed."""

from pathlib import Path
import json
import tomllib
import pandas as pd

HERE = Path(__file__).resolve().parent


def read_toml(path):
    return tomllib.loads(path.read_text(encoding="utf-8"))


records = []
seed = read_toml(HERE / "Q0_RESULT.toml")
records.append({
    "candidate_id": "eta_000_seed", "eta": 0.0, "result": seed,
    "design_path": HERE / "seed_uniform_875.toml", "output_dir": HERE,
})
for candidate in sorted((HERE / "m1_candidates").glob("eta_*.toml")):
    out = HERE / "m1_validations" / candidate.stem
    result_file = out / "Q0_RESULT.toml"
    if result_file.is_file():
        result = read_toml(result_file)
    else:
        result = None
    records.append({
        "candidate_id": candidate.stem,
        "eta": float(read_toml(candidate)["eta"]),
        "result": result,
        "design_path": candidate,
        "output_dir": out,
    })

rows = []
for rec in records:
    d = rec["result"]
    event_file = rec["output_dir"] / "Q0_EVENT_METRICS.csv"
    events = pd.read_csv(event_file) if event_file.is_file() else pd.DataFrame()
    feasible = bool(d is not None and d["all_five_events_pass"] and
                    len(events) == 5 and bool(events.complete.all()) and bool(events["pass"].all()) and
                    d["critical_real_part_s_inv"] <= -0.05 and
                    d["equilibrium_residual_inf"] < 1e-8)
    rows.append({
        "candidate_id": rec["candidate_id"], "eta": rec["eta"],
        "simulation_complete": bool(len(events)==5 and events.complete.all()),
        "all_five_events_pass": bool(len(events)==5 and events["pass"].all()),
        "fully_validated": feasible,
        "GFL_percent": d["GFL_percent"] if d else float("nan"),
        "GFL_MW": d["GFL_MW"] if d else float("nan"),
        "retained_SG_MW": d["retained_SG_MW"] if d else float("nan"),
        "critical_real_s_inv": d["critical_real_part_s_inv"] if d else float("nan"),
        "equilibrium_residual_inf": d["equilibrium_residual_inf"] if d else float("nan"),
        "max_frequency_deviation_Hz": events.F_peak_Hz.max() if len(events) else float("nan"),
        "max_RoCoF_Hz_s": events.RoCoF_peak_Hz_s.max() if len(events) else float("nan"),
        "min_actuator_slack": events.min_SG_actuator_fraction_slack.min() if len(events) else float("nan"),
        "worst_frequency_event": str(events.loc[events.F_peak_Hz.idxmax(), "event"]) if len(events) else "",
        "design_path": rec["design_path"].relative_to(HERE).as_posix(),
        "output_dir": rec["output_dir"].relative_to(HERE).as_posix() if rec["output_dir"]!=HERE else ".",
        "status": "FULLY_VALIDATED" if feasible else ("NOT_EXECUTED" if len(events)==0 else "EVENT_FAIL_OR_INCOMPLETE"),
    })

audit = pd.DataFrame(rows).sort_values("eta")
audit.to_csv(HERE / "M1_CANDIDATE_AUDIT.csv", index=False)
feasible_rows = audit[audit.fully_validated]
if feasible_rows.empty:
    raise SystemExit("M1 has no fully validated zero-delay witness")
winner = feasible_rows.sort_values(["GFL_MW","candidate_id"],ascending=[False,True]).iloc[0]
best = pd.DataFrame([winner.to_dict()])
best.to_csv(HERE / "M1_ZERO_DELAY_BEST_VALIDATED.csv", index=False)
chosen = next(rec for rec in records if rec["candidate_id"] == winner.candidate_id)
design = read_toml(chosen["design_path"])
lines = [f"{k} = {json.dumps(design[k])}" for k in ("rho","Kp","Ki")]
lines += [
    f"candidate_id = {json.dumps(winner.candidate_id)}",
    f"eta = {float(winner.eta)}",
    f"GFL_MW = {float(winner.GFL_MW)}",
    f"GFL_percent = {float(winner.GFL_percent)}",
    f"retained_SG_MW = {float(winner.retained_SG_MW)}",
    "status = \"BEST_FULLY_VALIDATED_FOUND_NOT_OPTIMUM\"",
    "global_or_local_optimality_certified = false",
]
(HERE / "M1_ZERO_DELAY_DESIGN.toml").write_text("\n".join(lines)+"\n",encoding="utf-8")
out = chosen["output_dir"]
pd.read_csv(out / "Q0_EVENT_METRICS.csv").to_csv(HERE / "M1_ZERO_DELAY_EVENTS.csv",index=False)
pd.read_csv(out / "Q0_FULL_PHYSICAL_SPECTRUM.csv").to_csv(HERE / "M1_ZERO_DELAY_POLES.csv",index=False)
print("M1_BEST",winner.candidate_id,"GFL_percent",winner.GFL_percent,"retained_SG_MW",winner.retained_SG_MW)
