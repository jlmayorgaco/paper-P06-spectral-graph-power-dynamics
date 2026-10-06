"""Summarize saved outputs; does not simulate or invent missing results."""
from pathlib import Path
import csv
import json
import hashlib
import tomllib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports/analytic_iteration_20261001/refined"


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def config(name):
    return tomllib.loads((OUT / name).read_text(encoding="utf-8"))


def main():
    data = rows(OUT / "evaluations.csv")
    bylabel = {r["label"]: r for r in data}
    candidates = {name: config(name + ".toml") for name in ("baseline", "joint_final", "fixed_gains_final")}
    histories = {name: rows(OUT / (name + "_history.csv")) for name in ("joint", "fixed_gains")}
    gate = config("gradient_gate_summary.toml")
    pd_path = OUT / "independent_pd/events.csv"
    independent = rows(pd_path) if pd_path.exists() else []
    refined = {"joint_final": bylabel.get("joint_refined"), "fixed_gains_final": bylabel.get("fixed_refined"),
               "baseline": bylabel["baseline"]}
    checks = []
    for row in independent:
        direct = refined[row["label"]]
        if direct is None:
            continue
        candidate_hash = hashlib.sha256((OUT / (row["label"] + ".toml")).read_bytes()).hexdigest()
        checks.append({"label": row["label"], "candidate_hash_matches": candidate_hash == row["sha256"],
                       "frequency_difference_Hz": abs(float(row["F"]) - float(direct["Fpeak_Hz"])),
                       "rocof_difference_Hz_s": abs(float(row["R"]) - float(direct["Rpeak_Hz_s"])),
                       "alpha_difference_s_inv": abs(float(row["alpha"]) - float(direct["alpha"])),
                       "frequency_pass": float(row["F"]) <= .5,
                       "rocof_pass": float(row["R"]) <= .5,
                       "modal_pass": float(row["alpha"]) <= -.05,
                       "voltage_pass": float(row["Vmin"]) >= .9 and float(row["Vmax"]) <= 1.1})
    holdouts = [r for r in data if r["label"].startswith("holdout_")]
    incomplete_holdouts = [{"file": p.name, **tomllib.loads(p.read_text(encoding="utf-8"))}
                           for p in sorted(OUT.glob("holdout_*_incomplete.toml"))]
    holdout_summary = [{"label": r["label"], "F_Hz": float(r["Fpeak_Hz"]),
                        "R_Hz_s": float(r["Rpeak_Hz_s"]), "max_normalized_constraint": float(r["max_constraint"]),
                        "frequency_pass": float(r["Fpeak_Hz"]) <= .5,
                        "rocof_pass": float(r["Rpeak_Hz_s"]) <= .5,
                        "actuator_margin_pass": float(r["limiter_fraction"]) >= .002,
                        "limiting_actuator_bus": int(r["limiter_bus"]), "limiting_actuator_kind": r["limiter_kind"],
                        "passes_declared_constraints": float(r["max_constraint"]) <= 0} for r in holdouts]
    summaries = {name: {"replacement_percent": c["replacement_percent"], "retained_SG_MW": c["retained_SG_MW"]}
                 for name, c in candidates.items()}
    summary = {"scope": "Finite-budget, finite-event, nonlinear experiment; not a capacity optimum or all-input guarantee",
               "candidates": summaries, "gradient_gate": gate, "independent_comparison": checks,
               "holdouts": holdout_summary, "incomplete_holdouts": incomplete_holdouts,
               "joint_extra_replacement_vs_start_MW": candidates["baseline"]["retained_SG_MW"] - candidates["joint_final"]["retained_SG_MW"],
               "joint_extra_replacement_vs_fixed_MW": candidates["fixed_gains_final"]["retained_SG_MW"] - candidates["joint_final"]["retained_SG_MW"],
               "max_explicit_step_residual": max(float(r["explicit_error"]) for h in histories.values() for r in h),
               "new_current_or_energy_certificate": False, "optimality_certified": False}
    summary["assessment"] = {
        "trajectory_gradient_gate_pass": gate["trajectory_gradient_relative_error"] < gate["threshold"],
        "explicit_step_algebra_pass": summary["max_explicit_step_residual"] < 1e-7,
        "all_accepted_design_steps_feasible": all(float(r["maxg"]) <= 0 for h in histories.values() for r in h),
        "joint_replacement_improved": summary["joint_extra_replacement_vs_start_MW"] > 0,
        "independent_validation_complete": len(checks) == 3,
        "independent_primary_case_pass": len(checks) == 3 and all(
            r["candidate_hash_matches"] and r["frequency_pass"] and r["rocof_pass"] and
            r["modal_pass"] and r["voltage_pass"] and r["frequency_difference_Hz"] < 1e-4 and
            r["rocof_difference_Hz_s"] < 1e-4 and r["alpha_difference_s_inv"] < 1e-5 for r in checks),
        "all_external_events_pass": len(holdout_summary) == 5 and not incomplete_holdouts and all(
            r["passes_declared_constraints"] for r in holdout_summary),
    }
    (OUT / "RESULT.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    operating = rows(ROOT / "reports/experiment_N/TABLE_N01_original_operating_point.csv")
    weights = {int(r["bus"]): float(r["P_gen_MW"]) for r in operating}
    parameter_rows = []
    for i, bus in enumerate(range(30, 40)):
        r = {"bus": bus, "original_dispatch_MW": weights[bus]}
        for label, candidate in candidates.items():
            for key in ("rho", "Kp", "Ki"):
                r[label + "_" + key] = candidate[key][i]
            r[label + "_SG_MW"] = weights[bus] * (1 - candidate["rho"][i])
        parameter_rows.append(r)
    with (OUT / "PARAMETERS_BY_BUS.csv").open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_rows[0].keys())
        writer.writeheader(); writer.writerows(parameter_rows)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.titleweight": "bold", "svg.fonttype": "none"})
    colors = {"baseline": "#687786", "joint": "#116d9b", "fixed_gains": "#b86425"}
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.3), layout="constrained")
    ax = axes[0, 0]
    for arm, label in (("joint", "Joint retention + PLL"), ("fixed_gains", "Retention, fixed PLL")):
        h = histories[arm]
        ax.plot([int(r["iteration"]) for r in h], [float(r["capacity_percent"]) for r in h],
                "o-", color=colors[arm], label=label)
    ax.set(xlabel="Accepted iteration", ylabel="Original dispatch replaced by GFL (%)", title="A  Explicit updates on the nonlinear model")
    ax.legend(frameon=False, fontsize=9);ax.grid(alpha=.15)

    for ax, metric, title, unit in ((axes[0, 1], "Fmax_Hz", "B  Frequency: all 39 buses", "Hz"),
                                   (axes[1, 0], "Rmax_Hz_s", "C  RoCoF: causal 0.5 s window", "Hz/s")):
        for filename, key, label in (("baseline_trajectory.csv", "baseline", "Start"),
                                     ("fixed_refined_trajectory.csv", "fixed_gains", "Fixed PLL"),
                                     ("joint_refined_trajectory.csv", "joint", "Joint")):
            path = OUT / filename
            if not path.exists():
                continue
            tr = rows(path)
            ax.plot([float(r["time_after_event_s"]) for r in tr], [float(r[metric]) for r in tr],
                    color=colors[key], label=label, linewidth=1.5)
        ax.axhline(.5, color="#b8343e", linestyle="--", linewidth=1, label="Declared limit")
        ax.set(xlabel="Time after event (s)", ylabel="Maximum absolute value (" + unit + ")", title=title)
        ax.legend(frameon=False, fontsize=9);ax.grid(alpha=.15)

    ax = axes[1, 1]
    if holdouts:
        idx = np.arange(len(holdouts))
        ax.bar(idx - .18, [float(r["Fpeak_Hz"])/.5 for r in holdouts], .36, label="Frequency", color="#116d9b")
        ax.bar(idx + .18, [float(r["Rpeak_Hz_s"])/.5 for r in holdouts], .36, label="RoCoF", color="#729bab")
        actuator_fail = [j for j, r in enumerate(holdouts) if float(r["limiter_fraction"]) < .002]
        if actuator_fail:
            ax.scatter(actuator_fail, np.full(len(actuator_fail), 1.10), marker="x", color="#b8343e", s=55,
                       label="SG actuator margin fails", zorder=5)
        ax.axhline(1, color="#b8343e", linestyle="--", linewidth=1)
        ax.set_xticks(idx, ["Bus " + r["bus"] + "\n" + ("+" if float(r["delta"]) > 0 else "") + str(int(float(r["delta"]))) + " MW" for r in holdouts])
        ax.set(ylabel="Peak / declared limit", title="D  Joint design: external-event checks", ylim=(0,1.32))
        ax.legend(frameon=False, fontsize=9);ax.grid(axis="y", alpha=.15)
    else:
        ax.text(.5, .5, "External-event checks not yet available", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
    fig.suptitle("Analytical iterative co-design: a bounded IEEE-39 experiment", fontsize=16, weight="bold")
    fig.supxlabel("Design event: +100 MW ZIP setpoint at bus 8 | Physical DC convention | 60 s validation | No global optimality claim", fontsize=9)
    for ext in ("png", "svg", "pdf"):
        fig.savefig(OUT / ("EXPERIMENT." + ext), dpi=180)
    plt.close(fig)

    lines = ["EXPERIMENTO DE ITERACIÓN ANALÍTICA NO LINEAL — IEEE-39", "Julia 1.11.9; proyecto y manifiesto preservados.", "",
             "Los resultados son de un experimento finito, no un certificado general ni un óptimo.",
             "La prueba inicial con rho uniforme 0.90 falló F<=0.5 Hz; se conserva en el directorio padre.",
             "La comparación usa rho=0.875 y PLL nominales como punto inicial común.", "",
             "El control de gradiente inicial con dt=0.025 s obtuvo 2.187% de error y se detuvo.",
             "Se refinó a dt=0.0125 s sin relajar el umbral de 2%; el error bajó a 0.5525%.",
             "La prueba inicial de malla se conserva en ../seed_0875/.", "",
             "RESULTADOS"]
    for name, c in summaries.items():
        lines.append(f"{name}: GFL {c['replacement_percent']:.8f}% del despacho; SG {c['retained_SG_MW']:.8f} MW.")
    lines.extend([f"Reemplazo adicional conjunto frente al inicio: {summary['joint_extra_replacement_vs_start_MW']:.8f} MW.",
                  f"Diferencia conjunta frente a PLL fijos: {summary['joint_extra_replacement_vs_fixed_MW']:.8f} MW.",
                  "La diferencia es a presupuesto de iteraciones declarado; no compara óptimos globales.", "",
                  "ESTADO DE LAS COMPROBACIONES", json.dumps(summary["assessment"], indent=2), "",
                  "VALIDACIÓN DE FÓRMULAS", json.dumps(gate, indent=2),
                  f"Máxima discrepancia del paso explícito con el QP: {summary['max_explicit_step_residual']:.3e}.", "",
                  "VALIDACIÓN INDEPENDIENTE POWERDYNAMICS", json.dumps(checks, indent=2), "",
                  "EVENTOS EXTERNOS", json.dumps(holdout_summary, indent=2),
                  "EVENTOS EXTERNOS INCOMPLETOS", json.dumps(incomplete_holdouts, indent=2), "",
                  "ALCANCE", "Sensibilidades exactas de la discretización SDIRK contrastadas con integración adaptativa Rodas5P.",
                  "La eliminación algebraica conserva la dinámica no lineal, con cargas Z y arquitectura fija.",
                  "Todos los pasos aceptados se verifican con la trayectoria no lineal durante 60 s.",
                  "El margen modal es de pequeña señal. La simulación finita no demuestra estabilidad regional general.",
                  "Se supervisan tensión y margen de actuadores SG. Corriente y DC se registran sin inventar límites de hardware.",
                  "Un evento externo que falla limita la generalización; una integración incompleta no prueba inestabilidad.",
                  "Los parámetros del caso y las fuentes están en protocol.toml y source_hashes_before/after.toml.", "",
                  "REPRODUCIR", "julia --startup-file=no --project=. experiments/analytic_iteration_20261001/run_experiment.jl",
                  "julia --startup-file=no --project=. experiments/analytic_iteration_20261001/validate_pd.jl",
                  "python experiments/analytic_iteration_20261001/report.py"])
    (OUT / "RESULTADOS_ES.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps({"candidates": summary["candidates"], "assessment": summary["assessment"],
                      "joint_extra_replacement_vs_start_MW": summary["joint_extra_replacement_vs_start_MW"],
                      "joint_extra_replacement_vs_fixed_MW": summary["joint_extra_replacement_vs_fixed_MW"]}, indent=2))


if __name__ == "__main__":
    main()
