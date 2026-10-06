"""Materialize the best strictly feasible saved analytical continuation point."""
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "reports" / "experiment_E" / "tables"
with (TABLES / "TABLE_E01_replaceable_generators.csv").open(newline="", encoding="utf-8") as f:
    generators = [r for r in csv.DictReader(f) if r["replaceable"] == "true"]
with (TABLES / "TABLE_E13_anchor38_gain_continuation.csv").open(newline="", encoding="utf-8") as f:
    path = list(csv.DictReader(f))

feasible = [r for r in path if r["gauge_detected"] == "true"
            and float(r["spectral_abscissa"]) <= -0.05]
if not feasible:
    raise RuntimeError("No strictly feasible saved continuation point")
point = max(feasible, key=lambda r: float(r["replacement_MW"]))
rho = [float(x) for x in point["rho"].split(";")]
kp = [float(x) for x in point["Kp"].split(";")]
ki = [float(x) for x in point["Ki"].split(";")]
if not (len(generators) == len(rho) == len(kp) == len(ki) == 10):
    raise RuntimeError("candidate dimension mismatch")

kp_nom = 5 * 2 * math.pi
ki_nom = kp_nom**2 / 4
rows = []
for j, (g, r, kpj, kij) in enumerate(zip(generators, rho, kp, ki)):
    p = float(g["SG_dispatch_initial_MW"])
    q = float(g["SG_reactive_initial_Mvar"])
    if not (0 <= r <= 1 and .25*kp_nom-1e-9 <= kpj <= 4*kp_nom+1e-9
            and .25*ki_nom-1e-9 <= kij <= 4*ki_nom+1e-9):
        raise RuntimeError(f"out of frozen domain at bus {g['bus']}")
    rows.append(dict(bus=int(g["bus"]), initial_SG_MW=p, initial_SG_Mvar=q,
                     rho=r, GFL_MW=r*p, retained_SG_MW=(1-r)*p,
                     GFL_Mvar=r*q, retained_SG_Mvar=(1-r)*q,
                     Kp=kpj, Ki=kij, Kp_over_nom=kpj/kp_nom,
                     Ki_over_nom=kij/ki_nom,
                     gain_bound_hit=any(abs(v-target)<1e-7 for v, target in (
                         (kpj,.25*kp_nom),(kpj,4*kp_nom),
                         (kij,.25*ki_nom),(kij,4*ki_nom)))))

mw = sum(r["GFL_MW"] for r in rows)
retained = sum(r["retained_SG_MW"] for r in rows)
if abs(mw - float(point["replacement_MW"])) > 1e-6:
    raise RuntimeError("objective mismatch")
out = TABLES / "TABLE_E14_provisional_per_generator_design.csv"
with out.open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
payload = dict(status="PROVISIONAL_ANALYTIC_ONLY", source="anchor38_gain_continuation",
               source_step=int(point["step"]), replaceable_buses=[r["bus"] for r in rows],
               initial_SG_MW=mw+retained, GFL_replacement_MW=mw,
               retained_SG_MW=retained, replacement_fraction=mw/(mw+retained),
               spectral_abscissa=float(point["spectral_abscissa"]),
               required_spectral_abscissa=-.05,
               global_optimality_certified=False, detailed_model_validated=False,
               time_domain_validated=False, candidate_frozen_for_validation=False)
with (ROOT / "reports" / "experiment_E" / "PROVISIONAL_ANALYTIC_RESULT.json").open("w", encoding="utf-8") as f:
    json.dump(payload, f, indent=2)
    f.write("\n")
print("PROVISIONAL_REPLACEMENT_MW:", mw)
print("PROVISIONAL_RETAINED_SG_MW:", retained)
print("PROVISIONAL_FRACTION:", payload["replacement_fraction"])
print("PROVISIONAL_SPECTRAL_ABSCISSA:", payload["spectral_abscissa"])
