"""Freeze the analytically certified fixed-support candidate before PD."""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports" / "experiment_N"


def rows(name: str) -> list[dict[str, str]]:
    with (OUT / name).open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, entries: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(entries[0]))
        writer.writeheader()
        writer.writerows(entries)


def main() -> None:
    model = json.loads((OUT / "MODEL_FREEZE.json").read_text(encoding="utf-8"))
    candidates = rows("TABLE_N10_all_KKT_candidates.csv")
    lattice = rows("TABLE_N09_support_lattice_budget_screen.csv")
    if len(lattice) != 2048:
        raise SystemExit("candidate freeze rejected: incomplete lattice screen")
    if any(r["feasible"].lower() == "true" for r in lattice):
        raise SystemExit("candidate freeze rejected: better sampled support")
    qualified = [r for r in candidates if r["KKT"].lower() == "true"]
    if not qualified:
        raise SystemExit("candidate freeze rejected: no KKT candidate")
    qualified.sort(key=lambda r: (float(r["retained_SG_MW"]),
                                  r["gainset"] != "lowKp_highKi"))
    c = qualified[0]
    if float(c["alpha"]) > -0.05 or c["LICQ"].lower() != "true" \
            or c["SOSC"].lower() != "true":
        raise SystemExit("candidate freeze rejected: analytical certificate")
    op = rows("TABLE_N01_original_operating_point.csv")
    ptotal = sum(float(r["P_gen_MW"]) for r in op)
    bus = int(c["bus"])
    rho = [1.0] * 10
    rho[bus - 30] = float(c["rho"])
    k0p = 2 * math.pi * 5
    k0i = k0p * k0p / 4
    kp = [0.25 * k0p] * 10
    ki = [4 * k0i] * 10
    retained = float(c["retained_SG_MW"])
    if abs(sum(float(op[i]["P_gen_MW"]) * (1 - rho[i]) for i in range(10))
           - retained) > 1e-8:
        raise SystemExit("candidate freeze rejected: dispatch mismatch")
    lines = [
        f'model_sha = "{model["MODEL_SHA"]}"',
        'status = "LOCAL_KKT_CERTIFIED_FIXED_SUPPORT"',
        'branch_completeness_certified = false',
        'global_certified = false',
        'sigma_req_per_s = 0.05',
        f'support_bus = {bus}',
        f'alpha_analytic_per_s = {float(c["alpha"]):.17g}',
        f'active_pole_real_per_s = {float(c["active_real"]):.17g}',
        f'active_pole_imag_per_s = {float(c["active_imag"]):.17g}',
        f'retained_SG_MW = {retained:.17g}',
        f'converted_GFL_MW = {ptotal-retained:.17g}',
        f'GFL_fraction = {(ptotal-retained)/ptotal:.17g}',
        f'spectral_multiplier = {float(c["mu"]):.17g}',
        f'KKT_stationarity_residual = {float(c["stationarity_residual"]):.17g}',
        f'KKT_complementarity_residual = {float(c["complementarity_residual"]):.17g}',
        'rho = [' + ', '.join(f'{v:.17g}' for v in rho) + ']',
        'Kp = [' + ', '.join(f'{v:.17g}' for v in kp) + ']',
        'Ki = [' + ', '.join(f'{v:.17g}' for v in ki) + ']',
    ]
    path = OUT / "Z_N_NOMINAL_FINAL.toml"
    content = "\n".join(lines) + "\n"
    if path.exists() and path.read_text(encoding="utf-8") != content:
        raise SystemExit("candidate freeze rejected: existing frozen candidate differs")
    path.write_text(content, encoding="utf-8")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    (OUT / (path.name + ".sha256")).write_text(digest + "\n", encoding="utf-8")
    write_csv(OUT / "TABLE_N11_nominal_best_candidate.csv", [
        {"support_bus": bus, "retained_SG_MW": retained,
         "converted_GFL_MW": ptotal-retained,
         "GFL_fraction": (ptotal-retained)/ptotal,
         "alpha_analytic": float(c["alpha"]), "KKT": True,
         "candidate_sha256": digest, "global_certified": False}])
    write_csv(OUT / "TABLE_N12_local_global_certificate.csv", [
        {"support_bus": bus, "primal_residual": c["primal_residual"],
         "stationarity_residual": c["stationarity_residual"],
         "complementarity_residual": c["complementarity_residual"],
         "spectral_multiplier": c["mu"], "LICQ": c["LICQ"],
         "SOSC": "STRICT_BOX_EMPTY_CRITICAL_CONE", "local_KKT_certified": True,
         "support_architectures_screened": 1024,
         "single_support_branches_corrected": len(candidates),
         "pair_budget_points": len(rows("TABLE_N09_pair_budget_audit.csv")),
         "branch_completeness_certified": False,
         "global_certified": False}])
    print("Z_NOMINAL_FROZEN", path.name, digest)
    print("RETAINED_SG_MW", retained)
    print("GFL_MW", ptotal-retained)


if __name__ == "__main__":
    main()
