from __future__ import annotations

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
E03 = ROOT / "experiments" / "tx3" / "E03_finite_externality"
E05B = ROOT / "experiments" / "tx3" / "E05B_dynamic_spectral_shift"
sys.path[:0] = [str(ROOT), str(E03), str(E05B)]

from campaign_core import action_vector, evaluate_operator, log_frequency_rule, operating_points, operator_logdet  # noqa: E402
from dynamic_spectral_shift import (  # noqa: E402
    all_coalitions,
    all_vertices,
    component_externality,
    integrated_pole_factor,
    mobius_vertices,
    parse_subset_key,
    pole_logdet,
    reduced_spectrum,
    subset_key,
)


ARTIFACT = ROOT / "artifacts" / "tx3" / "E05B_dynamic_spectral_shift"
NUMERICAL = ARTIFACT / "preregistration" / "E05B_NUMERICAL_FREEZE.json"
HOLDOUT_FREEZE = ARTIFACT / "preregistration" / "E05B_HOLDOUT_FREEZE.json"
CANDIDATE_FREEZE = ARTIFACT / "preregistration" / "E05B_CANDIDATE_CONTOUR_FREEZE.json"
E03_DIRECT = ROOT / "artifacts" / "tx3" / "E03_E04" / "direct"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def e03_reference(op_id: str) -> dict[tuple[int, ...], dict[str, Any]]:
    split = "discovery" if int(op_id.removeprefix("OP")) < 16 else "holdout"
    frame = pd.read_parquet(E03_DIRECT / f"{split}_vertices.parquet")
    frame = frame[frame["operating_point_id"] == op_id]
    return {parse_subset_key(str(row.subset)): row._asdict() for row in frame.itertuples(index=False)}


def worker(task: dict[str, Any]) -> dict[str, Any]:
    op = task["operating_point"]
    population = task["population"]
    numerical = task["numerical"]
    coarse_f, coarse_w = log_frequency_rule(
        *numerical["primary_frequency_band_hz"], int(numerical["coarse_frequency_nodes"])
    )
    fine_f, fine_w = log_frequency_rule(
        *numerical["primary_frequency_band_hz"], int(numerical["reference_frequency_nodes"])
    )
    sigma = float(numerical["sigma"])
    references = e03_reference(op["operating_point_id"]) if population == "development" else {}
    audit = numerical["pointwise_factorization_audit"]
    audit_subsets = {parse_subset_key(value) for value in audit["subsets"]}
    audit_frequencies = np.asarray(audit["frequencies_hz"], dtype=float)
    audit_this_op = op["operating_point_id"] in audit["operating_points"]
    started = time.perf_counter()
    records: dict[tuple[int, ...], dict[str, Any]] = {}
    spectra: list[np.ndarray] = []
    audit_rows: list[dict[str, Any]] = []
    for subset in all_vertices():
        point = evaluate_operator(op, action_vector({f"A{action}": 1.0 for action in subset}))
        record: dict[str, Any] = {
            "operating_point_id": op["operating_point_id"],
            "population": population,
            "subset_key": subset_key(subset),
            "order": len(subset),
            "status": point.status,
            **point.metadata,
        }
        if point.status != "SUCCESS":
            spectra.append(np.full(140, np.nan + 1j * np.nan, dtype=np.complex128))
            records[subset] = record
            continue
        poles, algebraic, mass = reduced_spectrum(point)
        spectra.append(poles)
        pole_coarse = integrated_pole_factor(poles, coarse_f, coarse_w, sigma)
        pole_fine = integrated_pole_factor(poles, fine_f, fine_w, sigma)
        record.update(
            {
                "algebraic_factor": algebraic,
                "mass_factor": mass,
                "pole_factor_coarse": pole_coarse,
                "pole_factor_fine": pole_fine,
                "full_factorized_coarse": algebraic + mass + pole_coarse,
                "full_factorized_fine": algebraic + mass + pole_fine,
                "finite_pole_count": len(poles),
            }
        )
        if references:
            reference = references[subset]
            record.update(
                {
                    "full_e03_coarse": float(reference["primary_coarse"]),
                    "full_e03_fine": float(reference["primary_fine"]),
                    "factorization_residual_coarse": algebraic + mass + pole_coarse - float(reference["primary_coarse"]),
                    "factorization_residual_fine": algebraic + mass + pole_fine - float(reference["primary_fine"]),
                }
            )
        if audit_this_op and subset in audit_subsets:
            direct, pivot = operator_logdet(point, audit_frequencies, sigma)
            factorized = algebraic + mass + pole_logdet(poles, audit_frequencies, sigma)
            for frequency, full, separated in zip(audit_frequencies, direct, factorized, strict=True):
                audit_rows.append(
                    {
                        "operating_point_id": op["operating_point_id"],
                        "population": population,
                        "subset_key": subset_key(subset),
                        "frequency_hz": float(frequency),
                        "full_sparse_logdet": float(full),
                        "factorized_logdet": float(separated),
                        "residual": float(full - separated),
                        "minimum_pivot_ratio": float(pivot),
                    }
                )
        records[subset] = record

    case_rows: list[dict[str, Any]] = []
    for coalition in all_coalitions():
        vertices = [records[subset] for subset in mobius_vertices(coalition)]
        valid = all(record["status"] == "SUCCESS" for record in vertices)
        row: dict[str, Any] = {
            "operating_point_id": op["operating_point_id"],
            "population": population,
            "coalition_key": subset_key(coalition),
            "order": len(coalition),
            "evaluation_status": "SUCCESS" if valid else "OPERATOR_FAIL",
            "resolved_pole_externality": False,
        }
        if valid:
            for field, output in (
                ("algebraic_factor", "delta_phi_algebraic"),
                ("mass_factor", "delta_phi_mass"),
                ("pole_factor_coarse", "delta_phi_poles_coarse"),
                ("pole_factor_fine", "delta_phi_poles"),
                ("full_factorized_fine", "delta_phi_full_factorized"),
            ):
                row[output] = component_externality(records, coalition, field)
            component_sum = row["delta_phi_algebraic"] + row["delta_phi_mass"] + row["delta_phi_poles"]
            row["component_closure_residual"] = row["delta_phi_full_factorized"] - component_sum
            if references:
                row["delta_phi_full_e03"] = component_externality(records, coalition, "full_e03_fine")
                row["full_externality_factorization_residual"] = row["delta_phi_full_factorized"] - row["delta_phi_full_e03"]
            epsilon_frequency = abs(row["delta_phi_poles"] - row["delta_phi_poles_coarse"])
            epsilon_eig = (2 ** len(coalition)) * 1e-10
            epsilon_pf = (2 ** len(coalition)) * max(float(record["initialization_residual"]) for record in vertices)
            uncertainty = epsilon_frequency + epsilon_eig + epsilon_pf
            row.update(
                {
                    "epsilon_pole_frequency": epsilon_frequency,
                    "epsilon_pole_eigensolution": epsilon_eig,
                    "epsilon_PF_Jacobian": epsilon_pf,
                    "pole_uncertainty": uncertainty,
                    "resolved_pole_externality": abs(row["delta_phi_poles"]) >= 5.0 * uncertainty,
                    "algebraic_fraction_absolute": abs(row["delta_phi_algebraic"])
                    / max(abs(row["delta_phi_algebraic"]) + abs(row["delta_phi_mass"]) + abs(row["delta_phi_poles"]), 1e-30),
                    "mass_fraction_absolute": abs(row["delta_phi_mass"])
                    / max(abs(row["delta_phi_algebraic"]) + abs(row["delta_phi_mass"]) + abs(row["delta_phi_poles"]), 1e-30),
                    "pole_fraction_absolute": abs(row["delta_phi_poles"])
                    / max(abs(row["delta_phi_algebraic"]) + abs(row["delta_phi_mass"]) + abs(row["delta_phi_poles"]), 1e-30),
                }
            )
        case_rows.append(row)
    return {
        "operating_point_id": op["operating_point_id"],
        "vertices": list(records.values()),
        "cases": case_rows,
        "audit": audit_rows,
        "spectra": np.stack(spectra),
        "wall_time_s": time.perf_counter() - started,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--population", choices=("development", "holdout"), required=True)
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    args = parser.parse_args()
    if not NUMERICAL.is_file() or not HOLDOUT_FREEZE.is_file():
        raise FileNotFoundError("E05B preregistration must exist before computation")
    if args.population == "holdout" and not CANDIDATE_FREEZE.is_file():
        raise RuntimeError("new E05B holdout is locked until candidate/contour freeze exists")
    numerical = load_json(NUMERICAL)
    if args.population == "development":
        selected_ops = operating_points()
    else:
        selected_ops = load_json(HOLDOUT_FREEZE)["operating_points"]
    tasks = [
        {"operating_point": op, "population": args.population, "numerical": numerical}
        for op in selected_ops
    ]
    started = time.perf_counter()
    results: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(worker, task): task["operating_point"]["operating_point_id"] for task in tasks}
        for count, future in enumerate(as_completed(futures), start=1):
            result = future.result()
            results.append(result)
            print(
                f"{args.population}: {count}/{len(tasks)} {result['operating_point_id']} "
                f"{result['wall_time_s']:.2f}s",
                flush=True,
            )
    results.sort(key=lambda value: value["operating_point_id"])
    vertices = pd.DataFrame([row for result in results for row in result["vertices"]])
    cases = pd.DataFrame([row for result in results for row in result["cases"]])
    audit = pd.DataFrame([row for result in results for row in result["audit"]])
    output = ARTIFACT / args.population
    vertices.to_parquet(output / "factor_vertices.parquet", index=False)
    cases.to_parquet(output / "component_externalities.parquet", index=False)
    if len(audit):
        audit.to_parquet(output / "pointwise_factorization_audit.parquet", index=False)
    spectra = np.stack([result["spectra"] for result in results])
    np.savez_compressed(
        ARTIFACT / "spectra" / f"{args.population}_spectra.npz",
        operating_point_ids=np.asarray([result["operating_point_id"] for result in results]),
        subset_keys=np.asarray([subset_key(subset) for subset in all_vertices()]),
        eigenvalues=spectra,
    )
    payload = {
        "population": args.population,
        "operating_points": len(results),
        "operator_builds": len(vertices),
        "pair_cases": int((cases["order"] == 2).sum()),
        "triple_cases": int((cases["order"] == 3).sum()),
        "successful_cases": int((cases["evaluation_status"] == "SUCCESS").sum()),
        "resolved_pole_cases": int(cases["resolved_pole_externality"].sum()),
        "pointwise_audit_rows": len(audit),
        "workers": args.workers,
        "wall_time_s": time.perf_counter() - started,
    }
    (ARTIFACT / "logs" / f"{args.population}_compute.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
