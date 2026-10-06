"""Run the single preregistered IAS26-050 nominal H4 DAE verification.

The intervention/configuration must already be frozen. This entry point runs
only H4 at g=0.25, saves the full transverse spectrum/eigenvectors and makes
no parameter sweep, reduced-root calculation, Monte Carlo, or figure.
"""

from __future__ import annotations

import csv
import hashlib
import json
import platform
import shutil
import subprocess
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy


CONFIG_REL = Path("reports/poster/ias2026/research/bnd_h4_mechanism/configs/IAS26-050_INTERVENTION_V1.json")
RESULTS_REL = Path("reports/poster/ias2026/research/bnd_h4_mechanism/results")
BASELINE_RUN_ID = "20260925T204017_549c3c07_h4_crossmode_v3"
F1_RUN_ID = "20260926T161116Z_d0fecb32_ias26_020_f1_audit_v1"
EXPECTED_INPUTS = {
    "f1_full_spectra_npz": "db6581d59fa0d042839e9a327c6d682d677ea4d6dd2209bcb2abb0b1e68f34db",
    "f1_final_audit_csv": "4ea8042a80bb1e68bee1f272b0daf7baa4f23a9611fa6b5975cd3239ac11f718",
    "baseline_sweep_csv": "2dd4559acdd847349f68faa26818687eb3a6194f5a2387bab17205955525429b",
    "baseline_runtime_json": "e99613e9b5848e21442f91fa33abc341e6bfc11c8db6c6c384507997c2f597b6",
}
EXPECTED_SOURCES = {
    "reports/poster/ias2026/research/experiments/tx4_contextual_return.py": "77d4ed0f76f8c4a6ec369bc30b306361cff9b789c6f432e287d5f180d9210b74",
    "reports/poster/ias2026/research/experiments/_f7_common.py": "2e6ac2b3490c061d5ee3c344043dbfa6dd663e50bfbe4aadab304ad0a8dfd44f",
    "reports/poster/ias2026/research/src/ibr_cycles/models/ieee39_case.py": "1f38a8e541780ab7205809dc1ba58bdafb6c923933f5772d46407064edc8c722",
    "reports/poster/ias2026/research/src/ibr_cycles/models/ieee39_devices.py": "bcf0c76be65b171b580684065cde4eb88d0433c6b4d7beefec966f8711036d90",
    "reports/poster/ias2026/research/src/ibr_cycles/models/ieee39_network.py": "f42144cd02a6dc779ed690fd866ad026f33eb7844584c7f10281fc7d98e3fe13",
    "reports/poster/ias2026/research/src/ibr_cycles/models/port_admittance.py": "b0fa694dd0111be158133a9cf8c4e36be670e50ef562dcd63f89d7ff87b8d4ef",
    "reports/poster/ias2026/research/src/ibr_cycles/models/port_core.py": "38c54b1a70f87fbbe2c301a717e6dd6b4112b3a9cf8deb663426ce08b729bb9b",
    "reports/poster/ias2026/research/src/ibr_cycles/certification/transverse.py": "1773203ef87389a1e0e661dd955a543961c671b05977908dc9ae30abb246c4d4",
    "reports/poster/ias2026/research/src/ibr_cycles/certification/symmetry.py": "053a3b1a28bf526b4b543aef23bb87f028dffe9fcf9c9e2b5f91bd84a9cab012",
    "reports/poster/ias2026/research/bnd_h4_mechanism/code/run_m1a_reconciliation.py": "97f6213f5de22b356cb509239878444c0833326ece81ecf09f8194be774cc231",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def array_sha256(value: np.ndarray) -> str:
    arr = np.ascontiguousarray(np.asarray(value))
    digest = hashlib.sha256()
    digest.update((arr.dtype.str + str(arr.shape)).encode("utf-8"))
    digest.update(arr.tobytes())
    return digest.hexdigest()


def write_json(path: Path, data: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True, default=float) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise RuntimeError(f"refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def mode_families(values: np.ndarray) -> list[list[int]]:
    remaining = set(range(values.size))
    families: list[list[int]] = []
    eps = np.finfo(np.float64).eps
    while remaining:
        i = max(remaining, key=lambda idx: values[idx].real)
        remaining.remove(i)
        value = values[i]
        tol = 1000.0 * eps * max(1.0, abs(value))
        if abs(value.imag) <= tol or not remaining:
            families.append([i])
            continue
        j = min(remaining, key=lambda idx: abs(values[idx] - value.conjugate()))
        if abs(values[j] - value.conjugate()) <= tol:
            remaining.remove(j)
            families.append([i, j])
        else:
            families.append([i])
    return families


def classify(alpha: float, tau: float) -> str:
    return "STABLE" if alpha < -tau else "UNSTABLE" if alpha > tau else "INDETERMINATE"


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("usage: run_ias26_050_nominal.py <repo-root> <run-id>")
    repo = Path(sys.argv[1]).resolve()
    run_id = sys.argv[2]
    if not run_id or "/" in run_id or "\\" in run_id:
        raise ValueError("RUN_ID must be a single path component")
    config_path = repo / CONFIG_REL
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if config.get("status") != "FROZEN_BEFORE_NOMINAL_VERIFICATION":
        raise RuntimeError("IAS26-050 config is not frozen before execution")
    if config.get("treatment", {}).get("primary_intervention", {}).get("g") != 0.25:
        raise RuntimeError("frozen primary treatment is not exactly g=0.25")
    if config.get("treatment", {}).get("control", {}).get("g") != 0.03625:
        raise RuntimeError("frozen control is not exactly g=0.03625")
    if config.get("system", {}).get("portfolio") != "30+33+35+37" or config.get("system", {}).get("port_order") != [30, 33, 35, 37]:
        raise RuntimeError("frozen H4 membership or port ordering changed")
    if config.get("system", {}).get("technology") != "all-GFL":
        raise RuntimeError("frozen system technology is not all-GFL")
    if config.get("system", {}).get("unchanged_parameters") != {
        "k": 1.425,
        "t": 1.5,
        "h": 1.0,
        "leak": 0.05,
        "network_loads_dispatch_locations_ratings": "identical to frozen P4 IEEE-39 TX4 source",
    }:
        raise RuntimeError("frozen non-g treatment parameters changed")
    if config.get("metrics_and_decisions", {}).get("tau_dec_s-1") != 1e-8:
        raise RuntimeError("frozen tau_dec changed")

    results = repo / RESULTS_REL
    run_root = results / run_id
    if run_root.exists():
        raise FileExistsError(f"immutable RUN_ID already exists: {run_root}")
    baseline = results / BASELINE_RUN_ID
    f1_run = results / F1_RUN_ID
    input_paths = {
        "f1_full_spectra_npz": f1_run / "raw/F1_FULL_SPECTRA.npz",
        "f1_final_audit_csv": f1_run / "tables/F1_FINAL_AUDIT.csv",
        "baseline_sweep_csv": baseline / "derived/TX4_CONTEXTUAL_RETURN_SWEEP.csv",
        "baseline_runtime_json": baseline / "environment/runtime.json",
    }
    input_hashes = {key: sha256_file(path) for key, path in input_paths.items()}
    if input_hashes != EXPECTED_INPUTS:
        raise RuntimeError("one or more frozen IAS26-050 input hashes changed")
    expected_config_hashes = config["existing_evidence"]["sha256"]
    for key, value in EXPECTED_INPUTS.items():
        config_key = {
            "f1_full_spectra_npz": "f1_full_spectra_npz",
            "f1_final_audit_csv": "f1_final_audit_csv",
            "baseline_sweep_csv": "baseline_sweep_csv",
        }.get(key)
        if config_key and expected_config_hashes[config_key] != value:
            raise RuntimeError(f"config source hash mismatch for {key}")

    source_hashes = {rel: sha256_file(repo / rel) for rel in EXPECTED_SOURCES}
    for rel, expected in EXPECTED_SOURCES.items():
        if source_hashes[rel] != expected:
            raise RuntimeError(f"frozen model source changed: {rel}")
    config_source_hashes = config["existing_evidence"]["sha256"]
    config_source_keys = {
        "reports/poster/ias2026/research/experiments/tx4_contextual_return.py": "tx4_contextual_return_source",
        "reports/poster/ias2026/research/experiments/_f7_common.py": "shared_parameter_source",
        "reports/poster/ias2026/research/src/ibr_cycles/models/ieee39_case.py": "ieee39_case_source",
        "reports/poster/ias2026/research/src/ibr_cycles/models/ieee39_devices.py": "ieee39_devices_source",
        "reports/poster/ias2026/research/src/ibr_cycles/models/ieee39_network.py": "ieee39_network_source",
        "reports/poster/ias2026/research/src/ibr_cycles/models/port_admittance.py": "port_admittance_source",
        "reports/poster/ias2026/research/src/ibr_cycles/models/port_core.py": "port_core_source",
        "reports/poster/ias2026/research/src/ibr_cycles/certification/transverse.py": "transverse_operator_source",
        "reports/poster/ias2026/research/src/ibr_cycles/certification/symmetry.py": "symmetry_source",
        "reports/poster/ias2026/research/bnd_h4_mechanism/code/run_m1a_reconciliation.py": "fingerprint_helper_source",
    }
    for rel, config_key in config_source_keys.items():
        if config_source_hashes[config_key] != source_hashes[rel]:
            raise RuntimeError(f"config source hash mismatch for {rel}")

    # Start the immutable run only after the treatment config and frozen inputs validate.
    for folder in ("config", "inputs", "raw", "tables", "claims", "report", "environment", "logs", "code"):
        (run_root / folder).mkdir(parents=True, exist_ok=True)
    shutil.copy2(config_path, run_root / "config" / config_path.name)
    shutil.copy2(Path(__file__).resolve(), run_root / "code" / "run_ias26_050_nominal.py")
    shutil.copy2(input_paths["f1_full_spectra_npz"], run_root / "inputs" / "F1_FULL_SPECTRA_H4_G0.npz")
    shutil.copy2(input_paths["f1_final_audit_csv"], run_root / "inputs" / "F1_FINAL_AUDIT.csv")
    shutil.copy2(input_paths["baseline_sweep_csv"], run_root / "inputs" / "TX4_CONTEXTUAL_RETURN_SWEEP.csv")
    shutil.copy2(input_paths["baseline_runtime_json"], run_root / "inputs" / "baseline_runtime.json")
    baseline_spectrum_path = input_paths["f1_full_spectra_npz"]
    npz = np.load(baseline_spectrum_path, allow_pickle=False)
    ids = npz["portfolio_ids"].astype(str).tolist()
    h4_index = ids.index("30+33+35+37")
    h4_dim = int(npz["transverse_dims"][h4_index])
    baseline_values = np.asarray(npz["eigenvalues"][h4_index, :h4_dim], dtype=complex)
    baseline_vectors = np.asarray(npz["eigenvectors"][h4_index, :h4_dim, :h4_dim], dtype=complex)
    audit_rows = read_csv(input_paths["f1_final_audit_csv"])
    h4_audit = next(row for row in audit_rows if row["portfolio_id"] == "30+33+35+37")
    expected_dims = (int(config["existing_evidence"]["f1_H4_state_count"]), int(config["existing_evidence"]["f1_H4_transverse_state_count"]))
    if (int(h4_audit["dae_state_dim"]), int(h4_audit["transverse_dim"])) != expected_dims:
        raise RuntimeError("frozen F1 H4 state dimensions disagree with the IAS26-050 config")
    if h4_audit["equilibrium_hash"] != config["existing_evidence"]["f1_H4_equilibrium_hash"]:
        raise RuntimeError("frozen F1 H4 equilibrium hash disagrees with the IAS26-050 config")
    if h4_audit["model_hash"] != config["existing_evidence"]["f1_H4_model_hash"]:
        raise RuntimeError("frozen F1 H4 model hash disagrees with the IAS26-050 config")
    if abs(float(np.max(baseline_values.real)) - float(h4_audit["alpha_perp_s-1"])) > 1e-12:
        raise RuntimeError("F1 H4 complete spectrum does not reproduce its audited alpha_perp")
    sweep_rows = read_csv(input_paths["baseline_sweep_csv"])
    nominal_row = next(row for row in sweep_rows if float(row["g"]) == 0.25)
    baseline_g0 = next(row for row in sweep_rows if float(row["g"]) == 0.03625)
    tau = float(config["metrics_and_decisions"]["tau_dec_s-1"])
    f_lo, f_hi = config["metrics_and_decisions"]["target_band_hz"]

    try:
        experiment = repo / "reports/poster/ias2026/research/experiments"
        source = repo / "reports/poster/ias2026/research/src"
        code = repo / "reports/poster/ias2026/research/bnd_h4_mechanism/code"
        sys.path[:0] = [str(experiment), str(source), str(code)]
        import tx4_contextual_return as tx  # type: ignore
        import run_m1a_reconciliation as m1a  # type: ignore
        from ibr_cycles.certification.symmetry import frequency_partner, rotation_generator
        from ibr_cycles.certification.transverse import transverse_operator

        if tuple(tx.CORE) != (30, 33, 35, 37):
            raise RuntimeError("H4 port order differs from the frozen configuration")
        theta = tx.Theta(g=0.25, k=1.425, t=1.5, h=1.0)
        case = tx.tx4_case(tx.CORE, theta)
        rx, _ = rotation_generator(case.dae, case.equilibrium.z)
        partner = frequency_partner(case.dae)
        transverse = transverse_operator(case.system.A, rx, partner.w)
        a_perp = np.asarray(transverse.a_perp)
        values, vectors = np.linalg.eig(a_perp)
        order = np.argsort(values.real)[::-1]
        values = values[order]
        vectors = vectors[:, order]
        finite = np.isfinite(values.real) & np.isfinite(values.imag)
        if not np.all(finite) or not np.all(np.isfinite(vectors.real)) or not np.all(np.isfinite(vectors.imag)):
            raise RuntimeError("candidate direct-DAE spectrum/eigenvectors contain non-finite values")
        alpha_perp = float(np.max(values.real))
        frequencies = np.abs(values.imag) / (2.0 * np.pi)
        band_mask = (frequencies >= f_lo - 1e-10) & (frequencies <= f_hi + 1e-10) & (values.imag >= -1e-10) & (np.abs(values) > 1e-3)
        if not np.any(band_mask):
            raise RuntimeError("candidate has no admissible transverse eigenvalue in the frozen target band")
        band_indices = np.flatnonzero(band_mask)
        band_index = int(band_indices[np.argmax(values[band_indices].real)])
        alpha_omega = float(values[band_index].real)
        eig_band = complex(values[band_index])
        full_indices = np.flatnonzero(np.abs(values.real - alpha_perp) <= 1e-12)
        full_index = int(max(full_indices, key=lambda idx: values[idx].imag))
        eig_full = complex(values[full_index])
        families = mode_families(values)
        family_alphas = sorted((float(max(values[index].real for index in family)) for family in families), reverse=True)
        next_family_alpha = family_alphas[1] if len(family_alphas) > 1 else float("nan")
        family_gap = alpha_perp - next_family_alpha
        pair_residuals = a_perp @ vectors - vectors * values[np.newaxis, :]
        denominators = max(float(np.linalg.norm(a_perp, 2)), 1e-300) * np.maximum(np.linalg.norm(vectors, axis=0), 1e-300)
        relative_residuals = np.linalg.norm(pair_residuals, axis=0) / denominators
        f_residual = float(np.max(np.abs(case.dae.f(case.equilibrium.x, case.equilibrium.z, {}))))
        g_residual = float(np.max(np.abs(case.dae.g(case.equilibrium.x, case.equilibrium.z, {}))))
        fingerprints = m1a.equilibrium_fingerprints(case)
        model_hash = m1a.model_fingerprint(case)
        eq_hash = fingerprints["xz_sha256"]
        tx_source_values = {
            "alpha_perp": float(np.max(baseline_values.real)),
            "alpha_Omega": float(h4_audit["alpha_Omega_s-1"]),
            "lambda_perp_real": float(h4_audit["dominant_lambda_real_s-1"]),
            "lambda_perp_imag": float(h4_audit["dominant_lambda_imag_rad_s-1"]),
            "frequency_perp_hz": float(h4_audit["dominant_frequency_hz"]),
            "alpha_Omega_source_g025": float(nominal_row["alpha"]),
            "frequency_Omega_source_g025_hz": float(nominal_row["frequency_hz"]),
        }
        delta_alpha_perp = alpha_perp - tx_source_values["alpha_perp"]
        delta_alpha_omega = alpha_omega - tx_source_values["alpha_Omega"]
        observed_nominal_delta_alpha_omega = delta_alpha_omega
        existing_sweep_reproduction_difference = alpha_omega - float(nominal_row["alpha"])
        preregistered_delta_alpha_omega = float(config["metrics_and_decisions"]["nominal_pre_run_prediction"]["delta_alpha_intervention_minus_original_s-1"])
        prediction_error = observed_nominal_delta_alpha_omega - preregistered_delta_alpha_omega

        observed_dims = (int(case.dae.n_x), int(a_perp.shape[0]))
        finite_and_residual_gate = (
            observed_dims == expected_dims
            and f_residual <= 1e-7
            and g_residual <= 1e-7
            and np.all(np.isfinite(values.real))
            and np.all(np.isfinite(values.imag))
            and np.all(np.isfinite(relative_residuals))
        )
        target_band_gate = classify(float(baseline_g0["alpha"]), tau) == "UNSTABLE" and classify(alpha_omega, tau) == "STABLE"
        full_spectrum_gate = classify(alpha_perp, tau) == "STABLE"
        overall_pass = bool(finite_and_residual_gate and target_band_gate and full_spectrum_gate)

        np.savez_compressed(
            run_root / "raw" / "H4_G025_FULL_TRANSVERSE_SPECTRUM.npz",
            eigenvalues=values,
            eigenvectors=vectors,
            eigenpair_relative_residuals=relative_residuals,
            a_perp=a_perp,
            state_dims=np.asarray(observed_dims, dtype=int),
            equilibrium_x=np.asarray(case.equilibrium.x),
            equilibrium_z=np.asarray(case.equilibrium.z),
            system_A=np.asarray(case.system.A),
        )
        write_csv(
            run_root / "tables" / "H4_G025_FULL_SPECTRUM.csv",
            [
                {
                    "eigenvalue_rank_real_desc": rank,
                    "eigenvalue_real_s-1": float(value.real),
                    "eigenvalue_imag_rad_s-1": float(value.imag),
                    "frequency_hz": float(abs(value.imag) / (2.0 * np.pi)),
                    "is_target_band_mode": bool(band_mask[rank]),
                    "is_full_alpha_perp_family": bool(abs(value.real - alpha_perp) <= 1e-12),
                    "eigenpair_relative_residual": float(relative_residuals[rank]),
                }
                for rank, value in enumerate(values)
            ],
        )
        write_csv(
            run_root / "tables" / "IAS26-050_NOMINAL_COMPARISON.csv",
            [
                {
                    "condition": "original_g_0.03625",
                    "g": 0.03625,
                    "alpha_perp_s-1": tx_source_values["alpha_perp"],
                    "alpha_Omega_s-1": tx_source_values["alpha_Omega"],
                    "target_lambda_real_s-1": tx_source_values["alpha_Omega"],
                    "target_lambda_imag_rad_s-1": tx_source_values["lambda_perp_imag"],
                    "target_frequency_hz": tx_source_values["frequency_perp_hz"],
                    "state_count": expected_dims[0],
                    "transverse_dimension": expected_dims[1],
                    "source": "IAS26-020 full H4 spectra at frozen P4",
                },
                {
                    "condition": "fixed_intervention_g_0.25",
                    "g": 0.25,
                    "alpha_perp_s-1": alpha_perp,
                    "alpha_Omega_s-1": alpha_omega,
                    "target_lambda_real_s-1": eig_band.real,
                    "target_lambda_imag_rad_s-1": eig_band.imag,
                    "target_frequency_hz": abs(eig_band.imag) / (2.0 * np.pi),
                    "state_count": observed_dims[0],
                    "transverse_dimension": observed_dims[1],
                    "source": "IAS26-050 one-point direct DAE verification",
                },
            ],
        )
        mode_row = {
            "alpha_perp_g025_s-1": alpha_perp,
            "alpha_Omega_g025_s-1": alpha_omega,
            "lambda_perp_g025_real_s-1": eig_full.real,
            "lambda_perp_g025_imag_rad_s-1": eig_full.imag,
            "frequency_perp_g025_hz": abs(eig_full.imag) / (2.0 * np.pi),
            "lambda_Omega_g025_real_s-1": eig_band.real,
            "lambda_Omega_g025_imag_rad_s-1": eig_band.imag,
            "frequency_Omega_g025_hz": abs(eig_band.imag) / (2.0 * np.pi),
            "next_family_alpha_g025_s-1": next_family_alpha,
            "gap_to_next_family_g025_s-1": family_gap,
            "family_id": "not_assigned; no cross-g eigenvector MAC performed",
        }
        write_json(run_root / "claims" / "IAS26-050_STATUS.json", {
            "ticket": "IAS26-050",
            "status": "PASS" if overall_pass else "NOMINAL_INTERVENTION_NOT_VALIDATED",
            "gates": {
                "frozen_input_and_source_hashes": "PASS",
                "equilibrium_finite_and_residuals": "PASS" if finite_and_residual_gate else "FAIL",
                "target_band_sign_change": "PASS" if target_band_gate else "FAIL",
                "complete_transverse_alpha_perp_stable_at_g025": "PASS" if full_spectrum_gate else "FAIL_OR_INDETERMINATE",
            },
            "tau_dec_s-1": tau,
            "nominal_metrics": {
                **mode_row,
                "alpha_perp_g003625_s-1": tx_source_values["alpha_perp"],
                "alpha_Omega_g003625_s-1": tx_source_values["alpha_Omega"],
                "delta_alpha_perp_g025_minus_g003625_s-1": delta_alpha_perp,
                "delta_alpha_Omega_g025_minus_g003625_s-1": delta_alpha_omega,
                "pre_run_predicted_delta_alpha_Omega_s-1": preregistered_delta_alpha_omega,
                "pre_run_prediction_error_s-1": prediction_error,
                "existing_sweep_g025_reproduction_difference_s-1": existing_sweep_reproduction_difference,
                "equilibrium_f_residual": f_residual,
                "equilibrium_g_residual": g_residual,
                "max_eigenpair_relative_residual": float(np.max(relative_residuals)),
                "coupling_residual": float(transverse.coupling_residual),
                "equilibrium_hash": eq_hash,
                "model_hash": model_hash,
                "state_count": observed_dims[0],
                "transverse_dimension": observed_dims[1],
            },
            "claims_safe": ["The fixed g=0.25 intervention moves the direct-DAE target-band eigenvalue left at the canonical nominal point."] if target_band_gate else [],
            "claims_unsafe": [
                "Do not claim universal rescue or operational robustness before IAS26-060.",
                "Do not claim the same mode family across g without eigenvector-based tracking.",
                "Do not compare reduced feedback poles with DAE poles; IAS26-010 remains BLOCKED_M1_STRICT.",
                "Do not claim positive local damping or optimality/cost superiority.",
            ],
            "mc_execution_authorized": False,
        })

        shutil.copy2(
            input_paths["f1_full_spectra_npz"],
            run_root / "raw" / "H4_G003625_FROZEN_FULL_TRANSVERSE_SPECTRUM.npz",
        )
        write_json(run_root / "environment" / "case_fingerprints.json", {
            "candidate_g025": {"model_hash": model_hash, "equilibrium_hash": eq_hash, "reduced_A_hash": fingerprints["reduced_A_sha256"]},
            "original_g003625_from_IAS26-020": {"model_hash": h4_audit["model_hash"], "equilibrium_hash": h4_audit["equilibrium_hash"], "reduced_A_hash": h4_audit["reduced_A_hash"]},
            "same_source_model_code": all(
                source_hashes[rel] == EXPECTED_SOURCES[rel]
                for rel in EXPECTED_SOURCES
                if rel != "reports/poster/ias2026/research/bnd_h4_mechanism/code/run_m1a_reconciliation.py"
            ),
            "parameterized_model_hash_matches_original": model_hash == h4_audit["model_hash"],
            "equilibrium_hash_changes_with_g": eq_hash != h4_audit["equilibrium_hash"],
        })
        (run_root / "report" / "IAS26-050_SUMMARY.md").write_text(
            "# IAS26-050 — fixed intervention nominal DAE check\n\n"
            f"Status: **{'PASS' if overall_pass else 'NOMINAL_INTERVENTION_NOT_VALIDATED'}**\n\n"
            f"Treatment was frozen before this one-point run in `{CONFIG_REL.as_posix()}`: H4 all-GFL, `g=0.25`; control `g=0.03625`. No sweep, MC, reduced-root calculation, eta, or figure was run.\n\n"
            f"- Complete-spectrum `alpha_perp`: `{tx_source_values['alpha_perp']:.12g}` → `{alpha_perp:.12g}` s^-1 (`Delta={delta_alpha_perp:.12g}`).\n"
            f"- Target-band `alpha_Omega`: `{tx_source_values['alpha_Omega']:.12g}` → `{alpha_omega:.12g}` s^-1 (`Delta={delta_alpha_omega:.12g}`).\n"
            f"- Frozen target-band Delta-alpha prediction: `{preregistered_delta_alpha_omega:.12g}` s^-1; direct-DAE comparison error: `{prediction_error:.6g}` s^-1.\n"
            f"- Candidate equilibrium residuals: f `{f_residual:.6g}`, g `{g_residual:.6g}`; state/transverse dimensions `{observed_dims[0]}/{observed_dims[1]}`.\n"
            f"- Candidate `alpha_perp` verdict at frozen tau_dec={tau:.0e}: **{classify(alpha_perp, tau)}**.\n\n"
            "Full spectra and eigenvectors for the candidate and the frozen original H4 point are retained. The direct-DAE result is not a validation of reduced poles. IAS26-060 remains prohibited until its operational distribution, correlations, limits, redispatch, seed, and 1,000 scenario IDs are frozen in its own pre-run config.\n",
            encoding="utf-8",
        )
    except Exception:
        (run_root / "logs" / "nominal_run_exception.txt").write_text(traceback.format_exc(), encoding="utf-8")
        write_json(run_root / "claims" / "IAS26-050_STATUS.json", {
            "ticket": "IAS26-050",
            "status": "NOMINAL_RUN_FAILED",
            "mc_execution_authorized": False,
            "config_sha256": sha256_file(config_path),
            "input_sha256": input_hashes,
        })
        raise

    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip()
    dirty = bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=repo, text=True).strip())
    provenance = {
        "run_id": run_id,
        "ticket": "IAS26-050",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "config_path": CONFIG_REL.as_posix(),
        "config_sha256": sha256_file(config_path),
        "current_git_commit": commit,
        "current_worktree_dirty": dirty,
        "baseline_commit": "549c3c07383c1ecc6e6a274cc7ca900f97bc47d3",
        "input_sha256": input_hashes,
        "source_sha256": source_hashes,
        "execution_script_sha256": sha256_file(Path(__file__).resolve()),
        "python": {"version": sys.version, "executable": sys.executable, "platform": platform.platform()},
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "run_scope": "one direct DAE nominal H4 case at g=0.25 only; no MC/sweep/reduced model/figure",
        "output_policy": "all outputs beneath immutable RUN_ID; no standalone .zip packages or root-level outputs; required numerical arrays use .npz within RUN_ID",
    }
    write_json(run_root / "environment" / "provenance.json", provenance)
    out_hashes = {}
    for path in sorted(item for item in run_root.rglob("*") if item.is_file()):
        out_hashes[path.relative_to(run_root).as_posix()] = sha256_file(path)
    write_json(run_root / "environment" / "OUTPUT_SHA256.json", out_hashes)
    print(json.dumps({
        "run_id": run_id,
        "status": "PASS" if overall_pass else "NOMINAL_INTERVENTION_NOT_VALIDATED",
        "alpha_perp_g025": alpha_perp,
        "alpha_Omega_g025": alpha_omega,
        "delta_alpha_perp": delta_alpha_perp,
        "delta_alpha_Omega": delta_alpha_omega,
        "output": str(run_root),
    }, indent=2))
    return 0 if overall_pass else 2


if __name__ == "__main__":
    raise SystemExit(main())
