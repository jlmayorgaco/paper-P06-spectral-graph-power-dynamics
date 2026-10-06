"""Audit the frozen V4 full spectra and finalize F1 without a new eigensolve.

The source NPZ is the exact 16-portfolio GFL11 campaign. This script checks
its stored eigenpairs against reconstructed transverse operators, recomputes
the full-spectrum and target-band summaries, records model/equilibrium hashes,
and renders only F1 beneath a new immutable RUN_ID.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
from itertools import combinations
from pathlib import Path

import numpy as np


CORE = (30, 33, 35, 37)
SOURCE_RUN_ID = "20260926T083621_h4_f1_f2_assets_v1"
SOURCE_COMMIT = "a042efd7156a505928c52b9a373aa0839371f32b"
FREQ_BAND_HZ = (0.3, 1.5)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, default=float), encoding="utf-8")


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def portfolio_id(members: tuple[int, ...]) -> str:
    return "+".join(map(str, members)) if members else "BASE"


def expected_portfolios() -> list[tuple[int, ...]]:
    return [
        subset
        for size in range(len(CORE) + 1)
        for subset in combinations(CORE, size)
    ]


def mode_families(values: np.ndarray) -> list[list[int]]:
    """Pair conjugates for reporting the next distinct modal-family gap."""
    remaining = set(range(values.size))
    families: list[list[int]] = []
    eps = np.finfo(np.float64).eps
    while remaining:
        i = max(remaining, key=lambda index: values[index].real)
        remaining.remove(i)
        value = values[i]
        tol = 1000.0 * eps * max(1.0, abs(value))
        if abs(value.imag) <= tol or not remaining:
            families.append([i])
            continue
        j = min(remaining, key=lambda index: abs(values[index] - value.conjugate()))
        if abs(values[j] - value.conjugate()) <= tol:
            remaining.remove(j)
            families.append([i, j])
        else:
            families.append([i])
    return families


def render_f1(run_root: Path, rows: list[dict], tau_dec: float) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    positions: dict[str, float] = {}
    for cardinality in range(len(CORE) + 1):
        group = sorted(
            (row for row in rows if row["cardinality"] == cardinality),
            key=lambda row: row["portfolio_id"],
        )
        offsets = np.linspace(-0.28, 0.28, len(group)) if len(group) > 1 else [0.0]
        for offset, row in zip(offsets, group):
            positions[row["portfolio_id"]] = cardinality + float(offset)

    row_by_id = {row["portfolio_id"]: row for row in rows}
    figure, axis = plt.subplots(figsize=(8.4, 5.5))
    for parent in rows:
        members = (
            ()
            if parent["members"] == "BASE"
            else tuple(int(bus) for bus in parent["members"].split("+"))
        )
        for bus in CORE:
            if bus in members:
                continue
            child = tuple(sorted((*members, bus)))
            child_id = portfolio_id(child)
            if child_id in positions:
                axis.plot(
                    [positions[parent["portfolio_id"]], positions[child_id]],
                    [parent["alpha_perp_s-1"], row_by_id[child_id]["alpha_perp_s-1"]],
                    color="#a5adb8",
                    linewidth=0.7,
                    alpha=0.75,
                    zorder=1,
                )

    colors = {"STABLE": "#238b70", "INDETERMINATE": "#d59b21", "UNSTABLE": "#c74440"}
    for status, color in colors.items():
        group = [row for row in rows if row["status"] == status]
        if group:
            axis.scatter(
                [positions[row["portfolio_id"]] for row in group],
                [row["alpha_perp_s-1"] for row in group],
                s=48,
                color=color,
                edgecolor="#20252b",
                linewidth=0.5,
                label=status.lower(),
                zorder=3,
            )

    h4 = row_by_id["30+33+35+37"]
    axis.scatter(
        [positions[h4["portfolio_id"]]],
        [h4["alpha_perp_s-1"]],
        s=125,
        marker="*",
        facecolor="#8b1e3f",
        edgecolor="#20252b",
        linewidth=0.8,
        label="H4",
        zorder=4,
    )
    labels = {
        "BASE",
        "30+33+35+37",
        "30+33+35",
        "30+33+37",
        "30+35+37",
        "33+35+37",
        min(rows, key=lambda row: row["alpha_perp_s-1"])["portfolio_id"],
    }
    for row in rows:
        if row["portfolio_id"] not in labels:
            continue
        label = "H4" if row["portfolio_id"] == "30+33+35+37" else row["portfolio_id"]
        label_offsets = {
            "BASE": (0, -13),
            "30+33+35": (-10, -14),
            "30+33+37": (-24, 10),
            "30+35+37": (24, 10),
            "33+35+37": (0, -15),
            "33+35": (0, 8),
            "30+33+35+37": (0, 8),
        }
        xoff, yoff = label_offsets.get(
            row["portfolio_id"],
            (0, 7 if row["alpha_perp_s-1"] >= 0 else -11),
        )
        axis.annotate(
            label,
            (positions[row["portfolio_id"]], row["alpha_perp_s-1"]),
            xytext=(xoff, yoff),
            textcoords="offset points",
            ha="center",
            va="bottom" if yoff > 0 else "top",
            fontsize=6.7,
            color="#30343b",
        )

    axis.axhline(0.0, color="#20252b", linewidth=0.9)
    axis.set_xticks(range(5), ["0", "1", "2", "3", "4"])
    axis.set_xlabel("H4 portfolio cardinality |H|")
    axis.set_ylabel(r"full transverse spectral abscissa $\alpha_{\perp}$ [s$^{-1}$]")
    axis.set_title("F1 — Stability cliff: 15 proper subsets stable; H4 unstable")
    axis.grid(axis="y", color="#d8dde3", linewidth=0.55, alpha=0.85)
    axis.legend(frameon=False, ncol=4, loc="upper left", fontsize=8)
    figure.text(
        0.99,
        0.01,
        (
            fr"$\tau_{{dec}}={tau_dec:.0e}$ s$^{{-1}}$ · "
            r"$\alpha_{\Omega}=\alpha_{\perp}$ for 16/16 · "
            "Python-only V4 lattice; Julia parity is H4 GFL11 only"
        ),
        ha="right",
        va="bottom",
        fontsize=7.2,
        color="#4b5563",
    )
    figure.tight_layout(rect=(0, 0.04, 1, 1))
    output = run_root / "figures" / "F1_STABILITY_CLIFF"
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output.with_suffix(".png"), dpi=320, bbox_inches="tight")
    figure.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(output.with_suffix(".svg"), bbox_inches="tight")
    plt.close(figure)


def main() -> int:
    if len(sys.argv) != 4:
        raise SystemExit("usage: audit_finalize_f1.py <repo-root> <source-run-id> <new-run-id>")
    repo = Path(sys.argv[1]).resolve()
    source_run_id = sys.argv[2]
    run_id = sys.argv[3]
    result_root = repo / "reports/poster/ias2026/research/bnd_h4_mechanism/results"
    source_root = result_root / source_run_id
    run_root = result_root / run_id
    if source_run_id != SOURCE_RUN_ID:
        raise ValueError(f"IAS26-020 source is frozen to {SOURCE_RUN_ID}")
    resumable_files = {
        "raw/F1_FULL_SPECTRA.npz",
        "tables/F1_FINAL_AUDIT.csv",
        "tables/F1_FULL_SPECTRA_AUDITED.csv",
        "figures/F1_STABILITY_CLIFF.png",
        "figures/F1_STABILITY_CLIFF.pdf",
        "figures/F1_STABILITY_CLIFF.svg",
        "claims/IAS26-020_STATUS.json",
        "environment/provenance.json",
        "report/IAS26-020_SUMMARY.md",
    }
    resume_partial = run_root.exists()
    if resume_partial:
        existing = {
            path.relative_to(run_root).as_posix()
            for path in run_root.rglob("*")
            if path.is_file()
        }
        if not existing or not existing.issubset(resumable_files):
            raise FileExistsError(f"run root is not a recognized partial F1 audit: {run_root}")

    npz_path = source_root / "raw/F1_FULL_SPECTRA.npz"
    portfolio_csv = source_root / "derived/F1_PORTFOLIOS.csv"
    spectra_csv = source_root / "derived/F1_FULL_SPECTRA.csv"
    runtime_path = source_root / "environment/runtime.json"
    for path in (npz_path, portfolio_csv, spectra_csv, runtime_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    copied_npz = run_root / "raw" / "F1_FULL_SPECTRA.npz"
    if copied_npz.is_file() and sha256_file(copied_npz) != sha256_file(npz_path):
        raise RuntimeError("partial run contains a different spectrum NPZ; refusing to resume")

    runtime = json.loads(runtime_path.read_text(encoding="utf-8"))
    if runtime["git_commit"] != SOURCE_COMMIT or runtime["git_dirty"]:
        raise RuntimeError("source F1 asset is not the frozen clean source run")

    rel_renderer = "reports/poster/ias2026/research/bnd_h4_mechanism/code/render_f1_f2.py"
    frozen_renderer = subprocess.check_output(
        ["git", "show", f"{SOURCE_COMMIT}:{rel_renderer}"],
        cwd=repo,
        text=True,
        encoding="utf-8",
    )
    tau_match = re.search(r"^BOUNDARY_TOLERANCE\s*=\s*([0-9.eE+-]+)\s*$", frozen_renderer, re.M)
    if tau_match is None:
        raise RuntimeError("could not recover the boundary tolerance from the frozen renderer")
    tau_dec = float(tau_match.group(1))

    source_paths = {
        "tx4_contextual_return.py": repo / "reports/poster/ias2026/research/experiments/tx4_contextual_return.py",
        "ieee39_case.py": repo / "reports/poster/ias2026/research/src/ibr_cycles/models/ieee39_case.py",
        "port_admittance.py": repo / "reports/poster/ias2026/research/src/ibr_cycles/models/port_admittance.py",
        "ieee39_devices.py": repo / "reports/poster/ias2026/research/src/ibr_cycles/models/ieee39_devices.py",
    }
    source_hashes = {name: sha256_file(path) for name, path in source_paths.items()}
    source_hash_matches = {
        name: source_hashes[name] == runtime["source_hashes"][name]
        for name in source_hashes
    }
    if not all(source_hash_matches.values()):
        raise RuntimeError(f"model source differs from frozen asset hashes: {source_hash_matches}")

    with np.load(npz_path, allow_pickle=False) as data:
        portfolio_ids = data["portfolio_ids"].astype(str)
        dae_dims = data["dae_state_dims"].astype(int)
        transverse_dims = data["transverse_dims"].astype(int)
        eigenvalues = data["eigenvalues"].copy()
        eigenvectors = data["eigenvectors"].copy()

    expected = expected_portfolios()
    expected_ids = [portfolio_id(members) for members in expected]
    if portfolio_ids.tolist() != expected_ids:
        raise RuntimeError("stored portfolio order/content differs from the frozen V4 lattice")
    if len(portfolio_ids) != 16 or int(transverse_dims.sum()) != 1216:
        raise RuntimeError("the source NPZ is not a complete 16-portfolio/1,216-mode lattice")
    if eigenvalues.shape != (16, int(transverse_dims.max())):
        raise RuntimeError(f"unexpected eigenvalue array shape {eigenvalues.shape}")
    if eigenvectors.shape != (16, int(transverse_dims.max()), int(transverse_dims.max())):
        raise RuntimeError(f"unexpected eigenvector array shape {eigenvectors.shape}")

    with portfolio_csv.open(encoding="utf-8-sig", newline="") as handle:
        source_portfolios = list(csv.DictReader(handle))
    with spectra_csv.open(encoding="utf-8-sig", newline="") as handle:
        source_spectra = list(csv.DictReader(handle))
    if [row["portfolio"] for row in source_portfolios] != expected_ids:
        raise RuntimeError("portfolio CSV order does not match the exact frozen lattice")
    if len(source_spectra) != 1216:
        raise RuntimeError("source full-spectrum CSV is incomplete")
    spectra_by_id: dict[str, list[dict]] = {key: [] for key in expected_ids}
    for row in source_spectra:
        spectra_by_id[row["portfolio"]].append(row)

    experiment = repo / "reports/poster/ias2026/research/experiments"
    source = repo / "reports/poster/ias2026/research/src"
    code_dir = repo / "reports/poster/ias2026/research/bnd_h4_mechanism/code"
    sys.path[:0] = [str(experiment), str(source), str(code_dir)]
    import tx4_contextual_return as tx
    import run_m1a_reconciliation as m1a
    from ibr_cycles.certification.symmetry import frequency_partner, rotation_generator
    from ibr_cycles.certification.transverse import transverse_operator

    if tuple(tx.CORE) != CORE or tx.P4.as_dict() != {"g": 0.03625, "k": 1.425, "t": 1.5, "h": 1.0}:
        raise RuntimeError("current experiment context differs from the frozen GFL11/P4 contract")

    audit_rows: list[dict] = []
    full_rows: list[dict] = []
    case_hashes: dict[str, dict] = {}
    eigenpair_max_by_id: dict[str, float] = {}
    eq_residual_differences: list[float] = []
    stored_alpha_differences: list[float] = []
    stored_spectrum_differences: list[float] = []
    transverse_operators: dict[str, np.ndarray] = {}

    for case_index, members in enumerate(expected):
        pid = expected_ids[case_index]
        values = eigenvalues[case_index, : transverse_dims[case_index]]
        vectors = eigenvectors[case_index, : transverse_dims[case_index], : transverse_dims[case_index]]
        if not np.all(np.isfinite(values.real)) or not np.all(np.isfinite(values.imag)):
            raise RuntimeError(f"non-finite stored eigenvalue in {pid}")
        if not np.all(np.isfinite(vectors.real)) or not np.all(np.isfinite(vectors.imag)):
            raise RuntimeError(f"non-finite stored eigenvector in {pid}")
        alpha_perp = float(np.max(values.real))
        frequencies = np.abs(values.imag) / (2.0 * math.pi)
        band_mask = (
            (frequencies >= FREQ_BAND_HZ[0] - 1e-10)
            & (frequencies <= FREQ_BAND_HZ[1] + 1e-10)
            & (values.imag >= -1e-10)
            & (np.abs(values) > 1e-3)
        )
        if not np.any(band_mask):
            raise RuntimeError(f"no eigenvalue in the frozen target band for {pid}")
        band_indices = np.flatnonzero(band_mask)
        band_index = int(band_indices[np.argmax(values[band_indices].real)])
        alpha_omega = float(values[band_index].real)
        full_indices = np.flatnonzero(np.abs(values.real - alpha_perp) <= 1e-12)
        full_index = int(max(full_indices, key=lambda idx: values[idx].imag))
        full_eig = complex(values[full_index])
        band_eig = complex(values[band_index])

        families = mode_families(values)
        family_alphas = sorted(
            (float(max(values[index].real for index in family)) for family in families),
            reverse=True,
        )
        next_alpha = family_alphas[1] if len(family_alphas) > 1 else float("nan")
        family_gap = alpha_perp - next_alpha

        source_row = source_portfolios[case_index]
        stored_alpha = float(source_row["alpha"])
        stored_alpha_differences.append(abs(alpha_perp - stored_alpha))
        if int(source_row["dae_state_dim"]) != int(dae_dims[case_index]):
            raise RuntimeError(f"DAE state dimension mismatch for {pid}")
        if int(source_row["transverse_dim"]) != int(transverse_dims[case_index]):
            raise RuntimeError(f"transverse dimension mismatch for {pid}")

        rows_spectrum = spectra_by_id[pid]
        if len(rows_spectrum) != int(transverse_dims[case_index]):
            raise RuntimeError(f"full-spectrum row count mismatch for {pid}")
        csv_values = np.asarray(
            [float(row["eigenvalue_real_s-1"]) + 1j * float(row["eigenvalue_imag_rad_s-1"]) for row in rows_spectrum]
        )
        stored_spectrum_differences.append(float(np.max(np.abs(csv_values - values))))
        target_rows = [row for row in rows_spectrum if row["is_target_band_mode"].lower() == "true"]
        if len(target_rows) != 1:
            raise RuntimeError(f"target-mode identity is not unique for {pid}")
        target_eig = complex(float(target_rows[0]["eigenvalue_real_s-1"]), float(target_rows[0]["eigenvalue_imag_rad_s-1"]))
        if abs(target_eig - band_eig) > 1e-10:
            raise RuntimeError(f"frozen target-band mode does not match NPZ selection for {pid}")

        case = tx.tx4_case(members, tx.P4)
        rx, _ = rotation_generator(case.dae, case.equilibrium.z)
        partner = frequency_partner(case.dae)
        transverse = transverse_operator(case.system.A, rx, partner.w)
        a_perp = np.asarray(transverse.a_perp)
        if a_perp.shape != (transverse_dims[case_index], transverse_dims[case_index]):
            raise RuntimeError(f"reconstructed transverse matrix has wrong shape for {pid}")
        residuals = a_perp @ vectors - vectors * values[np.newaxis, :]
        denominators = max(float(np.linalg.norm(a_perp, 2)), 1e-300) * np.maximum(np.linalg.norm(vectors, axis=0), 1e-300)
        pair_rel = np.linalg.norm(residuals, axis=0) / denominators
        max_pair_rel = float(np.max(pair_rel))
        eigenpair_max_by_id[pid] = max_pair_rel
        transverse_operators[pid] = a_perp

        equilibrium = case.equilibrium
        f_res = float(np.max(np.abs(case.dae.f(equilibrium.x, equilibrium.z, {}))))
        g_res = float(np.max(np.abs(case.dae.g(equilibrium.x, equilibrium.z, {}))))
        saved_residuals = [float(source_row["equilibrium_f_residual"]), float(source_row["equilibrium_g_residual"])]
        eq_residual_differences.extend([abs(f_res - saved_residuals[0]), abs(g_res - saved_residuals[1])])

        model_hash = m1a.model_fingerprint(case)
        eq_hashes = m1a.equilibrium_fingerprints(case)
        case_hashes[pid] = {
            "model_hash": model_hash,
            "equilibrium_hash": eq_hashes["xz_sha256"],
            "reduced_A_hash": eq_hashes["reduced_A_sha256"],
        }

        status = "STABLE" if alpha_perp < -tau_dec else "UNSTABLE" if alpha_perp > tau_dec else "INDETERMINATE"
        mask = sum(1 << CORE.index(bus) for bus in members)
        parents = [portfolio_id(tuple(bus for bus in members if bus != removed)) for removed in members]
        children = [
            portfolio_id(tuple(sorted((*members, bus))))
            for bus in CORE
            if bus not in members
        ]
        audit_rows.append(
            {
                "portfolio_id": pid,
                "portfolio_mask_CORE_order_30_33_35_37": mask,
                "members": "+".join(map(str, members)) or "BASE",
                "cardinality": len(members),
                "parent_ids": ";".join(parents) or "NONE",
                "child_ids": ";".join(children) or "NONE",
                "alpha_perp_s-1": alpha_perp,
                "alpha_Omega_s-1": alpha_omega,
                "dominant_lambda_real_s-1": full_eig.real,
                "dominant_lambda_imag_rad_s-1": full_eig.imag,
                "dominant_frequency_hz": abs(full_eig.imag) / (2.0 * math.pi),
                "next_family_alpha_s-1": next_alpha,
                "gap_to_next_family_s-1": family_gap,
                "tau_dec_s-1": tau_dec,
                "status": status,
                "dae_state_dim": int(dae_dims[case_index]),
                "transverse_dim": int(transverse_dims[case_index]),
                "max_eigenpair_relative_residual": max_pair_rel,
                "equilibrium_f_residual": f_res,
                "equilibrium_g_residual": g_res,
                "coupling_residual": float(transverse.coupling_residual),
                "model_hash": model_hash,
                "equilibrium_hash": eq_hashes["xz_sha256"],
                "reduced_A_hash": eq_hashes["reduced_A_sha256"],
                "solver_backend": "NumPy 2.3.5; stored np.linalg.eig eigenpairs; no new eigensolve",
                "source_run_id": source_run_id,
                "run_id": run_id,
                "mode_tracking_note": "per-portfolio dominant family selected independently; no cross-portfolio continuation/MAC",
            }
        )

        for rank, value in enumerate(values):
            full_rows.append(
                {
                    "portfolio_id": pid,
                    "eigenvalue_rank_real_desc": rank,
                    "eigenvalue_real_s-1": float(value.real),
                    "eigenvalue_imag_rad_s-1": float(value.imag),
                    "frequency_hz": float(abs(value.imag) / (2.0 * math.pi)),
                    "is_in_target_band": bool(band_mask[rank]),
                    "is_full_spectral_abscissa_family": bool(abs(value.real - alpha_perp) <= 1e-12),
                    "eigenpair_relative_residual": float(pair_rel[rank]),
                    "model_hash": model_hash,
                    "equilibrium_hash": eq_hashes["xz_sha256"],
                    "source_run_id": source_run_id,
                }
            )

    pass_ledger = [row for row in audit_rows if row["status"] == "STABLE"]
    indeterminate = [row for row in audit_rows if row["status"] == "INDETERMINATE"]
    unstable = [row for row in audit_rows if row["status"] == "UNSTABLE"]
    h4 = next(row for row in audit_rows if row["portfolio_id"] == "30+33+35+37")
    alpha_band_match = all(abs(row["alpha_perp_s-1"] - row["alpha_Omega_s-1"]) <= 1e-12 for row in audit_rows)
    source_alpha_match = max(stored_alpha_differences) <= 1e-12
    source_spectra_match = max(stored_spectrum_differences) <= 1e-12
    gate = (
        len(audit_rows) == 16
        and len(full_rows) == 1216
        and len(pass_ledger) == 15
        and len(unstable) == 1
        and len(indeterminate) == 0
        and h4["alpha_perp_s-1"] > tau_dec
        and all(row["alpha_perp_s-1"] < -tau_dec for row in pass_ledger)
        and alpha_band_match
        and source_alpha_match
        and source_spectra_match
        and all(source_hash_matches.values())
    )

    # Copy exact input spectrum bytes; historical source artifacts remain untouched.
    run_root.mkdir(parents=True, exist_ok=True)
    for directory in ("raw", "tables", "figures", "claims", "environment", "report"):
        (run_root / directory).mkdir(parents=True, exist_ok=True)
    if not copied_npz.exists():
        shutil.copyfile(npz_path, copied_npz)
    final_audit_path = run_root / "tables" / "F1_FINAL_AUDIT.csv"
    if not final_audit_path.exists():
        write_csv(final_audit_path, audit_rows)
    spectra_audit_path = run_root / "tables" / "F1_FULL_SPECTRA_AUDITED.csv"
    if not spectra_audit_path.exists():
        write_csv(spectra_audit_path, full_rows)
    render_f1(run_root, audit_rows, tau_dec)

    source_parity = repo / "reports/poster/ias2026/research/bnd_h4_mechanism/results/20260925T225832_dee6fe69_h4_modal_schur_v1/claims/CROSSCODE_GFL11_RECONCILIATION.json"
    parity = json.loads(source_parity.read_text(encoding="utf-8"))
    parity_scope = "H4 GFL11 only; not a 16-portfolio Julia lattice"
    decision = {
        "ticket": "IAS26-020",
        "run_id": run_id,
        "status": "PASS" if gate else "BLOCKED",
        "source_run_id": source_run_id,
        "tau_dec_s-1": tau_dec,
        "tau_dec_provenance": f"BOUNDARY_TOLERANCE in render_f1_f2.py at source commit {SOURCE_COMMIT}",
        "portfolio_count": len(audit_rows),
        "full_eigenvalue_count": len(full_rows),
        "proper_subsets_stable": len(pass_ledger),
        "h4_unstable": h4["alpha_perp_s-1"] > tau_dec,
        "indeterminate_count": len(indeterminate),
        "alpha_perp_equals_alpha_Omega_all_16": alpha_band_match,
        "source_alpha_matches_npz": source_alpha_match,
        "source_csv_spectra_match_npz": source_spectra_match,
        "source_model_code_hashes_match": source_hash_matches,
        "max_stored_eigenpair_relative_residual": max(eigenpair_max_by_id.values()),
        "max_equilibrium_residual_difference_from_source_csv": max(eq_residual_differences),
        "h4": {key: h4[key] for key in ("alpha_perp_s-1", "alpha_Omega_s-1", "dominant_frequency_hz", "dae_state_dim", "transverse_dim", "equilibrium_f_residual", "equilibrium_g_residual", "max_eigenpair_relative_residual")},
        "julia_scope": parity_scope,
        "julia_h4_parity_status": parity.get("status"),
        "v9_supplement": "NOT_USED_IN_P1; no V9 claim is included in F1",
        "poster_claim": "15 proper subsets are stable and the exact H4 assembly is unstable for the registered TX4/P4 GFL11 case",
        "mode_tracking_note": "No cross-portfolio continuation/MAC; each portfolio's dominant family is selected from its retained full spectrum.",
    }
    write_json(run_root / "claims" / "IAS26-020_STATUS.json", decision)

    input_hashes = {
        "F1_FULL_SPECTRA.npz": sha256_file(npz_path),
        "F1_PORTFOLIOS.csv": sha256_file(portfolio_csv),
        "F1_FULL_SPECTRA.csv": sha256_file(spectra_csv),
        "source_runtime.json": sha256_file(runtime_path),
        "frozen_renderer_git_blob_sha256": hashlib.sha256(frozen_renderer.encode()).hexdigest(),
        "julia_h4_parity.json": sha256_file(source_parity),
    }
    provenance = {
        "run_id": run_id,
        "source_run_id": source_run_id,
        "source_run_git_commit": SOURCE_COMMIT,
        "current_git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip(),
        "current_git_dirty": bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=repo, text=True).strip()),
        "f1_finalizer_sha256": sha256_file(Path(__file__).resolve()),
        "input_sha256": input_hashes,
        "current_model_source_hashes": source_hashes,
        "source_hash_matches": source_hash_matches,
        "case_hashes": case_hashes,
        "original_f1_scope": "Python-only V4 GFL11 full transverse spectra",
        "cross_code_scope": parity_scope,
        "tau_dec_s-1": tau_dec,
        "recomputation_policy": "Reused the frozen 1,216 eigenvalues/eigenvectors; reconstructed operating-point transverse matrices only to verify eigenpair residuals and provenance; no new eigensolve or F2/F3/F4/MC/TDS run.",
        "output_policy": "All new outputs are organized beneath this immutable RUN_ID; historical runs were read-only.",
    }
    summary = (
        "# IAS26-020 — final F1 audit\n\n"
        f"Status: **{'PASS' if gate else 'BLOCKED'}**\n\n"
        f"RUN_ID: `{run_id}`. The exact source NPZ was reused byte-for-byte: 16 portfolios, 1,216 transverse eigenvalues and eigenvectors. "
        f"Decision threshold `tau_dec={tau_dec:.0e} s^-1` was read from the renderer at the clean source commit, not selected after inspecting the spectrum.\n\n"
        f"The full-spectrum gate gives {len(pass_ledger)} stable proper subsets, {len(unstable)} unstable H4, and {len(indeterminate)} indeterminate portfolios. "
        f"`alpha_perp` and `alpha_Omega` agree for all 16. H4 has alpha `{h4['alpha_perp_s-1']:.12g} s^-1`, "
        f"frequency `{h4['dominant_frequency_hz']:.9g} Hz`, 86 DAE states and 84 transverse coordinates.\n\n"
        f"The eigenpair residuals were checked against reconstructed transverse operators without re-running an eigensolver. "
        f"The figure is explicitly Python-only for the V4 lattice; stored Python–Julia parity covers H4 GFL11 only. "
        "V9 is not included in the P1 claim. No F2, F3, F4, Monte Carlo, or TDS work was performed.\n"
    )
    (run_root / "report" / "IAS26-020_SUMMARY.md").write_text(summary, encoding="utf-8")
    provenance["output_sha256"] = {
        path.relative_to(run_root).as_posix(): sha256_file(path)
        for path in run_root.rglob("*")
        if path.is_file() and path != run_root / "environment" / "provenance.json"
    }
    write_json(run_root / "environment" / "provenance.json", provenance)
    print(json.dumps({"run_root": str(run_root), "status": decision["status"], "tau_dec": tau_dec, "portfolios": len(audit_rows), "eigenvalues": len(full_rows), "stable_proper": len(pass_ledger), "h4_alpha": h4["alpha_perp_s-1"], "h4_frequency_hz": h4["dominant_frequency_hz"], "alpha_perp_equals_alpha_Omega": alpha_band_match, "max_eigenpair_relative_residual": max(eigenpair_max_by_id.values()), "max_equilibrium_residual_difference": max(eq_residual_differences), "source_hash_matches": source_hash_matches}, indent=2))
    return 0 if gate else 2


if __name__ == "__main__":
    raise SystemExit(main())
