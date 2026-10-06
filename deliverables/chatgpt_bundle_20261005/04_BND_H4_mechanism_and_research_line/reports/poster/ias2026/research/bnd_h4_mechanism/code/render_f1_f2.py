"""Regenerate the poster-ready F1 and F2 assets from frozen GFL11 data.

This is an asset-generation run, not a new mechanism experiment. F1 solves
the exact 16-case inclusion lattice and stores every transverse eigenvalue and
eigenvector. F2 consumes only the previously validated physical local factors
I+M_ii and collective I+Q_H factors from the baseline campaign.
"""

from __future__ import annotations

import csv
import hashlib
import json
import platform
import subprocess
import sys
from itertools import combinations
from pathlib import Path
from shutil import copy2

import numpy as np


CORE = (30, 33, 35, 37)
BOUNDARY_TOLERANCE = 1e-8


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git(repo: Path, *args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=repo, text=True).strip()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, default=float), encoding="utf-8")


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def portfolio_id(members: tuple[int, ...]) -> str:
    return "+".join(str(bus) for bus in members) if members else "BASE"


def all_portfolios() -> list[tuple[int, ...]]:
    return [tuple(subset) for size in range(len(CORE) + 1) for subset in combinations(CORE, size)]


def save_figure(figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path.with_suffix(".png"), dpi=320, bbox_inches="tight")
    figure.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(path.with_suffix(".svg"), bbox_inches="tight")


def make_f1(repo: Path, run_root: Path) -> dict:
    experiment = repo / "reports/poster/ias2026/research/experiments"
    source = repo / "reports/poster/ias2026/research/src"
    sys.path.insert(0, str(experiment))
    sys.path.insert(0, str(source))
    import tx4_contextual_return as tx  # type: ignore  # noqa: PLC0415

    tx.OUT = run_root / "derived"
    tx.OUT.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    spectra: list[np.ndarray] = []
    vectors: list[np.ndarray] = []
    portfolio_rows: list[dict] = []

    for members in all_portfolios():
        case = tx.tx4_case(members, tx.P4)
        mode = tx.mode_of(case)
        values, eigenvectors = np.linalg.eig(mode.a_perp)
        order = np.argsort(values.real)[::-1]
        values = values[order]
        eigenvectors = eigenvectors[:, order]
        spectra.append(values)
        vectors.append(eigenvectors)
        critical_index = int(np.argmin(abs(values - mode.eig)))
        for rank, value in enumerate(values):
            rows.append(
                {
                    "portfolio": portfolio_id(members),
                    "members": "+".join(map(str, members)) or "BASE",
                    "cardinality": len(members),
                    "dae_state_dim": int(case.dae.n_x),
                    "transverse_dim": int(values.size),
                    "eigenvalue_rank_real_desc": rank,
                    "eigenvalue_real_s-1": float(value.real),
                    "eigenvalue_imag_rad_s-1": float(value.imag),
                    "frequency_hz": float(abs(value.imag) / (2.0 * np.pi)),
                    "is_target_band_mode": bool(rank == critical_index),
                }
            )
        status = (
            "BOUNDARY"
            if abs(mode.alpha) <= BOUNDARY_TOLERANCE
            else "UNSTABLE"
            if mode.alpha > 0.0
            else "STABLE"
        )
        portfolio_rows.append(
            {
                "portfolio": portfolio_id(members),
                "members": "+".join(map(str, members)) or "BASE",
                "cardinality": len(members),
                "dae_state_dim": int(case.dae.n_x),
                "transverse_dim": int(values.size),
                "alpha": mode.alpha,
                "dominant_lambda_real_s-1": mode.eig.real,
                "dominant_lambda_imag_rad_s-1": mode.eig.imag,
                "dominant_frequency_hz": mode.frequency_hz,
                "damping_ratio": mode.damping_ratio,
                "status": status,
                "equilibrium_f_residual": float(np.max(np.abs(case.dae.f(case.equilibrium.x, case.equilibrium.z, {})))),
                "equilibrium_g_residual": float(np.max(np.abs(case.dae.g(case.equilibrium.x, case.equilibrium.z, {})))),
                "coupling_residual": mode.coupling_residual,
            }
        )

    max_dim = max(values.size for values in spectra)
    eigenvalue_array = np.full((len(spectra), max_dim), np.nan + 1j * np.nan, dtype=complex)
    eigenvector_array = np.zeros((len(vectors), max_dim, max_dim), dtype=complex)
    state_dims = []
    for i, (values, basis) in enumerate(zip(spectra, vectors)):
        dim = values.size
        state_dims.append(dim)
        eigenvalue_array[i, :dim] = values
        eigenvector_array[i, :dim, :dim] = basis
    np.savez_compressed(
        run_root / "raw" / "F1_FULL_SPECTRA.npz",
        portfolio_ids=np.asarray([row["portfolio"] for row in portfolio_rows]),
        dae_state_dims=np.asarray([row["dae_state_dim"] for row in portfolio_rows], dtype=int),
        transverse_dims=np.asarray(state_dims, dtype=int),
        eigenvalues=eigenvalue_array,
        eigenvectors=eigenvector_array,
    )
    write_csv(run_root / "derived" / "F1_PORTFOLIOS.csv", portfolio_rows)
    write_csv(run_root / "derived" / "F1_FULL_SPECTRA.csv", rows)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    positions: dict[str, float] = {}
    for cardinality in range(5):
        group = [row for row in portfolio_rows if row["cardinality"] == cardinality]
        offsets = np.linspace(-0.28, 0.28, len(group)) if len(group) > 1 else np.asarray([0.0])
        for offset, row in zip(offsets, sorted(group, key=lambda item: item["portfolio"])):
            positions[row["portfolio"]] = cardinality + float(offset)

    figure, axis = plt.subplots(figsize=(7.7, 5.2))
    row_by_id = {row["portfolio"]: row for row in portfolio_rows}
    for parent in portfolio_rows:
        parent_members = () if parent["members"] == "BASE" else tuple(int(bus) for bus in parent["members"].split("+") if bus)
        for bus in CORE:
            if bus not in parent_members:
                child_members = tuple(sorted((*parent_members, bus)))
                child_id = portfolio_id(child_members)
                if child_id in positions:
                    axis.plot(
                        [positions[parent["portfolio"]], positions[child_id]],
                        [parent["alpha"], row_by_id[child_id]["alpha"]],
                        color="#a5adb8",
                        linewidth=0.65,
                        alpha=0.75,
                        zorder=1,
                    )

    stable = [row for row in portfolio_rows if row["status"] == "STABLE"]
    unstable = [row for row in portfolio_rows if row["status"] == "UNSTABLE"]
    boundary = [row for row in portfolio_rows if row["status"] == "BOUNDARY"]
    for group, color, label in (
        (stable, "#238b70", "stable"),
        (boundary, "#d59b21", "boundary"),
        (unstable, "#c74440", "unstable"),
    ):
        if group:
            axis.scatter(
                [positions[row["portfolio"]] for row in group],
                [row["alpha"] for row in group],
                s=42,
                color=color,
                edgecolor="#20252b",
                linewidth=0.45,
                label=label,
                zorder=3,
            )
    h4 = row_by_id["30+33+35+37"]
    axis.scatter(
        [positions["30+33+35+37"]],
        [h4["alpha"]],
        s=100,
        facecolor="#8b1e3f",
        edgecolor="#20252b",
        linewidth=0.8,
        marker="*",
        label="H4",
        zorder=4,
    )
    labels_to_show = {
        "BASE",
        "30+33+35+37",
        "30+33+35",
        "30+33+37",
        "30+35+37",
        "33+35+37",
        min(portfolio_rows, key=lambda item: item["alpha"])["portfolio"],
    }
    for row in portfolio_rows:
        if row["portfolio"] not in labels_to_show:
            continue
        label = "H4" if row["portfolio"] == "30+33+35+37" else row["portfolio"]
        if row["portfolio"] == "BASE":
            label = "BASE"
        axis.annotate(
            label,
            (positions[row["portfolio"]], row["alpha"]),
            xytext=(0, 7 if row["alpha"] >= 0 else -11),
            textcoords="offset points",
            ha="center",
            va="bottom" if row["alpha"] >= 0 else "top",
            fontsize=6.7,
            color="#30343b",
        )
    axis.axhline(0.0, color="#20252b", linewidth=0.8)
    axis.set_xticks(range(5), ["0", "1", "2", "3", "4"])
    axis.set_xlabel("H4 portfolio cardinality |H|")
    axis.set_ylabel(r"target-band spectral abscissa $\alpha$ [s$^{-1}$]")
    axis.set_title("F1 — Stability cliff over the H4 inclusion lattice")
    axis.grid(axis="y", color="#d8dde3", linewidth=0.5, alpha=0.8)
    axis.legend(frameon=False, ncol=4, loc="upper left", fontsize=8)
    figure.text(
        0.99,
        0.01,
        f"H4: α={h4['alpha']:.6f} s⁻¹, f={h4['dominant_frequency_hz']:.4f} Hz · full transverse spectra retained",
        ha="right",
        va="bottom",
        fontsize=7.2,
        color="#4b5563",
    )
    save_figure(figure, run_root / "figures" / "F1_stability_cliff")
    plt.close(figure)
    return {"portfolios": len(portfolio_rows), "spectra_rows": len(rows), "h4": h4}


def make_f2(repo: Path, run_root: Path, baseline: Path) -> dict:
    physical_path = baseline / "derived" / "TX4_PHYSICAL_LOCAL_FACTORS.csv"
    collective_path = baseline / "derived" / "TX4_CONTEXTUAL_RETURN_SWEEP.csv"
    summary_path = baseline / "derived" / "TX4_CONTEXTUAL_RETURN_CORE_SUMMARY.json"
    if not all(path.exists() for path in (physical_path, collective_path, summary_path)):
        raise FileNotFoundError("baseline F2 inputs are incomplete")
    with physical_path.open(encoding="utf-8") as handle:
        physical_rows = list(csv.DictReader(handle))
    with collective_path.open(encoding="utf-8") as handle:
        collective_rows = list(csv.DictReader(handle))
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if "physical_local_sigma_min" not in physical_rows[0]:
        raise RuntimeError("F2 requires physical pre-normalized local factors")
    source_copy = run_root / "raw" / "inputs"
    copy2(physical_path, source_copy / physical_path.name)
    copy2(collective_path, source_copy / collective_path.name)
    copy2(summary_path, source_copy / summary_path.name)

    by_g: dict[float, dict[int, float]] = {}
    for row in physical_rows:
        by_g.setdefault(float(row["g"]), {})[int(row["device_bus"])] = float(row["physical_local_sigma_min"])
    collective_by_g = {float(row["g"]): row for row in collective_rows}
    buses = list(CORE)
    rows: list[dict] = []
    for g in sorted(by_g):
        collective = collective_by_g.get(g)
        for bus in buses:
            rows.append(
                {
                    "g": g,
                    "device_bus": bus,
                    "physical_local_sigma_min": by_g[g][bus],
                    "collective_sigma_min": "" if collective is None else float(collective["collective_sigma_min"]),
                    "alpha": "" if collective is None else float(collective["alpha"]),
                    "boundary_g_star": float(summary["g_eigen_boundary"]),
                    "local_operator_convention": "I+M_ii physical pre-normalized",
                    "collective_operator_convention": "I+Q_H",
                }
            )
    write_csv(run_root / "derived" / "F2_LOCAL_COLLECTIVE.csv", rows)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(7.7, 5.2))
    colors = {30: "#386cb0", 33: "#f0027f", 35: "#1b9e77", 37: "#e6ab02"}
    for bus in buses:
        bus_rows = [row for row in rows if row["device_bus"] == bus]
        axis.semilogy(
            [row["g"] for row in bus_rows],
            [row["physical_local_sigma_min"] for row in bus_rows],
            "o-",
            color=colors[bus],
            linewidth=1.2,
            markersize=3.5,
            label=fr"physical local $I+M_{{{bus},{bus}}}$",
        )
    collective = [row for row in rows if row["collective_sigma_min"] != ""]
    axis.semilogy(
        [row["g"] for row in collective],
        [float(row["collective_sigma_min"]) for row in collective],
        "s--",
        color="#20252b",
        linewidth=1.8,
        markersize=3.8,
        label=r"collective $I+Q_H$",
    )
    g_star = float(summary["g_eigen_boundary"])
    axis.axvline(g_star, color="#8b1e3f", linewidth=1.0, linestyle=":")
    axis.text(g_star, 0.98, fr"$g^*={g_star:.6f}$", transform=axis.get_xaxis_transform(), ha="right", va="top", fontsize=8, color="#8b1e3f")
    axis.set_xlabel(r"feedback parameter $g$")
    axis.set_ylabel(r"minimum singular value $\sigma_{\min}$")
    axis.set_title(r"F2 — Physical local closure versus collective closure")
    axis.grid(which="both", color="#d8dde3", linewidth=0.5, alpha=0.8)
    axis.legend(frameon=False, fontsize=8, ncol=2, loc="best")
    figure.text(
        0.99,
        0.01,
        r"local: physical pre-normalized $I+M_{ii}$ · collective: $I+Q_H$",
        ha="right",
        va="bottom",
        fontsize=7.2,
        color="#4b5563",
    )
    save_figure(figure, run_root / "figures" / "F2_local_collective_closure")
    plt.close(figure)
    return {"physical_rows": len(physical_rows), "collective_rows": len(collective_rows), "g_star": g_star}


def main() -> int:
    if len(sys.argv) != 4:
        raise SystemExit("usage: render_f1_f2.py <repo-root> <run-root> <baseline-run-root>")
    repo = Path(sys.argv[1]).resolve()
    run_root = Path(sys.argv[2]).resolve()
    baseline = Path(sys.argv[3]).resolve()
    for name in ("raw", "raw/inputs", "derived", "figures", "claims", "report", "environment", "logs"):
        (run_root / name).mkdir(parents=True, exist_ok=True)
    if not baseline.is_dir():
        raise FileNotFoundError(baseline)

    f1 = make_f1(repo, run_root)
    f2 = make_f2(repo, run_root, baseline)
    source_files = {
        "tx4_contextual_return.py": repo / "reports/poster/ias2026/research/experiments/tx4_contextual_return.py",
        "ieee39_case.py": repo / "reports/poster/ias2026/research/src/ibr_cycles/models/ieee39_case.py",
        "port_admittance.py": repo / "reports/poster/ias2026/research/src/ibr_cycles/models/port_admittance.py",
        "ieee39_devices.py": repo / "reports/poster/ias2026/research/src/ibr_cycles/models/ieee39_devices.py",
    }
    manifest = {
        "run_root": str(run_root),
        "baseline_run": str(baseline),
        "git_commit": git(repo, "rev-parse", "HEAD"),
        "git_dirty": bool(git(repo, "status", "--porcelain")),
        "python": {"version": sys.version, "implementation": platform.python_implementation(), "executable": sys.executable},
        "packages": {"numpy": np.__version__},
        "source_hashes": {name: sha256(path) for name, path in source_files.items()},
        "output_policy": "all generated files are beneath this immutable run root; no zip is produced",
        "scientific_scope": "F1 exact 16-case spectra and F2 physical-local/collective asset regeneration only",
    }
    write_json(run_root / "environment" / "runtime.json", manifest)
    claims = {
        "status": "PASS",
        "F1": {"status": "READY", **f1, "full_spectra_retained": True, "inclusion_lattice_complete": f1["portfolios"] == 16},
        "F2": {"status": "READY", **f2, "physical_local_convention": "I+M_ii", "collective_convention": "I+Q_H"},
        "F3": {"status": "BLOCKED_M1_STRICT", "safe_to_run": False},
        "F4": {"status": "PENDING_DEDICATED_JULIA_LATTICE", "safe_to_run": False},
    }
    write_json(run_root / "claims" / "FIGURE_READINESS.json", claims)
    (run_root / "report" / "RUN_SUMMARY.md").write_text(
        "# F1/F2 poster asset regeneration\n\n"
        "Status: **PASS**\n\n"
        "This run regenerated the exact 16-case GFL11 target-band spectra for F1 and the physical-local/collective closure data for F2. It did not run M2, eta continuation, GFM substitutions, or TDS.\n\n"
        + json.dumps(claims, indent=2, default=float)
        + "\n",
        encoding="utf-8",
    )
    print(json.dumps(claims, indent=2, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
