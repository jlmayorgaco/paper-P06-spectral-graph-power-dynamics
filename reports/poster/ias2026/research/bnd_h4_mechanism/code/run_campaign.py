"""Run the frozen H4 contextual-return campaign into one immutable run root.

This adapter deliberately contains no model equations. It imports the frozen
TX4 implementation, redirects its output variable before any campaign call,
and writes gates/claims/manifests beside the generated tables.
"""

from __future__ import annotations

import csv
import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path


EXPECTED_HASHES = {
    "models/ieee39_case.py": "1f38a8e541780ab7205809dc1ba58bdafb6c923933f5772d46407064edc8c722",
    "models/ieee39_devices.py": "bcf0c76be65b171b580684065cde4eb88d0433c6b4d7beefec966f8711036d90",
    "models/ieee39_network.py": "f42144cd02a6dc779ed690fd866ad026f33eb7844584c7f10281fc7d98e3fe13",
    "models/port_admittance.py": "b0fa694dd0111be158133a9cf8c4e36be670e50ef562dcd63f89d7ff87b8d4ef",
    "models/port_core.py": "38c54b1a70f87fbbe2c301a717e6dd6b4112b3a9cf8deb663426ce08b729bb9b",
    "experiments/tx4_contextual_return.py": "77d4ed0f76f8c4a6ec369bc30b306361cff9b789c6f432e287d5f180d9210b74",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git(repo: Path, *args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=repo, text=True).strip()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def as_float(row: dict[str, str], key: str) -> float:
    return float(row[key])


def gate(name: str, passed: bool, observed: object, expected: object, tolerance: object) -> dict:
    return {
        "name": name,
        "status": "PASS" if passed else "BLOCKED",
        "observed": observed,
        "expected": expected,
        "tolerance": tolerance,
    }


def evaluate_gates(repo: Path, run_root: Path, audit: dict, core: dict) -> list[dict]:
    derived = run_root / "derived"
    truth = read_csv(derived / "TX4_CONTEXTUAL_RETURN_NUMERICAL_TRUTH.csv")
    proper = read_csv(derived / "TX4_PROPER_SUBSET_CLOSURE.csv")
    sweep = read_csv(derived / "TX4_CONTEXTUAL_RETURN_SWEEP.csv")
    boundary = read_csv(derived / "TX4_CONTEXTUAL_RETURN_AT_BOUNDARY.csv")
    derivative = read_csv(derived / "TX4_CONTEXTUAL_RETURN_DERIVATIVE.csv")

    h4 = [r for r in truth if r["portfolio"] == "30+33+35+37" and r["tolerance_factor"] == "1.0"]
    if len(h4) != 1:
        raise RuntimeError(f"expected one H4 truth row, found {len(h4)}")
    h4_row = h4[0]
    p4_alpha = as_float(h4_row, "alpha")
    p4_frequency = as_float(h4_row, "frequency_hz")
    root_row = min(sweep, key=lambda r: abs(as_float(r, "alpha")))
    source_hashes = {}
    src = repo / "reports/poster/ias2026/research"
    for rel, expected in EXPECTED_HASHES.items():
        path = src / ("src/ibr_cycles/" + rel if rel.startswith("models/") else rel)
        observed = sha256(path)
        source_hashes[rel] = {"observed": observed, "expected": expected, "match": observed == expected}

    max_eq_f = max(as_float(r, "equilibrium_f_residual") for r in truth)
    max_eq_g = max(as_float(r, "equilibrium_g_residual") for r in truth)
    max_eig_error = max(as_float(r, "numpy_scipy_eig_error") for r in truth)
    max_schur = max(as_float(r, "return_schur_residual_max") for r in sweep)
    max_port = max(as_float(r, "port_identity_residual") for r in sweep)
    min_local = min(as_float(r, "local_sigma_min") for r in sweep)
    root_collective = min(as_float(r, "collective_sigma_min") for r in sweep)
    proper_stable = all(r["stable"].lower() == "true" for r in proper)
    hash_pass = all(v["match"] for v in source_hashes.values())
    repo_status = git(repo, "status", "--porcelain")
    run_root_isolated = run_root.parent.name == "results" and run_root.parent.parent.name == "bnd_h4_mechanism"

    return [
        gate("G0 isolated immutable run root", run_root_isolated, str(run_root), ".../bnd_h4_mechanism/results/<RUN_ID>", True),
        gate("G1 frozen source hashes", hash_pass, source_hashes, "all registered hashes match", True),
        gate("G2 equilibrium residuals", max_eq_f <= 1e-7 and max_eq_g <= 1e-7, {"max_f": max_eq_f, "max_g": max_eq_g}, "<= 1e-7", 1e-7),
        gate("G3 H4 dynamic-state dimension", int(float(h4_row["state_dim"])) == 86, int(float(h4_row["state_dim"])), 86, 0),
        gate("G4 DAE/index-1 conditioning", all(as_float(r, "gz_condition_independent") < 1e14 for r in truth), max(as_float(r, "gz_condition_independent") for r in truth), "finite and < 1e14", 1e14),
        gate("G5 transverse target mode", 0.3 <= p4_frequency <= 1.5 and max_eig_error <= 1e-6, {"frequency_hz": p4_frequency, "numpy_scipy_error": max_eig_error}, "frequency in [0.3,1.5] Hz; eig error <= 1e-6", 1e-6),
        gate("G6 H4 P4 instability", abs(p4_alpha - 0.1270064671440848) <= 5e-6 and p4_alpha > 0.0 and abs(p4_frequency - 0.6222796695779256) <= 5e-6, {"alpha_s-1": p4_alpha, "frequency_hz": p4_frequency}, {"alpha_s-1": 0.1270064671440848, "frequency_hz": 0.6222796695779256}, 5e-6),
        gate("G7 proper-subset minimality", len(proper) == 15 and proper_stable, {"rows": len(proper), "all_stable": proper_stable}, "15 proper subsets stable", 0),
        gate("G8 eigenvalue boundary", abs(float(core["g_eigen_boundary"]) - 0.2076814051903784) <= 5e-5 and abs(as_float(root_row, "alpha")) <= 5e-5, {"g_root": core["g_eigen_boundary"], "nearest_grid_alpha": as_float(root_row, "alpha")}, 0.2076814051903784, 5e-5),
        gate("G9 port/Schur identities", max_schur <= 1e-8 and max_port <= 1e-8, {"max_schur": max_schur, "max_port": max_port}, "both <= 1e-8", 1e-8),
        gate("G10 collective not local", min_local >= 0.2973268809593455 - 1e-6 and root_collective <= 1e-6, {"min_local_sigma": min_local, "min_collective_sigma": root_collective}, {"min_local_sigma": ">= 0.2973259", "root_collective_sigma": "<= 1e-6"}, 1e-6),
        gate("G11 derivative and numerical audit", bool(audit["pass"]) and len(derivative) == 4 and float(core["max_return_derivative_error"]) <= 1e-4, {"audit": audit, "derivative_rows": len(derivative), "max_derivative_error": core["max_return_derivative_error"]}, "audit pass; 4 derivative rows; error <= 1e-4", 1e-4),
    ]


def write_report(run_root: Path, summary: dict, gates: list[dict], manifest: dict) -> None:
    claims = run_root / "claims"
    report = run_root / "report"
    claims.mkdir(parents=True, exist_ok=True)
    report.mkdir(parents=True, exist_ok=True)
    (claims / "gates.json").write_text(json.dumps(gates, indent=2), encoding="utf-8")
    status = "PASS" if all(item["status"] == "PASS" for item in gates) else "BLOCKED"
    (claims / "claims.md").write_text(
        "# Claims for isolated H4 run\n\n"
        f"Campaign status: **{status}**.\n\n"
        "## Safe claims\n\n"
        "- The reported values come from the frozen active Python TX4 path and are routed to this run root.\n"
        "- The H4 P4 mode, proper-subset minimality, collective boundary, local factors, Schur identity, and root-motion derivative are gated by G0--G11 below.\n\n"
        "## Unsafe claims\n\n"
        "- This run does not establish equivalence with the historical PowerDynamics TDS alternative model.\n"
        "- It does not justify a global stability claim outside the registered transverse frequency band or beyond the frozen IEEE-39/P4 contract.\n\n"
        "See `gates.json` for observed values and tolerances.\n",
        encoding="utf-8",
    )
    (report / "RUN_SUMMARY.md").write_text(
        "# H4 mechanism campaign run\n\n"
        f"Status: **{status}**\n\n"
        f"Run root: `{run_root}`\n\n"
        "The primary TX4 contextual-return campaign was executed through the frozen implementation with output redirection set before evaluation.\n\n"
        "## Gate status\n\n"
        + "\n".join(f"- {item['name']}: **{item['status']}** — observed `{item['observed']}`" for item in gates)
        + "\n\n## Core summary\n\n```json\n"
        + json.dumps(summary, indent=2)
        + "\n```\n",
        encoding="utf-8",
    )
    (run_root / "environment" / "runtime.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("usage: run_campaign.py <repo-root> <run-root>")
    repo = Path(sys.argv[1]).resolve()
    run_root = Path(sys.argv[2]).resolve()
    for name in ("raw/python", "raw/julia", "derived", "tables", "figures", "claims", "report", "environment", "logs"):
        (run_root / name).mkdir(parents=True, exist_ok=True)

    experiment = repo / "reports/poster/ias2026/research/experiments"
    source = repo / "reports/poster/ias2026/research/src"
    sys.path.insert(0, str(experiment))
    sys.path.insert(0, str(source))
    import tx4_contextual_return as tx  # type: ignore  # noqa: PLC0415

    tx.OUT = run_root / "derived"
    tx.OUT.mkdir(parents=True, exist_ok=True)
    audit = tx.audit_cases()
    core = tx.sweep_and_tables()
    summary = {"audit": audit, "core": core}
    (run_root / "claims" / "tx4_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    source_hashes = {}
    for rel in EXPECTED_HASHES:
        path = source / "ibr_cycles" / rel if rel.startswith("models/") else experiment / "tx4_contextual_return.py"
        source_hashes[rel] = sha256(path)
    manifest = {
        "run_root": str(run_root),
        "git_commit": git(repo, "rev-parse", "HEAD"),
        "git_shortsha": git(repo, "rev-parse", "--short=8", "HEAD"),
        "git_dirty": bool(git(repo, "status", "--porcelain")),
        "python": {"version": sys.version, "implementation": platform.python_implementation(), "executable": sys.executable},
        "packages": {
            "numpy": __import__("numpy").__version__,
            "scipy": __import__("scipy").__version__,
        },
        "command": [sys.executable, str(Path(__file__).resolve()), str(repo), str(run_root)],
        "source_hashes": source_hashes,
        "output_policy": "all generated files are beneath this immutable run root; no zip is produced",
    }
    gates = evaluate_gates(repo, run_root, audit, core)
    write_report(run_root, summary, gates, manifest)
    print(json.dumps({"status": "PASS" if all(g["status"] == "PASS" for g in gates) else "BLOCKED", "run_root": str(run_root), "gates": gates}, indent=2))
    return 0 if all(g["status"] == "PASS" for g in gates) else 2


if __name__ == "__main__":
    raise SystemExit(main())
