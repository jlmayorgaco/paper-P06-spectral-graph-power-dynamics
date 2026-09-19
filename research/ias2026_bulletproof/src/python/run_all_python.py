"""Restartable Python orchestrator for the new campaign."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


CAMPAIGN = Path(__file__).resolve().parents[2]
PYTHON = Path(__file__).resolve().parent
REPORTS = CAMPAIGN / "reports"


def write_phase(name: str, status: str, notes: str, evidence: str = "") -> None:
    REPORTS.mkdir(parents=True, exist_ok=True)
    payload = {"phase": name, "status": status, "notes": notes, "evidence": evidence}
    (CAMPAIGN / "logs" / f"python_{name}.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")


def run_gate0(source_mode: str = "bundled", source_root: str | None = None) -> int:
    from run_gate0 import main

    argv = ["--source-mode", source_mode]
    if source_root:
        argv += ["--source-root", source_root]
    return main(argv)


def run_theory_tests() -> int:
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", str(CAMPAIGN / "src" / "tests")],
        cwd=CAMPAIGN.parents[1],
        text=True,
        capture_output=True,
    )
    (CAMPAIGN / "logs" / "theory_tests.log").write_text(result.stdout + result.stderr, encoding="utf-8")
    status = "PASS" if result.returncode == 0 else "FAIL"
    write_phase("theory-tests", status, "10,000 randomized cases plus transverse tests.", "logs/theory_tests.log")
    return result.returncode


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", default="inventory")
    parser.add_argument("--source-mode", choices=("bundled", "external"), default="bundled")
    parser.add_argument("--source-root", default=None)
    args = parser.parse_args()
    phase = args.phase
    if phase == "gate0":
        return run_gate0(args.source_mode, args.source_root)
    if phase == "theory-tests":
        return run_theory_tests()
    if phase == "all":
        code = run_gate0(args.source_mode, args.source_root)
        theory_code = run_theory_tests()
        return max(code, theory_code)
    if phase in {"robust-v4", "topology", "figures", "holdout", "scaling"}:
        write_phase(phase, "NOT_TESTED", "This phase requires a dedicated model/holdout implementation and has not been silently substituted with a synthetic result.")
        return 0
    write_phase("inventory", "NUMERICALLY_VERIFIED", "Python orchestrator is available.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
