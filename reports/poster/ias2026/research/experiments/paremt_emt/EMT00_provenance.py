# ruff: noqa: E501  -- long diagnostic strings and code-patch literals kept on one line
"""EMT00a - ParaEMT provenance and environment manifest (run with .venv/xtool-paremt).

Writes results/EMT00/{environment_manifest.txt, pip_freeze.txt, upstream_provenance.json}.
Records the official upstream clone (external/ParaEMT_upstream, never modified) and the
working copy (external/ParaEMT_tx4) with its local diff SHA.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
REPO = HERE.parents[5]
UP = REPO / "external" / "ParaEMT_upstream"
TX4 = REPO / "external" / "ParaEMT_tx4"
OUT = RESEARCH / "results" / "EMT00"


def git(path: Path, *args: str) -> str:
    return subprocess.run(
        ["git", "-C", str(path), *args], capture_output=True, text=True, check=True
    ).stdout.strip()


def diff_sha(path: Path) -> tuple[str, str]:
    diff = git(path, "diff", "HEAD")
    untracked = git(path, "ls-files", "--others", "--exclude-standard")
    blob = (diff + "\n--untracked--\n" + untracked).encode()
    return hashlib.sha256(blob).hexdigest(), git(path, "diff", "--stat", "HEAD")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    freeze = subprocess.run(
        [sys.executable, "-m", "pip", "freeze"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    (OUT / "pip_freeze.txt").write_text(freeze, encoding="utf-8")
    cpu = platform.processor()
    try:
        cpu = (
            subprocess.run(
                [
                    "powershell",
                    "-NoProfile",
                    "-Command",
                    "(Get-CimInstance Win32_Processor).Name",
                ],
                capture_output=True,
                text=True,
            ).stdout.strip()
            or cpu
        )
    except OSError:
        pass
    lines = [
        f"os: {platform.platform()}",
        f"cpu: {cpu}",
        f"logical_cpus: {os.cpu_count()}",
        f"python: {sys.version.split()[0]} ({sys.executable})",
        "packages (pip freeze):",
        *["  " + x for x in freeze.splitlines()],
        f"requirements_upstream: {(UP / 'requirements.txt').read_text(encoding='utf-8').strip()!r}",
    ]
    (OUT / "environment_manifest.txt").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    up_sha = git(UP, "rev-parse", "HEAD")
    tx4_sha = git(TX4, "rev-parse", "HEAD")
    up_diff, _ = diff_sha(UP)
    tx4_diff, tx4_stat = diff_sha(TX4)
    prov = {
        "upstream_repository": git(UP, "remote", "get-url", "origin"),
        "upstream_commit": up_sha,
        "upstream_commit_date": git(UP, "log", "-1", "--format=%cI"),
        "upstream_commit_subject": git(UP, "log", "-1", "--format=%s"),
        "clone_date": "2026-09-11",
        "upstream_local_diff_sha256": up_diff,
        "upstream_is_pristine": git(UP, "status", "--porcelain") == "",
        "working_copy": "external/ParaEMT_tx4",
        "working_copy_base_commit": tx4_sha,
        "working_copy_diff_sha256": tx4_diff,
        "working_copy_diff_stat": tx4_stat,
        "license_file": "LICENSE.md",
        "license_first_line": (UP / "LICENSE.md")
        .read_text(encoding="utf-8")
        .splitlines()[0],
        "note": "vendor/paraemt (same upstream commit, locally modified by earlier TX3 work) is NOT used and not touched.",
    }
    (OUT / "upstream_provenance.json").write_text(
        json.dumps(prov, indent=1), encoding="utf-8"
    )
    print(json.dumps(prov, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
