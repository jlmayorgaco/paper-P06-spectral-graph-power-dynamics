"""Freeze the local software and IEEE-39 source inputs before model construction."""
from __future__ import annotations

import hashlib
import subprocess
import tomllib
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports" / "experiment_M"
PKG_ROOT = Path.home() / ".julia" / "packages"
PACKAGES = (
    "PowerDynamics", "NetworkDynamics", "ModelingToolkit",
    "ModelingToolkitBase", "SciMLBase", "ForwardDiff", "OrdinaryDiffEqRosenbrock",
    "IEEE39",
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def q(value: object) -> str:
    return '"' + str(value).replace("\\", "\\\\").replace('"', '\\"') + '"'


def cmd(*args: str) -> str:
    return subprocess.check_output(args, cwd=ROOT, text=True).strip()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    project = ROOT / "Project.toml"
    manifest = ROOT / "Manifest.toml"
    data = tomllib.loads(manifest.read_text(encoding="utf-8"))
    lines = [
        f"frozen_utc = {q(datetime.now(timezone.utc).isoformat())}",
        f"repository_commit = {q(cmd('git', 'rev-parse', 'HEAD'))}",
        f"repository_branch = {q(cmd('git', 'branch', '--show-current'))}",
        f"julia_version = {q(cmd('julia', '--version'))}",
        f"project_path = {q(project)}",
        f"project_sha256 = {q(sha(project))}",
        f"manifest_path = {q(manifest)}",
        f"manifest_sha256 = {q(sha(manifest))}",
        "",
    ]
    installed = data.get("deps", {})
    for name in PACKAGES:
        entry = installed.get(name, [])
        if isinstance(entry, list):
            entry = entry[0] if entry else {}
        lines += [f"[packages.{name}]", f"installed = {str(bool(entry)).lower()}"]
        for key in ("uuid", "version", "git-tree-sha1", "repo-rev", "repo-url", "path"):
            if key in entry:
                lines.append(f"{key.replace('-', '_')} = {q(entry[key])}")
        lines.append("")
    # Identify the exact installed source from the path used by this repository.
    pd_example = PKG_ROOT / "PowerDynamics" / "VzOiZ" / "docs" / "examples" / "ieee39_part1.jl"
    if pd_example.is_file():
        lines += [
            "[ieee39_source]",
            f"path = {q(pd_example)}",
            f"sha256 = {q(sha(pd_example))}",
            "origin = \"PowerDynamics installed example\"",
            "",
        ]
        for path in sorted((pd_example.parent / "ieee39data").glob("*.csv")):
            lines += [f"[ieee39_inputs.{path.stem}]", f"path = {q(path)}", f"sha256 = {q(sha(path))}", ""]
    (OUT / "SOFTWARE_PROVENANCE.toml").write_text("\n".join(lines), encoding="utf-8")
    print(OUT / "SOFTWARE_PROVENANCE.toml")


if __name__ == "__main__":
    main()
