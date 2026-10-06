"""Freeze and verify the active ExpN environment against Experiment M."""
from __future__ import annotations

import hashlib
import subprocess
import tomllib
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports" / "experiment_N"
M = ROOT / "reports" / "experiment_M" / "SOFTWARE_PROVENANCE.toml"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def quote(value: object) -> str:
    return '"' + str(value).replace('\\', '\\\\').replace('"', '\\"') + '"'


def command(*args: str) -> str:
    return subprocess.check_output(args, cwd=ROOT, text=True).strip()


def main() -> None:
    previous = tomllib.loads(M.read_text(encoding="utf-8"))
    current_project = sha(ROOT / "Project.toml")
    current_manifest = sha(ROOT / "Manifest.toml")
    current_source = sha(Path(previous["ieee39_source"]["path"]))
    input_hashes = {name: sha(Path(record["path"]))
                    for name, record in previous["ieee39_inputs"].items()}
    checks = {
        "project": current_project == previous["project_sha256"],
        "manifest": current_manifest == previous["manifest_sha256"],
        "ieee39_source": current_source == previous["ieee39_source"]["sha256"],
        **{f"ieee39_{name}": value == previous["ieee39_inputs"][name]["sha256"]
           for name, value in input_hashes.items()},
    }
    if not all(checks.values()):
        raise RuntimeError(f"active environment differs from frozen ExpM: {checks}")
    manifest = tomllib.loads((ROOT / "Manifest.toml").read_text(encoding="utf-8"))
    deps = manifest["deps"]
    def pkg(name: str) -> dict:
        entry = deps[name]
        return entry[0] if isinstance(entry, list) else entry
    lines = [
        f"frozen_utc = {quote(datetime.now(timezone.utc).isoformat())}",
        f"parent_commit = {quote(command('git', 'rev-parse', 'HEAD'))}",
        f"branch = {quote(command('git', 'branch', '--show-current'))}",
        f"julia_version = {quote(command('julia', '--version'))}",
        f"project_sha256 = {quote(current_project)}",
        f"manifest_sha256 = {quote(current_manifest)}",
        f"ieee39_source_path = {quote(previous['ieee39_source']['path'])}",
        f"ieee39_source_sha256 = {quote(current_source)}",
        f"same_environment_as_expM = {str(all(checks.values())).lower()}",
        "",
    ]
    for name in ("PowerDynamics", "NetworkDynamics"):
        entry = pkg(name)
        lines += [f"[packages.{name}]", f"version = {quote(entry['version'])}",
                  f"git_tree_sha1 = {quote(entry['git-tree-sha1'])}", ""]
    for name, value in sorted(input_hashes.items()):
        lines += [f"[ieee39_inputs.{name}]", f"sha256 = {quote(value)}", ""]
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "SOFTWARE_PROVENANCE.toml"
    path.write_text("\n".join(lines), encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
