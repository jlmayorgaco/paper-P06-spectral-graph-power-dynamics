"""Package the final independent validation bundle with an internal SHA256 manifest."""

from __future__ import annotations

import hashlib
import json
import shutil
import zipfile
from datetime import datetime
from pathlib import Path

from campaign_root import campaign_root_from_argv


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    campaign = campaign_root_from_argv()
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    name = f"IAS2026_LAST_INDEPENDENT_VALIDATION_{stamp}"
    bundle_root = campaign / "bundle" / name
    zip_path = campaign / "bundle" / f"{name}.zip"
    bundle_root.mkdir(parents=True, exist_ok=False)

    include = ["raw", "code", "env", "prereg", "reports", "docs", "derived", "figures", "logs"]
    for relative in include:
        source = campaign / relative
        if source.exists():
            shutil.copytree(source, bundle_root / relative, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    for relative in ("README.md", "run_all_last_validation.ps1"):
        source = campaign / relative
        if source.exists():
            shutil.copy2(source, bundle_root / relative)

    readme = bundle_root / "README.md"
    readme.write_text(
        f"# {name}\n\n"
        "This archive is the last independent IAS2026 validation bundle. It includes the frozen custom Python/Julia reconciliation, Julia collective mechanism, same-model TDS traces, corrected SimpleGFLDC scope, second-model policy search, negative holdout verdict, claim ledger, and standalone runners.\n\n"
        "The prohibited generated FINAL_PAPER.pdf and FINAL_POSTER_DRAFT.pdf artifacts are intentionally absent.\n",
        encoding="utf-8",
    )

    forbidden = {"FINAL_PAPER.pdf", "FINAL_POSTER_DRAFT.pdf"}
    removed = []
    for path in list(bundle_root.rglob("*")):
        if path.is_file() and path.name in forbidden:
            removed.append(str(path.relative_to(bundle_root)))
            path.unlink()
    files = sorted(path for path in bundle_root.rglob("*") if path.is_file() and path.name != "SHA256SUMS.txt")
    sums = "".join(f"{digest(path)}  {path.relative_to(bundle_root).as_posix()}\n" for path in files)
    (bundle_root / "SHA256SUMS.txt").write_text(sums, encoding="utf-8")
    manifest = {
        "bundle": name,
        "created_local": datetime.now().isoformat(),
        "files": len(files),
        "excluded_forbidden_artifacts": sorted(forbidden),
        "removed_if_present": removed,
        "required_components": include,
    }
    (bundle_root / "package_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    # Add the manifest itself to the checksum file after it is finalized.
    with (bundle_root / "SHA256SUMS.txt").open("a", encoding="utf-8") as handle:
        handle.write(f"{digest(bundle_root / 'package_manifest.json')}  package_manifest.json\n")
    files_for_zip = sorted(path for path in bundle_root.rglob("*") if path.is_file())
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for path in files_for_zip:
            archive.write(path, path.relative_to(bundle_root.parent).as_posix())
    print(zip_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
