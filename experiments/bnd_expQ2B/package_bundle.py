from __future__ import annotations

import hashlib
from pathlib import Path
import zipfile


ROOT = Path(__file__).resolve().parents[2]
REPORT_DIR = ROOT / "reports" / "experiment_Q2B"
ARCHIVE = ROOT / "reports" / "experiment_Q2B_bundle.zip"
CHECKSUM = ROOT / "reports" / "experiment_Q2B_bundle.zip.sha256"
README = REPORT_DIR / "BUNDLE_README.md"

README.write_text(
    "# Experiment Q2B bundle\n\n"
    "This archive contains the Q2B report, figures, CSV tables and sensor traces, "
    "the frozen candidate and validation records, plus the experiment scripts, "
    "source module and tests. Paths inside the archive are relative to the "
    "repository root.\n\n"
    "Primary status: `FAIL_OPTIMIZATION`. The best-found 100 MW candidate was "
    "independently validated in PowerDynamics for the tested event, but it is "
    "not a locally KKT-certified optimum. The normalized robustness interval "
    "bound is incomplete and global optimality is not established. See "
    "`reports/experiment_Q2B/REPORT_EXP_Q2B.md`, "
    "`reports/experiment_Q2B/FINAL_SUMMARY_EXP_Q2B.md`, and "
    "`reports/experiment_Q2B/TERMINAL_SUMMARY_EXP_Q2B.txt`.\n\n"
    "No dependency manifest was changed. No commit or push was made.\n",
    encoding="utf-8",
)

include_roots = [
    ROOT / "reports" / "experiment_Q2B",
    ROOT / "experiments" / "bnd_expQ2B",
    ROOT / "src" / "bnd_expQ2B",
    ROOT / "test" / "bnd_expQ2B",
]

files: list[Path] = []
for base in include_roots:
    if not base.exists():
        raise FileNotFoundError(f"Required bundle path does not exist: {base}")
    for path in base.rglob("*"):
        if not path.is_file():
            continue
        if "__pycache__" in path.parts or path.suffix in {".pyc", ".pyo"}:
            continue
        if path.resolve() in {ARCHIVE.resolve(), CHECKSUM.resolve()}:
            continue
        files.append(path)

files = sorted(set(files), key=lambda p: p.relative_to(ROOT).as_posix())
if not files:
    raise RuntimeError("No files found for Q2B bundle")

manifest_lines = ["sha256  repository-relative-path"]
for path in files:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest_lines.append(f"{digest}  {path.relative_to(ROOT).as_posix()}")
manifest = "\n".join(manifest_lines) + "\n"

ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
with zipfile.ZipFile(ARCHIVE, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
    for path in files:
        name = path.relative_to(ROOT).as_posix()
        info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
        info.compress_type = zipfile.ZIP_DEFLATED
        info.create_system = 3
        info.external_attr = (0o100644 & 0xFFFF) << 16
        zf.writestr(info, path.read_bytes(), compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
    info = zipfile.ZipInfo("Q2B_BUNDLE_SHA256_MANIFEST.txt", date_time=(1980, 1, 1, 0, 0, 0))
    info.compress_type = zipfile.ZIP_DEFLATED
    info.create_system = 3
    info.external_attr = (0o100644 & 0xFFFF) << 16
    zf.writestr(info, manifest.encode("utf-8"), compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)

archive_hash = hashlib.sha256(ARCHIVE.read_bytes()).hexdigest()
CHECKSUM.write_text(f"{archive_hash}  {ARCHIVE.name}\n", encoding="ascii")
print(f"files={len(files)}")
print(f"archive={ARCHIVE}")
print(f"archive_bytes={ARCHIVE.stat().st_size}")
print(f"sha256={archive_hash}")
