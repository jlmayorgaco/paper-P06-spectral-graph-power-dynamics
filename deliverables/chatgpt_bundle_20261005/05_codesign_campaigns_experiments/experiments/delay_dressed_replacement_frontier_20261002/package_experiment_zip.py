"""Refresh the artifact hashes and create a complete ZIP distribution."""
from pathlib import Path
import subprocess
import sys
import zipfile

HERE = Path(__file__).resolve().parent
ZIP_PATH = HERE / "delay_dressed_replacement_frontier_20261002.zip"
HASH_SCRIPT = HERE / "hash_experiment_artifacts.py"

subprocess.run([sys.executable, str(HASH_SCRIPT)], check=True)
files = sorted(p for p in HERE.rglob("*") if p.is_file() and p != ZIP_PATH)
with zipfile.ZipFile(ZIP_PATH, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=8) as archive:
    for path in files:
        archive.write(path, arcname=path.relative_to(HERE.parent).as_posix())

with zipfile.ZipFile(ZIP_PATH) as archive:
    bad = archive.testzip()
    if bad is not None:
        raise RuntimeError(f"ZIP integrity check failed at {bad}")
    count = len(archive.namelist())
print(f"created {ZIP_PATH} with {count} files; ZIP integrity passed")
