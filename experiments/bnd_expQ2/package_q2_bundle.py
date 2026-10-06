from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

root = Path(__file__).resolve().parents[2]
bundle = root / "reports" / "experiment_Q2_bundle.zip"
roots = [
    root / "experiments" / "bnd_expQ2",
    root / "src" / "bnd_expQ2",
    root / "test" / "bnd_expQ2",
    root / "reports" / "experiment_Q2",
]
files = sorted(
    p for base in roots for p in base.rglob("*")
    if p.is_file() and p != bundle and "__pycache__" not in p.parts
)
manifest = root / "reports" / "experiment_Q2" / "BUNDLE_CONTENTS.txt"
manifest.write_text("Experiment Q2 archive content\n\n" + "\n".join(
    p.relative_to(root).as_posix() for p in files if p != manifest
) + "\n", encoding="utf-8")
files = sorted(set(files + [manifest]))
with ZipFile(bundle, "w", ZIP_DEFLATED, compresslevel=9) as archive:
    for path in files:
        archive.write(path, path.relative_to(root).as_posix())
print(f"created {bundle} with {len(files)} files ({bundle.stat().st_size} bytes)")
