"""Package closure results and the repository inputs needed to rerun this stage."""
from __future__ import annotations

import hashlib
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
DEST = REPO / "deliverables" / "latency_robust_pll_closure_20261003_full_bundle.zip"

ROOT_FILES = ["Project.toml", "Manifest.toml"]
TREES = [
    "experiments/latency_robust_pll_closure_20261003",
    "experiments/latency_robust_pll_codesign_20261003",
    "experiments/delay_dressed_replacement_frontier_20261002",
    "experiments/nonlinear_codesign_20261001",
    "experiments/bnd_expQ2B",
    "src",
    "reports/experiment_A",
    "reports/experiment_D",
    "reports/experiment_N",
]


def main():
    DEST.parent.mkdir(parents=True, exist_ok=True)
    paths = [REPO / p for p in ROOT_FILES]
    for tree in TREES:
        paths.extend(p for p in (REPO / tree).rglob("*") if p.is_file())
    paths = sorted(set(paths))
    paths = [p for p in paths if "__pycache__" not in p.parts and p.suffix != ".pyc"]
    with zipfile.ZipFile(DEST, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6,
                         allowZip64=True) as archive:
        for path in paths:
            if not path.exists():
                raise FileNotFoundError(path)
            archive.write(path, path.relative_to(REPO).as_posix())
    with zipfile.ZipFile(DEST) as archive:
        bad = archive.testzip()
        if bad:
            raise RuntimeError(f"Corrupt ZIP member: {bad}")
        members = len(archive.namelist())
    digest = hashlib.sha256(DEST.read_bytes()).hexdigest()
    DEST.with_suffix(DEST.suffix + ".sha256").write_text(
        f"{digest}  {DEST.name}\n", encoding="utf-8")
    print(f"bundle={DEST}\nfiles={members}\nbytes={DEST.stat().st_size}\nsha256={digest}")


if __name__ == "__main__":
    main()
