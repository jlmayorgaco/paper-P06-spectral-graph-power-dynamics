"""Bundle this experiment without touching earlier experiments or Git state."""
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

HERE = Path(__file__).resolve().parent
target = HERE / "analytical_delay_codesign_mega_20261002.zip"
with ZipFile(target, "w", compression=ZIP_DEFLATED, compresslevel=6) as archive:
    for path in sorted(HERE.rglob("*")):
        if path.is_file() and path != target:
            archive.write(path, arcname=f"{HERE.name}/{path.relative_to(HERE).as_posix()}")
print(target, target.stat().st_size)
