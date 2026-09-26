# ruff: noqa: E501  -- sympy equation strings and parameter tables kept on one line
"""Create the isolated ANDES environment for the PCV Phase 8 reproduction (idempotent).

Spec: docs/20260911_GFL_REPRODUCTION_SPEC.md section 2. Run with any Python >= 3.11:

    python setup_xtool_env.py

Steps (each skipped when already done):
1. venv .venv/xtool-andes-gfl from the base CPython of .venv/tx3-andes;
2. pinned dependencies (the tx3-andes versions);
3. pristine ANDES copy: git archive of vendor/andes tag v2.0.0 -> src/andes;
4. editable install of the copy (SETUPTOOLS_SCM_PRETEND_VERSION=2.0.0);
5. copy pcv_models.py into the copy and register it in andes/models/__init__.py;
6. private home (USERPROFILE) and single-process code generation into it.
Never touches .venv/tx3-andes, vendor/andes or ~/.andes.
"""

from __future__ import annotations

import io
import os
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[6]
VENV = REPO / ".venv" / "xtool-andes-gfl"
SRC = VENV / "src" / "andes"
HOME = VENV / "home"
PY = VENV / "Scripts" / "python.exe"
PINS = [
    "numpy==2.5.2",
    "scipy==1.18.1",
    "sympy==1.14.0",
    "pandas==3.0.5",
    "kvxopt==1.3.3.1",
    "matplotlib==3.11.1",
    "openpyxl==3.1.5",
    "xlsxwriter==3.2.9",
    "dill==0.4.1",
    "pathos==0.3.5",
    "tqdm==4.70.0",
    "ipywidgets==8.1.9",
    "pyyaml==6.0.3",
    "chardet==7.6.0",
    "texttable==1.7.0",
    "setuptools>=64",
    "setuptools-scm>=8.0",
]
REGISTRATION = "    ('pcv_models', ['SG2AX', 'GFL11']),\n"
ANCHOR = "    ('static', ['PQ', 'PV', 'Slack']),\n"


def base_python() -> str:
    cfg = (REPO / ".venv" / "tx3-andes" / "pyvenv.cfg").read_text().splitlines()
    home = next(
        line.split("=", 1)[1].strip() for line in cfg if line.startswith("home")
    )
    return str(Path(home) / "python.exe")


def main() -> int:
    if not PY.exists():
        subprocess.run([base_python(), "-m", "venv", str(VENV)], check=True)
        subprocess.run([str(PY), "-m", "pip", "install", "-q", *PINS], check=True)
    if not (SRC / "pyproject.toml").exists():
        blob = subprocess.run(
            [
                "git",
                "-C",
                str(REPO / "vendor" / "andes"),
                "archive",
                "--format=zip",
                "v2.0.0",
            ],
            check=True,
            capture_output=True,
        ).stdout
        zipfile.ZipFile(io.BytesIO(blob)).extractall(SRC)
    probe = subprocess.run(
        [str(PY), "-c", "import andes; print(andes.__file__)"],
        capture_output=True,
        text=True,
    )
    if str(SRC) not in probe.stdout:
        env = dict(os.environ, SETUPTOOLS_SCM_PRETEND_VERSION="2.0.0")
        subprocess.run(
            [
                str(PY),
                "-m",
                "pip",
                "install",
                "-q",
                "--no-deps",
                "--no-build-isolation",
                "-e",
                str(SRC),
            ],
            check=True,
            env=env,
        )
    models = SRC / "andes" / "models"
    shutil.copyfile(HERE / "pcv_models.py", models / "pcv_models.py")
    init = (models / "__init__.py").read_text(encoding="utf-8")
    if REGISTRATION not in init:
        assert ANCHOR in init
        (models / "__init__.py").write_text(
            init.replace(ANCHOR, ANCHOR + REGISTRATION), encoding="utf-8"
        )
    HOME.mkdir(parents=True, exist_ok=True)
    # single-process code generation into the private pycode folder (Windows spawn +
    # parallel codegen deadlocked once; nomp avoids worker processes altogether)
    env = dict(os.environ, USERPROFILE=str(HOME), HOME=str(HOME))
    subprocess.run(
        [
            str(PY),
            "-c",
            f"import andes; andes.config_logger(stream_level=40); "
            f"andes.prepare(quick=True, nomp=True, pycode_path={str(HOME / 'pycode')!r})",
        ],
        check=True,
        env=env,
    )
    check = subprocess.run(
        [
            str(PY),
            "-c",
            "import importlib.util, andes; print(andes.__version__, importlib.util.find_spec('ibr_cycles'))",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    print("xtool-andes-gfl ready:", check.stdout.strip())
    return 0


if __name__ == "__main__":
    sys.exit(main())
