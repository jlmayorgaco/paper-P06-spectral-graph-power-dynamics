"""Provenance harness for the overnight Track-A final-validation run.

Every experiment writes into one timestamped run directory and records what it
would take to reproduce it: config, git commit, seed, environment, worker count,
accepted and rejected samples with reasons, raw and summary tables, a manifest
and a verdict. Nothing here writes to the frozen v1/v2A/v2B/v2C outputs.

BLAS threads are pinned to one inside workers. The linear algebra in this project
is small and dense, so process parallelism across cases is what pays; letting
each worker also spawn BLAS threads only oversubscribes 22 cores.
"""

from __future__ import annotations

import json
import os
import platform
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

RESEARCH = Path(__file__).resolve().parents[1]
REPO = RESEARCH.parents[3]
RUN_ROOT = RESEARCH / "outputs" / "ias2026"
STAMP_FILE = RUN_ROOT / "CURRENT_RUN"


def pin_blas_threads() -> None:
    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ):
        os.environ[name] = "1"


def git_commit() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO,
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
    except Exception:  # noqa: BLE001 - provenance must never break a run
        return "unknown"


def git_dirty() -> bool:
    try:
        out = subprocess.run(
            ["git", "status", "--porcelain"],
            cwd=REPO,
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        return bool(out.strip())
    except Exception:  # noqa: BLE001
        return True


def available_ram_gb() -> float:
    try:
        import ctypes

        class Status(ctypes.Structure):
            _fields_ = [
                ("dwLength", ctypes.c_ulong),
                ("dwMemoryLoad", ctypes.c_ulong),
                ("ullTotalPhys", ctypes.c_ulonglong),
                ("ullAvailPhys", ctypes.c_ulonglong),
                ("ullTotalPageFile", ctypes.c_ulonglong),
                ("ullAvailPageFile", ctypes.c_ulonglong),
                ("ullTotalVirtual", ctypes.c_ulonglong),
                ("ullAvailVirtual", ctypes.c_ulonglong),
                ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
            ]

        status = Status()
        status.dwLength = ctypes.sizeof(Status)
        ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status))
        return status.ullAvailPhys / 2**30
    except Exception:  # noqa: BLE001
        return float("nan")


def choose_workers(*, reserve: int = 4, ram_per_worker_gb: float = 0.6) -> int:
    """Conservative worker count: leave cores free and never exceed the RAM."""

    cores = os.cpu_count() or 1
    by_core = max(1, cores - reserve)
    ram = available_ram_gb()
    by_ram = int(ram / ram_per_worker_gb) if ram == ram else by_core
    return int(max(1, min(by_core, by_ram, 16)))


def environment() -> dict[str, Any]:
    import numpy
    import pandas
    import scipy

    return {
        "python": sys.version.split()[0],
        "executable": sys.executable,
        "platform": f"{platform.system()} {platform.release()} {platform.machine()}",
        "numpy": numpy.__version__,
        "scipy": scipy.__version__,
        "pandas": pandas.__version__,
        "cpu_count": os.cpu_count(),
        "available_ram_gb": round(available_ram_gb(), 2),
        "longdouble_mantissa_bits": int(numpy.finfo(numpy.longdouble).nmant) + 1,
    }


def run_directory(create: bool = True) -> Path:
    """The single run directory for tonight, stable across experiment scripts."""

    RUN_ROOT.mkdir(parents=True, exist_ok=True)
    if STAMP_FILE.exists():
        path = Path(STAMP_FILE.read_text(encoding="utf-8").strip())
        if path.exists() or not create:
            return path
    stamp = time.strftime("%Y%m%dT%H%M%S")
    path = RUN_ROOT / f"final_validation_overnight_{stamp}"
    if create:
        path.mkdir(parents=True, exist_ok=True)
        STAMP_FILE.write_text(str(path), encoding="utf-8")
    return path


@dataclass
class Experiment:
    """One overnight experiment: its own folder, manifest, verdict and log."""

    name: str
    question: str
    seed: int | None = None
    config: dict[str, Any] = field(default_factory=dict)
    workers: int = 1
    started: float = field(default_factory=time.time)
    notes: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.root = run_directory()
        self.directory = self.root / self.name
        self.directory.mkdir(parents=True, exist_ok=True)

    def path(self, filename: str) -> Path:
        return self.directory / filename

    def note(self, text: str) -> None:
        self.notes.append(text)
        print("    " + text)

    def save_table(self, frame, filename: str, *, parquet: bool = False) -> Path:
        target = self.path(filename)
        frame.to_csv(target, index=False)
        if parquet:
            try:
                frame.to_parquet(target.with_suffix(".parquet"), index=False)
            except Exception as error:  # noqa: BLE001
                self.note(f"parquet unavailable: {str(error)[:70]}")
        return target

    def finish(self, verdict: str, **results: Any) -> Path:
        payload = {
            "experiment": self.name,
            "question": self.question,
            "verdict": verdict,
            "seed": self.seed,
            "config": self.config,
            "workers": self.workers,
            "git_commit": git_commit(),
            "git_dirty": git_dirty(),
            "environment": environment(),
            "started_utc": time.strftime(
                "%Y-%m-%dT%H:%M:%SZ", time.gmtime(self.started)
            ),
            "elapsed_s": round(time.time() - self.started, 2),
            "notes": self.notes,
            "results": results,
        }
        target = self.path("manifest.json")
        target.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
        line = f"{self.name}: {verdict} ({payload['elapsed_s']:.0f}s)"
        with (self.root / "PROGRESS.log").open("a", encoding="utf-8") as handle:
            handle.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}  {line}\n")
        print(f"  -> {verdict}   manifest {target}")
        return target


def completed(name: str) -> bool:
    """True when this experiment already finished in the current run directory."""

    return (run_directory() / name / "manifest.json").exists()
