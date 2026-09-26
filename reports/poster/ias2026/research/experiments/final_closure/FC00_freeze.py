"""FC00: freeze record of the starting point of the final closure campaign."""

from __future__ import annotations

import os
import platform
import subprocess
import sys

import numpy
import pandas
import scipy
from _fc import ROOT, RUN, sha256, write_json

REPO = ROOT.parents[3]
MODEL_FILES = [
    "configs/ias2026/ieee39_network.json",
    "configs/kundur/kundur_network.json",
    "configs/ieee68/ieee68_network.json",
    "src/ibr_cycles/models/ieee39_devices.py",
    "src/ibr_cycles/models/ieee68_devices.py",
    "src/ibr_cycles/models/ieee39_case.py",
    "src/ibr_cycles/models/ieee39_network.py",
    "src/ibr_cycles/dynamics/linearize.py",
    "experiments/_f7_common.py",
    "experiments/F8_service_attribution.py",
    "experiments/F12_kundur.py",
    "experiments/G3_ieee68.py",
    "experiments/G2_tds.py",
]


def git(*args):
    return subprocess.run(
        ["git", *args], cwd=REPO, capture_output=True, text=True
    ).stdout.strip()


def main() -> int:
    andes = subprocess.run(
        [
            str(REPO / ".venv" / "tx3-andes" / "Scripts" / "python.exe"),
            "-c",
            "import andes; print(andes.__version__)",
        ],
        capture_output=True,
        text=True,
    ).stdout.strip()
    try:
        import ctypes

        class MS(ctypes.Structure):
            _fields_ = [
                ("dwLength", ctypes.c_ulong),
                ("dwMemoryLoad", ctypes.c_ulong),
                ("ullTotalPhys", ctypes.c_ulonglong),
                ("a", ctypes.c_ulonglong),
                ("b", ctypes.c_ulonglong),
                ("c", ctypes.c_ulonglong),
                ("d", ctypes.c_ulonglong),
                ("e", ctypes.c_ulonglong),
                ("f", ctypes.c_ulonglong),
            ]

        ms = MS()
        ms.dwLength = ctypes.sizeof(MS)
        ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(ms))
        ram = ms.ullTotalPhys
    except Exception:  # pragma: no cover
        ram = None
    record = {
        "git_commit": git("rev-parse", "HEAD"),
        "git_branch": git("branch", "--show-current"),
        "git_tag": "IAS2026_PRE_FINAL_VALIDATION",
        "git_status_research": git(
            "status", "--short", "--", "reports/poster/ias2026/research"
        ),
        "python": sys.version,
        "numpy": numpy.__version__,
        "scipy": scipy.__version__,
        "pandas": pandas.__version__,
        "andes": andes,
        "platform": platform.platform(),
        "cpu": platform.processor(),
        "logical_cpus": os.cpu_count(),
        "ram_bytes": ram,
        "seeds": {
            "campaign_base": 20260920,
            "note": "each FC script declares its own seed",
        },
        "solver_tolerances": {
            "power_flow_newton": 1e-12,
            "dae_equilibrium_newton": 1e-9,
            "tds_BDF_rtol": 1e-7,
            "tds_BDF_atol": 1e-9,
            "tds_network_newton": 1e-10,
            "jacobian": "central differences (dynamics.linearize), unit and 2x scale",
        },
        "model_hashes": {f: sha256(ROOT / f) for f in MODEL_FILES},
        "frozen_tags": [
            "IAS2026_TRACKA_F7_POLICY_HYPERGRAPH_FREEZE",
            "IAS2026_TRACKA_F8_F12_POST_F7_FREEZE",
        ],
        "run_directory": str(RUN),
    }
    write_json(RUN / "FC00_freeze_record.json", record)
    print(record["git_commit"], record["run_directory"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
