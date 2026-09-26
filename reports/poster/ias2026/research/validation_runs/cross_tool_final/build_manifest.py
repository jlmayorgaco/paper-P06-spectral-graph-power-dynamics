# ruff: noqa: E501  -- evidence labels, docstrings and verbatim source quotes kept on one line
"""Write results/CROSS_TOOL_MANIFEST.json: versions, hashes, tolerances, integrity.

Usage (any venv, research root): python build_manifest.py . <outer_handoff_zip>
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

repo = Path(sys.argv[1]).resolve()
outer_zip = Path(sys.argv[2])
HERE = Path(__file__).resolve().parent
HANDOFF = repo / "validation_inputs/cross_tool_handoff/claude_cross_tool_handoff"
PACK_SRC = HANDOFF / "ieee39_cross_tool_validation_pack/ieee39_cross_tool_validation"
FOREST = HANDOFF / "dynamic_forest_ieee39_validation/dynamic_forest_ieee39_validation"
ROOT = repo.parents[3]


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def rel(p: Path) -> str:
    try:
        return str(p.resolve().relative_to(ROOT)).replace("\\", "/")
    except ValueError:
        return p.name


def versions(name: str) -> dict:
    txt = (
        (HERE / f"supp/venvs/{name}.after.txt")
        .read_text(encoding="utf-8-sig")
        .splitlines()
    )
    pk = dict(line.split("==", 1) for line in txt[2:] if "==" in line)
    keep = (
        "andes",
        "pandapower",
        "numpy",
        "scipy",
        "pandas",
        "matplotlib",
        "kvxopt",
        "sympy",
        "pytest",
    )
    return {"python": txt[1].split()[0], **{k: pk[k] for k in keep if k in pk}}


def integrity(name: str) -> bool:
    b = HERE / f"supp/venvs/{name}.before.txt"
    a = HERE / f"supp/venvs/{name}.after.txt"
    return (
        b.read_text(encoding="utf-8-sig").splitlines()[2:]
        == a.read_text(encoding="utf-8-sig").splitlines()[2:]
    )


andes_static = json.loads((HERE / "supp/andes_static_supp.json").read_text())
inputs = [
    *(
        HANDOFF / n
        for n in (
            "ieee39_cross_tool_validation_pack.zip",
            "dynamic_forest_ieee39_validation.zip",
        )
    ),
    *sorted(PACK_SRC.glob("*.py")),
    PACK_SRC / "README.md",
    FOREST / "holdout_lines.json",
    FOREST / "data/F2c_holdout_reequilibrated_port_line_sensitivity.csv",
    FOREST / "dynamic_forest_line_port_holdout.py",
    FOREST / "CLAUDE_FINAL_ANDES_NETWORK_VALIDATION.md",
    repo / "configs/ias2026/ieee39_network.json",
    repo / "configs/ias2026/ieee39_andes_powerflow_reference.json",
    repo / "configs/ias2026/ieee39_ybus_andes.npy",
    repo / "results/F1_eigenvalue_reconciliation.csv",
    repo / "results/F1/F1_andes_spectra.npz",
    repo / "results/F1/andes_live_powerflow.json",
    repo / "experiments/F1_andes_equivalent_worker.py",
    repo / "experiments/F1_reconciliation.py",
    Path(andes_static["case_file"]),
]
scripts = sorted(HERE.glob("*.py"))
outputs = [
    repo / "docs/FINAL_CROSS_TOOL_VALIDATION.md",
    *(
        repo / "results" / n
        for n in (
            "CROSS_TOOL_EVIDENCE_TABLE.csv",
            "PANDAPOWER_STATIC_PARITY.csv",
            "ANDES_DYNAMIC_PARITY.csv",
            "ANDES_LINE_SENSITIVITY_HOLDOUT.csv",
            "CROSS_TOOL_VALIDATION_FIGURE.pdf",
            "CROSS_TOOL_VALIDATION_FIGURE.png",
            "CROSS_TOOL_VALIDATION_FIGURE_source.csv",
        )
    ),
    *sorted((HERE / "pack/results").glob("*")),
    *sorted((HERE / "pack_tapside/results").glob("*")),
    *sorted((HERE / "supp").glob("*.json")),
    *sorted((HERE / "supp").glob("*.csv")),
    *sorted((HERE / "phaseE").glob("*")),
]
pack_copies = {
    p.name: sha(p) == sha(PACK_SRC / p.name) for p in (HERE / "pack").glob("*.py")
}

manifest = {
    "campaign": "IEEE-39 final cross-tool validation (Phases A-G)",
    "date": "2026-09-11",
    "git_head_at_run": subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=repo, capture_output=True, text=True
    ).stdout.strip(),
    "handoff_outer_zip": {
        "name": outer_zip.name,
        "sha256": sha(outer_zip) if outer_zip.exists() else None,
    },
    "environments": {
        "tx3-analysis (Phase A, D-supp, E internal, assembly)": versions(
            "tx3-analysis"
        ),
        "tx3-andes (Phase C/D, E ANDES)": versions("tx3-andes"),
        "xtool-pandapower (Phase B; created for this campaign)": versions(
            "xtool-pandapower"
        ),
    },
    "frozen_venvs_package_lists_unchanged": {
        n: integrity(n) for n in ("bc-cert", "tx3-analysis", "tx3-andes", "tx3-paraemt")
    },
    "andes": {
        "version": andes_static["andes_version"],
        "case_file": rel(Path(andes_static["case_file"])),
        "case_sha256": andes_static["case_sha256"],
        "case_hash_matches_canonical_source": andes_static[
            "case_hash_matches_canonical_source"
        ],
        "vendor_git_tag": "v2.0.0 (eda5163c, clean tree)",
        "operating_point_sha256_8dp": andes_static["operating_point_sha256_8dp"],
        "stale_ybus_npy_used_for_any_verdict": False,
    },
    "pack_scripts_run_byte_identical_to_supplied": pack_copies,
    "pack_deviations": [
        "pandapower 3.4.0 has no pandapower.converter.from_ppc; bound to pandapower.converter.pypower.from_ppc "
        "(run_pack_pandapower_compat.py)",
        "canonical_ppc reads pv keys p/v/sn; canonical JSON has p0/v0/Sn; aliases added with identical values",
        "pandapower 3.4.0 from_ppc.py:303 indexes sn instead of sn_mva for zero-rated impedance branches; "
        "RATE_A=0 set to pandapower's own MAX_VAL=99999 fallback",
        "from_ppc places the tap on the HV winding when the ppc from-bus is the LV bus (branches 35, 37, 38); "
        "B-corrected reorients those rows HV->LV with tap 1/t and z*t^2 (two-port identical to 1.4e-14)",
    ],
    "tolerances": {
        "static": {
            "ybus_pu": 1e-9,
            "vm_pu": 1e-5,
            "va_rad": 1e-5,
            "gen_p_pu": 1e-8,
            "gen_q_pu": 1e-5,
            "branch_mva_preferred": 1e-3,
            "branch_mva_investigate_above": 1e-2,
        },
        "dynamic": {
            "alpha_1_per_s": 1e-4,
            "freq_hz": 1e-3,
            "band_nearest_mode_preferred": 1e-4,
            "rhp": "agree for every subset",
        },
        "theory": {"port_identity": 1e-8},
        "holdout_frozen_strong": {
            "spearman": 0.90,
            "signs": "11/12",
            "top5": "4/5",
            "median_rel": 0.10,
            "strong_opposite": 0,
            "preferred": "12/12 and spearman>=0.95",
        },
        "solver": {
            "python_pf_tol": 1e-12,
            "pandapower_tolerance_mva": 1e-10,
            "andes_pflow_tol": 1e-6,
            "internal_equilibrium_tol": 1e-9,
            "line_scaling_eps": 0.002,
        },
    },
    "holdout": json.loads((FOREST / "holdout_lines.json").read_text()),
    "inputs_sha256": {rel(p): sha(p) for p in inputs if p.exists()},
    "scripts_sha256": {rel(p): sha(p) for p in scripts},
    "outputs_sha256": {rel(p): sha(p) for p in outputs if p.exists() and p.is_file()},
}
(repo / "results/CROSS_TOOL_MANIFEST.json").write_text(json.dumps(manifest, indent=2))
print(
    json.dumps(
        {
            k: manifest[k]
            for k in (
                "git_head_at_run",
                "handoff_outer_zip",
                "environments",
                "frozen_venvs_package_lists_unchanged",
                "pack_scripts_run_byte_identical_to_supplied",
            )
        },
        indent=2,
    )
)
