"""FC17 (step 20, last): bundle and checksum the campaign deliverables.

Destination: results/FINAL_CLOSURE.

The full run directory stays under outputs/ (its per-experiment manifest.json files
are versioned). results/FINAL_CLOSURE holds the figures with their source CSVs, the
tables every final document cites, and MANIFEST_SHA256.json covering them together
with the named deliverables (results/*.csv, docs, theory notes, configs).
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

from _fc import CONFIGS, RESULTS, ROOT, RUN, sha256, write_json

BUNDLE = RESULTS / "FINAL_CLOSURE"
TABLES = [
    "FC00_freeze_record.json",
    "FC01_transverse_quotient/FC01_structure.csv",
    "FC01_transverse_quotient/FC01_summary.json",
    "FC02_port_holdout/FC02_summary.json",
    "FC02_port_holdout/FC02_kundur_holdout.csv",
    "FC02_port_holdout/FC02_synthetic.csv",
    "FC03_governed_replication/FC03_points.csv",
    "FC03_governed_replication/FC03_summary.json",
    "FC04_e14_n6_pg/FC04_summary.json",
    "FC05_nonlinear_thresholds/FC05_thresholds.csv",
    "FC06_resilience_complex/FC06_kappa_NL.csv",
    "FC06_resilience_complex/FC06_R_k.csv",
    "FC06_resilience_complex/FC06_first_event_complex.csv",
    "FC06_resilience_complex/FC06_summary.json",
    "FC07_second_order_curvature/FC07_errors.csv",
    "FC07_second_order_curvature/FC07_slopes.csv",
    "FC07_second_order_curvature/FC07_summary.json",
    "FC07_second_order_curvature/FC07_diagnostics.csv",
    "FC08_nonlinear_certificate/FC08_summary.json",
    "FC09_bc03_certificates/FC09_summary.json",
    "FC10_monotone_and_paths/FC10_summary.json",
    "FC10_monotone_and_paths/FC10_census_lattice.csv",
    "FC12_planning/FC12_plans.csv",
    "FC12_planning/FC12_abc.csv",
    "FC12_planning/FC12_P3_support.csv",
    "FC12_planning/FC12_summary.json",
    "FC13_hopf_l1/FC13_hopf.csv",
    "FC13_hopf_l1/FC13_summary.json",
    "FC14_andes_crosscheck/FC14_andes.json",
]
NAMED = [
    RESULTS / "FINAL_EVIDENCE_TABLE.csv",
    RESULTS / "TSQ_ieee39_reaudit.csv",
    RESULTS / "zero_frequency_port_validation.csv",
    ROOT / "docs" / "FINAL_TRANSACTION_THEORY_AND_EVIDENCE.md",
    ROOT / "docs" / "FINAL_TPWRS_CLAIMS.md",
    ROOT / "docs" / "FINAL_IAS_SAFE_CLAIMS.md",
    ROOT / "docs" / "FINAL_REJECTED_WORDING.md",
    ROOT / "docs" / "FINAL_NONLINEAR_MODEL_TRACEABILITY.md",
    ROOT / "docs" / "FINAL_UNIT_AUDIT.md",
    ROOT / "theory" / "TRANSVERSE_STABILITY_QUOTIENT.md",
    ROOT / "theory" / "SYMMETRY_DEFLECTED_PORT_CLOSURE.md",
    ROOT / "theory" / "PRINCIPAL_MINOR_PORTFOLIO_STRUCTURE.md",
    ROOT / "theory" / "FINAL_COMBINATORIAL_THEOREMS.md",
    ROOT / "configs" / "port_relocation_v1.yaml",
    CONFIGS / "final_nonlinear_composability_v1.yaml",
    CONFIGS / "ieee39_governed_documented_v1.json",
    CONFIGS / "ieee39_documented_limits_v1.json",
]


def main(argv) -> int:
    (BUNDLE / "figures").mkdir(parents=True, exist_ok=True)
    copied, missing = [], []
    for rel in TABLES:
        src = RUN / rel
        if not src.exists():
            missing.append(rel)
            continue
        copied.append(shutil.copy2(src, BUNDLE / src.name))
    for src in sorted((RUN / "FC15_figures").glob("FIG*")):
        copied.append(shutil.copy2(src, BUNDLE / "figures" / src.name))
    files = [Path(p) for p in copied] + [p for p in NAMED if p.exists()]
    missing += [str(p.relative_to(ROOT)) for p in NAMED if not p.exists()]
    manifest = {
        "run": str(RUN.relative_to(ROOT)),
        "files": {
            str(p.relative_to(ROOT)).replace("\\", "/"): sha256(p) for p in files
        },
        "missing": missing,
    }
    write_json(BUNDLE / "MANIFEST_SHA256.json", manifest)
    print(len(manifest["files"]), "files;", "missing:", missing)
    return 1 if missing else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
