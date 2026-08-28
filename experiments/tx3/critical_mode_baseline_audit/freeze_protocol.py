from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
AUDIT_SOURCE = ROOT / "experiments/tx3/critical_mode_baseline_audit"
ARTIFACT = ROOT / "artifacts/tx3/critical_mode_baseline_audit"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    protocol_dir = ARTIFACT / "protocol"; protocol_dir.mkdir(parents=True, exist_ok=True)
    if (ARTIFACT / "tables/critical_mode_externalities.parquet").exists():
        raise RuntimeError("cannot freeze protocol after audit coalition outcomes exist")
    source_paths = [
        AUDIT_SOURCE / "freeze_protocol.py", AUDIT_SOURCE / "run_audit.py",
        AUDIT_SOURCE / "finalize_audit.py", AUDIT_SOURCE / "validate_audit.py",
        ROOT / "tests/unit/test_critical_mode_baseline_audit.py",
    ]
    input_paths = [
        ROOT / "systems/ieee39_tx3_gfl3/andes_case.xlsx",
        ROOT / "artifacts/tx3/E03_E04/preregistration/ACTION_SET_FREEZE.json",
        ROOT / "artifacts/tx3/E05B_dynamic_spectral_shift/tables/holdout_dynamic_reproducibility.parquet",
        ROOT / "artifacts/tx3/E05D_weak_grid/preregistration/E05D_CANDIDATE_FREEZE.json",
        ROOT / "artifacts/tx3/E05D_weak_grid/preregistration/E05D_HOLDOUT_FREEZE.json",
        ROOT / "artifacts/tx3/E05D_weak_grid/E05D_FINAL_DECISION.json",
        ROOT / "artifacts/tx3/E05C_margin_conditioned/preregistration/e05c_mode_family_freeze.parquet",
        ROOT / "experiments/tx3/E05D_weak_grid/weak_grid.py",
        ROOT / "experiments/tx3/E05C_margin_conditioned/margin_conditioned.py",
    ]
    for path in source_paths + input_paths:
        if not path.is_file(): raise FileNotFoundError(path)
    payload = {
        "freeze_id": "TX3-CRITICAL-MODE-BASELINE-AUDIT-1.0",
        "created_utc": datetime.now(UTC).isoformat(),
        "git_sha_before_freeze": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "scope": "descriptive post-mortem at kappa_grid=1 only; no new claim gate",
        "historical_statuses_immutable": {"C3a": "REJECTED", "C3b-A": "REJECTED", "C3b-B": "REJECTED", "C3c": "UNRESOLVED"},
        "population": {"E05D_frozen_seed_count": 8, "E05B_frozen_coalition_count": 15, "pair_count": 10, "triple_count": 5},
        "stress": {"kappa_grid": 1.0, "continuation": False, "new_stress_coordinate": False},
        "mode_rule": "minimum-damping valid positive-imaginary mode in 0.1-30 Hz, selected independently per seed from EMPTY only",
        "classification": {
            "participation": "abs(conj(left_k)*right_k), normalized over dynamic states",
            "descriptor_endogeneity": "conj(left_k)*M_k*right_k/(left^H M right); absolute magnitudes normalized for group comparison",
            "groups": {
                "synchronous_electromechanical": ["GENROU", "GAST", "SEXS"],
                "PLL_dominated": ["PLL2", "BusFreq"],
                "current_control_dominated": ["REGCP1", "REECB1", "REPCA1"],
            },
            "dominance_threshold": 0.60,
            "hybrid_rule": "no non-other group reaches 0.60 and the two largest groups each reach 0.20",
            "otherwise": "other",
        },
        "tracking": {
            "algorithm": "global one-to-one Hungarian biorthogonal MAC at every action-homotopy step",
            "BMAC_minimum": 0.90, "action_homotopy_nodes": [0.0, 0.25, 0.5, 0.75, 1.0],
            "refinement_node_counts": [9, 17], "manual_mode_replacement": False,
        },
        "BMAC_minimum": 0.90, "eigenpair_residual_maximum": 1e-10,
        "outputs": ["Delta_S lambda_crit", "Delta_S zeta_crit", "lambda_<q", "rho_comp", "rho_baseline"],
        "E06_authorized": False, "new_ParaEMT_authorized": False, "mechanism_surgery_authorized": False,
        "source_hashes": {str(path.relative_to(ROOT)).replace("\\", "/"): sha256(path) for path in source_paths},
        "input_hashes": {str(path.relative_to(ROOT)).replace("\\", "/"): sha256(path) for path in input_paths},
    }
    path = protocol_dir / "CRITICAL_MODE_BASELINE_AUDIT_FREEZE.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "FROZEN", "path": str(path), "source_count": len(source_paths), "input_count": len(input_paths)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
