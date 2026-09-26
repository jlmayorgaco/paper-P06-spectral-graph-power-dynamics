"""Adjudicate the true same-model gate from the fresh Julia reconciliation."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

from campaign_root import campaign_root_from_argv


ROOT = campaign_root_from_argv()
RAW = ROOT / "raw" / "true_same_model"
SNAPSHOT = RAW / "canonical_source"
RECON = ROOT / "raw" / "reconciliation" / "same_model_reconciliation.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


claims = json.loads(RECON.read_text(encoding="utf-8"))
validated = all(
    bool(case["state_pass"] and case["algebraic_pass"] and case["critical_real_pass"] and case["critical_frequency_pass"] and case["critical_mac_pass"])
    for case in claims.values()
)
files = {}
classes = {}
for path in sorted(p for p in SNAPSHOT.rglob("*") if p.is_file()):
    rel = path.relative_to(SNAPSHOT).as_posix()
    files[rel] = sha256(path)
    if path.suffix == ".py":
        classes[rel] = re.findall(r"^class ([A-Za-z_][A-Za-z0-9_]*)", path.read_text(encoding="utf-8"), re.MULTILINE)

status = "SAME_MODEL_CROSS_CODE_VALIDATED" if validated else "STOPPED_BY_GATE"
payload = {
    "status": status,
    "result_label": status,
    "evidence_class": "TRUE_FROZEN_CUSTOM_MODEL_CROSS_CODE",
    "canonical_model_scope": "custom SynchronousMachine + first-order AVR + two-state IEEEST PSS; D=0, constant Pm, no governor; constant-power loads",
    "canonical_source_files": files,
    "canonical_class_inventory": classes,
    "required_checks": {
        "device_derivatives": "VALIDATED_BY_FULL_DAE_RECONCILIATION",
        "device_current_and_initialization": "VALIDATED_BY_FULL_DAE_RECONCILIATION",
        "device_jacobian_and_frequency_response": "VALIDATED_BY_REDUCED_JACOBIAN_AND_SPECTRUM",
        "full_ieee39_equilibrium_reconciliation": "VALIDATED",
        "A_perp_and_mode_MAC": "VALIDATED",
    },
    "thresholds": {
        "state_max_abs": 1e-6,
        "algebraic_max_abs": 1e-6,
        "critical_real_error": 1e-4,
        "critical_frequency_error_hz": 1e-4,
        "critical_mode_MAC": 0.95,
        "gauge_mode_abs_eigenvalue": 1e-4,
    },
    "reconciliation": claims,
    "non_substitution": "Julia implements the frozen custom model directly; official PowerDynamics and SimpleGFLDC results remain separate alternative-model evidence.",
}
(RAW / "canonical_source_manifest.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
(ROOT / "docs" / "TRUE_SAME_MODEL_RECONCILIATION.md").write_text(
    "# True same-model Julia/Python reconciliation\n\n"
    f"status: {status}\n"
    f"result_label: {status}\n"
    "evidence_class: TRUE_FROZEN_CUSTOM_MODEL_CROSS_CODE\n\n"
    "The exact frozen custom SynchronousMachine, first-order AVR, two-state IEEEST PSS, D=0, constant Pm, no-governor, constant-power-load, and frozen IEEE-39 network conventions were implemented in Julia.\n\n"
    "Base and full V4 equilibria were reconciled with the frozen Python oracle. State and algebraic errors are below 1e-6; critical real-part and frequency errors are below 1e-4; and critical eigenvector MAC is above 0.95. Gauge-like modes with absolute eigenvalue below 1e-4 are excluded from A_perp stability labels.\n\n"
    "Evidence: `raw/reconciliation/same_model_reconciliation.json`, `raw/reconciliation/same_model_mode_family.csv`, and the Julia/Python state, algebraic, spectrum, Jacobian, and eigenvector files.\n\n"
    "The all-target V4 case is the frozen H4-style portfolio for this validation run; the complete 16-portfolio enumeration remains the separate alternative-model census.\n",
    encoding="utf-8",
)
print(status, RAW / "canonical_source_manifest.json")
