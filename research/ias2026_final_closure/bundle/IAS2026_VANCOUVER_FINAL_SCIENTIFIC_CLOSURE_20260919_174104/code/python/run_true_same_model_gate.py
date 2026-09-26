"""Freeze the canonical IEEE-39 model and record the true cross-code gate."""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
from pathlib import Path

from campaign_root import campaign_root_from_argv


CAMPAIGN = campaign_root_from_argv()
RAW = CAMPAIGN / "raw" / "true_same_model"
SNAPSHOT = RAW / "canonical_source"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    files = {}
    classes = {}
    for path in sorted(p for p in SNAPSHOT.rglob("*") if p.is_file()):
        rel = path.relative_to(SNAPSHOT).as_posix()
        files[rel] = sha256(path)
        if path.suffix == ".py":
            text = path.read_text(encoding="utf-8")
            classes[rel] = re.findall(r"^class ([A-Za-z_][A-Za-z0-9_]*)", text, re.MULTILINE)
    try:
        canonical_commit = subprocess.check_output(
            ["git", "-C", r"C:\w\tx4modalaudit", "rev-parse", "HEAD"], text=True
        ).strip()
    except Exception:
        canonical_commit = "UNAVAILABLE"
    payload = {
        "status": "STOPPED_BY_GATE",
        "result_label": "STOPPED_BY_GATE",
        "evidence_class": "STOPPED_BY_GATE",
        "canonical_commit": canonical_commit,
        "canonical_model_scope": "custom SynchronousMachine + first-order AVR + two-state IEEEST PSS; D=0, constant Pm, no governor",
        "canonical_source_files": files,
        "canonical_class_inventory": classes,
        "required_checks": {
            "device_derivatives": "NOT_RUN",
            "device_current_and_initialization": "NOT_RUN",
            "device_jacobian_and_frequency_response": "NOT_RUN",
            "full_ieee39_equilibrium_reconciliation": "NOT_RUN",
            "A_perp_and_mode_MAC": "NOT_RUN",
        },
        "stop_reason": "The fresh P2/P5 runs use PowerDynamics machine/controller implementations, not the frozen custom no-governor SynchronousMachine/AVR/PSS model. No Julia port of that exact full model was completed, so the same-model cross-code gate is not promoted.",
        "non_substitution": "P2 is a fresh alternative dynamic-model negative holdout; P5 is a second-model negative holdout. Neither is evidence for TRUE_SAME_MODEL_CROSS_CODE_PASS.",
    }
    (RAW / "canonical_source_manifest.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    report = CAMPAIGN / "docs" / "TRUE_SAME_MODEL_RECONCILIATION.md"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(
        "# True same-model Julia/Python reconciliation\n\n"
        "status: STOPPED_BY_GATE\n"
        "result_label: STOPPED_BY_GATE\n"
        "evidence_class: STOPPED_BY_GATE\n\n"
        "The canonical frozen model is the custom SynchronousMachine with first-order AVR and two-state IEEEST PSS, D=0, constant mechanical power, and no governor. Its source and network snapshots are frozen under `raw/true_same_model/canonical_source/`.\n\n"
        "The fresh P2 implementation uses the official PowerDynamics IEEE-39 machines/AVR/governors around the GFL replacement, and P5 uses official SimpleGFLDC. Those are valid alternative-model negative results but are not the same-model cross-code gate. The exact Julia port required for device derivatives/current/initialization/Jacobian/frequency response, full equilibrium reconciliation, A-perp, and mode MAC was not completed; those checks remain NOT_RUN. No retuning or model substitution is used to upgrade this gate.\n\n"
        "Manifest: `raw/true_same_model/canonical_source_manifest.json`.\n",
        encoding="utf-8",
    )
    print("TRUE_SAME_MODEL_STOPPED_BY_GATE", RAW / "canonical_source_manifest.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
