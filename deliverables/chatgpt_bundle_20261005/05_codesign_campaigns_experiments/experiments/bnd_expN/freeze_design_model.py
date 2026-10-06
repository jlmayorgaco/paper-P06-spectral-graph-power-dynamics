"""Gate and freeze the PD-free ExpN design model before Z optimization."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports" / "experiment_N"


def rows(name: str) -> list[dict[str, str]]:
    with (OUT / name).open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise SystemExit("MODEL_FREEZE_REJECTED: " + reason)


def main() -> None:
    trim = rows("TABLE_N04_trim_validation.csv")
    identity = rows("TABLE_N05_parametric_identity.csv")
    ports = rows("TABLE_N06_withheld_port_identity.csv")
    poles = rows("TABLE_N07_withheld_pole_identity.csv")
    modes = rows("TABLE_N03_near_zero_mode_classification.csv")
    gradients = rows("TABLE_N08_derivative_validation.csv")
    require(len(trim) == 60 and all(r["status"] == "PASS" for r in trim),
            "equilibrium trim grid")
    require(len(identity) == 69 and all(r["matrix_gate"] == "PASS" for r in identity),
            "withheld reduced matrices")
    require(len(ports) >= 60 and all(r["port_gate"] == "PASS" for r in ports),
            "withheld GFL ports")
    require(len(poles) == 69 and all(r["pole_gate"] == "PASS" for r in poles),
            "withheld finite poles")
    require(len(modes) >= 70 and all(r["classification"] != "UNRESOLVED" for r in modes),
            "near-zero classification")
    require(len(gradients) == 90 and all(r["status"] == "PASS" for r in gradients),
            "total derivatives")

    paths = [
        "src/bnd_model_expN/PDExactDesignN.jl",
        "src/bnd_model_expN/PDReferenceN.jl",
        "src/bnd_design_e/CollectiveModel.jl",
        "src/bnd_design/AnalyticSG.jl",
        "src/bnd_design/AnalyticGFLPLL.jl",
        "src/pd39/PD39.jl",
        "src/pd39/model.jl",
        "src/pd39/equilibrium.jl",
        "reports/experiment_N/TABLE_N01_original_operating_point.csv",
        "reports/experiment_N/SOFTWARE_PROVENANCE.toml",
        "reports/experiment_N/TABLE_N02_rho_semantics.csv",
        "reports/experiment_N/TABLE_N03_near_zero_mode_classification.csv",
        "reports/experiment_N/TABLE_N04_trim_validation.csv",
        "reports/experiment_N/TABLE_N05_parametric_identity.csv",
        "reports/experiment_N/TABLE_N06_withheld_port_identity.csv",
        "reports/experiment_N/TABLE_N07_withheld_pole_identity.csv",
        "reports/experiment_N/TABLE_N08_derivative_validation.csv",
    ]
    for p in paths:
        require((ROOT / p).is_file(), f"missing {p}")
    hashes = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}
    manifest = {"algorithm": "SHA-256", "inputs": hashes,
                "gates": {"trim": len(trim), "identity": len(identity),
                          "ports": len(ports), "poles": len(poles),
                          "near_zero": len(modes), "gradients": len(gradients)}}
    serialized = json.dumps(manifest, sort_keys=True, separators=(",", ":")).encode()
    manifest["MODEL_SHA"] = hashlib.sha256(serialized).hexdigest()
    path = OUT / "MODEL_FREEZE.json"
    if path.exists():
        require(json.loads(path.read_text(encoding="utf-8")) == manifest,
                "existing freeze differs")
    else:
        path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
    print("MODEL_STAGE_PASSED YES")
    print("MODEL_SHA", manifest["MODEL_SHA"])


if __name__ == "__main__":
    main()
