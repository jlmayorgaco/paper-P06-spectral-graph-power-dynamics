"""Check scientific gates, hash frozen inputs and outputs, make a portable ZIP.

The ZIP contains this experiment directory. Executing the Julia model from the
ZIP still requires the adjacent repository source and its frozen Julia project.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
ZIP = HERE / "optimal_latency_margin_20261002.zip"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def table(name: str) -> list[dict[str, str]]:
    with (HERE / name).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    frozen = json.loads((HERE / "FROZEN_PROTOCOL.json").read_text(encoding="utf-8"))
    for source in frozen["source"].values():
        target = ROOT / source["path"]
        assert target.is_file(), f"missing frozen source {target}"
        assert digest(target) == source["sha256"], f"source drift {target}"

    t01 = table("T01_DELAY_TRANSITION_REPRODUCTION.csv")
    assert len(t01) == 16 and all(r["status"] in {"SAFE", "UNSAFE"} for r in t01)
    counts = {(r["design_id"], float(r["tau_ms"])): int(r["roots_violating_margin"]) for r in t01}
    assert counts[("seed_875", 38.0)] == 0
    assert counts[("seed_875", 40.0)] == 6
    assert counts[("best_zero_delay_88455", 38.0)] == 4
    assert counts[("best_zero_delay_88455", 40.0)] == 12
    assert len(table("T01_RIGHTMOST_EXCLUSION.csv")) == 13
    assert all(r["status"] == "RIGHTMOST_WITHIN_0P001_S_INV_NUMERICALLY" for r in table("T01_RIGHTMOST_EXCLUSION.csv"))

    crossings = table("T02_PRECISE_CROSSINGS.csv")
    assert len(crossings) == 2
    assert 39.37 < float(crossings[0]["local_root_crossing_ms"]) < 39.39
    assert 37.38 < float(crossings[1]["local_root_crossing_ms"]) < 37.40
    assert all(float(r["local_bracket_width_ms"]) < 0.0001 for r in crossings)
    assert [r["status"] for r in table("T02_BOUNDARY_VERIFICATION.csv")] == ["SAFE", "UNSAFE", "SAFE", "UNSAFE"]
    fixed = table("T04_PRECISE_CROSSINGS.csv")
    assert len(fixed) == 5
    assert all(float(r["root_residual"]) < 1e-9 for r in fixed)
    assert [r["status"] for r in table("T04_ENDPOINT_BOUNDARY_VERIFICATION_0P01MS.csv")] == ["SAFE", "UNSAFE", "SAFE", "UNSAFE"]
    assert all(r["status"] == "INDETERMINATE" for r in table("T04_ENDPOINT_BOUNDARY_VERIFICATION.csv"))
    observed_change = float(fixed[-1]["fixed_gain_tau_crit_ms"]) - float(fixed[0]["fixed_gain_tau_crit_ms"])
    assert 0 < observed_change < 0.01
    assert max(float(r["relative_error"]) for r in table("T03_TAU_DERIVATIVE_VALIDATION.csv")) < 0.01
    assert math.isclose(float(table("MATERIALITY_GATE.csv")[1]["delay_loss_ms"]), -observed_change, abs_tol=1e-10)
    for filename in (
        "T05_GAIN_OPTIMIZATION_TRACE.csv",
        "T06_OPTIMAL_RHO_TAU_FRONTIER.csv",
        "T07_MAX_SPECTRAL_REPLACEMENT_VS_DELAY.csv",
        "T09_HETEROGENEOUS_DELAY_RESULTS.csv",
        "T10_SPATIAL_PREDICTOR_VALIDATION.csv",
    ):
        assert table(filename)[0]["status"] == "BLOCKED_WEAK_EFFECT_STOP_RULE"

    for filename in (
        "F01_ROOT_LOCUS_30_40MS.png",
        "F02_FIXED_GAIN_RHO_TAU_FRONTIER.png",
        "F05_FAST_MODE_ROOT_LOCUS.png",
        "F06_RETUNING_DELAY_MARGIN_RECOVERY.png",
        "README.md",
        "STATUS.md",
        "METHOD.md",
        "FINAL_REPORT.md",
        "POSTER_CLAIMS_LATENCY_FRONTIER.md",
        "REPRODUCE.md",
    ):
        assert (HERE / filename).is_file(), f"missing {filename}"
    assert not (HERE / "F03_OPTIMAL_RHO_TAU_FRONTIER.png").exists()
    assert not (HERE / "F04_MAX_REPLACEMENT_VS_DELAY.png").exists()

    files = sorted(p for p in HERE.rglob("*") if p.is_file() and p not in {ZIP, HERE / "RESULT_HASHES.json"})
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "only new targeted experiment files; source repository needed to execute Julia scripts",
        "gate_status": "WEAK_EFFECT_STOP_RULE",
        "file_count_excluding_manifest_and_zip": len(files),
        "files_sha256": {p.relative_to(HERE).as_posix(): digest(p) for p in files},
    }
    (HERE / "RESULT_HASHES.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    with ZipFile(ZIP, "w", compression=ZIP_DEFLATED, compresslevel=9) as archive:
        for path in files + [HERE / "RESULT_HASHES.json"]:
            archive.write(path, arcname=f"optimal_latency_margin_20261002/{path.relative_to(HERE).as_posix()}")
    with ZipFile(ZIP) as archive:
        assert archive.testzip() is None
        assert len(archive.namelist()) == len(files) + 1
    print(f"VALIDATED; files={len(files)+1}; zip={ZIP}; bytes={ZIP.stat().st_size}")


if __name__ == "__main__":
    main()
