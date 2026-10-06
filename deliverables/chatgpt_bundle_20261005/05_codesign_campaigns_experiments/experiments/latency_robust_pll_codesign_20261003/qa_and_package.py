"""Validate the frozen provenance and create one immutable experiment archive."""

from __future__ import annotations

import csv
import hashlib
import json
import tomllib
import zipfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
ZIP = HERE.parent / (HERE.name + "_bundle.zip")


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read(name: str) -> list[dict[str,str]]:
    with (HERE / name).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def main() -> None:
    manifest = json.loads((HERE / "EXPERIMENT_MANIFEST.json").read_text())
    mismatches = []
    for rel, expected in manifest["model_sources_sha256"].items():
        path = ROOT / rel
        if not path.exists() or sha(path) != expected:
            mismatches.append(rel)
    assert not mismatches, f"model source hashes changed: {mismatches}"
    snapshot = json.loads((HERE / "POSTRUN_SOURCE_SNAPSHOT.json").read_text())
    for rel, rec in snapshot["files"].items():
        assert sha(HERE / "source_snapshot" / rel) == rec["sha256"]
        assert sha(ROOT / rel) == rec["sha256"]
    inputs = json.loads((HERE / "POSTRUN_INPUT_SNAPSHOT.json").read_text())
    for rel, rec in inputs["files"].items():
        assert sha(HERE / "input_snapshot" / rel) == rec["sha256"]
        assert sha(ROOT / rel) == rec["sha256"]
    summary = json.loads((HERE / "RESULT_SUMMARY.json").read_text())
    best = summary["best_fully_validated_found_id"]
    event = tomllib.loads((HERE / "event_validations" / best / "Q0_RESULT.toml").read_text())
    assert event["all_five_events_pass"] and len(event["events"]) == 5
    assert all(e["pass"] and e["complete"] for e in event["events"])
    assert abs(summary["fixed_GFL_percent"]-manifest["physical_replacement_percent"]) < 1e-10
    action = read("L1_ACTION_SPACE_VALIDATION.csv")
    assert len(action) >= 48 and all(int(r["numerical_difference_rank"]) == 10 for r in action)
    assert max(float(r["relative_reconstruction_error"]) for r in action) < 1e-12
    g = read("L2_TAU_MARGIN_GRADIENT_VALIDATION.csv")
    assert len(g) >= 9 and all(r["sign_agreement"] == "true" for r in g)
    trace = read("L3_OPTIMIZATION_TRACE.csv")
    accepted = [r for r in trace if r["accepted"] == "True"]
    assert best in {r["candidate_id"] for r in accepted}
    assert all(float(r["tau_crit_ms"]) <= summary["best_local_numerical_uniform_latency_threshold_ms"] + 1e-8 for r in accepted)
    counts = read(f"evaluations/{best}/ROOT_COUNTS.csv")
    assert any(r["status"] == "SAFE" for r in counts)
    assert any(r["status"] == "UNSAFE" for r in counts)
    roots = read(f"evaluations/{best}/ROOTS.csv")
    earliest = min(float(r["local_crossing_ms"]) for r in roots)
    assert abs(earliest-summary["best_local_numerical_uniform_latency_threshold_ms"]) < 1e-5
    unsafe = min((r for r in counts if r["status"]=="UNSAFE"),key=lambda r:float(r["tau_ms"]))
    assert 2*sum(float(r["local_crossing_ms"])<=float(unsafe["tau_ms"]) for r in roots)==int(unsafe["root_count"])
    required = ["README.md","STATUS.md","PROVENANCE.md","EXPERIMENT_MANIFEST.json",
                "THEORY_LATENCY_ROBUST_CODESIGN.md","POSTER_CLAIM_LEDGER.md","REPRODUCE.md",
                "FINAL_REPORT.md","PENDING_CANDIDATES.md","BEST_FOUND_DESIGN.toml","RESULT_SUMMARY.json",
                "L0_REPRODUCTION.csv","L1_ACTION_SPACE_VALIDATION.csv",
                "L2_TAU_MARGIN_GRADIENT_VALIDATION.csv","L3_OPTIMIZATION_TRACE.csv",
                "L4_MODAL_FAMILY_MAP.csv"]
    assert all((HERE / n).is_file() for n in required)
    figures = list((HERE / "figures").glob("*.png"))
    assert len(figures) >= 6 and all(p.stat().st_size > 10_000 for p in figures)
    report = dict(status="PASS_NUMERICAL_EVIDENCE_AND_ARCHIVE_COMPLETENESS",
                  best_design_id=best, frozen_source_hashes_match=True,
                  zero_delay_frozen_events_pass=5, action_checks=len(action),
                  gradient_checks=len(g), figure_count=len(figures),
                  postrun_source_snapshot_files=len(snapshot["files"]),
                  postrun_input_snapshot_files=len(inputs["files"]),
                  global_optimality_certified=False,
                  positive_delay_nonlinear_events_validated=False,
                  model_requires_referenced_repository_sources=True)
    (HERE / "QA_REPORT.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    assert not ZIP.exists(), f"Archive exists; refusing to overwrite {ZIP}"
    files = sorted(p for p in HERE.rglob("*") if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc")
    with zipfile.ZipFile(ZIP,"w",compression=zipfile.ZIP_DEFLATED,compresslevel=7) as z:
        for path in files:
            z.write(path,path.relative_to(ROOT))
    with zipfile.ZipFile(ZIP) as z:
        assert z.testzip() is None
        assert len(z.namelist()) == len(files)
    digest = sha(ZIP)
    (ZIP.parent / (ZIP.name+".sha256")).write_text(f"{digest}  {ZIP.name}\n",encoding="ascii")
    print("QA_PASS",best,"archive",ZIP,"files",len(files),"bytes",ZIP.stat().st_size,"sha256",digest)


if __name__ == "__main__":
    main()
