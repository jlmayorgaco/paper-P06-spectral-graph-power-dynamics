from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
E03_E04 = ROOT / "artifacts" / "tx3" / "E03_E04"
ARTIFACT = ROOT / "artifacts" / "tx3" / "E05_mechanism_consequence"
PAPER = ROOT / "reports" / "papers" / "tx3_connected_intervention_calculus"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def candidate_statistics(frame: pd.DataFrame, order: int) -> list[dict[str, Any]]:
    records = []
    selected = frame[(frame["split"] == "holdout") & (frame["order"] == order) & frame["resolved_flag"]]
    for coalition, group in selected.groupby("coalition_key"):
        same_sign = max((group["direct_externality"] > 0).mean(), (group["direct_externality"] < 0).mean())
        if len(group) < 4 or same_sign < 0.75:
            continue
        records.append(
            {
                "coalition_key": coalition,
                "actions": coalition.split("-"),
                "resolved_holdout_points": len(group),
                "same_sign_fraction": float(same_sign),
                "median_externality": float(group["direct_externality"].median()),
                "median_absolute_externality": float(group["direct_externality"].abs().median()),
                "median_eta_feedback": float(group["eta_feedback"].median()),
            }
        )
    return records


def main() -> int:
    for name in (
        "preregistration",
        "candidates",
        "discovery",
        "holdout",
        "density",
        "figures",
        "tables",
        "logs",
        "manifests",
        "reports",
    ):
        (ARTIFACT / name).mkdir(parents=True, exist_ok=True)
    pair = pd.read_parquet(E03_E04 / "tables" / "pair_externalities.parquet")
    triple = pd.read_parquet(E03_E04 / "tables" / "triple_externalities.parquet")
    combined = pd.concat([pair, triple], ignore_index=True)
    pair_candidates = candidate_statistics(combined, 2)
    triple_candidates = candidate_statistics(combined, 3)
    curvature = max(
        (record for record in pair_candidates if record["median_eta_feedback"] <= 0.5),
        key=lambda record: record["median_absolute_externality"],
    )
    feedback = max(
        (record for record in pair_candidates if record["median_eta_feedback"] >= 0.9),
        key=lambda record: record["median_absolute_externality"],
    )
    third_order = max(triple_candidates, key=lambda record: record["median_absolute_externality"])
    candidates = []
    for candidate_id, label, selection_rule, record in (
        (
            "C_R",
            "curvature_dominated",
            "largest absolute median reproducible holdout pair externality with median eta_feedback <= 0.5",
            curvature,
        ),
        (
            "C_F",
            "feedback_dominated",
            "largest absolute median reproducible holdout pair externality with median eta_feedback >= 0.9",
            feedback,
        ),
        (
            "C_3",
            "genuine_third_order",
            "largest absolute median reproducible holdout triple externality",
            third_order,
        ),
    ):
        candidates.append({"candidate_id": candidate_id, "class": label, "selection_rule": selection_rule, **record})
    git_sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    candidate_payload = {
        "freeze_id": "TX3-E05-CANDIDATES-1.0",
        "created_utc": datetime.now(UTC).isoformat(),
        "source_stage": "accepted E03/E04 holdout",
        "selection_population_rule": "resolved on >=4/8 holdout points with same-sign fraction >=0.75",
        "selection_used_e05_consequence_results": False,
        "candidates": candidates,
    }
    numerical = {
        "freeze_id": "TX3-E05-NUMERICS-1.0",
        "git_sha_before_e05": git_sha,
        "reference_operating_point": "OP00",
        "reference_mode_band_hz": [0.1, 30.0],
        "tracking_target_band_hz": [0.05, 45.0],
        "tracking_method": "one-to-one Hungarian assignment maximizing scale-invariant left/right biorthogonal MAC",
        "minimum_tracking_bmac": 0.70,
        "discovery_mode_family_minimum_valid_points": 12,
        "holdout_minimum_valid_points": 6,
        "density_frequency_band_hz": [0.1, 30.0],
        "density_frequency_nodes": 256,
        "density_grid": "uniform in log frequency including endpoints",
        "density_integral": "trapezoid in log frequency divided by log(omega_max/omega_min)",
        "sigma": 0.0,
        "modal_reduction": "exact algebraic Schur complement followed by descriptor-mass scaling",
        "primary_consequence_metric": "finite Mobius externality of the BMAC-tracked modal damping ratio",
        "secondary_metrics": [
            "tracked-mode real-part externality per second",
            "tracked-mode frequency externality Hz",
            "tracked-mode log pole-resolvent externality",
            "joint endpoint changes",
            "system spectral-abscissa externality",
        ],
        "source_hashes": {
            str(path.relative_to(ROOT)).replace("\\", "/"): sha256(path)
            for path in (
                ROOT / "experiments" / "tx3" / "E05_discovery_holdout" / "mechanism_consequence.py",
                ROOT / "experiments" / "tx3" / "E05_discovery_holdout" / "run_e05_campaign.py",
                ROOT / "tests" / "unit" / "test_e05_mode_tracking.py",
            )
        },
        "holdout_consequences_inspected_before_freeze": False,
        "new_paraemt_runs_permitted": False,
        "E06_or_later_permitted": False,
    }
    gate = {
        "freeze_id": "TX3-E05-GATE-1.0",
        "candidate_supported_rule": {
            "valid_holdout_mode_tracking_points_minimum": 6,
            "absolute_damping_ratio_externality_material_threshold": 0.0025,
            "material_holdout_points_minimum": 4,
            "minimum_same_sign_fraction": 0.75,
        },
        "E05_supported_rule": "at least one of C_R, C_F, or C_3 satisfies the candidate-supported rule",
        "E05_unresolved_rule": "fewer than 6/8 valid tracked holdout points for every candidate",
        "cycle_or_SCC_claim_assessed": False,
        "harmful_label_rule": "harmful only when damping externality is negative or real-part externality is positive; DeltaPhi sign alone is insufficient",
    }
    write_json(ARTIFACT / "preregistration" / "E05_CANDIDATE_FREEZE.json", candidate_payload)
    write_json(ARTIFACT / "preregistration" / "E05_NUMERICAL_FREEZE.json", numerical)
    write_json(ARTIFACT / "preregistration" / "E05_GATE_FREEZE.json", gate)
    prereg = f"""# TX3 E05 mechanism--consequence preregistration

**Freeze IDs:** TX3-E05-CANDIDATES-1.0 / TX3-E05-NUMERICS-1.0 / TX3-E05-GATE-1.0  
**Parent evidence:** C1=SUPPORTED and C2=SUPPORTED.  
**Scope:** E05 only. No ParaEMT, TDS surgery, negative control, repair, E06, or later stage is permitted.

E05 does not search for a dangerous cycle. It maps three objectively selected, reproducible E03 coalitions to externality-density spectra and BMAC-tracked small-signal modal consequences. Candidate selection is frozen from accepted E03/E04 holdout statistics before any E05 consequence is evaluated: `{curvature['coalition_key']}` is curvature-dominated, `{feedback['coalition_key']}` is feedback-dominated, and `{third_order['coalition_key']}` is the largest reproducible third-order coalition.

The OP00 positive-imaginary oscillatory modes in 0.1--30 Hz define reference families. After discovery only, each candidate locks the reference family having the largest median absolute damping-ratio Mobius externality among families with valid BMAC tracking on at least 12/16 discovery points. Holdout remains untouched until `E05_MODE_FAMILY_FREEZE.json` is created and committed.

A candidate is dynamically consequential only if its locked family is valid at at least 6/8 holdout points, at least 4 holdout points satisfy `|Delta_S zeta| >= 0.0025`, and the same-sign fraction is at least 0.75. The sign of the E03 logdet externality never labels an action harmful. A negative damping externality or positive real-part externality is required for a harmful interpretation.

The original cycle/SCC C3 claim remains unassessed. E05 instead tests the narrower preregistered claim `C3a`: at least one accepted finite spectral externality maps reproducibly to a material tracked-mode consequence. Negative results are retained and stop progression.
"""
    (ARTIFACT / "preregistration" / "E05_PREREGISTRATION.md").write_text(prereg, encoding="utf-8")
    theory = """# TX3 theory freeze for E05 mechanism--consequence mapping

**Freeze ID:** `TX3-TF-E05-1.0`  
**Parent freezes:** `TX3-TF-0.2` and `TX3-TF-C1-C2-1.0`.

E05 preserves the C1/C2 action coordinates, fully re-equilibrated descriptor operator, primary logdet functional, frequency band, and direct finite Mobius convention. The frequency-resolved externality density is the signed Mobius sum of `log|det T(j omega)|` at the registered coalition vertices. Its normalized log-frequency integral must recover the E03 direct externality within the registered numerical comparison tolerance.

Dynamic modes are obtained from the exact algebraic Schur complement `A = M^{-1}(f_x-f_y g_y^{-1}g_x)`. Positive-imaginary OP00 modes in 0.1--30 Hz define immutable reference families. Families are mapped across operating points and action vertices by one-to-one Hungarian assignment maximizing scale-invariant left/right biorthogonal MAC. No `min(zeta)` mode switching is used.

For each tracked modal scalar `q`, its finite consequence externality uses the same Mobius operator as C1: `Delta_S q = sum_U (-1)^(|S|-|U|) q(1_U)`. The primary scalar is damping ratio; real part, frequency, log pole-resolvent factor, joint endpoint changes, and system spectral abscissa are secondary. A large negative logdet externality is not defined as harmful. Harm requires a negative damping-ratio interaction or positive real-part interaction on the locked family.

E05 tests `C3a`, not the original cycle-localization C3: at least one objectively selected E03 coalition has a reproducible material tracked-mode consequence under the frozen holdout gate. No quiver cycle, SCC, surgery target, correction, causality, TDS consequence, or EMT transfer is claimed.
"""
    (PAPER / "TX3_THEORY_FREEZE_E05.md").write_text(theory, encoding="utf-8")
    print(json.dumps(candidate_payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
