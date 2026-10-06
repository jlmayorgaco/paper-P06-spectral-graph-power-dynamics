from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import qmc


ROOT = Path(__file__).resolve().parents[3]
ARTIFACT = ROOT / "artifacts" / "tx3" / "E05B_dynamic_spectral_shift"
PAPER = ROOT / "reports" / "papers" / "tx3_connected_intervention_calculus"
SEED = 20260828


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    for name in (
        "preregistration", "development", "holdout", "spectra", "contours", "tables",
        "figures", "logs", "manifests", "reports",
    ):
        (ARTIFACT / name).mkdir(parents=True, exist_ok=True)
    git_sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    source_paths = [
        ROOT / "experiments/tx3/E05B_dynamic_spectral_shift/dynamic_spectral_shift.py",
        ROOT / "experiments/tx3/E05B_dynamic_spectral_shift/run_e05b_decomposition.py",
        ROOT / "experiments/tx3/E05B_dynamic_spectral_shift/freeze_e05b_candidates.py",
        ROOT / "experiments/tx3/E05B_dynamic_spectral_shift/run_e05b_contours.py",
        ROOT / "tests/unit/test_e05b_dynamic_spectral_shift.py",
    ]
    missing = [str(path) for path in source_paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"frozen E05B sources missing: {missing}")

    samples = qmc.Sobol(d=4, scramble=True, seed=SEED).random_base2(m=4)
    lower = np.asarray([0.92, 0.95, 0.95, -15.0])
    upper = np.asarray([1.08, 1.05, 1.05, 15.0])
    scaled = qmc.scale(samples, lower, upper)
    holdout = [
        {
            "operating_point_id": f"B{index:02d}",
            "split": "e05b_holdout",
            "load_scale": float(row[0]),
            "gfl_dispatch_scale": float(row[1]),
            "reactive_load_scale": float(row[2]),
            "redispatch_mw": float(row[3]),
            "role": "independent_scrambled_sobol_holdout",
        }
        for index, row in enumerate(scaled)
    ]
    pd.DataFrame(holdout).to_parquet(ARTIFACT / "preregistration" / "E05B_HOLDOUT_MANIFEST.parquet", index=False)
    write_json(
        ARTIFACT / "preregistration" / "E05B_HOLDOUT_FREEZE.json",
        {
            "freeze_id": "TX3-E05B-HOLDOUT-1.0",
            "created_utc": datetime.now(UTC).isoformat(),
            "seed": SEED,
            "generator": "independent scrambled Sobol base-2 draw of 16 points",
            "ranges": {
                "load_scale": [0.92, 1.08],
                "gfl_dispatch_scale": [0.95, 1.05],
                "reactive_load_scale": [0.95, 1.05],
                "redispatch_mw": [-15.0, 15.0],
            },
            "operating_points": holdout,
            "E03_or_E05_holdout_reused": False,
            "outcomes_inspected_before_freeze": False,
        },
    )
    numerical = {
        "freeze_id": "TX3-E05B-NUMERICS-1.0",
        "git_sha_before_preregistration_commit": git_sha,
        "primary_frequency_band_hz": [0.1, 30.0],
        "coarse_frequency_nodes": 48,
        "reference_frequency_nodes": 96,
        "frequency_rule": "normalized Gauss-Legendre integration in log frequency",
        "sigma": 0.0,
        "finite_pole_factor": "eigenvalues of M^-1(fx-fy gy^-1 gx), evaluated independently of full descriptor logdet",
        "algebraic_factor": "sparse LU logabsdet(-gy)",
        "mass_factor": "sum log(M_ii), with M_ii strictly positive",
        "full_descriptor_reference": "frozen E03 48/96-node full logdet vertices for all 24 development operating points",
        "pointwise_factorization_audit": {
            "operating_points": ["OP00", "B00", "B15"],
            "subsets": ["EMPTY", "A3", "A6", "A3-A6", "A7", "A8", "A7-A8", "A2", "A2-A7", "A2-A8", "A2-A7-A8"],
            "frequencies_hz": [0.1, 0.3, 1.0, 3.0, 10.0, 30.0],
            "maximum_absolute_logdet_residual": 1e-7,
        },
        "contour_gap_fraction": 0.45,
        "contour_nodes": 128,
        "contour_reference_band_hz": [0.1, 30.0],
        "contour_numerical_tolerance": 1e-8,
        "source_hashes": {
            str(path.relative_to(ROOT)).replace("\\", "/"): sha256(path) for path in source_paths
        },
        "new_paraemt_runs_permitted": False,
        "E05C_or_E06_permitted": False,
    }
    gates = {
        "freeze_id": "TX3-E05B-GATES-1.0",
        "significance_multiplier_K": 5.0,
        "D1_factorization_supported_rule": {
            "all_24_development_points_and_all_pair_triple_cases_present": True,
            "maximum_pointwise_factorization_residual": 1e-7,
            "maximum_development_integrated_full_factorization_residual": 1e-7,
            "maximum_component_closure_residual": 1e-10,
        },
        "dynamic_resolution_rule": "abs(DeltaPhi_poles_96) >= 5*(abs(DeltaPhi_poles_96-DeltaPhi_poles_48) + 2^order*1e-10 + 2^order*max_initialization_residual)",
        "development_selection_rule": {
            "eligible": "resolved at >=12/24 development points and same-sign fraction >=0.75",
            "rank": "descending median absolute DeltaPhi_poles_96; select top two pairs and top two triples",
            "mandatory_diagnostics": ["A7-A8", "A3-A6", "A2-A7-A8"],
            "modal_family": "at OP00 choose the isolated positive-imaginary pole in 0.1-30 Hz with single-pole occupancy at every coalition vertex and largest absolute connected mu1",
        },
        "D2_holdout_dynamic_supported_rule": {
            "minimum_valid_smooth_coverage": 0.90,
            "minimum_distinct_pairs": 2,
            "minimum_distinct_triples": 2,
            "resolved_holdout_points_per_coalition": 8,
            "minimum_same_sign_fraction": 0.75,
        },
        "D3_contour_localization_supported_rule": {
            "minimum_reproducible_pairs": 1,
            "minimum_reproducible_triples": 1,
            "valid_single_pole_holdout_points_per_candidate": 8,
            "resolved_mu1_holdout_points_per_candidate": 8,
            "minimum_same_sign_fraction_of_mu1_real_or_imag_dominant_component": 0.75,
            "mu0_required": 0,
            "numerical_contour_error_maximum": 1e-8,
        },
        "material_damping_gate_changed": False,
        "C3a_reopened": False,
        "progression_requires_independent_review": True,
    }
    write_json(ARTIFACT / "preregistration" / "E05B_NUMERICAL_FREEZE.json", numerical)
    write_json(ARTIFACT / "preregistration" / "E05B_GATE_FREEZE.json", gates)

    prereg = f"""# E05B — Dynamic spectral-shift decomposition preregistration

**Freeze:** TX3-E05B-HOLDOUT/NUMERICS/GATES-1.0. **Git SHA before freeze:** `{git_sha}`.

E05B tests a new claim and does not rescue E05. `C3a=REJECTED`, its damping threshold
`|Delta_S zeta| >= 0.0025`, and the accepted C1/C2 decisions are immutable. E05C, E06,
TDS surgery, and new ParaEMT science are outside scope.

For every E03 pair and triple at all 24 prior operating points, and subsequently at
16 newly generated independent Sobol points, the descriptor determinant is separated as

`log|det T| = log|det(-gy)| + log|det M| + log|det(sI-Ad)|`,

where `Ad=M^-1(fx-fy gy^-1 gx)`. The three finite Mobius components are evaluated at
identical coalition vertices. Finite poles are computed independently by dense eigensolution;
the frozen E03 full-descriptor vertices and a pointwise sparse-LU audit test the factorization.

Only the 24 prior points are development data. They select two pair and two triple dynamic
candidates by the frozen reproducibility/rank rule; A7-A8, A3-A6, and A2-A7-A8 remain mandatory
diagnostics. Candidate identities and OP00 contour families are hashed and committed before
the new holdout may be evaluated.

The connected dynamic determinant is represented by
`Xi_S(s)=sum_U (-1)^(|S|-|U|) tr[(sI-Ad(U))^-1]`. Circular contours are determined solely by
the local baseline pole gap. Exact residue moments and a 128-node numerical Xi integral must
agree. D1 is exact-factor closure, D2 is reproducible finite-pole externality, and D3 is
localized connected pole motion. These are not an engineering damping-materiality gate.
"""
    (ARTIFACT / "preregistration" / "E05B_PREREGISTRATION.md").write_text(prereg, encoding="utf-8")
    theory = """# TX3 theory freeze for E05B dynamic spectral shift

**Freeze ID:** `TX3-TF-E05B-1.0`

For the index-one finite DAE operator
`T(s)=[[sM-fx,-fy],[-gx,-gy]]`, nonsingular `gy` and strictly positive diagonal
`M` give the exact Schur identity

`det T(s)=det(-gy) det(M) det(sI-Ad)`,

with `Ad=M^-1(fx-fy gy^-1 gx)`. Therefore every finite Mobius externality separates
linearly into algebraic, mass-scaling, and finite-pole terms. Algebraic and mass factors
are independent of `s`; they vanish from the logarithmic derivative.

For coalition S, define the meromorphic connected logarithmic derivative directly as
`Xi_S(s)=sum_{U subset S} (-1)^(|S|-|U|) tr[(sI-Ad(U))^-1]`. No branch of a complex
logarithm is needed. In this finite-dimensional, generally nonnormal setting the identity
follows from Jacobi's determinant formula. Perturbation-determinant/spectral-shift literature
is conceptual background, not a claim that self-adjoint operator theorems apply unchanged.

For a contour avoiding all poles,
`mu_{S,k}=(2 pi i)^-1 integral s^k Xi_S(s) ds` is the Mobius sum of enclosed pole moments.
`mu_0` is connected pole count and `mu_1` is connected enclosed pole sum. A single isolated
positive-imaginary pole is used; its conjugate carries the conjugate moment. Stable one-pole
occupancy is required at every coalition vertex. Contour results establish localized pole
motion, not material damping-margin consequence.
"""
    (PAPER / "TX3_THEORY_FREEZE_E05B.md").write_text(theory, encoding="utf-8")
    print(json.dumps({"status": "FROZEN", "holdout_points": 16, "seed": SEED}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
