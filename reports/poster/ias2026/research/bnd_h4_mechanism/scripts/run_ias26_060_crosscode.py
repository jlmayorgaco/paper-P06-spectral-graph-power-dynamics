"""Run the preregistered same-model Python cross-language holdout."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
RESEARCH = HERE.parents[2]
CONTRACT = HERE.parents[0].parent
sys.path.insert(0, str(RESEARCH / "src"))
sys.path.insert(0, str(RESEARCH / "experiments"))

from ibr_cycles.certification.symmetry import frequency_partner, rotation_generator  # noqa: E402
from ibr_cycles.certification.transverse import transverse_operator  # noqa: E402
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case  # noqa: E402
from ibr_cycles.models.ieee39_devices import ConverterParameters  # noqa: E402
from ibr_cycles.models.ieee39_network import load_network  # noqa: E402


TAU = 1e-8
GAUGE_TOL = 1e-3
H4 = (30, 33, 35, 37)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def scenario_network(base, row):
    active = {int(k): float(v) for k, v in json.loads(row.P_L_s_by_bus_pu_json).items()}
    reactive = {int(k): float(v) for k, v in json.loads(row.Q_L_s_by_bus_pu_json).items()}
    loads = {bus: complex(active[bus], reactive[bus]) for bus in active}
    scheduled = {int(k): float(v) for k, v in json.loads(row.P_G_scheduled_by_bus_pu_json).items()}
    pv = {bus: dict(spec) for bus, spec in base.pv.items()}
    for bus, p in scheduled.items():
        pv[bus]["p"] = p
    return replace(base, loads=loads, pv=pv)


def verdict(alpha: float) -> str:
    if not np.isfinite(alpha):
        return "INVALID"
    if alpha < -TAU:
        return "STABLE"
    if alpha > TAU:
        return "UNSTABLE"
    return "INDETERMINATE"


def solve_python(network, g: float):
    converter = ConverterParameters(voltage_control=True, voltage_gain=g, voltage_leak=0.05)
    case = solve_case(
        ReplacementPlan.of({bus: 1.0 for bus in H4}), network=network,
        converter=converter, machine_scaling={"ka": 1.425, "ta": 1.5}, tol=1e-9,
    )
    rx, _ = rotation_generator(case.dae, case.equilibrium.z)
    w = frequency_partner(case.dae).w
    transverse = transverse_operator(case.system.A, rx, w)
    transverse_values = np.linalg.eigvals(transverse.a_perp)
    keep = np.abs(transverse_values) >= GAUGE_TOL
    if not np.any(keep):
        raise RuntimeError("empty complete transverse spectrum")
    candidates = np.flatnonzero(keep)
    critical_index = candidates[np.argmax(transverse_values[candidates].real)]
    critical = complex(transverse_values[critical_index])
    alpha = float(transverse_values[candidates].real.max())
    band = transverse_values[keep & (transverse_values.imag >= -1e-10)]
    band_f = np.abs(band.imag) / (2 * np.pi)
    band = band[(band_f >= 0.3 - 1e-10) & (band_f <= 1.5 + 1e-10) & (np.abs(band) > 1e-3)]
    alpha_omega = float(band.real.max()) if band.size else math.nan
    # The full-state realization is needed for a direct cross-code eigenvector MAC.
    full_values, full_vectors = np.linalg.eig(case.system.A)
    full_keep = np.flatnonzero(np.abs(full_values) >= GAUGE_TOL)
    full_critical_idx = full_keep[np.argmax(full_values[full_keep].real)]
    if abs(float(full_values[full_critical_idx].real) - alpha) > 1e-7:
        raise RuntimeError("full-matrix and quotient alpha disagree beyond 1e-7")
    return case, alpha, critical, alpha_omega, full_values, full_vectors, int(full_critical_idx)


def normalized_mac(a: np.ndarray, b: np.ndarray) -> float:
    den = float(np.vdot(a, a).real * np.vdot(b, b).real)
    return float(abs(np.vdot(a, b)) ** 2 / den) if den > 0 else math.nan


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    args = parser.parse_args()
    run = args.run_root.resolve()
    cfg = json.loads((run / "config/IAS26-060_MC_OPERATING_V1.json").read_text(encoding="utf-8"))
    expected_hashes = cfg["source_provenance"]["sha256"]
    repo = Path(__file__).resolve().parents[6]
    for rel in (
        "reports/poster/ias2026/research/configs/ias2026/ieee39_network.json",
        "reports/poster/ias2026/research/src/ibr_cycles/models/ieee39_network.py",
        "reports/poster/ias2026/research/src/ibr_cycles/models/ieee39_case.py",
        "reports/poster/ias2026/research/src/ibr_cycles/models/ieee39_devices.py",
    ):
        p = repo / rel
        expected = expected_hashes.get(rel)
        if expected and sha256(p) != expected:
            raise RuntimeError(f"frozen Python source hash changed: {rel}")

    selected = pd.read_csv(run / "tables/CROSSCODE_HOLDOUT_SELECTION.csv")
    selected = selected[selected.selection_status != "EMPTY_STRATUM_NO_REPLACEMENT"]
    if selected.scenario_id.nunique() != 50:
        raise RuntimeError(f"expected 50 unique selected IDs, got {selected.scenario_id.nunique()}")
    case_table = pd.read_csv(run / "derived/MC_CASES.csv")
    manifest = pd.read_csv(run / "inputs/IAS26-060_SCENARIOS_V1.csv").set_index("scenario_id")
    vectors = pd.read_csv(run / "raw/crosscode/JULIA_CRITICAL_EIGENVECTORS.csv")
    base = load_network()
    rows = []
    for n, sel in enumerate(selected.itertuples(index=False), 1):
        sid = sel.scenario_id
        row = manifest.loc[sid]
        network = scenario_network(base, row)
        for treatment, g in (("ORIGINAL", 0.03625), ("INTERVENTION", 0.25)):
            jul = case_table[(case_table.scenario_id == sid) & (case_table.portfolio_mask == 15) & (case_table.treatment == treatment)]
            if len(jul) != 1 or jul.iloc[0].status != "PASS":
                raise RuntimeError(f"missing Julia reference for {sid}/{treatment}")
            jr = jul.iloc[0]
            case, py_alpha, py_lam, py_omega, py_vals, py_full_vecs, pycrit_idx = solve_python(network, g)
            jvec_rows = vectors[(vectors.scenario_id == sid) & (vectors.treatment == treatment)].sort_values("vector_index")
            jvec = jvec_rows.vector_real.to_numpy() + 1j * jvec_rows.vector_imag.to_numpy()
            jlam = complex(float(jvec_rows.lambda_real.iloc[0]), float(jvec_rows.lambda_imag.iloc[0]))
            py_vec = py_full_vecs[:, pycrit_idx]
            if jvec.size != py_vec.size:
                raise RuntimeError(f"state-order dimension mismatch {sid}: Julia {jvec.size} Python {py_vec.size}")
            # Compare the Julia critical family against every Python full-system eigenvector,
            # allowing conjugate representatives. The coordinate order is the frozen TX4 order.
            py_keep = np.flatnonzero(np.abs(py_vals) >= GAUGE_TOL)
            macs = np.array([max(normalized_mac(jvec, py_full_vecs[:, j]), normalized_mac(np.conj(jvec), py_full_vecs[:, j])) for j in py_keep])
            mac_order = np.argsort(macs)[::-1]
            best_full_idx = int(py_keep[mac_order[0]])
            best_mac = float(macs[mac_order[0]])
            margin = float(best_mac - macs[mac_order[1]]) if macs.size > 1 else math.nan
            nearest_idx = int(np.argmin(np.abs(py_vals[py_keep] - jlam)))
            nearest_gap = float(abs(py_vals[py_keep[nearest_idx]] - jlam))
            alpha_gap = float(abs(py_alpha - float(jr.alpha_perp)))
            freq_gap = float(abs(abs(py_lam.imag) / (2*np.pi) - float(jr.critical_frequency_hz)))
            resolved = best_mac >= 0.70 and margin >= 0.10
            status = "RESOLVED_SAME_FAMILY" if resolved and best_full_idx == pycrit_idx else ("RESOLVED_DIFFERENT_DOMINANT_FAMILY" if resolved else "IDENTITY_UNRESOLVED")
            rows.append({
                "scenario_id": sid, "stratum": sel.stratum, "rank_within_stratum": sel.rank_within_stratum,
                "treatment": treatment, "g": g, "status": "PASS",
                "julia_alpha_perp": float(jr.alpha_perp), "python_alpha_perp": py_alpha,
                "abs_delta_alpha": alpha_gap, "julia_critical_frequency_hz": float(jr.critical_frequency_hz),
                "python_critical_frequency_hz": abs(py_lam.imag)/(2*np.pi), "abs_delta_frequency_hz": freq_gap,
                "julia_verdict": verdict(float(jr.alpha_perp)), "python_verdict": verdict(py_alpha),
                "verdict_agrees": verdict(float(jr.alpha_perp)) == verdict(py_alpha),
                "critical_root_julia_real": jlam.real, "critical_root_julia_imag": jlam.imag,
                "critical_root_python_real": py_lam.real, "critical_root_python_imag": py_lam.imag,
                "nearest_root_distance": nearest_gap, "critical_mode_mac": best_mac,
                "mac_score_margin": margin, "mode_identity_status": status,
                "state_count": case.system.n, "python_pf_residual": case.dae.power_flow.max_mismatch,
                "python_dae_f_residual": float(np.linalg.norm(case.dae.f(case.equilibrium.x, case.equilibrium.z, {}), np.inf)),
                "python_dae_g_residual": float(np.linalg.norm(case.dae.g(case.equilibrium.x, case.equilibrium.z, {}), np.inf)),
                "python_alpha_omega": py_omega, "julia_alpha_omega": float(jr.alpha_omega) if pd.notna(jr.alpha_omega) else math.nan,
            })
        print(f"IAS26_060_CROSSCODE_COMPLETE scenarios={n}/50 id={sid}", flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(run / "tables/CROSSCODE_JULIA_PYTHON_HOLDOUT.csv", index=False)
    usable = out[out.status == "PASS"]
    summary = {
        "N_scenarios": int(out.scenario_id.nunique()), "N_case_comparisons": int(len(out)),
        "max_abs_delta_alpha": float(usable.abs_delta_alpha.max()),
        "median_abs_delta_alpha": float(usable.abs_delta_alpha.median()),
        "max_abs_delta_f_hz": float(usable.abs_delta_frequency_hz.max()),
        "verdict_agreement_count": int(usable.verdict_agrees.sum()),
        "verdict_agreement_rate": float(usable.verdict_agrees.mean()),
        "mode_identity_resolved_count": int(usable.mode_identity_status.str.startswith("RESOLVED").sum()),
        "mode_identity_status_counts": usable.mode_identity_status.value_counts().to_dict(),
        "crosscode_scope": "same-model cross-language reproduction; not independent physical validation",
        "powerdynamics_exact_model_available": False,
        "thresholds_unchanged": True,
    }
    (run / "derived/CROSSCODE_SUMMARY.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
