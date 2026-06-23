"""
Phase-0B controller-aware rational bridge via Schur complement.

This script tests whether the full ANDES state matrix can be reduced to an
electromechanical subspace by condensing controller/exciter/PSS/IBR states:

    A_eff(s) = A_ee + A_ec (s I - A_cc)^(-1) A_ce.

The validation is non-circular. The ANDES critical pole is never used as an
evaluation point for declaring success. Instead, independent initial poles
(A_ee modes, scalar bridge modes where available, and a fixed generic grid) are
iterated/refined through the fixed-point condition eig(A_eff(s)) = s. The
converged bridge poles are then compared against the full ANDES eigensolve.

The manuscript must not be edited from this script. It writes only artifacts
and a report for user review.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "phase0b_rational_bridge"
OMEGA_S = 2 * np.pi * 60.0


def damping_ratio(s: complex) -> float:
    if abs(s) < 1e-12:
        return float("inf")
    return float(-s.real / abs(s))


def freq_hz(s: complex) -> float:
    return float(abs(s.imag) / (2 * np.pi))


def write_json(path: Path, data: Any) -> None:
    def default(obj):
        if isinstance(obj, complex):
            return {"real": float(obj.real), "imag": float(obj.imag)}
        if isinstance(obj, np.generic):
            return obj.item()
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return str(obj)

    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, default=default)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def load_andes_case(path: str):
    import andes  # type: ignore

    ss = andes.load(path, setup=False, no_output=True)
    ss.setup()
    ss.PFlow.run()
    ss.EIG.run()
    return ss


def case_paths() -> list[tuple[str, str]]:
    import andes  # type: ignore

    return [
        ("no_pss", str(ROOT / "outputs" / "ieee39_mode_diagnostics" / "ieee39_no_ieeest.xlsx")),
        ("base", str(andes.get_case("ieee39/ieee39_full.xlsx"))),
        ("mix60", str(ROOT / "validation" / "ieee39_rational_filter" / "cases" / "ieee39_ibr_mix60.xlsx")),
    ]


def classify_partition(names: list[str]) -> tuple[list[int], list[int], str]:
    """Return electromechanical and condensed-state indices.

    Primary electromechanical states are:
    - synchronous machine angle/speed: delta GENROU, omega GENROU;
    - grid-forming inverter angle state: delta REGF1.

    PLL, BusFreq, REGCP1/REECA1/REPCA1, and all other REGF1 internal states are
    treated as controller/filter states and are condensed.
    """
    e_idx = []
    for i, name in enumerate(names):
        if name.startswith("delta GENROU") or name.startswith("omega GENROU"):
            e_idx.append(i)
        elif name.startswith("delta REGF1"):
            e_idx.append(i)
    c_idx = [i for i in range(len(names)) if i not in set(e_idx)]
    rationale = (
        "x_e contains physical synchronous-machine angle/speed states "
        "(delta/omega GENROU) and the explicit grid-forming angle state "
        "(delta REGF1). PLL, BusFreq, REGCP1/REECA1/REPCA1, exciters, PSS, "
        "governors, GENROU internal voltage states, and REGF1 internal control "
        "states are condensed into x_c."
    )
    return e_idx, c_idx, rationale


def block_partition(A: np.ndarray, e_idx: list[int], c_idx: list[int]) -> dict[str, np.ndarray]:
    return {
        "Aee": A[np.ix_(e_idx, e_idx)],
        "Aec": A[np.ix_(e_idx, c_idx)],
        "Ace": A[np.ix_(c_idx, e_idx)],
        "Acc": A[np.ix_(c_idx, c_idx)],
    }


def aeff(blocks: dict[str, np.ndarray], s: complex) -> np.ndarray:
    Aee, Aec, Ace, Acc = blocks["Aee"], blocks["Aec"], blocks["Ace"], blocks["Acc"]
    n_c = Acc.shape[0]
    try:
        X = np.linalg.solve(s * np.eye(n_c, dtype=complex) - Acc, Ace)
    except np.linalg.LinAlgError:
        X = np.linalg.pinv(s * np.eye(n_c, dtype=complex) - Acc) @ Ace
    return Aee + Aec @ X


def one_sided_modes(poles: np.ndarray, f_min: float = 0.05, f_max: float = 25.0) -> list[complex]:
    modes = [
        complex(s) for s in poles
        if s.imag > 1e-7 and f_min <= freq_hz(complex(s)) <= f_max and abs(s) > 1e-9
    ]
    return sorted(modes, key=lambda s: (damping_ratio(s), freq_hz(s)))


def full_andes_modes(ss, max_modes: int = 12) -> list[dict[str, Any]]:
    mu = np.asarray(ss.EIG.mu, dtype=complex)
    pf = np.asarray(ss.EIG.pfactors)
    names = [str(x) for x in ss.EIG.x_name]
    modes = one_sided_modes(mu)
    rows = []
    for s in modes[:max_modes]:
        mode_idx = int(np.argmin(np.abs(mu - s)))
        vals = np.abs(pf[mode_idx, :]) if pf.shape[0] == len(mu) else np.abs(pf[:, mode_idx])
        top = np.argsort(vals)[-8:][::-1]
        rows.append(
            {
                "mode_index": mode_idx,
                "pole": complex(s),
                "real": float(s.real),
                "imag": float(s.imag),
                "freq_hz": freq_hz(s),
                "zeta": damping_ratio(s),
                "top_states": "; ".join(f"{names[i]}:{vals[i]:.4g}" for i in top),
            }
        )
    return rows


def scalar_bridge_initials(case_name: str) -> list[complex]:
    """Independent initial poles from the failed scalar bridge, when available."""
    path = ROOT / "outputs" / "phase0_physical_bridge" / f"{case_name}_matrices.json"
    if not path.exists():
        return []
    data = json.loads(path.read_text(encoding="utf-8"))
    M = np.asarray(data["M_eff"], dtype=float)
    D = np.diag(np.asarray(data["D_eff"], dtype=float))
    L = np.asarray(data["L"], dtype=float)
    n = len(M)
    A = np.block(
        [
            [np.zeros((n, n)), np.eye(n)],
            [-np.diag(1.0 / M) @ L, -np.diag(1.0 / M) @ D],
        ]
    )
    return one_sided_modes(np.linalg.eigvals(A), f_min=0.05, f_max=5.0)


def generic_initials() -> list[complex]:
    freqs = [0.10, 0.25, 0.50, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 5.0, 10.0, 14.0]
    zetas = [0.03, 0.15, 0.35]
    vals = []
    for f in freqs:
        w = 2 * np.pi * f
        for zeta in zetas:
            vals.append(complex(-zeta * w, w * np.sqrt(max(1 - zeta * zeta, 1e-9))))
    return vals


def initial_poles(blocks: dict[str, np.ndarray], case_name: str) -> list[complex]:
    starts = []
    starts.extend(scalar_bridge_initials(case_name))
    starts.extend(one_sided_modes(np.linalg.eigvals(blocks["Aee"]), f_min=0.05, f_max=25.0))
    starts.extend(generic_initials())
    # Deduplicate loosely.
    unique = []
    for s in starts:
        if all(abs(s - u) > 1e-3 for u in unique):
            unique.append(s)
    return unique


def fixed_point_refine(blocks: dict[str, np.ndarray], s0: complex, max_iter: int = 80) -> dict[str, Any]:
    s = complex(s0)
    history = []
    converged = False
    residual = float("inf")
    for _ in range(max_iter):
        vals = np.linalg.eigvals(aeff(blocks, s))
        lam = complex(vals[np.argmin(np.abs(vals - s))])
        residual = abs(lam - s) / max(abs(s), 1.0)
        history.append({"s": [float(s.real), float(s.imag)], "lambda": [float(lam.real), float(lam.imag)], "residual": residual})
        if residual < 1e-9:
            s = lam
            converged = True
            break
        # Relaxation avoids two-cycle behavior in strongly frequency-dependent cases.
        s = 0.45 * s + 0.55 * lam

    # Optional branch-preserving root polish. This still starts from the
    # independent fixed-point result, not from an ANDES pole.
    try:
        from scipy.optimize import root  # type: ignore

        def fun(y):
            z = complex(y[0], y[1])
            vals = np.linalg.eigvals(aeff(blocks, z))
            lam = complex(vals[np.argmin(np.abs(vals - z))])
            r = lam - z
            return [r.real, r.imag]

        sol = root(fun, [s.real, s.imag], method="hybr", options={"maxfev": 200})
        if sol.success:
            s_polished = complex(float(sol.x[0]), float(sol.x[1]))
            vals = np.linalg.eigvals(aeff(blocks, s_polished))
            lam = complex(vals[np.argmin(np.abs(vals - s_polished))])
            polished_resid = abs(lam - s_polished) / max(abs(s_polished), 1.0)
            if polished_resid < residual:
                s = s_polished
                residual = polished_resid
                converged = polished_resid < 1e-7
    except Exception:
        pass

    return {
        "initial": complex(s0),
        "pole": s,
        "converged": converged,
        "residual": float(residual),
        "iterations": len(history),
        "history": history[-8:],
    }


def cluster_poles(results: list[dict[str, Any]], tol: float = 1e-4) -> list[dict[str, Any]]:
    accepted = [
        r for r in results
        if r["converged"] and r["pole"].imag > 1e-6 and freq_hz(r["pole"]) <= 25.0
    ]
    clusters: list[list[dict[str, Any]]] = []
    for r in accepted:
        placed = False
        for cl in clusters:
            if abs(r["pole"] - cl[0]["pole"]) <= tol * max(abs(cl[0]["pole"]), 1.0):
                cl.append(r)
                placed = True
                break
        if not placed:
            clusters.append([r])
    rows = []
    for cl in clusters:
        best = min(cl, key=lambda r: r["residual"])
        p = best["pole"]
        rows.append(
            {
                "pole": p,
                "real": float(p.real),
                "imag": float(p.imag),
                "freq_hz": freq_hz(p),
                "zeta": damping_ratio(p),
                "residual": float(best["residual"]),
                "n_initials": len(cl),
                "initial_real": float(best["initial"].real),
                "initial_imag": float(best["initial"].imag),
            }
        )
    return sorted(rows, key=lambda r: (r["zeta"], r["freq_hz"]))


def match_modes(full_modes: list[dict[str, Any]], bridge_modes: list[dict[str, Any]], n: int = 5) -> list[dict[str, Any]]:
    rows = []
    for full in full_modes[:n]:
        if not bridge_modes:
            break
        best = min(bridge_modes, key=lambda b: abs(b["freq_hz"] - full["freq_hz"]))
        freq_err = abs(best["freq_hz"] - full["freq_hz"]) / max(full["freq_hz"], 1e-12)
        zeta_err = abs(best["zeta"] - full["zeta"]) / max(abs(full["zeta"]), 1e-12)
        rows.append(
            {
                "andes_mode_index": full["mode_index"],
                "andes_real": full["real"],
                "andes_imag": full["imag"],
                "andes_freq_hz": full["freq_hz"],
                "andes_zeta": full["zeta"],
                "bridge_real": best["real"],
                "bridge_imag": best["imag"],
                "bridge_freq_hz": best["freq_hz"],
                "bridge_zeta": best["zeta"],
                "relative_freq_error": float(freq_err),
                "relative_zeta_error": float(zeta_err),
                "damping_sign_correct": bool(full["zeta"] > 0 and best["zeta"] > 0),
                "bridge_residual": best["residual"],
                "andes_top_states": full["top_states"],
            }
        )
    return rows


def effective_damping_summary(blocks: dict[str, np.ndarray], e_names: list[str], s_eval: complex, ss) -> dict[str, Any]:
    Ae = aeff(blocks, s_eval)
    delta = [i for i, name in enumerate(e_names) if name.startswith("delta GENROU")]
    omega = [i for i, name in enumerate(e_names) if name.startswith("omega GENROU")]
    # Match GENROU ids appearing in both delta and omega names.
    delta_by_id = {e_names[i].split()[-1]: i for i in delta}
    omega_by_id = {e_names[i].split()[-1]: i for i in omega}
    ids = sorted(set(delta_by_id) & set(omega_by_id), key=lambda x: int(x) if x.isdigit() else x)
    dpos = [delta_by_id[i] for i in ids]
    wpos = [omega_by_id[i] for i in ids]
    if not dpos or not wpos:
        return {"available": False, "reason": "No paired GENROU delta/omega states in x_e."}

    D_block = -Ae[np.ix_(wpos, wpos)]
    L_block = -Ae[np.ix_(wpos, dpos)]
    scalar_D = np.asarray(ss.GENROU.D.v, dtype=float) / OMEGA_S if hasattr(ss, "GENROU") and hasattr(ss.GENROU, "D") else np.array([])
    return {
        "available": True,
        "eval_pole": {"real": float(s_eval.real), "imag": float(s_eval.imag)},
        "paired_genrou_ids": ids,
        "scalar_D_raw": scalar_D.tolist(),
        "scalar_D_norm": float(np.linalg.norm(scalar_D)),
        "D_eff_first_order_block_real_diag": [float(x) for x in np.real(np.diag(D_block))],
        "D_eff_first_order_block_imag_diag": [float(x) for x in np.imag(np.diag(D_block))],
        "D_eff_first_order_block_real_fro": float(np.linalg.norm(np.real(D_block))),
        "D_eff_first_order_block_imag_fro": float(np.linalg.norm(np.imag(D_block))),
        "L_eff_first_order_block_real_fro": float(np.linalg.norm(np.real(L_block))),
        "note": (
            "D_eff is reported as the local first-order omega-omega damping block "
            "-A_eff_omega,omega(s), not as a globally valid scalar nodal D."
        ),
    }


def run_case(case_name: str, path: str) -> dict[str, Any]:
    ss = load_andes_case(path)
    names = [str(x) for x in ss.EIG.x_name]
    A = np.asarray(ss.EIG.As, dtype=float)
    e_idx, c_idx, rationale = classify_partition(names)
    blocks = block_partition(A, e_idx, c_idx)
    e_names = [names[i] for i in e_idx]
    c_names = [names[i] for i in c_idx]

    starts = initial_poles(blocks, case_name)
    refined = [fixed_point_refine(blocks, s) for s in starts]
    bridge_modes = cluster_poles(refined)
    full_modes = full_andes_modes(ss, max_modes=12)
    matches = match_modes(full_modes, bridge_modes, n=5)

    critical = matches[0] if matches else None
    n_good_modes = sum(
        1 for r in matches
        if r["relative_freq_error"] <= 0.10 and r["damping_sign_correct"]
    )
    accepted = bool(
        critical
        and critical["relative_freq_error"] <= 0.05
        and critical["relative_zeta_error"] <= 0.15
        and critical["damping_sign_correct"]
        and n_good_modes >= min(3, len(matches))
    )
    eval_pole = complex(bridge_modes[0]["real"], bridge_modes[0]["imag"]) if bridge_modes else starts[0]
    damping = effective_damping_summary(blocks, e_names, eval_pole, ss)

    partition_rows = []
    for role, idxs in [("electromechanical", e_idx), ("condensed_control", c_idx)]:
        for i in idxs:
            partition_rows.append({"case": case_name, "role": role, "state_index": i, "state_name": names[i]})
    write_csv(OUT / f"{case_name}_partition.csv", partition_rows)
    write_csv(OUT / f"{case_name}_mode_validation.csv", matches)

    bridge_rows = []
    for r in bridge_modes:
        bridge_rows.append({k: v for k, v in r.items() if k != "pole"})
    write_csv(OUT / f"{case_name}_converged_bridge_poles.csv", bridge_rows)

    result = {
        "case": case_name,
        "path": path,
        "n_states": len(names),
        "n_electromechanical_states": len(e_idx),
        "n_condensed_states": len(c_idx),
        "partition_rationale": rationale,
        "electromechanical_state_names": e_names,
        "condensed_state_name_sample": c_names[:40],
        "n_initials": len(starts),
        "n_converged_bridge_poles": len(bridge_modes),
        "accepted": accepted,
        "acceptance": {
            "critical_freq_error_le_5pct": bool(critical and critical["relative_freq_error"] <= 0.05),
            "critical_zeta_error_le_15pct": bool(critical and critical["relative_zeta_error"] <= 0.15),
            "critical_damping_sign_correct": bool(critical and critical["damping_sign_correct"]),
            "n_good_modes_freq10pct_sign": n_good_modes,
        },
        "critical_match": critical,
        "full_modes": full_modes,
        "bridge_modes": bridge_modes[:12],
        "effective_damping": damping,
        "outputs": {
            "partition_csv": str((OUT / f"{case_name}_partition.csv").relative_to(ROOT)),
            "mode_validation_csv": str((OUT / f"{case_name}_mode_validation.csv").relative_to(ROOT)),
            "converged_bridge_poles_csv": str((OUT / f"{case_name}_converged_bridge_poles.csv").relative_to(ROOT)),
        },
    }
    return result


def write_report(status: dict[str, Any]) -> None:
    lines = [
        "# Phase-0B Controller-Aware Rational Bridge",
        "",
        f"Gate status: **{status['status']}**",
        "",
        "This audit uses the Schur complement `A_eff(s)=A_ee + A_ec (sI-A_cc)^(-1) A_ce`.",
        "Success is non-circular: poles are obtained by fixed-point/root refinement from independent initial poles, not by evaluating the bridge at the known ANDES pole.",
        "",
        "## Cases",
    ]
    for case in status["cases"]:
        cm = case.get("critical_match") or {}
        ed = case.get("effective_damping") or {}
        lines.extend(
            [
                f"### {case['case']}",
                f"- Accepted: {case['accepted']}",
                f"- States: x_e={case['n_electromechanical_states']}, x_c={case['n_condensed_states']}",
                f"- Converged bridge poles: {case['n_converged_bridge_poles']} from {case['n_initials']} independent starts",
            ]
        )
        if cm:
            lines.extend(
                [
                    f"- Critical ANDES pole: {cm['andes_real']:.6g} + j{cm['andes_imag']:.6g}, f={cm['andes_freq_hz']:.4g} Hz, zeta={cm['andes_zeta']:.4g}",
                    f"- Critical bridge pole: {cm['bridge_real']:.6g} + j{cm['bridge_imag']:.6g}, f={cm['bridge_freq_hz']:.4g} Hz, zeta={cm['bridge_zeta']:.4g}",
                    f"- Errors: frequency={100*cm['relative_freq_error']:.3g}%, damping={100*cm['relative_zeta_error']:.3g}%, sign={cm['damping_sign_correct']}",
                ]
            )
        if ed.get("available"):
            lines.extend(
                [
                    f"- Scalar D norm: {ed['scalar_D_norm']:.4g}",
                    f"- Controller-aware D_eff real Frobenius norm: {ed['D_eff_first_order_block_real_fro']:.4g}",
                    f"- Controller-aware D_eff imaginary Frobenius norm: {ed['D_eff_first_order_block_imag_fro']:.4g}",
                ]
            )
        lines.extend(
            [
                f"- Partition CSV: `{case['outputs']['partition_csv']}`",
                f"- Validation CSV: `{case['outputs']['mode_validation_csv']}`",
                "",
            ]
        )
    lines.extend(
        [
            "## Interpretation",
            status["interpretation"],
            "",
        ]
    )
    (OUT / "phase0b_rational_bridge_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    results = []
    for case_name, path in case_paths():
        if not Path(path).exists():
            results.append({"case": case_name, "path": path, "accepted": False, "error": "case path missing"})
            continue
        results.append(run_case(case_name, path))

    accepted = all(r.get("accepted", False) for r in results)
    status = {
        "gate": "G0.5B controller-aware rational Schur bridge",
        "status": "PASS" if accepted else "BLOCKED",
        "acceptance_rule": (
            "critical frequency error <= 5%, critical damping-ratio error <= 15%, "
            "correct damping sign, and at least three dominant modes within 10% frequency with correct sign"
        ),
        "cases": results,
        "proceed_to_phase2": bool(accepted),
        "interpretation": (
            "The controller-aware Schur bridge passed the non-circular gate and may be used to unlock Phase 2."
            if accepted
            else "The controller-aware Schur bridge did not pass the non-circular gate for all required cases. "
            "Do not advance to Phase 2-4 as full-ANDES evidence; inspect the validation tables and participation "
            "factors to identify which modes are not captured by the selected electromechanical partition."
        ),
    }
    write_json(OUT / "phase0b_rational_bridge_status.json", status)
    write_report(status)
    print(json.dumps(status, indent=2, default=str))
    return 0 if accepted else 2


if __name__ == "__main__":
    raise SystemExit(main())
