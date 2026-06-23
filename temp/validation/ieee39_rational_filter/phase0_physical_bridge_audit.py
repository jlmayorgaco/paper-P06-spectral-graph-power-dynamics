"""
Phase-0 physical ANDES -> (L, M, D) bridge audit.

This script implements the strict gate requested by the validation roadmap:
extract a physically justified scalar swing surrogate from an ANDES IEEE39
linearization and compare its electromechanical modes against the full ANDES
eigensolve. It intentionally does not fit damping to make the surrogate look
good. D is read from the machine damping field where available.

If the extracted surrogate does not meet the modal matching tolerance, the
output records a BLOCKED gate and the dependent estimator phases should not be
run as full-ANDES evidence.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from run_ieee39_surrogate_experiments import (
    build_laplacian_from_andes,
    damping_ratio,
    kron_reduce,
    qep_poles,
)


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "phase0_physical_bridge"
OMEGA_S = 2 * np.pi * 60.0


def write_json(path: Path, data) -> None:
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def load_case(path: str):
    import andes  # type: ignore

    ss = andes.load(path, setup=False, no_output=True)
    ss.setup()
    ss.PFlow.run()
    ss.EIG.run()
    return ss


def full_andes_modes(ss, f_min: float = 0.1, f_max: float = 3.0, zeta_max: float = 0.8) -> list[dict]:
    """Return one side of stable oscillatory full-ANDES modes in the EM band."""
    mu = np.asarray(ss.EIG.mu, dtype=complex)
    pf = np.asarray(ss.EIG.pfactors)
    names = [str(x) for x in ss.EIG.x_name]
    em_idx = [
        i for i, name in enumerate(names)
        if name.startswith("delta GENROU") or name.startswith("omega GENROU")
    ]

    rows = []
    for mode_idx, s in enumerate(mu):
        freq = abs(s.imag) / (2 * np.pi)
        if not (s.real < 0 and s.imag > 1e-7 and f_min <= freq <= f_max):
            continue
        zeta = damping_ratio(s)
        if zeta > zeta_max:
            continue
        if pf.shape[0] == len(mu):
            em_participation = float(np.sum(np.abs(pf[mode_idx, em_idx])))
            top_vals = np.abs(pf[mode_idx, :])
        else:
            em_participation = float(np.sum(np.abs(pf[em_idx, mode_idx])))
            top_vals = np.abs(pf[:, mode_idx])
        top = np.argsort(top_vals)[-5:][::-1]
        rows.append(
            {
                "mode_index": int(mode_idx),
                "real": float(s.real),
                "imag": float(s.imag),
                "freq_hz": float(freq),
                "zeta": float(zeta),
                "em_participation": em_participation,
                "top_states": "; ".join(f"{names[i]}:{top_vals[i]:.4g}" for i in top),
            }
        )
    return sorted(rows, key=lambda r: r["zeta"])


def extracted_surrogate(ss) -> dict:
    """Build the physically extracted scalar swing surrogate."""
    Lbus, bus_ids, _ = build_laplacian_from_andes(ss)
    gen_buses = [int(x) for x in ss.GENROU.bus.v]
    bus_pos = {bus: i for i, bus in enumerate(bus_ids)}
    keep = [bus_pos[b] for b in gen_buses]
    Lred = kron_reduce(Lbus, keep)
    Lred = (Lred + Lred.T) / 2
    Lred -= np.diag(np.sum(Lred, axis=1))
    Lred = (Lred + Lred.T) / 2

    # ANDES uses delta_dot = omega_s * speed_deviation. For a delta-coordinate
    # QEP M_eff delta_ddot + D_eff delta_dot + L delta = 0, this gives
    # M_eff = M / omega_s and D_eff = D / omega_s.
    M_raw = np.asarray(ss.GENROU.M.v, dtype=float)
    D_raw = np.asarray(ss.GENROU.D.v, dtype=float) if hasattr(ss.GENROU, "D") else np.zeros_like(M_raw)
    M_eff = M_raw / OMEGA_S
    D_eff = D_raw / OMEGA_S
    return {
        "bus_ids": bus_ids,
        "generator_buses": gen_buses,
        "L": Lred,
        "M_eff": M_eff,
        "D_eff": D_eff,
        "M_raw": M_raw,
        "D_raw": D_raw,
    }


def surrogate_modes(surr: dict) -> list[dict]:
    poles = qep_poles(surr["M_eff"], np.diag(surr["D_eff"]), surr["L"])
    rows = []
    for s in poles:
        if s.imag <= 1e-7:
            continue
        rows.append(
            {
                "real": float(s.real),
                "imag": float(s.imag),
                "freq_hz": float(abs(s.imag) / (2 * np.pi)),
                "zeta": float(damping_ratio(s)),
            }
        )
    return sorted(rows, key=lambda r: r["freq_hz"])


def match_modes(full_modes: list[dict], surrogate: list[dict], n_modes: int = 6) -> list[dict]:
    rows = []
    for full in full_modes[:n_modes]:
        if not surrogate:
            break
        best = min(surrogate, key=lambda s: abs(s["freq_hz"] - full["freq_hz"]))
        freq_error = abs(best["freq_hz"] - full["freq_hz"]) / max(full["freq_hz"], 1e-12)
        zeta_error = abs(best["zeta"] - full["zeta"]) / max(full["zeta"], 1e-12)
        rows.append(
            {
                "andes_mode_index": full["mode_index"],
                "andes_freq_hz": full["freq_hz"],
                "surrogate_freq_hz": best["freq_hz"],
                "relative_freq_error": float(freq_error),
                "andes_zeta": full["zeta"],
                "surrogate_zeta": best["zeta"],
                "relative_zeta_error": float(zeta_error),
                "damping_sign_correct": bool(best["zeta"] > 1e-5 and full["zeta"] > 0),
                "em_participation": full["em_participation"],
                "andes_top_states": full["top_states"],
            }
        )
    return rows


def audit_one(case_name: str, path: str) -> dict:
    ss = load_case(path)
    full = full_andes_modes(ss)
    surr = extracted_surrogate(ss)
    smodes = surrogate_modes(surr)
    matches = match_modes(full, smodes)

    critical = matches[0] if matches else None
    critical_freq_pass = bool(critical and critical["relative_freq_error"] <= 0.05)
    critical_damp_pass = bool(critical and critical["damping_sign_correct"])
    first_six_freq_pass = bool(matches and np.mean([r["relative_freq_error"] for r in matches]) <= 0.05)
    accepted = bool(critical_freq_pass and critical_damp_pass and first_six_freq_pass)

    write_csv(OUT / f"{case_name}_mode_match.csv", matches)
    write_json(
        OUT / f"{case_name}_matrices.json",
        {
            "generator_buses": surr["generator_buses"],
            "M_raw": surr["M_raw"].tolist(),
            "D_raw": surr["D_raw"].tolist(),
            "M_eff": surr["M_eff"].tolist(),
            "D_eff": surr["D_eff"].tolist(),
            "L": surr["L"].tolist(),
        },
    )

    return {
        "case": case_name,
        "path": path,
        "n_full_modes_in_band": len(full),
        "n_surrogate_modes": len(smodes),
        "critical_freq_pass_5pct": critical_freq_pass,
        "critical_damping_sign_pass": critical_damp_pass,
        "mean_first_six_freq_pass_5pct": first_six_freq_pass,
        "accepted": accepted,
        "blocking_reason": (
            "PASS"
            if accepted
            else "Physical L and M place several frequencies in the electromechanical band, "
            "but GENROU.D is zero in the packaged case, so the extracted scalar D cannot "
            "reproduce the positive damping of the full ANDES modes. Controller/exciter/PSS "
            "states provide damping that is not captured by a scalar nodal D without a "
            "validated reduction."
        ),
        "critical_match": critical,
        "full_modes": full[:10],
        "surrogate_modes": smodes,
        "outputs": {
            "mode_match_csv": str((OUT / f"{case_name}_mode_match.csv").relative_to(ROOT)),
            "matrices_json": str((OUT / f"{case_name}_matrices.json").relative_to(ROOT)),
        },
    }


def main() -> int:
    import andes  # type: ignore

    OUT.mkdir(parents=True, exist_ok=True)
    base = str(andes.get_case("ieee39/ieee39_full.xlsx"))
    no_pss = ROOT / "outputs" / "ieee39_mode_diagnostics" / "ieee39_no_ieeest.xlsx"
    cases = [("base", base)]
    if no_pss.exists():
        cases.append(("no_pss", str(no_pss)))

    results = [audit_one(name, path) for name, path in cases]
    accepted = all(r["accepted"] for r in results)
    status = {
        "gate": "G0.5 physical ANDES-to-(L,M,D) bridge",
        "status": "PASS" if accepted else "BLOCKED",
        "acceptance_rule": (
            "critical mode frequency error <= 5%, first-six mean frequency error <= 5%, "
            "and positive damping sign reproduced by physically extracted scalar D"
        ),
        "results": results,
        "proceed_to_phase2": bool(accepted),
    }
    write_json(OUT / "phase0_physical_bridge_status.json", status)

    lines = [
        "# Phase-0 Physical ANDES-to-(L,M,D) Bridge Audit",
        "",
        f"Gate status: **{status['status']}**",
        "",
        f"Acceptance rule: {status['acceptance_rule']}.",
        "",
    ]
    for r in results:
        lines.extend(
            [
                f"## {r['case']}",
                f"- Accepted: {r['accepted']}",
                f"- Full modes in 0.1--3 Hz band: {r['n_full_modes_in_band']}",
                f"- Surrogate modes: {r['n_surrogate_modes']}",
                f"- Blocking reason: {r['blocking_reason']}",
                f"- Mode match CSV: `{r['outputs']['mode_match_csv']}`",
                f"- Matrices JSON: `{r['outputs']['matrices_json']}`",
                "",
            ]
        )
    if not accepted:
        lines.append(
            "Dependent full-ANDES estimator validation phases must remain blocked. "
            "The result is a useful negative audit: the scalar network surrogate needs "
            "a validated controller-state reduction before it can be used as full-ANDES ground-truth evidence."
        )
    (OUT / "phase0_physical_bridge_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(json.dumps(status, indent=2))
    return 0 if accepted else 2


if __name__ == "__main__":
    raise SystemExit(main())
