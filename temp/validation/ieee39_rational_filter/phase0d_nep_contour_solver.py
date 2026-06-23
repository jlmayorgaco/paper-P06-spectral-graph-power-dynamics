"""
Phase-0D rational NEP validation by Beyn contour integration.

This script replaces the Phase-0B fixed-point Schur bridge with a genuine
contour-integral eigensolver.  The reference remains the full ANDES eigensolve.

The exact object solved here is the first-order Schur nonlinear eigenproblem

    S_e(s) x_e = 0,
    S_e(s) = s I - A_ee - A_ec (s I - A_cc)^(-1) A_ce,

where x_e contains the retained network/interface states and x_c contains all
controller, exciter, governor, PLL, and auxiliary states.  For cases with clean
delta/omega GENROU pairs, the report also extracts the local omega-omega
damping block used by Phase-0B as the Sigma/D_eff diagnostic.

No manuscript files are edited by this script.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "phase0d_nep_contour_solver"
OMEGA_S = 2 * np.pi * 60.0
RNG_SEED = 20260608


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

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, default=default)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
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
    """Retain the same interface states as Phase-0B.

    Keeping PLL and REGF1 internal states in x_c is intentional: Phase-0D tests
    whether a contour solver can recover controller-family zeros through the
    rational Schur term instead of inserting those states into x_e.
    """
    e_idx = []
    for i, name in enumerate(names):
        if name.startswith("delta GENROU") or name.startswith("omega GENROU"):
            e_idx.append(i)
        elif name.startswith("delta REGF1"):
            e_idx.append(i)
    e_set = set(e_idx)
    c_idx = [i for i in range(len(names)) if i not in e_set]
    rationale = (
        "x_e matches Phase-0B: delta/omega GENROU states plus explicit "
        "grid-forming angle states (delta REGF1). PLL, BusFreq, REGCP1/REECA1/"
        "REPCA1, exciters, PSS, governors, GENROU internal voltage states, and "
        "REGF1 internal controls are kept in x_c and enter through the rational "
        "Schur term."
    )
    return e_idx, c_idx, rationale


def block_partition(A: np.ndarray, e_idx: list[int], c_idx: list[int]) -> dict[str, np.ndarray]:
    return {
        "Aee": A[np.ix_(e_idx, e_idx)].astype(complex),
        "Aec": A[np.ix_(e_idx, c_idx)].astype(complex),
        "Ace": A[np.ix_(c_idx, e_idx)].astype(complex),
        "Acc": A[np.ix_(c_idx, c_idx)].astype(complex),
    }


def schur_nep(blocks: dict[str, np.ndarray], s: complex) -> np.ndarray:
    Aee, Aec, Ace, Acc = blocks["Aee"], blocks["Aec"], blocks["Ace"], blocks["Acc"]
    n_c = Acc.shape[0]
    R_Ace = np.linalg.solve(s * np.eye(n_c, dtype=complex) - Acc, Ace)
    return s * np.eye(Aee.shape[0], dtype=complex) - Aee - Aec @ R_Ace


def schur_nep_derivative(blocks: dict[str, np.ndarray], s: complex) -> np.ndarray:
    Aee, Aec, Ace, Acc = blocks["Aee"], blocks["Aec"], blocks["Ace"], blocks["Acc"]
    n_c = Acc.shape[0]
    eye_c = np.eye(n_c, dtype=complex)
    R_Ace = np.linalg.solve(s * eye_c - Acc, Ace)
    R2_Ace = np.linalg.solve(s * eye_c - Acc, R_Ace)
    return np.eye(Aee.shape[0], dtype=complex) + Aec @ R2_Ace


def schur_aeff(blocks: dict[str, np.ndarray], s: complex) -> np.ndarray:
    """A_eff(s) used only for the paired delta/omega damping diagnostic."""
    Aee, Aec, Ace, Acc = blocks["Aee"], blocks["Aec"], blocks["Ace"], blocks["Acc"]
    n_c = Acc.shape[0]
    mat = s * np.eye(n_c, dtype=complex) - Acc
    try:
        R_Ace = np.linalg.solve(mat, Ace)
    except np.linalg.LinAlgError:
        # Only used for the static Sigma(0)-style diagnostic.  The actual
        # contour solver never accepts a pseudoinverse result.
        R_Ace = np.linalg.pinv(mat) @ Ace
    return Aee + Aec @ R_Ace


def one_sided_modes(poles: np.ndarray, f_min: float = 0.03, f_max: float = 25.0) -> list[complex]:
    modes = [
        complex(s)
        for s in poles
        if s.imag > 1e-7 and f_min <= freq_hz(complex(s)) <= f_max and abs(s) > 1e-9
    ]
    return sorted(modes, key=lambda s: (damping_ratio(s), freq_hz(s)))


def full_andes_modes(ss, max_modes: int = 16) -> list[dict[str, Any]]:
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
                "pole": s,
                "real": float(s.real),
                "imag": float(s.imag),
                "freq_hz": freq_hz(s),
                "zeta": damping_ratio(s),
                "top_states": "; ".join(f"{names[i]}:{vals[i]:.4g}" for i in top),
            }
        )
    return rows


@dataclass(frozen=True)
class RectangleContour:
    name: str
    real_min: float
    real_max: float
    imag_min: float
    imag_max: float
    nodes_per_side: int

    def contains(self, z: complex, margin: float = 1e-9) -> bool:
        return (
            self.real_min - margin <= z.real <= self.real_max + margin
            and self.imag_min - margin <= z.imag <= self.imag_max + margin
        )

    def nodes_weights(self) -> tuple[np.ndarray, np.ndarray]:
        """Midpoint trapezoid nodes on a counter-clockwise rectangle."""
        n = self.nodes_per_side
        pts: list[complex] = []
        dzs: list[complex] = []

        dx = (self.real_max - self.real_min) / n
        dy = (self.imag_max - self.imag_min) / n

        for k in range(n):
            pts.append(complex(self.real_min + (k + 0.5) * dx, self.imag_min))
            dzs.append(complex(dx, 0.0))
        for k in range(n):
            pts.append(complex(self.real_max, self.imag_min + (k + 0.5) * dy))
            dzs.append(complex(0.0, dy))
        for k in range(n):
            pts.append(complex(self.real_max - (k + 0.5) * dx, self.imag_max))
            dzs.append(complex(-dx, 0.0))
        for k in range(n):
            pts.append(complex(self.real_min, self.imag_max - (k + 0.5) * dy))
            dzs.append(complex(0.0, -dy))

        return np.asarray(pts, dtype=complex), np.asarray(dzs, dtype=complex)

    def descriptor(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "real_min": self.real_min,
            "real_max": self.real_max,
            "imag_min": self.imag_min,
            "imag_max": self.imag_max,
            "nodes_per_side": self.nodes_per_side,
            "quadrature_nodes": 4 * self.nodes_per_side,
            "frequency_band_hz": [self.imag_min / (2 * np.pi), self.imag_max / (2 * np.pi)],
        }


def default_contours() -> list[RectangleContour]:
    return [
        # The network-margin contour is intentionally not the entire
        # low-frequency half-plane.  A very wide low-frequency rectangle has
        # more zeros/poles than the retained interface dimension can extract
        # with the two-moment Beyn form.  Splitting the search by physical
        # family makes the argument-count check well conditioned and still
        # covers the critical damping modes in the three validation cases.
        RectangleContour("network_margin", -4.0, 0.8, 3.5, 12.0, 128),
        RectangleContour("pll_control", -25.0, 0.8, 60.0, 95.0, 128),
    ]


def beyn_contour(
    blocks: dict[str, np.ndarray],
    contour: RectangleContour,
    m_probe: int,
    rank_tol: float = 1e-9,
    residual_tol: float = 1e-6,
    seed: int = RNG_SEED,
) -> dict[str, Any]:
    n = blocks["Aee"].shape[0]
    rng = np.random.default_rng(seed)
    Vhat = rng.standard_normal((n, m_probe)) + 1j * rng.standard_normal((n, m_probe))
    nodes, dzs = contour.nodes_weights()

    A0 = np.zeros((n, m_probe), dtype=complex)
    A1 = np.zeros((n, m_probe), dtype=complex)
    trace_count = 0.0 + 0.0j
    failures = 0

    for z, dz in zip(nodes, dzs):
        try:
            Tz = schur_nep(blocks, z)
            X = np.linalg.solve(Tz, Vhat)
            A0 += X * dz
            A1 += z * X * dz

            Tp = schur_nep_derivative(blocks, z)
            Y = np.linalg.solve(Tz, Tp)
            trace_count += np.trace(Y) * dz
        except np.linalg.LinAlgError:
            failures += 1

    factor = 1.0 / (2j * np.pi)
    A0 *= factor
    A1 *= factor
    trace_count *= factor

    U, sigma, Vh = np.linalg.svd(A0, full_matrices=False)
    if sigma.size == 0:
        return {"zeros": [], "sigma": [], "rank": 0, "trace_count": trace_count, "failures": failures}

    threshold = max(sigma[0] * rank_tol, 1e-13)
    K = int(np.sum(sigma > threshold))
    # A contour with weakly projected control modes can create a soft tail.
    # Keep at least one singular direction if the argument count indicates
    # content, but never exceed the probe dimension.
    arg_count_est = int(round(float(trace_count.real)))
    if arg_count_est > 0:
        K = max(K, min(arg_count_est, len(sigma), n))
    K = min(K, len(sigma), n)

    if K <= 0:
        return {
            "zeros": [],
            "sigma": [float(x) for x in sigma],
            "rank": 0,
            "trace_count": trace_count,
            "failures": failures,
        }

    U_K = U[:, :K]
    V_K = Vh.conj().T[:, :K]
    S_inv = np.diag(1.0 / sigma[:K])
    B_red = U_K.conj().T @ A1 @ V_K @ S_inv
    eigvals, eigvecs = np.linalg.eig(B_red)

    zeros = []
    for j, lam0 in enumerate(eigvals):
        lam = polish_nep_zero(blocks, complex(lam0), contour)
        x = null_vector(blocks, lam)
        try:
            Tlam = schur_nep(blocks, lam)
            res = np.linalg.norm(Tlam @ x) / max(np.linalg.norm(x), 1e-15)
            rel_res = res / max(np.linalg.norm(Tlam, ord="fro"), 1e-15)
        except np.linalg.LinAlgError:
            res = float("inf")
            rel_res = float("inf")
        accepted = bool(
            contour.contains(complex(lam), margin=1e-6)
            and lam.imag > 1e-7
            and rel_res <= residual_tol
        )
        zeros.append(
            {
                "pole": complex(lam),
                "real": float(lam.real),
                "imag": float(lam.imag),
                "freq_hz": freq_hz(complex(lam)),
                "zeta": damping_ratio(complex(lam)),
                "abs_residual": float(res),
                "relative_residual": float(rel_res),
                "accepted": accepted,
                "x_e": x,
            }
        )

    return {
        "zeros": zeros,
        "sigma": [float(x) for x in sigma],
        "rank": K,
        "trace_count": trace_count,
        "failures": failures,
    }


def null_vector(blocks: dict[str, np.ndarray], s: complex) -> np.ndarray:
    """Right singular vector associated with the smallest singular value."""
    T = schur_nep(blocks, s)
    _, _, vh = np.linalg.svd(T)
    x = vh.conj().T[:, -1]
    norm = np.linalg.norm(x)
    return x / norm if norm > 0 else x


def polish_nep_zero(blocks: dict[str, np.ndarray], s0: complex, contour: RectangleContour) -> complex:
    """Local root polish for a contour-discovered NEP zero.

    This is not a fixed-point iteration and it does not use ANDES poles.  It
    starts only from the Beyn Ritz value and solves for the smallest eigenvalue
    of S_e(s) to vanish.
    """
    try:
        from scipy.optimize import root  # type: ignore
    except Exception:
        return s0

    def fun(y):
        z = complex(float(y[0]), float(y[1]))
        try:
            vals = np.linalg.eigvals(schur_nep(blocks, z))
        except np.linalg.LinAlgError:
            return [1e6, 1e6]
        val = complex(vals[int(np.argmin(np.abs(vals)))])
        return [float(val.real), float(val.imag)]

    try:
        sol = root(fun, [s0.real, s0.imag], method="hybr", options={"maxfev": 120})
    except Exception:
        return s0

    if not sol.success:
        return s0
    z = complex(float(sol.x[0]), float(sol.x[1]))
    # Keep the polish local to the contour-discovered feature.  This prevents a
    # bad Ritz value from jumping to an unrelated low-frequency mode.
    contour_span = max(contour.real_max - contour.real_min, contour.imag_max - contour.imag_min)
    if abs(z - s0) > 0.25 * contour_span:
        return s0
    return z


def cluster_zeros(zeros: list[dict[str, Any]], tol: float = 5e-5) -> list[dict[str, Any]]:
    accepted = [z for z in zeros if z.get("accepted") and z["pole"].imag > 1e-7]
    accepted = sorted(accepted, key=lambda z: z["relative_residual"])
    clusters: list[list[dict[str, Any]]] = []
    for z in accepted:
        placed = False
        for cl in clusters:
            if abs(z["pole"] - cl[0]["pole"]) <= tol * max(abs(cl[0]["pole"]), 1.0):
                cl.append(z)
                placed = True
                break
        if not placed:
            clusters.append([z])

    rows = []
    for cl in clusters:
        best = min(cl, key=lambda z: z["relative_residual"])
        row = dict(best)
        row["n_contour_duplicates"] = len(cl)
        rows.append(row)
    return sorted(rows, key=lambda z: (z["zeta"], z["freq_hz"]))


def control_participation(blocks: dict[str, np.ndarray], s: complex, x_e: np.ndarray) -> tuple[float, float, np.ndarray]:
    Acc, Ace = blocks["Acc"], blocks["Ace"]
    n_c = Acc.shape[0]
    try:
        x_c = np.linalg.solve(s * np.eye(n_c, dtype=complex) - Acc, Ace @ x_e)
    except np.linalg.LinAlgError:
        x_c = np.linalg.pinv(s * np.eye(n_c, dtype=complex) - Acc) @ (Ace @ x_e)
    norm_e = float(np.linalg.norm(x_e))
    norm_c = float(np.linalg.norm(x_c))
    raw_denom = max(norm_e + norm_c, 1e-15)
    rms_e = norm_e / math.sqrt(max(len(x_e), 1))
    rms_c = norm_c / math.sqrt(max(len(x_c), 1))
    rms_denom = max(rms_e + rms_c, 1e-15)
    return rms_c / rms_denom, norm_c / raw_denom, x_c


def annotate_zeros(
    zeros: list[dict[str, Any]],
    blocks: dict[str, np.ndarray],
    c_names: list[str],
    tau: float = 0.80,
) -> list[dict[str, Any]]:
    annotated = []
    acc_eig = np.linalg.eigvals(blocks["Acc"])
    for z in zeros:
        pi_c, pi_c_raw, x_c = control_participation(blocks, z["pole"], z["x_e"])
        family = "II_control" if pi_c >= tau else "I_network"
        nearest_acc = complex(acc_eig[int(np.argmin(np.abs(acc_eig - z["pole"])))])
        vals = np.abs(x_c)
        top = np.argsort(vals)[-6:][::-1] if len(vals) else []
        top_states = "; ".join(f"{c_names[i]}:{vals[i]:.4g}" for i in top)
        row = {
            k: v
            for k, v in z.items()
            if k not in {"x_e", "accepted"}
        }
        row.update(
            {
                "pi_c": float(pi_c),
                "pi_c_raw_norm": float(pi_c_raw),
                "participation_note": "pi_c is RMS-normalized by block size; pi_c_raw_norm is the unscaled norm ratio.",
                "family": family,
                "nearest_acc_real": float(nearest_acc.real),
                "nearest_acc_imag": float(nearest_acc.imag),
                "nearest_acc_freq_hz": freq_hz(nearest_acc),
                "distance_to_nearest_acc": float(abs(nearest_acc - z["pole"])),
                "top_condensed_states": top_states,
            }
        )
        annotated.append(row)
    return sorted(annotated, key=lambda z: (z["zeta"], z["freq_hz"]))


def match_critical(full_modes: list[dict[str, Any]], zeros: list[dict[str, Any]]) -> dict[str, Any] | None:
    if not full_modes or not zeros:
        return None
    full = full_modes[0]
    s_ref = complex(full["real"], full["imag"])
    best = min(zeros, key=lambda z: abs(z["pole"] - s_ref))
    freq_err = abs(best["freq_hz"] - full["freq_hz"]) / max(full["freq_hz"], 1e-12)
    zeta_err = abs(best["zeta"] - full["zeta"]) / max(abs(full["zeta"]), 1e-12)
    return {
        "andes_mode_index": full["mode_index"],
        "andes_real": full["real"],
        "andes_imag": full["imag"],
        "andes_freq_hz": full["freq_hz"],
        "andes_zeta": full["zeta"],
        "andes_top_states": full["top_states"],
        "nep_real": best["real"],
        "nep_imag": best["imag"],
        "nep_freq_hz": best["freq_hz"],
        "nep_zeta": best["zeta"],
        "nep_family": best["family"],
        "nep_pi_c": best["pi_c"],
        "nep_relative_residual": best["relative_residual"],
        "relative_freq_error": float(freq_err),
        "relative_zeta_error": float(zeta_err),
        "damping_sign_correct": bool(full["zeta"] > 0 and best["zeta"] > 0),
        "pole_distance": float(abs(best["pole"] - s_ref)),
        "nearest_acc_real": best["nearest_acc_real"],
        "nearest_acc_imag": best["nearest_acc_imag"],
        "nearest_acc_freq_hz": best["nearest_acc_freq_hz"],
        "distance_to_nearest_acc": best["distance_to_nearest_acc"],
        "top_condensed_states": best["top_condensed_states"],
    }


def paired_delta_omega_positions(e_names: list[str]) -> tuple[list[int], list[int], list[str]]:
    delta = [i for i, name in enumerate(e_names) if name.startswith("delta GENROU")]
    omega = [i for i, name in enumerate(e_names) if name.startswith("omega GENROU")]
    delta_by_id = {e_names[i].split()[-1]: i for i in delta}
    omega_by_id = {e_names[i].split()[-1]: i for i in omega}
    ids = sorted(set(delta_by_id) & set(omega_by_id), key=lambda x: int(x) if x.isdigit() else x)
    return [delta_by_id[i] for i in ids], [omega_by_id[i] for i in ids], ids


def effective_damping_summary(blocks: dict[str, np.ndarray], e_names: list[str], s_eval: complex, ss) -> dict[str, Any]:
    dpos, wpos, ids = paired_delta_omega_positions(e_names)
    if not dpos or not wpos:
        return {"available": False, "reason": "No paired GENROU delta/omega states in x_e."}

    Acc = blocks["Acc"]
    mat = s_eval * np.eye(Acc.shape[0], dtype=complex) - Acc
    solve_mode = "solve"
    cond_est = float("inf")
    try:
        cond_est = float(np.linalg.cond(mat))
        np.linalg.solve(mat, blocks["Ace"])
    except np.linalg.LinAlgError:
        solve_mode = "pinv"
    Ae = schur_aeff(blocks, s_eval)
    D_block = -Ae[np.ix_(wpos, wpos)]
    L_block = -Ae[np.ix_(wpos, dpos)]
    scalar_D = np.asarray(ss.GENROU.D.v, dtype=float) / OMEGA_S if hasattr(ss, "GENROU") and hasattr(ss.GENROU, "D") else np.array([])
    return {
        "available": True,
        "eval_pole": {"real": float(s_eval.real), "imag": float(s_eval.imag)},
        "cond_sI_minus_Acc": cond_est,
        "linear_solve_mode": solve_mode,
        "paired_genrou_ids": ids,
        "scalar_D_norm": float(np.linalg.norm(scalar_D)),
        "D_eff_real_fro": float(np.linalg.norm(np.real(D_block))),
        "D_eff_imag_fro": float(np.linalg.norm(np.imag(D_block))),
        "D_eff_real_diag": [float(x) for x in np.real(np.diag(D_block))],
        "D_eff_imag_diag": [float(x) for x in np.imag(np.diag(D_block))],
        "L_eff_real_fro": float(np.linalg.norm(np.real(L_block))),
        "note": (
            "D_eff is the local omega-omega block -A_eff_ww(s_eval). "
            "At s=0 this is the Sigma(0)-style static Schur damping diagnostic; "
            "at the critical pole it matches the Phase-0B reported diagnostic."
        ),
    }


def load_phase0b_reference() -> dict[str, Any]:
    path = ROOT / "outputs" / "phase0b_rational_bridge" / "phase0b_rational_bridge_status.json"
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def phase0b_case_reference(phase0b: dict[str, Any], case_name: str) -> dict[str, Any] | None:
    for case in phase0b.get("cases", []):
        if case.get("case") == case_name:
            return case
    return None


def run_case(case_name: str, path: str, phase0b: dict[str, Any]) -> dict[str, Any]:
    ss = load_andes_case(path)
    names = [str(x) for x in ss.EIG.x_name]
    A = np.asarray(ss.EIG.As, dtype=float)
    e_idx, c_idx, rationale = classify_partition(names)
    blocks = block_partition(A, e_idx, c_idx)
    e_names = [names[i] for i in e_idx]
    c_names = [names[i] for i in c_idx]

    full_modes = full_andes_modes(ss, max_modes=16)
    full_eigs = np.asarray(ss.EIG.mu, dtype=complex)
    contours = default_contours()

    raw_zeros = []
    contour_rows = []
    acc_eigs = np.linalg.eigvals(blocks["Acc"])
    for ci, contour in enumerate(contours):
        result = beyn_contour(
            blocks,
            contour,
            m_probe=min(blocks["Aee"].shape[0], max(12, blocks["Aee"].shape[0])),
            rank_tol=1e-10,
            residual_tol=1e-5,
            seed=RNG_SEED + ci,
        )
        poles_in_contour = sum(1 for z in acc_eigs if contour.contains(complex(z)) and z.imag > 0)
        andes_in_contour = sum(1 for z in full_eigs if contour.contains(complex(z)) and z.imag > 0)
        accepted_count = sum(1 for z in result["zeros"] if z.get("accepted"))
        zero_count_estimate = float(np.real(result["trace_count"])) + float(poles_in_contour)
        count_error = abs(float(accepted_count) - zero_count_estimate)
        contour_rows.append(
            {
                "case": case_name,
                **contour.descriptor(),
                "beyn_rank": result["rank"],
                "accepted_zeros": accepted_count,
                "trace_count_real": float(np.real(result["trace_count"])),
                "trace_count_imag": float(np.imag(result["trace_count"])),
                "trace_count_rounded": int(round(float(np.real(result["trace_count"])))),
                "acc_poles_in_contour_upper": int(poles_in_contour),
                "andes_poles_in_contour_upper": int(andes_in_contour),
                "argument_count_plus_acc_poles": zero_count_estimate,
                "count_error_abs": count_error,
                "count_check_pass": bool(count_error <= 0.10),
                "linear_solve_failures": result["failures"],
                "leading_singular_values": "; ".join(f"{x:.3e}" for x in result["sigma"][:8]),
            }
        )
        for z in result["zeros"]:
            z["contour"] = contour.name
        raw_zeros.extend(result["zeros"])

    clustered = cluster_zeros(raw_zeros)
    annotated = annotate_zeros(clustered, blocks, c_names, tau=0.70)
    critical = match_critical(full_modes, annotated)

    count_checks_pass = all(row["count_check_pass"] for row in contour_rows)
    accepted = bool(
        critical
        and critical["relative_freq_error"] <= 0.05
        and critical["relative_zeta_error"] <= 0.10
        and critical["damping_sign_correct"]
        and count_checks_pass
    )

    # Sigma/D_eff diagnostics at s=0 and at the matched critical pole.
    damping_zero = effective_damping_summary(blocks, e_names, 0.0 + 0.0j, ss)
    eval_pole = complex(critical["nep_real"], critical["nep_imag"]) if critical else 0.0 + 0.0j
    damping_critical = effective_damping_summary(blocks, e_names, eval_pole, ss)
    ref_case = phase0b_case_reference(phase0b, case_name)
    if ref_case and ref_case.get("effective_damping", {}).get("available"):
        ed_ref = ref_case["effective_damping"]
        damping_critical["phase0b_D_eff_real_fro"] = ed_ref.get("D_eff_first_order_block_real_fro")
        damping_critical["phase0b_D_eff_imag_fro"] = ed_ref.get("D_eff_first_order_block_imag_fro")
        damping_critical["phase0b_eval_pole"] = ed_ref.get("eval_pole")

    partition_rows = []
    for role, idxs in [("retained_interface", e_idx), ("condensed_control", c_idx)]:
        for i in idxs:
            partition_rows.append({"case": case_name, "role": role, "state_index": i, "state_name": names[i]})
    write_csv(OUT / f"{case_name}_partition.csv", partition_rows)

    zero_rows = []
    for z in annotated:
        zero_rows.append({k: v for k, v in z.items() if k != "pole"})
    write_csv(OUT / f"{case_name}_zeros.csv", zero_rows)
    write_csv(OUT / f"{case_name}_contour_counts.csv", contour_rows)
    if critical:
        write_csv(OUT / f"{case_name}_critical_validation.csv", [critical])

    return {
        "case": case_name,
        "path": path,
        "n_states": len(names),
        "n_retained_interface_states": len(e_idx),
        "n_condensed_states": len(c_idx),
        "partition_rationale": rationale,
        "retained_state_names": e_names,
        "condensed_state_name_sample": c_names[:50],
        "contours": contour_rows,
        "n_beyn_zeros": len(annotated),
        "accepted": accepted,
        "count_checks_pass": bool(count_checks_pass),
        "critical_match": critical,
        "full_modes": full_modes,
        "beyn_zeros": annotated,
        "effective_damping": {
            "s0_static": damping_zero,
            "critical_pole": damping_critical,
        },
        "outputs": {
            "partition_csv": str((OUT / f"{case_name}_partition.csv").relative_to(ROOT)),
            "zeros_csv": str((OUT / f"{case_name}_zeros.csv").relative_to(ROOT)),
            "contour_counts_csv": str((OUT / f"{case_name}_contour_counts.csv").relative_to(ROOT)),
            "critical_validation_csv": str((OUT / f"{case_name}_critical_validation.csv").relative_to(ROOT)),
        },
    }


def plot_results(status: dict[str, Any]) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return

    n_cases = len(status["cases"])
    fig, axes = plt.subplots(1, n_cases, figsize=(5.2 * n_cases, 4.6), squeeze=False)
    contours = default_contours()
    for ax, case in zip(axes[0], status["cases"]):
        full = case.get("full_modes", [])
        zeros = case.get("beyn_zeros", [])
        if full:
            ax.scatter(
                [m["real"] for m in full],
                [m["imag"] for m in full],
                marker="x",
                color="black",
                label="ANDES",
                s=36,
            )
        if zeros:
            net = [z for z in zeros if z["family"] == "I_network"]
            ctrl = [z for z in zeros if z["family"] == "II_control"]
            ax.scatter([z["real"] for z in net], [z["imag"] for z in net], color="#1f77b4", label="Beyn Family I", s=30)
            ax.scatter([z["real"] for z in ctrl], [z["imag"] for z in ctrl], color="#d62728", label="Beyn Family II", s=30)

        # Acc poles are summarized, not all plotted, to keep the figure readable.
        for contour in contours:
            xs = [contour.real_min, contour.real_max, contour.real_max, contour.real_min, contour.real_min]
            ys = [contour.imag_min, contour.imag_min, contour.imag_max, contour.imag_max, contour.imag_min]
            ax.plot(xs, ys, color="#777777", linewidth=0.8, linestyle="--")
        ax.axvline(0, color="#999999", linewidth=0.8)
        ax.set_title(f"{case['case']} ({case['status_label']})")
        ax.set_xlabel("Re(s)")
        ax.set_ylabel("Im(s) [rad/s]")
        ax.set_ylim(0, 100)
        ax.grid(True, alpha=0.25)
    handles, labels = axes[0][0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=4, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "fig_phase0d_s_plane.png", dpi=220)
    fig.savefig(OUT / "fig_phase0d_s_plane.pdf")
    plt.close(fig)


def write_report(status: dict[str, Any]) -> None:
    lines = [
        "# Phase-0D Rational NEP Contour Solver",
        "",
        f"Gate status: **{status['status']}**",
        "",
        "This audit solves the exact first-order Schur nonlinear eigenproblem",
        "`S_e(s)=sI-A_ee-A_ec(sI-A_cc)^(-1)A_ce` with Beyn contour integration.",
        "No fixed-point iteration is used for acceptance.  The full ANDES eigensolve",
        "is the reference.",
        "",
        "## Non-Circularity Audit",
        "- Beyn zeros are first extracted from the contour moments `A0` and `A1` and the reduced eigenproblem.",
        "- The local polish starts only from each Beyn Ritz value and solves for the smallest eigenvalue of `S_e(s)` to vanish.",
        "- The polish function receives no ANDES pole, mode index, frequency, or damping target. ANDES eigenvalues are used only after extraction for validation and count diagnostics.",
        "- Therefore a zero can be produced from the Jacobian/Schur NEP even if the corresponding ANDES pole is not supplied to the solver.",
        "",
        "## Contours",
    ]
    for contour in default_contours():
        d = contour.descriptor()
        lines.append(
            f"- `{d['name']}`: Re=[{d['real_min']}, {d['real_max']}], "
            f"Im=[{d['imag_min']}, {d['imag_max']}] rad/s "
            f"({d['frequency_band_hz'][0]:.3g}-{d['frequency_band_hz'][1]:.3g} Hz), "
            f"{d['quadrature_nodes']} quadrature nodes."
        )
    lines.extend(["", "## Required Gate Checks"])
    for case in status["cases"]:
        cm = case.get("critical_match") or {}
        ed0 = case.get("effective_damping", {}).get("s0_static", {})
        edc = case.get("effective_damping", {}).get("critical_pole", {})
        lines.extend(
            [
                f"### {case['case']}",
                f"- Accepted: **{case['accepted']}**",
                f"- Argument-count checks pass: **{case['count_checks_pass']}**",
                f"- States: retained={case['n_retained_interface_states']}, condensed={case['n_condensed_states']}",
                f"- Beyn zeros accepted: {case['n_beyn_zeros']}",
            ]
        )
        if ed0.get("available"):
            lines.append(
                f"- Sigma/D_eff at s=0: scalar D norm={ed0['scalar_D_norm']:.4g}, "
                f"real Fro={ed0['D_eff_real_fro']:.4g}, imag Fro={ed0['D_eff_imag_fro']:.4g}, "
                f"solve={ed0['linear_solve_mode']}, cond={ed0['cond_sI_minus_Acc']:.3g}."
            )
        if edc.get("available"):
            ref = ""
            if "phase0b_D_eff_real_fro" in edc:
                ref = (
                    f" Phase-0B critical-pole reference: real Fro="
                    f"{edc['phase0b_D_eff_real_fro']:.4g}, imag Fro={edc['phase0b_D_eff_imag_fro']:.4g}."
                )
            lines.append(
                f"- Sigma/D_eff at matched critical pole: real Fro={edc['D_eff_real_fro']:.4g}, "
                f"imag Fro={edc['D_eff_imag_fro']:.4g}.{ref}"
            )
        if cm:
            lines.extend(
                [
                    f"- ANDES critical pole: {cm['andes_real']:.6g} + j{cm['andes_imag']:.6g}, "
                    f"f={cm['andes_freq_hz']:.5g} Hz, zeta={cm['andes_zeta']:.5g}.",
                    f"- NEP-Beyn matched pole: {cm['nep_real']:.6g} + j{cm['nep_imag']:.6g}, "
                    f"f={cm['nep_freq_hz']:.5g} Hz, zeta={cm['nep_zeta']:.5g}, "
                    f"family={cm['nep_family']}, pi_c={cm['nep_pi_c']:.3f}.",
                    f"- Critical errors: frequency={100*cm['relative_freq_error']:.4g}%, "
                    f"damping={100*cm['relative_zeta_error']:.4g}%, "
                    f"relative residual={cm['nep_relative_residual']:.3e}.",
                    f"- Nearest Acc pole to matched NEP zero: {cm['nearest_acc_real']:.6g} + "
                    f"j{cm['nearest_acc_imag']:.6g}, distance={cm['distance_to_nearest_acc']:.3g}.",
                    f"- Top condensed-state reconstruction: {cm['top_condensed_states']}",
                ]
            )
        lines.append("- Contour counts:")
        for row in case["contours"]:
            lines.append(
                f"  - {row['name']}: extracted={row['accepted_zeros']}, rank={row['beyn_rank']}, "
                f"argument-count={row['trace_count_real']:.3f}+j{row['trace_count_imag']:.2e}, "
                f"Acc poles inside={row['acc_poles_in_contour_upper']}, "
                f"argument+Acc={row['argument_count_plus_acc_poles']:.3f}, "
                f"ANDES poles inside={row['andes_poles_in_contour_upper']}, "
                f"count error={row['count_error_abs']:.3g}, pass={row['count_check_pass']}."
            )
        lines.extend(
            [
                f"- Zeros table: `{case['outputs']['zeros_csv']}`",
                f"- Contour diagnostics: `{case['outputs']['contour_counts_csv']}`",
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
    (OUT / "phase0d_nep_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    phase0b = load_phase0b_reference()
    results = []
    for case_name, path in case_paths():
        if not Path(path).exists():
            results.append({"case": case_name, "path": path, "accepted": False, "error": "case path missing"})
            continue
        print(f"[phase0d] running {case_name}: {path}")
        results.append(run_case(case_name, path, phase0b))

    for r in results:
        r["status_label"] = "PASS" if r.get("accepted") else "BLOCKED"

    accepted = all(r.get("accepted", False) for r in results)
    status = {
        "gate": "G0.5D rational Schur NEP by Beyn contour integration",
        "status": "PASS" if accepted else "BLOCKED",
        "random_seed": RNG_SEED,
        "acceptance_rule": (
            "For each of no_pss, base, and mix60: critical mode frequency error <= 5%, "
            "critical damping-ratio error <= 10%, correct damping sign, and accepted Beyn "
            "zero from contour with residual <= 1e-5."
        ),
        "solver": {
            "type": "Beyn contour integral",
            "nep": "S_e(s)=sI-Aee-Aec(sI-Acc)^(-1)Ace",
            "fixed_point_used_for_acceptance": False,
            "contours": [c.descriptor() for c in default_contours()],
        },
        "cases": results,
        "proceed_to_phase2": bool(accepted),
        "interpretation": (
            "Phase-0D passed for all three cases. The rational Schur NEP contour solver can be used "
            "as the bridge for later validation phases."
            if accepted
            else "Phase-0D did not pass the gate for all required cases. Do not advance to Phase 2-4 "
            "as full-ANDES evidence. Inspect the zero tables, control participation, and contour counts; "
            "the likely failure mode is insufficient projection of a controller-family mode onto the retained "
            "interface states or a contour/rank issue."
        ),
        "outputs": {
            "report": str((OUT / "phase0d_nep_report.md").relative_to(ROOT)),
            "status_json": str((OUT / "phase0d_nep_status.json").relative_to(ROOT)),
            "s_plane_png": str((OUT / "fig_phase0d_s_plane.png").relative_to(ROOT)),
            "s_plane_pdf": str((OUT / "fig_phase0d_s_plane.pdf").relative_to(ROOT)),
        },
    }
    write_json(OUT / "phase0d_nep_status.json", status)
    write_report(status)
    plot_results(status)
    print(json.dumps(status, indent=2, default=str))
    return 0 if accepted else 2


if __name__ == "__main__":
    raise SystemExit(main())
