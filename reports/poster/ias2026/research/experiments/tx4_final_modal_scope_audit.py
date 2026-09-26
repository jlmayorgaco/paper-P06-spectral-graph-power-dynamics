"""TX4 final correction and V9 modal-scope audit.

This script is deliberately a finite correction campaign.  It does not alter
the frozen TX4 policy or rerun a parameter search.  It recomputes the H4
physical local factors and the complete transverse spectra of the frozen V9
census so that the reduced 0.3--1.5 Hz target is not confused with global
safety.
"""

from __future__ import annotations

import json
import math
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
EXP = HERE.parent
SRC = EXP.parent / "src"
ROOT = HERE.parents[5]
sys.path.insert(0, str(EXP))
sys.path.insert(0, str(SRC))

from _f7_common import LEAK, Theta  # noqa: E402
from ibr_cycles.certification.symmetry import frequency_partner, rotation_generator  # noqa: E402
from ibr_cycles.certification.transverse import transverse_operator  # noqa: E402
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case  # noqa: E402
from ibr_cycles.models.ieee39_devices import ConverterParameters  # noqa: E402
from ibr_cycles.models.port_admittance import build_action_space  # noqa: E402

OUT = ROOT / "results"
P4 = Theta(g=0.03625, k=1.425, t=1.5, h=1.0)
H4 = (30, 33, 35, 37)
V9 = tuple(range(30, 39))
BAND = (0.3, 1.5)
ZERO_EIG = 1e-3
REAL_FREQ_HZ = 1e-6
G_SWEEP = (0.03625, 0.075, 0.125, 0.175, 0.2, 0.205, 0.2075, 0.20768, 0.2076814045, 0.2077, 0.21, 0.225, 0.25)


def make_case(members: tuple[int, ...], g: float = P4.g):
    return solve_case(
        ReplacementPlan.of({b: 1.0 for b in members}),
        converter=ConverterParameters(voltage_control=True, voltage_gain=g, voltage_leak=LEAK),
        machine_scaling={"ka": P4.k, "ta": P4.t},
    )


def transverse(case):
    rx, _ = rotation_generator(case.dae, case.equilibrium.z)
    fp = frequency_partner(case.dae)
    tr = transverse_operator(case.system.A, rx, fp.w)
    vals, vecs = np.linalg.eig(tr.a_perp)
    return tr, vals, vecs


def admissible_values(vals: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    keep = np.isfinite(vals.real) & np.isfinite(vals.imag) & (np.abs(vals) > ZERO_EIG)
    return vals[keep], np.flatnonzero(keep)


def mode_summary(case) -> dict:
    tr, vals, vecs = transverse(case)
    finite, indices = admissible_values(vals)
    if len(finite) == 0:
        raise RuntimeError("no finite transverse eigenvalues")
    local_idx = int(np.argmax(finite.real))
    global_idx = int(indices[local_idx])
    lam_global = complex(vals[global_idx])
    f_global = abs(lam_global.imag) / (2.0 * math.pi)
    band_mask = (abs(finite.imag) / (2.0 * math.pi) >= BAND[0]) & (abs(finite.imag) / (2.0 * math.pi) <= BAND[1])
    if band_mask.any():
        band_values = finite[band_mask]
        em_idx = int(np.argmax(band_values.real))
        lam_em = complex(band_values[em_idx])
        alpha_em = float(lam_em.real)
        f_em = float(abs(lam_em.imag) / (2.0 * math.pi))
    else:
        lam_em = complex(np.nan, np.nan)
        alpha_em = float("nan")
        f_em = float("nan")

    full_vec = tr.z @ vecs[:, global_idx]
    full_vec = full_vec / max(float(np.linalg.norm(full_vec)), 1e-300)
    sg_part = 0.0
    gfl_part = 0.0
    labels = []
    for slot in case.dae.slots:
        weight = float(np.sum(np.abs(full_vec[slot.start : slot.stop]) ** 2))
        labels.append((weight, slot.kind, int(slot.bus)))
        if slot.kind == "sg":
            sg_part += weight
        elif slot.kind == "gfl":
            gfl_part += weight
    labels.sort(reverse=True)
    if f_global <= REAL_FREQ_HZ:
        mode_class = "APERIODIC_REAL"
    elif f_global < BAND[0]:
        mode_class = "SLOW_OUTSIDE_BAND"
    elif f_global <= BAND[1]:
        mode_class = "TARGET_EM_BAND"
    else:
        mode_class = "FAST_OSCILLATORY"
    return {
        "alpha_global": float(lam_global.real),
        "lambda_global_real": float(lam_global.real),
        "lambda_global_imag": float(lam_global.imag),
        "f_global_hz": float(f_global),
        "alpha_EM": alpha_em,
        "lambda_EM_real": float(lam_em.real),
        "lambda_EM_imag": float(lam_em.imag),
        "f_EM_hz": f_em,
        "global_mode_class": mode_class,
        "global_mode_in_EM_band": bool(mode_class == "TARGET_EM_BAND"),
        "global_mode_real_aperiodic": bool(mode_class == "APERIODIC_REAL"),
        "global_mode_fast_oscillatory": bool(mode_class == "FAST_OSCILLATORY"),
        "gfl_participation_sum": gfl_part,
        "sg_participation_sum": sg_part,
        "top_participating_slots": ";".join(f"{kind}:{bus}:{weight:.4g}" for weight, kind, bus in labels[:6]),
        "transverse_dimension": int(tr.a_perp.shape[0]),
        "equilibrium_g_residual": float(np.max(np.abs(case.dae.g(case.equilibrium.x, case.equilibrium.z, {})))),
    }


def contextual_return(q: np.ndarray, i: int) -> tuple[np.ndarray, float]:
    n = q.shape[0] // 2
    ridx = sum(([2 * k, 2 * k + 1] for k in range(n) if k != i), [])
    ir = np.ix_([2 * i, 2 * i + 1], ridx)
    ri = np.ix_(ridx, [2 * i, 2 * i + 1])
    rr = np.ix_(ridx, ridx)
    ret = q[ir] @ np.linalg.solve(np.eye(len(ridx), dtype=complex) + q[rr], q[ri])
    ch = np.eye(2 * n, dtype=complex) + q
    cr = np.eye(len(ridx), dtype=complex) + q[rr]
    residual = abs(np.linalg.det(ch) - np.linalg.det(cr) * np.linalg.det(np.eye(2) - ret))
    return np.linalg.eigvals(ret), float(residual)


def q_from_m(m: np.ndarray) -> np.ndarray:
    n = m.shape[0] // 2
    total = np.eye(2 * n, dtype=complex) + m
    local = np.zeros_like(total)
    for i in range(n):
        local[2 * i : 2 * i + 2, 2 * i : 2 * i + 2] = total[2 * i : 2 * i + 2, 2 * i : 2 * i + 2]
    return np.linalg.solve(local, total) - np.eye(2 * n, dtype=complex)


def sign_audit() -> pd.DataFrame:
    raw = pd.read_csv(OUT / "TX4_CONTEXTUAL_RETURN_BOUNDARY.csv")
    g_star = float(raw.g_root.iloc[0])
    omega_star = float(raw.s_imag.iloc[0])
    base = make_case(())
    case = make_case(H4, g_star)
    space = build_action_space(base, case, H4)
    m = space.m(1j * omega_star)
    q = q_from_m(m)
    q_lam = np.linalg.eigvals(q)
    q_near = complex(q_lam[np.argmin(abs(q_lam + 1.0))])
    rows = []
    for i in range(4):
        mus, residual = contextual_return(q, i)
        mu = complex(mus[np.argmin(abs(mus - 1.0))])
        rows.append(
            {
                "device_i": H4[i],
                "q_eigenvalue_real": q_near.real,
                "q_eigenvalue_imag": q_near.imag,
                "abs(q_eigenvalue_plus_one)": abs(q_near + 1.0),
                "return_eigenvalue_real": mu.real,
                "return_eigenvalue_imag": mu.imag,
                "abs(return_eigenvalue_minus_one)": abs(mu - 1.0),
                "schur_residual": residual,
            }
        )
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "TX4_Q_VS_RETURN_SIGN_AUDIT.csv", index=False)
    return frame


def physical_local_audit() -> tuple[pd.DataFrame, pd.DataFrame]:
    base = make_case(())
    local_rows = []
    trajectories = []
    for g in G_SWEEP:
        case = make_case(H4, float(g))
        info = mode_summary(case)
        s = 1j * 2.0 * math.pi * info["f_EM_hz"]
        space = build_action_space(base, case, H4)
        m = space.m(s)
        q = q_from_m(m)
        collective = float(np.linalg.svd(np.eye(2 * len(H4)) + q, compute_uv=False)[-1])
        q_near = complex(np.linalg.eigvals(q)[np.argmin(abs(np.linalg.eigvals(q) + 1.0))])
        for i, bus in enumerate(H4):
            block = np.eye(2, dtype=complex) + space.block(m, i, i)
            mus, _ = contextual_return(q, i)
            mu = complex(mus[np.argmin(abs(mus - 1.0))])
            local_rows.append(
                {
                    "portfolio": "30+33+35+37",
                    "device_i": bus,
                    "g": float(g),
                    "frequency_hz": info["f_EM_hz"],
                    "alpha_EM": info["alpha_EM"],
                    "local_sigma_min_I_plus_Mii": float(np.linalg.svd(block, compute_uv=False)[-1]),
                    "local_factor_det_abs": float(abs(np.linalg.det(block))),
                    "collective_sigma_min_I_plus_QH": collective,
                    "q_eigenvalue_real": q_near.real,
                    "q_eigenvalue_imag": q_near.imag,
                    "return_eigenvalue_real": mu.real,
                    "return_eigenvalue_imag": mu.imag,
                }
            )
        trajectories.append(
            {
                "g": float(g),
                "alpha_EM": info["alpha_EM"],
                "frequency_hz": info["f_EM_hz"],
                "q_real": q_near.real,
                "q_imag": q_near.imag,
                "collective_sigma_min": collective,
            }
        )
    local_frame = pd.DataFrame(local_rows)
    local_frame.to_csv(OUT / "TX4_TRUE_LOCAL_VS_COLLECTIVE.csv", index=False)

    raw = pd.read_csv(OUT / "TX4_CONTEXTUAL_RETURN_BOUNDARY.csv")
    g_star = float(raw.g_root.iloc[0])
    omega_star = float(raw.s_imag.iloc[0])
    proper_rows = []
    for k in range(len(H4) + 1):
        for subset in combinations(H4, k):
            label = "+".join(map(str, subset)) or "BASE"
            case = make_case(tuple(subset))
            if subset:
                space = build_action_space(base, case, tuple(subset))
                split = space.split(1j * omega_star)
                sigma_q = split["sigma_min_i_plus_q"]
                det_q = abs(split["collective"])
            else:
                sigma_q = 1.0
                det_q = 1.0
            proper_rows.append(
                {
                    "portfolio": label,
                    "cardinality": len(subset),
                    "frequency_hz": omega_star / (2.0 * math.pi),
                    "sigma_min_I_plus_Q": sigma_q,
                    "det_I_plus_Q_abs": det_q,
                    "is_H4": label == "30+33+35+37",
                }
            )
    proper_frame = pd.DataFrame(proper_rows)
    proper_frame.to_csv(OUT / "TX4_TRUE_COLLECTIVE_SUBSET_FACTORS.csv", index=False)
    return local_frame, proper_frame


def v9_modal_census() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    pred = pd.read_csv(OUT / "TX4_V9_BLIND_VS_FULL.csv")[["portfolio", "predicted_verdict", "predicted_root_real", "predicted_root_frequency_hz", "predicted_local_factor_min_sigma", "predicted_local_factor_min_det_abs"]]
    rows = []
    for n in range(len(V9) + 1):
        for subset in combinations(V9, n):
            label = "+".join(map(str, subset)) or "BASE"
            case = make_case(tuple(subset))
            info = mode_summary(case)
            row = {"portfolio": label, "cardinality": n, **info}
            rows.append(row)
            if len(rows) == 1 or len(rows) % 64 == 0 or len(rows) == 512:
                print(f"V9 modal spectrum {len(rows)}/512", flush=True)
    truth = pd.DataFrame(rows)
    joined = truth.merge(pred, on="portfolio", how="left", validate="one_to_one")
    joined["global_status"] = np.where(joined.alpha_global >= 0.0, "UNSTABLE", "STABLE")
    joined["EM_status"] = np.where(joined.alpha_EM >= 0.0, "UNSTABLE", "STABLE")
    joined["global_false_safe"] = (joined.predicted_verdict == "STABLE_PREDICTED") & (joined.global_status == "UNSTABLE")
    joined["EM_false_safe"] = (joined.predicted_verdict == "STABLE_PREDICTED") & (joined.EM_status == "UNSTABLE")
    joined["EM_false_unstable"] = (joined.predicted_verdict == "UNSTABLE_PREDICTED") & (joined.EM_status == "STABLE")
    joined.to_csv(OUT / "TX4_V9_GLOBAL_VS_EM_TRUTH.csv", index=False)

    truth_em = joined[joined.alpha_EM.notna()].copy()
    tp = int(((truth_em.predicted_verdict == "UNSTABLE_PREDICTED") & (truth_em.EM_status == "UNSTABLE")).sum())
    tn = int(((truth_em.predicted_verdict == "STABLE_PREDICTED") & (truth_em.EM_status == "STABLE")).sum())
    fs = int(truth_em.EM_false_safe.sum())
    fu = int(truth_em.EM_false_unstable.sum())
    precision = tp / max(tp + fu, 1)
    recall = tp / max(tp + fs, 1)
    specificity = tn / max(tn + fu, 1)
    confusion = pd.DataFrame(
        [
            {
                "n_portfolios": len(truth_em),
                "accuracy": (tp + tn) / max(len(truth_em), 1),
                "true_positive": tp,
                "true_negative": tn,
                "false_safe": fs,
                "false_unstable": fu,
                "precision": precision,
                "recall": recall,
                "specificity": specificity,
                "unresolved_EM": int(joined.alpha_EM.isna().sum()),
                "band_hz_low": BAND[0],
                "band_hz_high": BAND[1],
            }
        ]
    )
    confusion.to_csv(OUT / "TX4_V9_EM_BAND_CONFUSION.csv", index=False)

    false_safe = joined[joined.global_false_safe].copy()
    taxonomy = false_safe[["portfolio", "cardinality", "alpha_global", "f_global_hz", "global_mode_class", "lambda_global_real", "lambda_global_imag", "alpha_EM", "f_EM_hz", "predicted_root_real", "predicted_root_frequency_hz"]].copy()
    taxonomy["taxonomy"] = taxonomy.global_mode_class
    taxonomy.to_csv(OUT / "TX4_V9_FALSE_SAFE_TAXONOMY.csv", index=False)

    slow = joined[joined.global_false_safe & (joined.f_global_hz > REAL_FREQ_HZ) & (joined.f_global_hz < BAND[0])].copy()
    slow["slow_audit_selection"] = "global false-safe with 0 < f_global < 0.3 Hz"
    slow["mode_mac_to_reduced"] = "NOT_AVAILABLE: reduced closure has no common full-state mode vector"
    slow["closure_root"] = slow.predicted_root_real
    slow["local_factor_regular"] = slow.predicted_local_factor_min_sigma > 1e-7
    slow.to_csv(OUT / "TX4_V9_SLOW_FALSE_SAFE_AUDIT.csv", index=False)
    summary = {
        "global_n": int(len(joined)),
        "global_false_safe": int(joined.global_false_safe.sum()),
        "em_n": int(len(truth_em)),
        "em_false_safe": fs,
        "em_false_unstable": fu,
        "taxonomy_counts": taxonomy.taxonomy.value_counts().to_dict(),
        "slow_false_safe_count": int(len(slow)),
        "slow_selection_rule": "all global false-safe cases with 0 < f_global < 0.3 Hz",
    }
    (OUT / "TX4_FINAL_MODAL_SCOPE_V9_SUMMARY.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return joined, taxonomy, slow


def main() -> None:
    sign = sign_audit()
    local, proper = physical_local_audit()
    joined, taxonomy, slow = v9_modal_census()
    print("SIGN AUDIT")
    print(sign.to_string(index=False))
    print("PHYSICAL LOCAL MIN", local.local_sigma_min_I_plus_Mii.min())
    print("COLLECTIVE MIN", local.collective_sigma_min_I_plus_QH.min())
    print("V9 TAXONOMY")
    print(taxonomy.taxonomy.value_counts().to_string())
    print("SLOW CASES", len(slow))
    print("PROPER MIN", proper.loc[~proper.is_H4, "sigma_min_I_plus_Q"].min())


if __name__ == "__main__":
    main()
