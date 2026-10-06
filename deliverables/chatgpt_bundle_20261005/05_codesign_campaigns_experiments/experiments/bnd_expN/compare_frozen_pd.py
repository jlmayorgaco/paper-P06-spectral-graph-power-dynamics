"""Independent complete-spectrum audit of the frozen ExpN candidate."""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy.linalg import eig, solve
from scipy.optimize import linear_sum_assignment

from compare_identity_probe import pd_key

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports" / "experiment_N"
P = OUT / "postfreeze_pd"


def read(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def matrix(name: str) -> np.ndarray:
    return np.atleast_2d(np.loadtxt(P / name, delimiter=","))


def poles(name: str) -> np.ndarray:
    return np.array([complex(float(r["real"]), float(r["imag"]))
                     for r in read(P / name)])


def main() -> None:
    raw = read(OUT / "TABLE_N13_powerdynamics_raw.csv")[0]
    pd_map = read(P / "PD_state_map.csv")
    an_map = read(P / "AN_state_map.csv")
    a = matrix("PD_A.csv")
    m = matrix("PD_M.csv")
    d = np.flatnonzero(np.diag(m) == 1)
    z = np.flatnonzero(np.diag(m) == 0)
    reduced = a[np.ix_(d, d)] - a[np.ix_(d, z)] @ solve(
        a[np.ix_(z, z)], a[np.ix_(z, d)])
    lookup = {pd_key(pd_map[j]["state_name"]): i for i, j in enumerate(d)}
    keys = [(int(r["bus"]), r["kind"], r["state_name"]) for r in an_map]
    perm = [lookup[k] for k in keys]
    aligned = reduced[np.ix_(perm, perm)]
    an_a = matrix("AN_Ared.csv")
    matrix_rel = np.linalg.norm(aligned - an_a) / np.linalg.norm(an_a)
    pd_poles = poles("PD_poles.csv")
    pd_poles = np.delete(pd_poles, np.argmin(abs(pd_poles)))
    an_poles = poles("AN_poles.csv")
    if len(pd_poles) != len(an_poles):
        raise RuntimeError("finite pole count mismatch")
    row, col = linear_sum_assignment(abs(pd_poles[:, None] - an_poles[None, :]))
    errors = abs(pd_poles[row] - an_poles[col])
    j = int(np.argmax(pd_poles.real))
    pd_right = pd_poles[j]
    matched = an_poles[col[np.where(row == j)[0][0]]]
    eig_pd, vec_pd = eig(aligned)
    eig_an, vec_an = eig(an_a)
    ip = int(np.argmin(abs(eig_pd - pd_right)))
    ia = int(np.argmin(abs(eig_an - matched)))
    overlap = abs(np.vdot(vec_pd[:, ip], vec_an[:, ia])) / (
        np.linalg.norm(vec_pd[:, ip]) * np.linalg.norm(vec_an[:, ia]))
    alpha_pd = float(pd_poles.real.max())
    alpha_an = float(an_poles.real.max())
    valid = (float(raw["trim_residual"]) < 1e-10 and
             float(raw["max_P_error_pu"]) < 1e-8 and
             float(raw["max_Q_error_pu"]) < 1e-8 and
             raw["physical_bounds_pass"].lower() == "true" and
             matrix_rel < 1e-8 and errors.max() < 1e-5 and
             abs(alpha_pd-alpha_an) < 1e-4 and alpha_pd <= -0.05 and
             overlap > 1-1e-6)
    result = {"candidate_sha": raw["candidate_sha"],
              "trim_residual": raw["trim_residual"],
              "P_error_pu": raw["max_P_error_pu"],
              "Q_error_pu": raw["max_Q_error_pu"],
              "physical_bounds_pass": raw["physical_bounds_pass"],
              "PD_alpha": alpha_pd, "AN_alpha": alpha_an,
              "alpha_error": abs(alpha_pd-alpha_an),
              "PD_rightmost_real": pd_right.real,
              "PD_rightmost_imag": pd_right.imag,
              "AN_matched_real": matched.real,
              "AN_matched_imag": matched.imag,
              "rightmost_pole_error": abs(pd_right-matched),
              "complete_finite_pole_max_error": float(errors.max()),
              "reduced_A_relative_error": matrix_rel,
              "critical_right_eigenvector_overlap": overlap,
              "PD_finite_count": len(pd_poles),
              "AN_finite_count": len(an_poles),
              "PD_validation": "PASS" if valid else "FAIL"}
    path = OUT / "TABLE_N13_powerdynamics_postfreeze_validation.csv"
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(result))
        writer.writeheader()
        writer.writerow(result)
    print(result)
    if not valid:
        raise SystemExit("FAIL_POWERDYNAMICS")


if __name__ == "__main__":
    main()
