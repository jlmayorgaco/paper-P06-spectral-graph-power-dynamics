"""Evaluate ExpN withheld identity gates against independently compiled PD cases."""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy.linalg import solve
from scipy.optimize import linear_sum_assignment

from compare_identity_probe import pd_key

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports" / "experiment_N"
CASES = OUT / "validation_cases"
FREQUENCIES = np.geomspace(1e-4, 1e3, 31)


def read(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def write(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def matrix(path: Path) -> np.ndarray:
    return np.atleast_2d(np.loadtxt(path, delimiter=","))


def poles(path: Path) -> np.ndarray:
    return np.array([complex(float(r["real"]), float(r["imag"])) for r in read(path)])


def transfer(m: np.ndarray, a: np.ndarray, b: np.ndarray,
             c: np.ndarray, d: np.ndarray, s: complex) -> np.ndarray:
    return c @ solve(s * m - a, b) + d


def load_admittance(bus: int, original: dict[int, dict]) -> np.ndarray:
    if bus not in (31, 39):
        return np.zeros((2, 2))
    row = original[bus]
    v = complex(float(row["V_real_pu"]), float(row["V_imag_pu"]))
    s = complex(float(row["P_load_MW"]), float(row["Q_load_Mvar"])) / 100
    y = np.conj(s / v) / v
    return np.array([[y.real, -y.imag], [y.imag, y.real]])


def main() -> None:
    exported = read(OUT / "TABLE_N05_validation_case_export.csv")
    original = {int(r["bus"]): r for r in read(OUT / "TABLE_N01_original_operating_point.csv")}
    matrix_rows, port_rows, pole_rows = [], [], []
    for case in exported:
        label = case["case"]
        path = CASES / label
        pds = read(path / "PD_state_map.csv")
        ans = read(path / "AN_state_map.csv")
        pd_a = matrix(path / "PD_A.csv")
        pd_m = matrix(path / "PD_M.csv")
        dynamic = np.flatnonzero(np.diag(pd_m) == 1)
        algebraic = np.flatnonzero(np.diag(pd_m) == 0)
        pd_red = (pd_a[np.ix_(dynamic, dynamic)] -
                  pd_a[np.ix_(dynamic, algebraic)] @
                  solve(pd_a[np.ix_(algebraic, algebraic)],
                        pd_a[np.ix_(algebraic, dynamic)]))
        lookup = {pd_key(pds[j]["state_name"]): k for k, j in enumerate(dynamic)}
        keys = [(int(r["bus"]), r["kind"], r["state_name"]) for r in ans]
        if set(keys) != set(lookup):
            raise RuntimeError(f"state inventory mismatch at {label}: "
                               f"missing={set(keys)-set(lookup)}, extra={set(lookup)-set(keys)}")
        permutation = [lookup[k] for k in keys]
        an_red = matrix(path / "AN_Ared.csv")
        difference = pd_red[np.ix_(permutation, permutation)] - an_red
        rel = np.linalg.norm(difference) / max(np.linalg.norm(an_red), 1e-300)
        absmax = np.max(np.abs(difference))
        matrix_rows.append(dict(case=label, development=case["development"],
            PD_dynamic=len(dynamic), AN_dynamic=len(keys),
            reduced_A_relative_error=float(rel), reduced_A_max_absolute=float(absmax),
            trim_residual_max=case["residual_max"],
            max_P_error_pu=case["max_P_error_pu"], max_Q_error_pu=case["max_Q_error_pu"],
            bounds_ok=case["bounds_ok"],
            matrix_gate="PASS" if rel < 1e-8 else "FAIL_MODEL"))

        pd_lambda = poles(path / "PD_poles.csv")
        an_lambda = poles(path / "AN_poles.csv")
        gauge_pd = int(np.argmin(np.abs(pd_lambda)))
        pd_finite = np.delete(pd_lambda, gauge_pd)
        if len(pd_finite) != len(an_lambda):
            raise RuntimeError(f"pole count mismatch at {label}")
        rr, cc = linear_sum_assignment(np.abs(pd_finite[:, None] - an_lambda[None, :]))
        errors = np.abs(pd_finite[rr] - an_lambda[cc])
        pole_rows.append(dict(case=label, development=case["development"],
            PD_gauge_real=float(pd_lambda[gauge_pd].real),
            PD_gauge_imag=float(pd_lambda[gauge_pd].imag),
            PD_finite_count=len(pd_finite), AN_finite_count=len(an_lambda),
            max_matched_pole_error=float(max(errors)),
            median_matched_pole_error=float(np.median(errors)),
            PD_alpha=float(max(pd_finite.real)), AN_alpha=float(max(an_lambda.real)),
            alpha_difference=float(max(pd_finite.real)-max(an_lambda.real)),
            pole_gate="PASS" if max(errors) < 1e-6 else "FAIL_MODEL"))

        aan = matrix(path / "AN_A_local.csv")
        ban = matrix(path / "AN_B_port.csv")
        can = matrix(path / "AN_C_port.csv")
        dan = matrix(path / "AN_D_port.csv")
        active = [int(r["bus"]) for r in read(path / "case_definition.csv")
                  if float(r["rho"]) > 0]
        for bus in active:
            ids = np.array([int(r["index"]) - 1 for r in ans if int(r["bus"]) == bus])
            yi = np.array([2 * bus - 2, 2 * bus - 1])
            mp = matrix(path / f"PD_bus{bus}_M.csv")
            ap = matrix(path / f"PD_bus{bus}_A.csv")
            bp = matrix(path / f"PD_bus{bus}_B.csv")
            cp = matrix(path / f"PD_bus{bus}_C.csv")
            dp = matrix(path / f"PD_bus{bus}_D.csv")
            adm = load_admittance(bus, original)
            port_errors = []
            for omega in FREQUENCIES:
                s = 1j * omega
                yan = (transfer(np.eye(len(ids)), aan[np.ix_(ids, ids)],
                                ban[np.ix_(ids, yi)], can[np.ix_(yi, ids)],
                                dan[np.ix_(yi, yi)], s) + adm)
                zpd = transfer(mp, ap, bp, cp, dp, s)
                ypd = -np.linalg.inv(zpd)  # installed positive-feedback convention
                port_errors.append(np.linalg.norm(ypd - yan) /
                                   max(np.linalg.norm(yan), 1e-15))
            port_rows.append(dict(case=label, bus=bus,
                development=case["development"], samples=len(FREQUENCIES),
                max_relative_error=float(max(port_errors)),
                median_relative_error=float(np.median(port_errors)),
                omega_max=float(FREQUENCIES[int(np.argmax(port_errors))]),
                port_gate="PASS" if max(port_errors) < 1e-8 else "FAIL_MODEL"))
        print(label, "A", rel, "pole", max(errors),
              "port", max((r["max_relative_error"] for r in port_rows if r["case"] == label), default="NA"))
    write(OUT / "TABLE_N05_parametric_identity.csv", matrix_rows)
    write(OUT / "TABLE_N06_withheld_port_identity.csv", port_rows)
    write(OUT / "TABLE_N07_withheld_pole_identity.csv", pole_rows)
    withheld = [r for r in matrix_rows if r["development"] == "false"]
    print("N_IDENTITY_SUMMARY cases", len(exported), "withheld", len(withheld),
          "max_A", max(r["reduced_A_relative_error"] for r in withheld),
          "max_port", max(r["max_relative_error"] for r in port_rows if r["development"] == "false"),
          "max_pole", max(r["max_matched_pole_error"] for r in pole_rows if r["development"] == "false"))


if __name__ == "__main__":
    main()
