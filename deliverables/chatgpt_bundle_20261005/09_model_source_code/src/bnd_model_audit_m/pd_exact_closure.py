"""Exact local linear reconstruction from installed NetworkDynamics port operators.

This is a fixed-operating-point reconstruction. It does not define a dispatch
continuation or a replacement optimizer.
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

import numpy as np
from scipy.linalg import eigvals, solve
from scipy.optimize import linear_sum_assignment


def load(path: Path) -> np.ndarray:
    return np.atleast_2d(np.loadtxt(path, delimiter=","))


def _bus_state_order(state_map: Path) -> tuple[list[int], list[dict]]:
    with state_map.open(newline="", encoding="utf-8") as fh:
        states = list(csv.DictReader(fh))
    def bus(index: int) -> int:
        match = re.match(r"VIndex\((\d+),", states[index]["state_name"])
        if not match:
            raise ValueError(f"not a bus state: {states[index]['state_name']}")
        return int(match[1])
    order = sorted(range(len(states)), key=bus)
    return order, [states[i] for i in order]


def reconstruct(case_dir: Path) -> dict:
    """Positive-feedback closure for static IEEE-39 lines, with exact port IO."""
    m = load(case_dir / "PD_Zbus_M.csv")
    a = load(case_dir / "PD_Zbus_A.csv")
    b = load(case_dir / "PD_Zbus_B.csv")
    c = load(case_dir / "PD_Zbus_C.csv")
    d = load(case_dir / "PD_Zbus_D.csv")
    y = load(case_dir / "PD_Ynw_D.csv")
    yin = case_dir / "PD_Yinj_D.csv"
    if yin.exists():
        y += load(yin)
    # Installed NetworkDynamics tests call feedback(Zbus, Ynw + Yinj; pos=true).
    w = np.eye(y.shape[0]) - y @ d
    wyc = solve(w, y @ c)
    wr = solve(w, np.eye(w.shape[0]))
    return dict(M=m, A=a + b @ wyc, B=b @ wr,
                C=c + d @ wyc, D=d @ wr)


def finite_poles(a: np.ndarray, m: np.ndarray) -> np.ndarray:
    dynamic = np.flatnonzero(np.diag(m) == 1)
    algebraic = np.flatnonzero(np.diag(m) == 0)
    reduced = a[np.ix_(dynamic, dynamic)] - a[np.ix_(dynamic, algebraic)] @ solve(
        a[np.ix_(algebraic, algebraic)], a[np.ix_(algebraic, dynamic)])
    return eigvals(reduced)


def export_and_validate(case_dir: Path) -> dict:
    sys = reconstruct(case_dir)
    for key, matrix in sys.items():
        np.savetxt(case_dir / f"ANALYTICAL_PD_EXACT_{key}.csv", matrix, delimiter=",", fmt="%.17g")
    order, states = _bus_state_order(case_dir / "PD_state_map.csv")
    with (case_dir / "ANALYTICAL_PD_EXACT_state_map.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["closure_index", "PD_index", "state_name", "differential"])
        writer.writeheader()
        for i, state in enumerate(states, 1):
            writer.writerow(dict(closure_index=i, PD_index=order[i-1]+1,
                                 state_name=state["state_name"], differential=state["differential"]))
    pd_a = load(case_dir / "PD_A.csv")[np.ix_(order, order)]
    pd_m = load(case_dir / "PD_M.csv")[np.ix_(order, order)]
    da = sys["A"] - pd_a
    dm = sys["M"] - pd_m
    poles = finite_poles(sys["A"], sys["M"])
    with (case_dir / "ANALYTICAL_PD_EXACT_poles.csv").open("w", newline="", encoding="utf-8") as fh:
        writer=csv.writer(fh);writer.writerow(("real","imag"));writer.writerows(zip(poles.real,poles.imag))
    pdp=np.loadtxt(case_dir/"PD_poles.csv",delimiter=",",skiprows=1)
    reference=pdp[:,0]+1j*pdp[:,1]
    if len(reference)==len(poles):
        row,col=linear_sum_assignment(np.abs(reference[:,None]-poles[None,:]))
        pole_error=float(np.max(np.abs(reference[row]-poles[col])))
    else:pole_error=float("nan")
    return dict(case=case_dir.name, state_count=len(order),
        M_max_abs_error=float(np.max(np.abs(dm))),
        A_max_abs_error=float(np.max(np.abs(da))),
        A_relative_frobenius=float(np.linalg.norm(da)/np.linalg.norm(pd_a)),
        max_finite_pole_error_s_inv=pole_error,
        matrix_identity_pass=bool(np.linalg.norm(da)/np.linalg.norm(pd_a)<1e-9),
        pole_identity_pass=bool(pole_error<1e-8))
