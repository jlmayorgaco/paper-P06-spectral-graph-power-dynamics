"""Independent matrix, pole, and bus-port comparison for the first ExpN probe."""
from __future__ import annotations

import csv
import re
from pathlib import Path

import numpy as np
from scipy.linalg import solve
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports" / "experiment_N"
PROBE = OUT / "identity_probe_bus38_half"


def table(name: str) -> list[dict]:
    with (PROBE / name).open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def matrix(name: str) -> np.ndarray:
    return np.atleast_2d(np.loadtxt(PROBE / name, delimiter=","))


def pd_key(raw: str) -> tuple[int, str, str]:
    match = re.fullmatch(r"VIndex\((\d+), :(.*)\)", raw)
    if not match:
        raise ValueError(f"unmapped PD state {raw}")
    bus, name = int(match[1]), match[2]
    tail = name.split("₊")[-1]
    if "₊gov₊" in name:
        return bus, "SG", "gov_" + tail
    if "₊avr₊" in name:
        return bus, "SG", "avr_" + tail
    if "machine₊" in name:
        machine = {"ψ″_q": "psi2q", "ψ″_d": "psi2d", "E′_d": "Epd",
                   "E′_q": "Epq", "ω": "omega", "δ": "delta"}
        return bus, "SG", "machine_" + machine[tail]
    if "₊cc1₊" in name:
        return bus, "GFL", {"γ_q": "gamma_q", "γ_d": "gamma_d"}[tail]
    if "₊pll₊" in name:
        return bus, "GFL", {"θ": "theta", "Δω_rad_s": "delta_omega_rad_s",
                            "Δω_i_rad_s": "delta_omega_i_rad_s"}[tail]
    if "₊filter₊" in name or "₊v_dc_" in name:
        return bus, "GFL", tail
    raise ValueError(f"unmapped PD state {raw}")


def transfer(m: np.ndarray, a: np.ndarray, b: np.ndarray,
             c: np.ndarray, d: np.ndarray, s: complex) -> np.ndarray:
    return c @ solve(s * m - a, b) + d


def main() -> None:
    pds = table("PD_state_map.csv")
    ans = table("AN_state_map.csv")
    pd_a = matrix("PD_A.csv")
    pd_m = matrix("PD_M.csv")
    dynamic = np.flatnonzero(np.diag(pd_m) == 1)
    algebraic = np.flatnonzero(np.diag(pd_m) == 0)
    pd_red = (pd_a[np.ix_(dynamic, dynamic)] -
              pd_a[np.ix_(dynamic, algebraic)] @
              solve(pd_a[np.ix_(algebraic, algebraic)],
                    pd_a[np.ix_(algebraic, dynamic)]))
    pd_lookup = {pd_key(pds[j]["state_name"]): k for k, j in enumerate(dynamic)}
    keys = [(int(r["bus"]), r["kind"], r["state_name"]) for r in ans]
    permutation = [pd_lookup[k] for k in keys]
    an_red = matrix("AN_Ared.csv")
    delta = pd_red[np.ix_(permutation, permutation)] - an_red
    matrix_rel = np.linalg.norm(delta) / np.linalg.norm(an_red)
    pd_poles = np.array([complex(float(r["real"]), float(r["imag"]))
                         for r in table("PD_poles.csv")])
    an_poles = np.array([complex(float(r["real"]), float(r["imag"]))
                         for r in table("AN_poles.csv")])
    pd_poles = np.delete(pd_poles, np.argmin(np.abs(pd_poles)))
    if len(pd_poles) != len(an_poles):
        raise RuntimeError(f"pole count {len(pd_poles)} vs {len(an_poles)}")
    rows, cols = linear_sum_assignment(np.abs(pd_poles[:, None] - an_poles[None, :]))
    pole_error = np.abs(pd_poles[rows] - an_poles[cols])

    states38 = np.array([int(r["index"]) - 1 for r in ans if int(r["bus"]) == 38])
    port38 = np.array([2 * 38 - 2, 2 * 38 - 1])
    aan = matrix("AN_A_local.csv")
    ban = matrix("AN_B_port.csv")
    can = matrix("AN_C_port.csv")
    dan = matrix("AN_D_port.csv")
    mp = matrix("PD_bus38_M.csv")
    ap = matrix("PD_bus38_A.csv")
    bp = matrix("PD_bus38_B.csv")
    cp = matrix("PD_bus38_C.csv")
    dp = matrix("PD_bus38_D.csv")
    errors = []
    signs = []
    for omega in np.geomspace(1e-4, 1e3, 51):
        s = 1j * omega
        yan = transfer(np.eye(len(states38)), aan[np.ix_(states38, states38)],
                       ban[np.ix_(states38, port38)],
                       can[np.ix_(port38, states38)],
                       dan[np.ix_(port38, port38)], s)
        zpd = transfer(mp, ap, bp, cp, dp, s)
        ypd = np.linalg.inv(zpd)
        sign = 1 if np.linalg.norm(ypd - yan) < np.linalg.norm(-ypd - yan) else -1
        signs.append(sign)
        errors.append(np.linalg.norm(sign * ypd - yan) /
                      max(np.linalg.norm(yan), 1e-15))
    result = dict(case="bus38_rho050_nominal", PD_dynamic=len(dynamic),
                  AN_dynamic=len(keys), reduced_A_relative_error=matrix_rel,
                  reduced_A_max_absolute_error=float(np.max(np.abs(delta))),
                  finite_pole_max_error=float(max(pole_error)),
                  PD_alpha=float(max(pd_poles.real)), AN_alpha=float(max(an_poles.real)),
                  port_max_relative_error=float(max(errors)),
                  port_sign_consistent=len(set(signs)) == 1,
                  status="PASS" if matrix_rel < 1e-8 and max(pole_error) < 1e-6 and
                  max(errors) < 1e-8 else "FAIL_MODEL")
    path = OUT / "TABLE_N05_identity_probe.csv"
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(result))
        writer.writeheader()
        writer.writerow(result)
    print(result)


if __name__ == "__main__":
    main()
