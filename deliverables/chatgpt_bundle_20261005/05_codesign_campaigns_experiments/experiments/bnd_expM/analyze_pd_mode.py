"""Locate and classify the failing ExpK pole in the direct PD Jacobian."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eig, solve


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports" / "experiment_M"
CASE = sys.argv[1] if len(sys.argv) > 1 else "ExpK_nominal"
MATRICES = OUT / "matrices" / CASE


def main() -> None:
    a = np.loadtxt(MATRICES / "PD_A.csv", delimiter=",")
    m = np.loadtxt(MATRICES / "PD_M.csv", delimiter=",")
    with (MATRICES / "PD_state_map.csv").open(newline="", encoding="utf-8") as fh:
        states = list(csv.DictReader(fh))
    d = np.flatnonzero(np.diag(m) == 1)
    z = np.flatnonzero(np.diag(m) == 0)
    reduced = a[np.ix_(d, d)] - a[np.ix_(d, z)] @ solve(a[np.ix_(z, z)], a[np.ix_(z, d)])
    eigenvalues, left, right = eig(reduced, left=True, right=True)
    if CASE == "ExpG_candidate":
        # The old stability_audit labels both near-zero values as gauge. Keep
        # the second one visible and inspect its state structure explicitly.
        j = int(np.argsort(np.abs(eigenvalues))[1])
    else:
        physical = np.where(np.abs(eigenvalues) > 1e-6, eigenvalues.real, -np.inf)
        j = int(np.argmax(physical))
    contributions = np.abs(left[:, j].conj() * right[:, j])
    contributions /= np.sum(contributions)
    rows = []
    for k in np.argsort(contributions)[::-1]:
        name = states[int(d[k])]["state_name"]
        family = ("GFL_PLL" if "pll" in name else "GFL_CURRENT" if "cc1" in name else
                  "GFL_FILTER" if "filter" in name else "GFL_DC" if "v_dc" in name else
                  "SG_GOV" if "gov" in name else "SG_AVR" if "avr" in name else
                  "SG_MACHINE" if "machine" in name else "OTHER")
        rows.append(dict(case=CASE, pole_real=eigenvalues[j].real, pole_imag=eigenvalues[j].imag,
                         rank=len(rows) + 1, state_index=int(d[k]) + 1, state_name=name,
                         family=family, participation=float(contributions[k])))
    suffix = "" if CASE == "ExpK_nominal" else "_" + CASE
    table = OUT / "tables" / ("TABLE_M17_critical_pole_participation" + suffix + ".csv")
    table.parent.mkdir(parents=True, exist_ok=True)
    with table.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print("critical pole", eigenvalues[j], "states", len(d), "algebraic", len(z))
    for row in rows[:20]:
        print(f'{row["participation"]:.6f}', row["family"], row["state_name"].encode("ascii", "backslashreplace").decode())


if __name__ == "__main__":
    main()
