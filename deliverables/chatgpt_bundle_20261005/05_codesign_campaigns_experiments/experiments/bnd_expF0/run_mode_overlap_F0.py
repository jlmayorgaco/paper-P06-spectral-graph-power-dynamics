"""Post-preregistration projection of frozen ExpC q-modes onto the F0 graph basis."""
from __future__ import annotations

from pathlib import Path
import re

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[2]
F0 = ROOT / "reports" / "experiment_F0"
C0 = ROOT / "reports" / "experiment_C" / "matrices"
C0_TABLE = ROOT / "reports" / "experiment_C" / "tables" / "TABLE_C07_physical_pole_graph_mapping.csv"
TABLES = F0 / "tables"


def matrix(path: Path) -> np.ndarray:
    return pd.read_csv(path).iloc[:, 1:].to_numpy(float)


def main() -> None:
    M = matrix(C0 / "C0_M.csv")
    A = matrix(C0 / "C0_Ared.csv")
    states = pd.read_csv(C0 / "C0_states.csv")
    q_by_bus: dict[int, int] = {}
    for row in states.itertuples():
        name = str(row.state_name)
        if "machine" not in name or "δ" not in name:
            continue
        found = re.search(r"VIndex\((\d+),", name)
        if found:
            q_by_bus[int(found.group(1))] = int(row.state_index) - 1
    q_buses = sorted(q_by_bus)
    Lc = matrix(F0 / "matrices" / "Lc_conductance.csv")
    ports = pd.read_csv(TABLES / "TABLE_F02_graph_spectrum.csv").shape[0]
    expected_buses = list(range(30, 40))
    if q_buses != expected_buses or M.shape != (ports, ports) or A.shape[0] != len(states):
        status = (f"STOP: q-state map or mass matrix does not match F0 ports; "
                  f"q_buses={q_buses}, F0_ports={expected_buses}, M={M.shape}")
        (F0 / "F0_MODE_OVERLAP_STATUS.md").write_text(status + "\n", encoding="utf-8")
        print(status)
        return

    with threadpool_limits(limits=1):
        mval, mv = np.linalg.eigh((M + M.T) / 2)
        if mval.min() <= 0:
            status = f"STOP: ExpC mass matrix is not SPD; min eigenvalue={mval.min()}"
            (F0 / "F0_MODE_OVERLAP_STATUS.md").write_text(status + "\n", encoding="utf-8")
            print(status)
            return
        Mh = (mv * np.sqrt(mval)) @ mv.T
        Mih = (mv * (1 / np.sqrt(mval))) @ mv.T
        gval, Ug = np.linalg.eigh((Mih @ Lc @ Mih + Mih @ Lc.T @ Mih) / 2)
        poles = pd.read_csv(C0_TABLE)
        poles = poles[(poles.case == "C0_all_SG") &
                      poles.role.isin(["spectral_abscissa_pair", "lowest_damping_band_pair",
                                       "high_retained_participation"]) &
                      (poles.lambda_imag > 0) & (poles.q_norm > 0.05)].drop_duplicates("pole_id")
        lam, V = np.linalg.eig(A)
        rows = []
        qcols = []
        qidx = [q_by_bus[b] for b in q_buses]
        for pole in poles.itertuples():
            target = complex(float(pole.lambda_real), float(pole.lambda_imag))
            j = int(np.argmin(np.abs(lam - target)))
            err = float(abs(lam[j] - target))
            if err > 1e-6:
                continue
            q = Mh @ V[qidx, j]
            q /= np.linalg.norm(q)
            qcols.append(q)
            energy = np.abs(Ug.T.conj() @ q) ** 2
            order = np.argsort(energy)[::-1]
            rows.append({"pole_id": pole.pole_id, "lambda_real": target.real,
                         "lambda_imag": target.imag, "frequency_hz": pole.frequency_hz,
                         "q_norm_from_ExpC": pole.q_norm, "eigenvalue_match_error": err,
                         "top_graph_mode": int(order[0] + 1),
                         "top_mode_energy": float(energy[order[0]]),
                         "top3_graph_mode_energy": float(energy[order[:3]].sum())})
        if not rows:
            status = "STOP: no frozen ExpC critical poles matched the saved Ared matrix"
            (F0 / "F0_MODE_OVERLAP_STATUS.md").write_text(status + "\n", encoding="utf-8")
            print(status)
            return
        pd.DataFrame(rows).to_csv(TABLES / "TABLE_F05_graph_vs_real_modes.csv", index=False)
        Q = np.column_stack(qcols)
        Uq, singular, _ = np.linalg.svd(Q, full_matrices=False)
        rank = int(np.sum(singular > 1e-10 * singular.max()))
        Qorth = Uq[:, :rank]
        bands = [("low_modes_1_3", np.arange(0, min(3, ports))),
                 ("middle_modes_4_7", np.arange(3, min(7, ports))),
                 ("high_modes_8_10", np.arange(max(0, ports - 3), ports))]
        angle_rows = []
        for label, indices in bands:
            sigma = np.clip(np.linalg.svd(Ug[:, indices].T.conj() @ Qorth,
                                         compute_uv=False), 0.0, 1.0)
            for j, value in enumerate(sigma, start=1):
                angle_rows.append({"graph_band": label, "principal_angle_index": j,
                                   "cosine_overlap": value,
                                   "angle_degrees": np.degrees(np.arccos(value)),
                                   "critical_q_subspace_rank": rank})
        pd.DataFrame(angle_rows).to_csv(TABLES / "TABLE_F06_principal_angles.csv", index=False)

    status = (f"Computed post-preregistration M-weighted projections for {len(rows)} frozen ExpC critical poles; "
              f"critical q-subspace rank={rank}. Lc was not selected or modified using these overlaps.")
    (F0 / "F0_MODE_OVERLAP_STATUS.md").write_text(status + "\n", encoding="utf-8")
    report = F0 / "REPORT_EXP_F0.md"
    with report.open("a", encoding="utf-8") as stream:
        stream.write("\n## F0.4 Graph-mode relevance diagnostic\n\n")
        stream.write(status + " See TABLE_F05 for top-mode energy and TABLE_F06 for preregistered low/middle/high-band principal angles.\n")
    print(status)


if __name__ == "__main__":
    main()
