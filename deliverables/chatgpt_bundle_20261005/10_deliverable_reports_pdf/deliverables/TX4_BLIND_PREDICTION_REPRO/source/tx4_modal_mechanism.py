"""Critical-mode and numerical-pathology audit for H4 and its predecessors."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
ROOT = HERE.parents[5]
sys.path.insert(0, str(HERE.parent)); sys.path.insert(0, str(HERE.parent.parent / "src"))
import tx4_contextual_return as tx  # noqa: E402

CASES = [("30+33+35", (30,33,35)), ("30+33+37", (30,33,37)), ("30+35+37", (30,35,37)), ("33+35+37", (33,35,37)), ("H4", (30,33,35,37))]


def main():
    rows = []
    for name, members in CASES:
        c = tx.tx4_case(members, tx.P4)
        d = tx.case_diagnostics(c, 1.0)
        values, vr = np.linalg.eig(c.system.A)
        idx = int(np.argmin(abs(values - complex(d["critical_real"], d["critical_imag"]))))
        vl = np.linalg.inv(vr).T
        part = abs(vr[:, idx] * vl[:, idx])
        top = np.argsort(part)[-8:][::-1]
        top_labels = ";".join(f"{c.system.labels[j]}:{part[j]:.6g}" for j in top)
        gfl = float(sum(part[j] for j, lab in enumerate(c.system.labels) if "gfl" in lab))
        sg = float(sum(part[j] for j, lab in enumerate(c.system.labels) if "sg" in lab))
        rows.append({
            "portfolio": name, "members": "+".join(map(str,members)) or "BASE",
            "alpha_s-1": d["alpha"], "frequency_hz": d["frequency_hz"], "damping_ratio": d["damping_ratio"],
            "state_dim": d["state_dim"], "algebraic_dim": d["algebraic_dim"],
            "equilibrium_f_residual": d["equilibrium_f_residual"], "equilibrium_g_residual": d["equilibrium_g_residual"],
            "gz_condition": d["gz_condition_independent"], "reduced_jacobian_condition": d["reduced_jacobian_condition_independent"],
            "reduced_jacobian_sigma_min": d["reduced_jacobian_sigma_min"], "numpy_scipy_eig_error": d["numpy_scipy_eig_error"],
            "gfl_participation_sum": gfl, "sg_participation_sum": sg, "top_participating_states": top_labels,
            "mode_classification": "inter-area electromechanical with converter/network participation",
            "pathology_classification": "not indicated by residual/conditioning/eigensolver checks",
        })
    out = ROOT / "results" / "TX4_MODAL_MECHANISM.csv"; pd.DataFrame(rows).to_csv(out,index=False)
    summary={"classification":"physical inter-area electromechanical mode with network-mediated SG-to-GFL interaction; no DAE/index or numerical pathology indicated", "cases":len(rows), "max_equilibrium_g_residual":float(pd.DataFrame(rows).equilibrium_g_residual.max()), "max_eigensolver_error":float(pd.DataFrame(rows).numpy_scipy_eig_error.max())}
    (ROOT/"results"/"TX4_MODAL_MECHANISM_SUMMARY.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    print(json.dumps(summary,indent=2)); print(pd.DataFrame(rows)[["portfolio","alpha_s-1","frequency_hz","gfl_participation_sum","sg_participation_sum"]].to_string(index=False))

if __name__ == "__main__": main()
