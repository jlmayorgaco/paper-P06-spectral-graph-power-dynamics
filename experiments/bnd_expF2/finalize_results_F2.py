"""Rebuild the compact F2 modal-structure table from frozen numeric artifacts."""
from pathlib import Path
import tomllib

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
F0=ROOT/"reports"/"experiment_F0"
F2=ROOT/"reports"/"experiment_F2"
result=tomllib.loads((F2/"RESULTS_EXP_F2.toml").read_text(encoding="utf-8"))
U=pd.read_csv(F0/"matrices"/"Lc_eigenvectors.csv").iloc[:,1:].to_numpy(float)
G=pd.read_csv(F0/"matrices"/"Yport_real.csv").iloc[:,1:].to_numpy(float)
B=pd.read_csv(F0/"matrices"/"Yport_imag.csv").iloc[:,1:].to_numpy(float)
Y=G+1j*B
Gm=U.T@G@U
Bm=U.T@B@U
diag=lambda A: A-np.diag(np.diag(A))
metrics={
    "conductance_offdiag_ratio":np.linalg.norm(diag(Gm))/np.linalg.norm(Y),
    "susceptance_offdiag_ratio":np.linalg.norm(diag(Bm))/np.linalg.norm(Y),
    "commutator_ratio":np.linalg.norm(G@B-B@G)/(np.linalg.norm(G)*np.linalg.norm(B)),
    "network_admittance_offdiag_ratio":result["homogeneous_network_offdiag_ratio"],
    "homogeneous_closed_loop_port_offdiag_ratio":result["homogeneous_closed_loop_port_offdiag_ratio"],
    "homogeneous_descriptor_offdiag_ratio":result["homogeneous_descriptor_offdiag_ratio"],
    "rank_dA_dKp":result["gain_update_ranks"]["Kp"],
    "rank_dA_dKi":result["gain_update_ranks"]["Ki"],
    "rank_joint_gain_update":result["gain_update_ranks"]["joint"],
}
pd.DataFrame({"metric":list(metrics),"value":list(metrics.values())}).to_csv(
    F2/"tables"/"TABLE_F01_modal_structure.csv",index=False)
print("F2 modal structure table rebuilt from F0 network matrices and frozen F2 metrics")
