
from __future__ import annotations
import sys, json, subprocess
from pathlib import Path
import numpy as np
import pandas as pd

repo = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(repo/"src"))
from ibr_cycles.models.ieee39_network import load_network, solve_power_flow
from canonical_ppc import branch_flows_from_voltages

out = Path(__file__).resolve().parent / "results"
out.mkdir(exist_ok=True)

cfg = repo/"configs/ias2026"
net = load_network(cfg/"ieee39_network.json")
pf = solve_power_flow(net)
ref = json.loads((cfg/"ieee39_andes_powerflow_reference.json").read_text())

V = pf.voltages
vref = np.asarray(ref["v"])
aref = np.asarray(ref["a"])
a = np.angle(V)
offset = np.median(a-aref)

ylines = np.load(cfg/"ieee39_ybus_andes.npy")
yandes = ylines.copy()
for bus,b in ((4,1.0),(5,2.0)):
    i=net.position(bus); yandes[i,i] += 1j*b

metrics = {
    "tool":"pure_python_vs_frozen_ANDES_"+ref["tool"],
    "pf_converged":bool(pf.converged),
    "pf_max_mismatch":float(pf.max_mismatch),
    "ybus_max_abs_error":float(np.max(np.abs(net.ybus-yandes))),
    "vm_max_abs_error_pu":float(np.max(np.abs(np.abs(V)-vref))),
    "va_max_gauge_aligned_error_rad":float(np.max(np.abs((a-offset)-aref))),
}
flows = branch_flows_from_voltages(repo, V)
pd.DataFrame(flows, columns=["idx","from_bus","to_bus","P_from_MW","Q_from_Mvar","P_to_MW","Q_to_Mvar"]).to_csv(out/"python_branch_flows.csv",index=False)

# project tests
cp = subprocess.run([sys.executable,"-m","pytest","-q",
                     "tests/test_ieee39_network.py","tests/test_ieee39_baseline.py"],
                    cwd=repo, text=True, capture_output=True)
metrics["pytest_returncode"]=cp.returncode
metrics["pytest_stdout"]=cp.stdout.strip()

# stored equation-equivalent ANDES reconciliation
r = pd.read_csv(repo/"results/F1_eigenvalue_reconciliation.csv")
r3=r[r.stage.eq("R3 first-order AVR")]
metrics.update({
    "R3_max_alpha_error":float(r3.alpha_band_error.max()),
    "R3_max_frequency_error_hz":float(r3.freq_band_error_hz.max()),
    "R3_RHP_counts_all_agree":bool(r3.rhp_agrees.all()),
})
(out/"python_reference.json").write_text(json.dumps(metrics,indent=2))
print(json.dumps(metrics,indent=2))
