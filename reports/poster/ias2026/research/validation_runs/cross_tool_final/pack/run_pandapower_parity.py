
from __future__ import annotations
import sys, json
from pathlib import Path
import numpy as np
import pandas as pd

try:
    import pandapower as pp
    from pandapower.converter import from_ppc
except Exception as e:
    raise SystemExit(f"pandapower unavailable: {e}")

repo=Path(sys.argv[1]).resolve()
sys.path.insert(0,str(repo/"src"))
from ibr_cycles.models.ieee39_network import load_network, solve_power_flow
from canonical_ppc import canonical_to_ppc, branch_flows_from_voltages

out=Path(__file__).resolve().parent/"results"
out.mkdir(exist_ok=True)

canonical=canonical_to_ppc(repo)
net=from_ppc(canonical, f_hz=60, validate_conversion=False)
pp.runpp(net, algorithm="nr", calculate_voltage_angles=True, init="flat",
         enforce_q_lims=False, tolerance_mva=1e-10, max_iteration=50)

# from_ppc preserves ppc BUS_I as pandapower bus index
ids=np.asarray(sorted(net.bus.index),int)
res=net.res_bus.loc[ids]
Vpp=res.vm_pu.to_numpy()*np.exp(1j*np.radians(res.va_degree.to_numpy()))

pynet=load_network(repo/"configs/ias2026/ieee39_network.json")
py=solve_power_flow(pynet)
ids_ref=np.asarray(pynet.bus_idx,int)
if not np.array_equal(ids,ids_ref):
    order=[list(ids).index(b) for b in ids_ref]
    Vpp=Vpp[order]

# global angle gauge
off=np.median(np.angle(Vpp)-np.angle(py.voltages))
vm_err=float(np.max(np.abs(np.abs(Vpp)-np.abs(py.voltages))))
va_err=float(np.max(np.abs((np.angle(Vpp)-off)-np.angle(py.voltages))))

fpp=branch_flows_from_voltages(repo,Vpp)
fpy=branch_flows_from_voltages(repo,py.voltages)
flow_err=float(np.max(np.abs(fpp[:,3:7]-fpy[:,3:7])))

metrics={
    "pandapower_version":getattr(pp,"__version__","unknown"),
    "converged":bool(net.converged),
    "vm_max_abs_error_pu":vm_err,
    "va_max_gauge_aligned_error_rad":va_err,
    "canonical_branch_terminal_flow_max_error_MVA":flow_err,
}
pd.DataFrame(fpp,columns=["idx","from_bus","to_bus","P_from_MW","Q_from_Mvar","P_to_MW","Q_to_Mvar"]).to_csv(out/"pandapower_branch_flows.csv",index=False)
(out/"pandapower_parity.json").write_text(json.dumps(metrics,indent=2))
print(json.dumps(metrics,indent=2))

# Informative only: stock case39 is a different solved case, so do not use it as parity evidence.
try:
    from pandapower.networks import case39
    stock=case39()
    pp.runpp(stock, calculate_voltage_angles=True, init="results", enforce_q_lims=False)
    stock_metrics={
        "note":"INFORMATIVE_ONLY_stock_case39_is_not_the_frozen_ANDES_case",
        "stock_bus1_vm":float(stock.res_bus.loc[0,"vm_pu"]) if 0 in stock.res_bus.index else None,
        "canonical_bus1_vm":float(abs(py.voltages[0])),
    }
    (out/"pandapower_stock_case39_difference.json").write_text(json.dumps(stock_metrics,indent=2))
except Exception:
    pass
