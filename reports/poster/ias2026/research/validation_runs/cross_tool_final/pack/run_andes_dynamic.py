
from __future__ import annotations
import sys, json, subprocess, tempfile
from pathlib import Path
import pandas as pd

try:
    import andes
except Exception as e:
    raise SystemExit(f"ANDES unavailable: {e}")

repo=Path(sys.argv[1]).resolve()
out=Path(__file__).resolve().parent/"results"
out.mkdir(exist_ok=True)

# 1) live static reference
andes.config_logger(stream_level=40)
case=Path(andes.__file__).resolve().parent/"cases/ieee39/ieee39_full.xlsx"
ss=andes.load(str(case), setup=True, no_output=True)
ss.PFlow.run()
static={
    "andes_version":getattr(andes,"__version__","unknown"),
    "case":str(case),
    "pf_converged":bool(ss.PFlow.converged),
    "bus_v":[float(x) for x in ss.Bus.v.v],
    "bus_a":[float(x) for x in ss.Bus.a.v],
}
(out/"andes_live_static.json").write_text(json.dumps(static,indent=2))

# 2) independent equation-equivalent dynamic worker already maintained by the project
worker=repo/"experiments/F1_andes_equivalent_worker.py"
target=out/"andes_equivalent.csv"
with tempfile.TemporaryDirectory(prefix="andes_eq_") as work:
    cp=subprocess.run([sys.executable,str(worker),str(target),work],
                      cwd=repo,text=True,capture_output=True)
if cp.returncode:
    print(cp.stdout); print(cp.stderr,file=sys.stderr)
    raise SystemExit(cp.returncode)

# compare new ANDES result with current internal Python columns from the frozen reconciliation
internal=pd.read_csv(repo/"results/F1_eigenvalue_reconciliation.csv")
fresh=pd.read_csv(target)
merged=internal.merge(fresh[["stage","members","alpha_band","freq_band_hz","band_rhp_count"]],
                      on=["stage","members"],suffixes=("_internal","_andes_new"))
merged["alpha_error"]=abs(merged.internal_alpha_band-merged.alpha_band)
merged["freq_error_hz"]=abs(merged.internal_freq_band_hz-merged.freq_band_hz)
merged["rhp_agrees"]=merged.internal_band_rhp==merged.band_rhp_count
merged.to_csv(out/"andes_dynamic_reconciliation.csv",index=False)

r3=merged[merged.stage.eq("R3 first-order AVR")]
summary={
    "andes_version":getattr(andes,"__version__","unknown"),
    "R3_max_alpha_error":float(r3.alpha_error.max()),
    "R3_max_frequency_error_hz":float(r3.freq_error_hz.max()),
    "R3_all_RHP_counts_agree":bool(r3.rhp_agrees.all()),
}
(out/"andes_dynamic_summary.json").write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
