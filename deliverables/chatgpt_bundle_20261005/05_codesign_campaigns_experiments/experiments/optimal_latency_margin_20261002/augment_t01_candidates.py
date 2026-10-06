"""Join exact local root continuations to the complete margin-root counts."""
from pathlib import Path
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
counts=pd.read_csv(HERE/"T01_DELAY_TRANSITION_REPRODUCTION.csv")
fast=pd.read_csv(HERE/"T08_FAST_MODE_PROVENANCE.csv")
slow=pd.read_csv(HERE/"T01_SLOW_MODE_CANDIDATES.csv")
rows=[]
for rec in counts.itertuples():
    f=fast[(fast.design_id==rec.design_id)&np.isclose(fast.tau_ms,rec.tau_ms,atol=1e-7)]
    f=f[f.status=="LOCAL_EXACT_CHARACTERISTIC_CONTINUATION"]
    s=slow[(slow.design_id==rec.design_id)&np.isclose(slow.tau_ms,rec.tau_ms,atol=1e-7)]
    candidates=[]
    for root in f.itertuples():
        candidates.append((root.root_real,root.frequency_hz,f"fast_{int(root.family_id)}",
                           root.small_factor_sigma_min))
    for root in s.itertuples():
        if root.converged:
            candidates.append((root.root_real,root.frequency_hz,"slow_ODE",
                               root.residual))
    candidates.sort(key=lambda x:x[0],reverse=True)
    best=candidates[0] if candidates else (np.nan,np.nan,"UNKNOWN",np.nan)
    refined_violating=2*int((f.root_real>-0.05).sum())
    coverage=("COMPLETE_MARGIN_VIOLATIONS_NUMERICALLY" if
              rec.status=="UNSAFE" and refined_violating==rec.roots_violating_margin else
              "SAFE_NEEDS_RIGHTMOST_EXCLUSION" if rec.status=="SAFE" else
              "INDETERMINATE_ROOT_COVERAGE")
    rows.append({**rec._asdict(),
                 "rightmost_real_part_candidate":best[0],
                 "rightmost_frequency_hz_candidate":best[1],
                 "critical_family_id_candidate":best[2],
                 "critical_root_residual":best[3],
                 "refined_violating_root_count":refined_violating,
                 "rightmost_coverage_status":coverage})
pd.DataFrame(rows).to_csv(HERE/"T01_DELAY_TRANSITION_REPRODUCTION.csv",index=False)
print("AUGMENTED_T01",len(rows),"cases")
