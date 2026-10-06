"""T11: does the transport rule enlarge the feasible (inverter share, delay) region?  Spectral screen only.
For uniform share rho: count roots right of the -0.05 margin (gauge excluded) for baseline gains and for the rule-selected transport."""
import os, sys, pathlib
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np, pandas as pd
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'experiments' / 'interaction_decision_20261004'))
from model import Model, P0, TAU
from t3_full_spectrum_count import count
from t8_taylor_and_second_design import select_node, gains
rows = []
for rho in [0.80, 0.85, 0.875, 0.90, 0.925, 0.95]:
    p0 = P0.copy(); p0[:10] = rho
    base = count(Model(p0, TAU), -0.05)[0] - 1
    row = dict(rho=rho, base_40ms=base); sel = None
    if base == 0:
        try: sel = select_node(p0)[0]
        except Exception as e: row['error'] = str(e)[:60]
    row['node_Hz'] = sel[0].imag / 2 / np.pi if sel else np.nan
    for t in [44, 48]:
        h = (t - 40) * 1e-3
        row[f'none_{t}'] = count(Model(p0, TAU + h), -0.05)[0] - 1 if base == 0 else np.nan
        row[f'rule_{t}'] = count(Model(gains(p0, 'exact', sel[0], h), TAU + h), -0.05)[0] - 1 if sel else np.nan
        row[f'Ki_{t}'] = gains(p0, 'exact', sel[0], h)[20] if sel else np.nan
    rows.append(row); print(row, flush=True); pd.DataFrame(rows).to_csv('T11_share_delay_map.csv', index=False)
