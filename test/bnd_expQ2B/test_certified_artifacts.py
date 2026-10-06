"""Integration acceptance gate over independently produced design/PD evidence."""
from pathlib import Path
import hashlib, json, tomllib
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'reports/experiment_Q2B/CERTIFIED_SEARCH'
def read(p):return tomllib.loads((OUT/p).read_text())
c=read('Z_LOCAL_SECURE_FINAL.toml');k=read('LOCAL_CERTIFICATE_FEASIBLE.toml')
assert hashlib.sha256((OUT/'Z_LOCAL_SECURE_FINAL.toml').read_bytes()).hexdigest()==(OUT/'Z_LOCAL_SECURE_FINAL.toml.sha256').read_text().strip()
assert k['local_KKT_certificate'] and k['primal']<2e-8 and k['stationarity']<1e-6
assert k['complementarity']<1e-8 and k['LICQ_rank']==k['active_count']
assert min(k['multipliers'])>0
assert min(k['projected_Hessian_eigenvalues'])>10*k['Hessian_step_variation']
for p in ('RICCATI_FULL_BAND_CERTIFICATE.toml','PD/RICCATI_FULL_BAND_CERTIFICATE.toml'):
    r=read(p);assert r['X_positive'] and r['strict_negative_LMI']
    assert r['beta_lower_tested']>=c['beta_requirement']
    assert r['Riccati_residual_big']<r['negative_LMI_min_eigenvalue_Float64']/1e6
assert read('ALL_TIME_CERTIFICATE.toml')['pass']
p=read('PD/SPECTRUM.toml');assert p['margin_pass'] and p['alpha_error']<1e-5
assert p['max_P_error_pu']<1e-9 and p['max_Q_error_pu']<1e-9
assert pd.read_csv(OUT/'PD/BIJECTIVE_POLE_MATCHING.csv').error.max()<1e-5
events=pd.read_csv(OUT/'PD/EVENTS.csv');target=events[events.event=='bus16_100MW'].iloc[0]
assert target.frequency_pass and target.rocof_pass and target.Pdc_change<1e-12
assert target.power_balance_error_MW<1e-6
gap=json.loads((OUT/'GLOBAL_GAP.json').read_text());assert gap['lower_MW']==0 and not gap['global_certified']
assert pd.read_csv(OUT/'IMMUTABLE_INPUT_HASH_AUDIT.csv').passed.all()
print('PASS: frozen hash, numerical KKT/LICQ/SOSC, all-frequency robustness, all-time bounds, independent PD spectrum/trim/event, global-claim discipline')
