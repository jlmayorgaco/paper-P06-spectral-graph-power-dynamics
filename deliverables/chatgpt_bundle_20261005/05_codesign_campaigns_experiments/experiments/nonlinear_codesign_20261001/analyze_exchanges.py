import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'reports/nonlinear_codesign_20261001/exchange_checkpoint'
scope=json.loads((OUT/'scope.json').read_text());active=np.array(scope['active_rows_zero_based'])
rows=list(csv.DictReader((OUT/'validation.csv').open()));c0=np.fromstring(rows[0]['constraints'],sep=';')
result=[]
for row in rows[1:]:
    c=np.fromstring(row['constraints'],sep=';');shift=c[active]-c0[active]
    result.append({k:row[k] for k in ['id','increase_bus','decrease_bus','MW','tuned','feasible']}|
        {'actual_signature_change_norm':float(np.linalg.norm(shift)),
         'maximum_constraint':float(max(c)),'limiting_constraint_zero_based':int(np.argmax(c))})
with (OUT/'finite_exchange_results.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=result[0]);w.writeheader();w.writerows(result)
fig,ax=plt.subplots(figsize=(8,4.3));xx=np.arange(4)
fixed=np.array([r['actual_signature_change_norm'] for r in result if r['tuned']=='false'])
tuned=np.array([r['actual_signature_change_norm'] for r in result if r['tuned']=='true'])
ax.bar(xx-.18,fixed,.35,label='Fixed PLL gains',color='#C99A20')
ax.bar(xx+.18,tuned,.35,label='Analytical-gradient gain compensation',color='#006B4A')
ax.set_yscale('log');ax.set_xticks(xx,['31 / 32\n0.1 MW','31 / 32\n1 MW','35 / 36\n0.1 MW','35 / 36\n1 MW'])
ax.set_ylabel('Measured change in selected constraint vector')
ax.set_xlabel('Increase GFL / decrease GFL, equal MW exchanged')
ax.set_title('Node exchangeability depends on control authority')
ax.legend(frameon=False,fontsize=8);ax.spines[['top','right']].set_visible(False)
fig.tight_layout();fig.savefig(OUT/'finite_exchange_comparison.png',dpi=180)
fig.savefig(OUT/'finite_exchange_comparison.pdf')
print(json.dumps(result,indent=2))
