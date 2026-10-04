"""Figures use saved measurements only; no new simulations."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
OUT=Path(__file__).resolve().parent
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
d=pd.read_csv(OUT/'TABLE_01_WALK_REMAINDER.csv');x=np.arange(len(d))
fig,ax=plt.subplots(figsize=(7.2,4.3),layout='constrained')
ax.bar(x-.23,d.scalar_norm_remainder,.23,label='Scalar norm tail',color='#8295a8')
ax.bar(x,d.walk_remainder_order1,.23,label='Closed-walk tail',color='#147d92')
ax.bar(x+.23,d.actual_error_order1,.23,label='Observed first-order error',color='#dc8c34')
ax.set(yscale='log',xticks=x,xticklabels=[str(v) for v in d.frequency_Hz],xlabel='Frequency [Hz]',ylabel='Trace-log magnitude (dimensionless)',title='Tighter finite-change bounds at three frozen IEEE-39 points')
ax.legend(frameon=False,loc='upper right',fontsize=8);ax.grid(axis='y',alpha=.15)
fig.text(.5,-.025,'Pointwise algebraic checks; these are not replacement-capacity bounds.',ha='center',fontsize=9)
fig.savefig(OUT/'FIG_01_WALK_BOUNDS.png',dpi=240,bbox_inches='tight');fig.savefig(OUT/'FIG_01_WALK_BOUNDS.svg',bbox_inches='tight');plt.close(fig)

r=json.loads((OUT/'CONTOUR_RESULT.json').read_text());panels=json.loads((OUT/'PANELS.json').read_text())
c=-.12229688588491999+.32408378548166134j
fig,ax=plt.subplots(figsize=(6.7,4.3),layout='constrained')
ax.add_patch(Rectangle((c.real-.04,c.imag-.04),.08,.08,facecolor='#d8eef1',edgecolor='#147d92',lw=2))
for p in panels:
    a,b=np.array(p['a']),np.array(p['b']);mid=(a+b)/2
    ax.plot(c.real+mid[0],c.imag+mid[1],'.',color='#147d92',ms=1.6)
ax.plot(c.real,c.imag,'x',color='#1a3248',label='Preserved numerical root estimate')
ax.axvline(-.05,color='#b04040',ls='--',label='Required margin')
ax.set(xlim=(-.18,-.035),ylim=(.275,.375),xlabel='Real part [1/s]',ylabel='Imaginary part [rad/s]',title='One continuous contour and its full parameter box')
ax.text(-.174,.369,f"{r['panel_count']} panels; verified q < {r['q_max_upper']+1e-6:.6f}",fontsize=9,va='top')
label='Exactly one root retained in this square' if r['count_one_proved'] else 'Boundary enclosed; root count unresolved'
ax.text(-.174,.280,label,fontsize=9);ax.legend(frameon=False,fontsize=8,loc='upper right')
fig.text(.5,-.025,'Exported model, uniform 40 ms delay; not a complete-spectrum certificate.',ha='center',fontsize=9)
fig.savefig(OUT/'FIG_02_CONTINUOUS_CONTOUR.png',dpi=240,bbox_inches='tight');fig.savefig(OUT/'FIG_02_CONTINUOUS_CONTOUR.svg',bbox_inches='tight');plt.close(fig)
print('Generated two PNG/SVG figures from saved data.')
