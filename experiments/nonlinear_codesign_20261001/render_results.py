"""Render only generated physical-model search data; never insert claimed optima."""
import argparse,csv,json,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/'reports/nonlinear_codesign_20261001'
pa=argparse.ArgumentParser();pa.add_argument('candidate_json');args=pa.parse_args()
path=Path(args.candidate_json);d=json.loads(path.read_text());folder=path.parent
weights=np.array([float(r['P_gen_MW']) for r in csv.DictReader((ROOT/'reports/experiment_N/TABLE_N01_original_operating_point.csv').open())])
matches=[]
for f in folder.glob('*_eval.csv'):
    rows=list(csv.DictReader(f.open()))
    if len(rows)==6 and all(r['complete']=='true' for r in rows) and all(abs(float(r['retained_MW'])-d['retained_MW'])<1e-6 for r in rows):matches.append((f,rows))
if not matches:raise RuntimeError('No complete matching physical evaluation for the saved candidate')
case_file,rows=sorted(matches)[-1]
history_paths=[folder/'history.json']
if folder.name=='search_full_horizon':
    history_paths=[BASE/'search_physical/history.json',BASE/'search_physical_fast/history.json',folder/'history.json']
history=[h for file in history_paths if file.is_file() for h in json.loads(file.read_text())]
rho=np.array(d['rho']);kp=np.array(d['Kp']);ki=np.array(d['Ki']);buses=np.arange(30,40)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                     'axes.labelcolor':'#17372D','text.color':'#17372D','axes.titleweight':'bold'})
fig,axs=plt.subplots(2,2,figsize=(9,10.8),gridspec_kw={'hspace':.42,'wspace':.31})
fig.patch.set_facecolor('#FAFAF5')
fig.suptitle('BEYOND NODAL DAMPING',x=.09,y=.966,ha='left',fontsize=21,color='#003B2D',fontweight='bold')
fig.text(.09,.926,'Nonlinear SG-to-GFL co-design — preliminary IEEE-39 instance',fontsize=12)
fig.text(.09,.886,f'{d["replacement_percent"]:.3f}% GFL replacement   |   {d["retained_MW"]:.2f} MW SG retained',fontsize=15,color='#006B4A',weight='bold')
ax=axs[0,0]
for j,h in enumerate(history,1):
    value=100*(1-h['retained_MW']/weights.sum())
    ax.scatter(j,value,c='#006B4A' if h['accepted'] else '#D6292E',marker='o' if h['accepted'] else 'x',s=42)
ax.set(xlabel='Nonlinear trial',ylabel='Replacement (%)',title='A  Accept the trajectory, not just the poles')
ax.plot([],[],color='#006B4A',marker='o',ls='',label='Accepted');ax.plot([],[],color='#D6292E',marker='x',ls='',label='Rejected')
ax.legend(frameon=False,fontsize=9);ax.grid(alpha=.16)
baseline_file=BASE/'uniform_baseline/summary.csv'
if baseline_file.is_file():
    uniform=max(float(r['rho']) for r in csv.DictReader(baseline_file.open()) if r['feasible']=='true')
    ax.axhline(100*uniform,color='#777777',ls=':',lw=1)
    ax.text(.98,.08,f'Uniform / nominal gains: {100*uniform:.2f}%',transform=ax.transAxes,
            ha='right',fontsize=8,color='#555555')
ax=axs[0,1]
f=np.array([float(r['Fpeak_Hz']) for r in rows]);r=np.array([float(r['Rpeak_Hz_s']) for r in rows])
xx=np.arange(6);cfg=d.get('problem_instance',{})
ax.bar(xx-.18,f/cfg.get('frequency_limit_Hz',.5),.35,color='#1768AC',label='Frequency')
ax.bar(xx+.18,r/cfg.get('rocof_limit_Hz_s',.5),.35,color='#C99A20',label='RoCoF')
ax.axhline(1,color='#D6292E',ls='--',lw=1)
ax.set_xticks(xx,[f'{a["bus"]}\n{"+" if float(a["delta"])>0 else "−"}100' for a in rows])
ax.set(xlabel='Load bus / commanded step (MW)',ylabel='Peak / declared limit',title='B  Six nonlinear disturbance cases',ylim=(0,max(1.25,float(max(f/.5))+.15)))
ax.legend(frameon=False,fontsize=9)
ax=axs[1,0]
ax.bar(buses,weights*(1-rho),color='#006B4A');ax.set_xticks(buses)
ax.set(xlabel='Generator bus',ylabel='Retained synchronous generation (MW)',title='C  Location is a design variable')
ax.grid(axis='y',alpha=.16)
ax=axs[1,1]
ax.plot(buses,kp/(10*np.pi),'o-',color='#1768AC',label=r'$K_p/K_{p,0}$')
ax.plot(buses,ki/((10*np.pi)**2/4),'s-',color='#C99A20',label=r'$K_i/K_{i,0}$')
ax.set_xticks(buses);ax.set(xlabel='Converter bus',ylabel='PLL gain / initial gain',title='D  Joint local gain tuning')
ax.set_ylim(0,1.12*max(1.,max(kp/(10*np.pi)),max(ki/((10*np.pi)**2/4))))
ax.ticklabel_format(axis='y',style='plain',useOffset=False)
ax.legend(frameon=False);ax.grid(alpha=.16)
fig.subplots_adjust(left=.09,right=.97,top=.81,bottom=.17)
fig.text(.09,.113,'Physical DC-supply convention • all 39 buses monitored • 60 s • 0.5 s measurement window',fontsize=9)
fig.text(.09,.083,'Best feasible found; global optimality and a nonlinear region-of-attraction certificate are not established.',fontsize=9)
fig.text(.09,.055,'Voltage band and actuator margin are declared experiment settings. Equipment current / energy limits are pending.',fontsize=8.5)
for ext in ('png','pdf','svg'):fig.savefig(BASE/f'nonlinear_codesign_result.{ext}',dpi=180,facecolor=fig.get_facecolor())
table=[]
for i,b in enumerate(buses):table.append({'bus':int(b),'rho':rho[i],'retained_SG_MW':weights[i]*(1-rho[i]),'converted_GFL_MW':weights[i]*rho[i],'Kp':kp[i],'Ki':ki[i]})
with (BASE/'physical_candidate_parameters.csv').open('w',newline='') as fp:
    w=csv.DictWriter(fp,fieldnames=table[0]);w.writeheader();w.writerows(table)
manifest={'candidate_json':str(path.resolve()),'candidate_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
          'history_files':[str(f.resolve()) for f in history_paths],
          'evaluation_csv':str(case_file.resolve()),'evaluation_sha256':hashlib.sha256(case_file.read_bytes()).hexdigest(),
          'replacement_percent':d['replacement_percent'],'retained_MW':d['retained_MW'],
          'max_frequency_Hz':float(max(f)),'max_RoCoF_Hz_s':float(max(r)),
          'voltage_min_pu':min(float(z['Vmin']) for z in rows),'voltage_max_pu':max(float(z['Vmax']) for z in rows),
          'maximum_current_to_initial_ratio':max(float(z['current_ratio_max']) for z in rows),
          'dc_voltage_min':min(float(z['DC_min']) for z in rows),'dc_voltage_max':max(float(z['DC_max']) for z in rows),
          'minimum_actuator_fraction_slack':min(float(z['limiter_fraction']) for z in rows),
          'nominal_alpha':.05*(d['constraints'][0]-1),'global_optimality_certified':False}
(BASE/'physical_candidate_summary.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(manifest,indent=2))
