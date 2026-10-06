"""Aggregate only executed cases; preserve failures and finite-horizon scope."""
from pathlib import Path
import json, hashlib, tomllib
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent; ROOT=OUT.parents[1]
OLD=ROOT/'reports/experiment_Q2B/CERTIFIED_SEARCH/PD'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def cumtrap(y,t): return np.r_[0,np.cumsum(np.diff(t)*(y[1:]+y[:-1])/2)]
def trap(y,t): return float(cumtrap(np.asarray(y),np.asarray(t))[-1])
protocol=json.loads((OUT/'PREREGISTRATION.json').read_text())
candidate=tomllib.loads((ROOT/protocol['candidate']).read_text())
rho=np.array(candidate['rho']); machine=pd.read_csv(ROOT/'reports/experiment_D/inputs/machine.csv').set_index('bus')
design=pd.read_csv(OUT/'PHYSICAL_TRIM_PQ.csv')[['bus','rho','SG_P_MW','GFL_P_MW','SG_Q_Mvar','GFL_Q_Mvar']]
design['Kp']=candidate['Kp'];design['Ki']=candidate['Ki'];design['pure_PLL_delay_s']=0.0
design.to_csv(OUT/'FIXED_DESIGN.csv',index=False)
gov=pd.read_csv(ROOT/'reports/experiment_D/inputs/gov.csv').set_index('bus')
results=pd.read_csv(OUT/'EVENT_RESULTS.csv')
parity=tomllib.loads((OUT/'PARITY.toml').read_text())
op=pd.read_csv(OUT/'LEGACY_POLES.csv');pp=pd.read_csv(OUT/'PHYSICAL_POLES.csv')
z=op.real.to_numpy()+1j*op.imag.to_numpy();w=pp.real.to_numpy()+1j*pp.imag.to_numpy()
i,j=linear_sum_assignment(abs(z[:,None]-w[None,:]))
pd.DataFrame({'legacy_real':z[i].real,'legacy_imag':z[i].imag,'physical_real':w[j].real,'physical_imag':w[j].imag,'error':abs(z[i]-w[j])}).to_csv(OUT/'POLE_MATCHING.csv',index=False)
parity['bijective_pole_matching_absolute']=float(max(abs(z[i]-w[j])))
comparisons=[];balances=[];quads=[];traces={};govrows=[]
for row in results.itertuples():
    if not row.complete: continue
    path=OUT/f'{row.model}_{row.event}_SENSORS.csv';d=pd.read_csv(path);t=d.time_s.to_numpy();traces[(row.model,row.event)]=d
    e=d.aggregate_energy_MJ.to_numpy();pin=d.net_input_MW.to_numpy();E=e-e[0]
    # Discontinuous power at t=1: use left limit on the interval ending at the jump.
    def integrate_jump(y,t):
        y=np.asarray(y);dt=np.diff(t);left=y[:-1].copy();right=y[1:].copy()
        right[np.isclose(t[1:],1)]=y[0]
        return np.r_[0,np.cumsum(.5*dt*(left+right))]
    integ=integrate_jump(pin,t);err=E-integ
    for stride in [1,2,4]:
        tt=t[::stride];yy=pin[::stride];ee=E[::stride]-integrate_jump(yy,tt)
        quads.append(dict(model=row.model,event=row.event,dt_s=.005*stride,final_error_MJ=ee[-1],max_error_MJ=max(abs(ee))))
    L=(d.load_MW+d.network_losses_MW+d.filter_losses_MW).to_numpy();deltaL=L-L[0]
    freqarea=np.zeros(len(t));boundary=np.zeros(len(t));awpower=np.zeros(len(t));mechint=np.zeros(len(t))
    for bus in range(30,40):
        if rho[bus-30]==1:continue
        nu=d[f'SG_Hz_bus{bus}'].to_numpy()/60
        pm=d[f'SG_mechanical_MW_bus{bus}'].to_numpy();mechint+=cumtrap(pm-pm[0],t)
        if bus==39:
            freqarea-=pm[0]*cumtrap(nu,t);continue
        sn=machine.loc[bus,'Sn']*(1-rho[bus-30]);g=gov.loc[bus]
        x1=d[f'gov_xg1_bus{bus}'].to_numpy();x2=d[f'gov_xg2_bus{bus}'].to_numpy()
        ref=x1[0]-nu/g.R
        limited=((x1>g.V_max)&(ref>x1))|((x1<g.V_min)&(ref<x1))
        aw=np.where(limited,x1-ref,0)
        assert abs(aw[0]) < 1e-12, 'Frozen initial governor is not interior'
        freqarea+=sn*(1/g.R+g.DT)*cumtrap(nu,t)
        boundary+=sn*((g.T2-g.T1)*(x1-x1[0])-g.T3*(x2-x2[0]));awpower+=sn*aw
        govrows.append(dict(model=row.model,event=row.event,bus=bus,limited_samples=int(sum(limited)),limited_sample_duration_s=float(sum(limited)*.005),max_valve=x1.max(),min_valve=x1.min(),antiwindup_integral_MJ=trap(sn*aw,t),power_output_identity_max_MW=float(max(abs(pm-sn*(x2-g.DT*nu))))))
    awarea=cumtrap(awpower,t);loadarea=integrate_jump(deltaL,t)
    pdc=sum(100*rho[b-30]*d[f'Pdc_pu_bus{b}'].to_numpy() for b in range(30,40))
    pdcarea=cumtrap(pdc-pdc[0],t)
    energy_area_rhs=boundary+awarea+pdcarea-loadarea-E
    areaerr=freqarea-energy_area_rhs
    gov_err=mechint-(-freqarea+boundary+awarea)
    output=pd.DataFrame({'time_s':t,'delta_energy_MJ':E,'integrated_net_input_MJ':integ,'energy_integral_error_MJ':err,'weighted_rotor_area_MJ':freqarea,'governor_terminal_MJ':boundary,'antiwindup_MJ':awarea,'delta_load_losses_integral_MJ':loadarea,'frequency_area_RHS_MJ':energy_area_rhs,'frequency_area_error_MJ':areaerr,'governor_integral_error_MJ':gov_err})
    output.to_csv(OUT/f'{row.model}_{row.event}_BALANCE.csv',index=False)
    balances.append(dict(model=row.model,event=row.event,energy_RHS_error_MW=row.energy_RHS_error_MW,energy_integral_final_error_MJ=err[-1],energy_integral_max_error_MJ=max(abs(err)),frequency_area_MJ=freqarea[-1],frequency_area_final_error_MJ=areaerr[-1],frequency_area_max_error_MJ=max(abs(areaerr)),governor_final_error_MJ=gov_err[-1],antiwindup_MJ=awarea[-1],terminal_governor_MJ=boundary[-1],delta_energy_MJ=E[-1],load_losses_integral_MJ=loadarea[-1],realized_tail_delta_load_MW=float(d.load_MW.iloc[-1]-d.load_MW.iloc[0])))
    ref=pd.read_csv(OLD/f'{row.event}_SENSORS.csv')
    assert len(ref)==len(d) and np.max(abs(ref.time_s-d.time_s))<1e-10
    def maxdiff(stem): return max(np.max(abs(d[f'{stem}_bus{b}']-ref[f'{stem}_bus{b}'])) for b in range(30,40))
    oldf=max(np.max(abs(ref[f'grid_Hz_bus{b}'])) for b in range(30,40))
    oldr=max(np.max(abs(ref[f'rocof_Hz_s_bus{b}'])) for b in range(30,40))
    comparisons.append(dict(model=row.model,event=row.event,legacy_archived_Fpeak_Hz=oldf,new_Fpeak_Hz=row.Fpeak,legacy_archived_Rpeak_Hz_s=oldr,new_Rpeak_Hz_s=row.Rpeak,max_grid_difference_Hz=maxdiff('grid_Hz'),max_RoCoF_difference_Hz_s=maxdiff('rocof_Hz_s'),max_GFL_power_difference_MW=maxdiff('GFL_MW'),max_vdc_reflection_error_pu=max(np.max(abs(d[f'Vdc_pu_bus{b}']+ref[f'Vdc_pu_bus{b}']-5)) for b in range(30,40)) if row.model=='physical' else np.nan))
pd.DataFrame(comparisons).to_csv(OUT/'PAIRED_COMPARISON.csv',index=False)
pd.DataFrame(balances).to_csv(OUT/'ENERGY_AND_FREQUENCY_AREA.csv',index=False)
pd.DataFrame(quads).to_csv(OUT/'QUADRATURE_SENSITIVITY.csv',index=False)
pd.DataFrame(govrows).to_csv(OUT/'GOVERNOR_BRANCHES.csv',index=False)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
phys=results[results.model=='physical'];fig,ax=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
x=np.arange(len(phys));labels=[v.replace('_',' + ').replace('MW',' MW') for v in phys.event]
ax[0].bar(x-.17,phys.Fpeak/.5,width=.32,label='Frequency / 0.5 Hz',color='#087e8b');ax[0].bar(x+.17,phys.Rpeak/.5,width=.32,label='RoCoF / (0.5 Hz/s)',color='#dd7b27')
ax[0].axhline(1,color='#ab2430',linestyle='--',label='Declared limit')
ax[0].set_xticks(x,labels,rotation=17,ha='right');ax[0].legend(fontsize=9);ax[0].set_title('Fixed design: no retuning');ax[0].set_ylabel('Fraction of each declared limit')
for ev,col in [('bus16_100MW','#087e8b'),('bus8_100MW','#dd7b27'),('bus29_100MW','#ab2430')]:
    d=traces.get(('physical',ev))
    if d is None: continue
    b=pd.read_csv(OUT/f'physical_{ev}_BALANCE.csv')
    ax[1].plot(b.time_s,b.delta_energy_MJ,color=col,label=ev.replace('_',' + '))
    ax[1].plot(b.time_s[::100],b.integrated_net_input_MJ[::100],color=col,linestyle=':',linewidth=3,alpha=.6)
ax[1].set(xlabel='Time [s]',ylabel='Energy change [MJ]',title='Positive storage = integrated net power')
ax[1].legend(fontsize=9);fig.suptitle('IEEE-39: physical energy closes; security remains event dependent',fontweight='bold',fontsize=13)
fig.savefig(OUT/'FIG_PHYSICAL_BRIDGE.png',dpi=220);fig.savefig(OUT/'FIG_PHYSICAL_BRIDGE.svg');plt.close(fig)
fig,ax=plt.subplots(1,2,figsize=(11,4.3),layout='constrained')
for ev,col in [('bus16_100MW','#087e8b'),('bus8_100MW','#dd7b27'),('bus29_100MW','#ab2430')]:
    if ('physical',ev) not in traces:continue
    b=pd.read_csv(OUT/f'physical_{ev}_BALANCE.csv');ax[0].plot(b.time_s,b.weighted_rotor_area_MJ,color=col,label=ev)
    ax[0].plot(b.time_s[::100],b.frequency_area_RHS_MJ[::100],color=col,linestyle=':',linewidth=3,alpha=.6)
    ax[1].plot(b.time_s,b.frequency_area_error_MJ,color=col,label=ev)
ax[0].set(xlabel='Time [s]',ylabel='Weighted rotor-frequency area [MJ]',title='Exact finite-time identity, numerical quadrature')
ax[1].set(xlabel='Time [s]',ylabel='Area identity residual [MJ]',title='Residuals retained, not called a certificate')
ax[0].legend(fontsize=9);fig.savefig(OUT/'FIG_FREQUENCY_AREA_BRIDGE.png',dpi=220);fig.savefig(OUT/'FIG_FREQUENCY_AREA_BRIDGE.svg');plt.close(fig)
all_physical_complete=len(phys)==len(protocol['events']) and phys.complete.all()
status={'status':'EXECUTED_LIMITED_VALIDATION' if all_physical_complete else 'INCOMPLETE_CAMPAIGN','parity':parity,'physical_event_count':int(len(phys)),'physical_frequency_and_rocof_pass_count':int(sum(phys.complete & phys.frequency_pass & phys.rocof_pass)), 'max_energy_RHS_error_MW':float(phys.energy_RHS_error_MW.max()),'positive_energy_gate_pass':bool(all_physical_complete and np.isfinite(phys.energy_RHS_error_MW).all() and (phys.energy_RHS_error_MW<1e-6).all()),'fixed_GFL_MW':candidate['converted_GFL_MW'],'fixed_GFL_percent':100*candidate['GFL_fraction'],'retained_SG_MW':candidate['J_MW'],'new_optimization':False,'replacement_optimum_claim':False,'lossless_graph_floor_transfer':'BLOCKED_UNBOUNDED_FULL_MODEL_DEFECTS','delay_validation':False,'poster_main_claim_ready':False,'all_frozen_inputs_unchanged':all(sha(ROOT/p)==h for p,h in protocol['inputs'].items()),'generated_model_unchanged':sha(OUT/'PhysicalBridge.jl')==protocol['generated_model_sha256']}
(OUT/'RESULTS.json').write_text(json.dumps(status,indent=2),encoding='utf-8')
print(json.dumps(status,indent=2));print(pd.DataFrame(comparisons).to_string(index=False));print(pd.DataFrame(balances).to_string(index=False))
