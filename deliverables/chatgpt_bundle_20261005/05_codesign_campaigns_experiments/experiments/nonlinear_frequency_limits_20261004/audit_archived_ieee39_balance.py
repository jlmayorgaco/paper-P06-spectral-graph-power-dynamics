"""Read-only audit of seven frozen, zero-pure-delay Q2B trajectories.

Run from any directory: python path/to/audit_archived_ieee39_balance.py
Writes only beside this script. No simulator or model package is imported.
The as-executed DC sign is deliberately retained and contrasted with positive
physical capacitor storage. Numerical differentiation/integration is diagnostic,
not a new ODE solution or a validated-arithmetic certificate.
"""
from pathlib import Path
import hashlib
import json
import platform
import tomllib
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
Q = ROOT / "reports/experiment_Q2B/CERTIFIED_SEARCH"
PD = Q / "PD"
BUSES = np.arange(30, 40)
SBASE, FBASE, RF, XF, CDC, WINDOW = 100., 60., .01, .03, 1.25, .5
WBASE = 2*np.pi*FBASE


def readcsv(path):
    return pd.read_csv(path, float_precision="round_trip")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fields(df, stem):
    return df[[f"{stem}_bus{b}" for b in BUSES]].to_numpy(float)


def integral(t, y, jumps=()):
    """Trapezoids, replacing the right sample by a left extrapolation before jumps.

    The stored event sample is right-consistent. It must not be used as the
    left-side endpoint of the interval preceding the event. At turn-off, use
    the last two preceding samples for a one-sided linear extrapolation.
    This approximation is separately checked by coarsening the stored samples.
    """
    y = np.asarray(y)
    dt = np.diff(t)
    inc = .5*(y[1:]+y[:-1])*dt.reshape((-1,)+(1,)*(y.ndim-1))
    for jump in jumps:
        ix = np.flatnonzero(np.isclose(t, jump, rtol=0, atol=1e-10))
        if len(ix) and ix[0] >= 2:
            i = ix[0]
            left = y[i-1]+(y[i-1]-y[i-2])*(t[i]-t[i-1])/(t[i-1]-t[i-2])
            inc[i-1] = .5*(y[i-1]+left)*(t[i]-t[i-1])
    return np.concatenate((np.zeros_like(y[:1]), np.cumsum(inc, axis=0)))


def slope(y, x, mask):
    xx, yy = np.asarray(x)[mask].ravel(), np.asarray(y)[mask].ravel()
    return float(xx@yy/(xx@xx))


def rms(x):
    return float(np.sqrt(np.mean(np.asarray(x)**2)))


def main():
    candidate = Q / "Z_LOCAL_SECURE_FINAL.toml"
    expected = (Q / "Z_LOCAL_SECURE_FINAL.toml.sha256").read_text().strip()
    assert sha(candidate) == expected
    design = tomllib.loads(candidate.read_text())
    rho = np.array(design["rho"])
    op = readcsv(ROOT / "reports/experiment_N/TABLE_N01_original_operating_point.csv").set_index("bus").loc[BUSES]
    machine = readcsv(ROOT / "reports/experiment_D/inputs/machine.csv").set_index("bus").loc[BUSES]
    governors = readcsv(ROOT / "reports/experiment_D/inputs/gov.csv").set_index("bus")
    sn = (1-rho)*op.Sn_original_MVA.to_numpy(float)
    h = op.H_seconds.to_numpy(float)
    damping = machine.D.to_numpy(float)
    rs = machine.R_s.to_numpy(float)
    m = 2*h*sn
    active = sn > 0
    assert np.all(rho > 0) and np.all(damping == 0) and np.all(rs == 0)
    assert np.array_equal(BUSES[active], np.array(design["support"]))
    frozen = json.loads((ROOT / "reports/experiment_N/MODEL_FREEZE.json").read_text())
    paths = ["src/pd39/model.jl", "src/bnd_model_expN/PDExactDesignN.jl",
             "src/bnd_model_expN/PDReferenceN.jl", "src/bnd_design/AnalyticSG.jl"]
    source_hashes = []
    for rel in paths:
        actual = sha(ROOT/rel)
        source_hashes.append(dict(path=rel, expected=frozen["inputs"][rel], actual=actual,
                                  pass_exact=actual == frozen["inputs"][rel]))
    assert all(x["pass_exact"] for x in source_hashes)
    assert "(p_ac-P_dc)/v_dc_state" in (ROOT/paths[0]).read_text()
    assert "(pac-p.Pdc)/(p.Cdc*vdc)" in (ROOT/paths[1]).read_text()

    summaries, node_rows, quadrature_rows, traces, input_files, gov_rows = [], [], [], [], [], []
    for path in sorted(PD.glob("*_SENSORS.csv")):
        event = path.name.removesuffix("_SENSORS.csv")
        df = readcsv(path)
        t = df.time_s.to_numpy(float)
        dt = float(np.median(np.diff(t)))
        assert abs(dt-.005) < 1e-10
        jumps = [1.,1.1] if event == "bus16_pulse" else [1.]
        input_files.append(dict(path=str(path.relative_to(ROOT)).replace("\\","/"),
                                sha256=sha(path), rows=len(df), time_start=float(t[0]),
                                time_end=float(t[-1]), sample_period=dt))
        v = fields(df, "V_pu")
        pg, qg = fields(df,"GFL_MW"), fields(df,"GFL_Mvar")
        ps, qs, pm = fields(df,"SG_MW"), fields(df,"SG_Mvar"), fields(df,"SG_mechanical_MW")
        dc, pdc_pu = fields(df,"Vdc_pu"), fields(df,"Pdc_pu")
        sf, gf = fields(df,"SG_Hz"), fields(df,"grid_Hz")
        omega = 1+np.nan_to_num(sf, nan=0.)/FBASE
        # Rotations are orthogonal: P,Q,V recover the squared internal current.
        current2 = (pg**2+qg**2)/((SBASE*rho)**2*v**2)
        ef_node = SBASE*rho*XF/(2*WBASE)*current2
        edc_node = SBASE*rho*CDC/2*dc**2
        lf_node = SBASE*rho*RF*current2
        pdc_node = SBASE*rho*pdc_pu
        copper_node = np.zeros_like(ps)
        copper_node[:,active] = rs[active]/sn[active]*(ps[:,active]**2+qs[:,active]**2)/v[:,active]**2
        damping_node = sn*damping*omega*(omega-1)
        kin_node = .5*m*omega**2
        kin, ef, edc = kin_node.sum(1), ef_node.sum(1), edc_node.sum(1)
        losses_filter, pdc = lf_node.sum(1), pdc_node.sum(1)
        p_mech, p_sg, p_gfl = pm.sum(1), ps.sum(1), pg.sum(1)
        load, lossnet = df.load_MW.to_numpy(float), df.network_losses_MW.to_numpy(float)
        sg_rhs = p_mech-p_sg-copper_node.sum(1)-damping_node.sum(1)
        all_rhs = p_mech+pdc-load-lossnet-losses_filter-copper_node.sum(1)-damping_node.sum(1)
        executed_energy = kin+ef-edc
        positive_storage = kin+ef+edc
        d_kin = np.gradient(kin,t,edge_order=2)
        d_ef_node = np.gradient(ef_node,t,axis=0,edge_order=2)
        d_dc_node = np.gradient(edc_node,t,axis=0,edge_order=2)
        dc_rhs_node = pg+lf_node+d_ef_node-pdc_node
        d_executed = np.gradient(executed_energy,t,edge_order=2)
        # Torque/MW-at-nominal-speed balance: no coherent-frequency approximation.
        torque_rhs = ((pm-ps-copper_node)/omega-sn*damping*(omega-1)).sum(1)
        coi = (omega*m).sum(1)/m.sum()
        coi_derivative = np.gradient(coi,t,edge_order=2)*m.sum()
        inc_kin = (.5*m*(omega-1)**2).sum(1)
        dispersion = (.5*m*(omega-coi[:,None])**2).sum(1)
        kinetic_identity_error = kin-kin[0]-m.sum()*(coi-coi[0])-(inc_kin-inc_kin[0])
        metric_phase = np.zeros_like(gf)
        lag = int(round(WINDOW/dt))
        for j in range(len(t)):
            metric_phase[j] = 2*np.pi*WINDOW*gf[j]+(metric_phase[j-lag] if j>=lag else 0)
        f_area = integral(t,gf)
        boundary_area = np.trapezoid(metric_phase[-lag-1:],t[-lag-1:],axis=0)/(2*np.pi*WINDOW)
        smooth = (t>=.02)&(t<=t[-1]-.02)
        dcfit = (t>=1.05)&(t<=t[-1]-.02)
        for j in jumps:
            smooth &= np.abs(t-j) > .025
            dcfit &= np.abs(t-j) > .05
        post = t>=1
        tail = t>=51
        int_all = integral(t,all_rhs,jumps)
        int_sg = integral(t,sg_rhs,jumps)
        int_gfl_rhs = integral(t,pdc-p_gfl-losses_filter,jumps)
        residual_exec = executed_energy-executed_energy[0]-int_all
        residual_positive = positive_storage-positive_storage[0]-int_all
        residual_sg = kin-kin[0]-int_sg
        residual_gfl = (ef-ef[0])-(edc-edc[0])-int_gfl_rhs
        # Endpoint identity for interior TGOV1, with hidden xg1 reconstructed
        # from xg2=(Pm/Sn)+DT*nu and the second governor ODE. This is a
        # data-derived reconstruction, not a new simulated trajectory.
        gov_boundary = 0.
        gov_antiwindup = 0.
        weighted_frequency_area = 0.
        gov_ode_max_rms = 0.
        inferred_clipped_samples=0
        for k,b in enumerate(BUSES):
            if not active[k]:
                continue
            nu = omega[:,k]-1
            area_nu = float(integral(t,nu)[-1])
            if b == 39:
                weighted_frequency_area -= pm[0,k]*area_nu
                continue
            g = governors.loc[b]
            r,t1,t2,t3,dd = (float(g[z]) for z in ["R","T1","T2","T3","DT"])
            pref = float(op.loc[b,"governor_p_ref_pu"])
            ref = (pref-nu)/r
            x2 = pm[:,k]/sn[k]+dd*nu
            dx2 = np.gradient(x2,t,edge_order=2)
            x1_interior = (t3*dx2+x2-(t2/t1)*ref)/(1-t2/t1)
            # Invalidate the smooth branch using the implemented derivative
            # limiter. On a clipped branch dxg1=0 => xg1=T3*dxg2+xg2.
            # Branch classification from differentiated sensors is diagnostic.
            x1_sat=t3*dx2+x2
            # The inverse of the output equation is not globally one-to-one.
            # Select only samples visibly on the valve boundary (1e-4 pu
            # diagnostic tolerance); ambiguous switch samples remain interior.
            # This does not constitute direct observation of the hidden state.
            hi=(x1_interior>float(g.V_max))&(ref>x1_sat)&(np.abs(x1_sat-float(g.V_max))<1e-4)
            lo=(x1_interior<float(g.V_min))&(ref<x1_sat)&(np.abs(x1_sat-float(g.V_min))<1e-4)
            clipped=hi|lo
            x1=np.where(clipped,x1_sat,x1_interior)
            rhs1 = (ref-x1)/t1
            rhs1=np.where(clipped,0.,rhs1)
            x1_residual = np.gradient(x1,t,edge_order=2)-rhs1
            aw=np.where(clipped,x1-ref,0.)
            aw_integral=sn[k]*float(integral(t,aw)[-1])
            endpoint = sn[k]*((t2-t1)*(x1[-1]-x1[0])-t3*(x2[-1]-x2[0]))
            droop_area = sn[k]*(1/r+dd)*area_nu
            pm_integral = float(integral(t,pm[:,k]-pm[0,k])[-1])
            gov_boundary += endpoint
            gov_antiwindup += aw_integral
            inferred_clipped_samples+=int(clipped.sum())
            weighted_frequency_area += droop_area
            gov_ode_max_rms=max(gov_ode_max_rms,rms(x1_residual[smooth]))
            gov_rows.append(dict(event=event,bus=int(b),frequency_weight_MW_per_Hz=sn[k]*(1/r+dd)/FBASE,
                reconstructed_xg1_min=float(x1.min()),reconstructed_xg1_max=float(x1.max()),
                interior_formula_xg1_max=float(x1_interior.max()),
                inferred_clipped_samples=int(clipped.sum()),
                first_inferred_clipped_s=float(t[clipped][0]) if clipped.any() else None,
                last_inferred_clipped_s=float(t[clipped][-1]) if clipped.any() else None,
                valve_limit_min=float(g.V_min),valve_limit_max=float(g.V_max),
                reconstructed_first_ODE_rms_pu_s=rms(x1_residual[smooth]),
                integrated_mechanical_delta_MJ=pm_integral,
                governor_endpoint_MJ=endpoint,negative_weighted_frequency_area_MJ=-droop_area,
                antiwindup_integral_MJ=aw_integral,
                interior_governor_integral_identity_error_MJ=pm_integral+droop_area-endpoint,
                hybrid_governor_integral_identity_error_MJ=pm_integral+droop_area-endpoint-aw_integral))
        load_loss = load+lossnet+losses_filter+copper_node.sum(1)+damping_node.sum(1)
        integral_load_loss = float(integral(t,load_loss-load_loss[0],jumps)[-1])
        integrated_frequency_prediction = gov_boundary+gov_antiwindup-(executed_energy[-1]-executed_energy[0])-integral_load_loss
        restart=np.flatnonzero(t>=1.15)[0]
        restart_residual=(executed_energy[-1]-executed_energy[restart])-float(integral(t[restart:],all_rhs[restart:])[-1])
        fit_plus = slope(d_dc_node,dc_rhs_node,dcfit)
        dc_error_plus = rms((d_dc_node-dc_rhs_node)[dcfit])
        dc_error_minus = rms((d_dc_node+dc_rhs_node)[dcfit])
        summary = dict(event=event, n=len(df), end_time_s=float(t[-1]),
            M_total_MW_s=float(m.sum()), kinetic_initial_MJ=float(kin[0]),
            filter_initial_MJ=float(ef[0]), dc_initial_MJ=float(edc[0]),
            kinetic_change_final_MJ=float(kin[-1]-kin[0]),
            filter_change_final_MJ=float(ef[-1]-ef[0]), dc_change_final_MJ=float(edc[-1]-edc[0]),
            dc_change_max_abs_MJ=float(np.max(np.abs(edc-edc[0]))),
            external_net_energy_final_MJ=float(int_all[-1]),
            executed_energy_integral_error_final_MJ=float(residual_exec[-1]),
            executed_energy_integral_error_max_MJ=float(np.max(np.abs(residual_exec))),
            positive_energy_integral_error_max_MJ=float(np.max(np.abs(residual_positive))),
            positive_energy_integral_error_final_MJ=float(residual_positive[-1]),
            sg_energy_integral_error_final_MJ=float(residual_sg[-1]),
            gfl_signed_energy_integral_error_final_MJ=float(residual_gfl[-1]),
            sg_power_derivative_error_rms_MW=rms((d_kin-sg_rhs)[smooth]),
            total_derivative_error_rms_MW=rms((d_executed-all_rhs)[smooth]),
            torque_derivative_error_rms_MW=rms((coi_derivative-torque_rhs)[smooth]),
            naive_coi_power_vs_exact_torque_max_MW=float(np.max(np.abs(sg_rhs-torque_rhs))),
            dc_sign_fit_slope_plus=fit_plus, dc_sign_plus_error_rms_MW=dc_error_plus,
            dc_sign_minus_error_rms_MW=dc_error_minus,
            dc_sign_minus_to_plus_error_ratio=dc_error_minus/dc_error_plus,
            electrical_KCL_error_max_MW=float(np.max(np.abs(p_sg+p_gfl-load-lossnet))),
            pdc_change_max_MW=float(np.max(np.abs(pdc-pdc[0]))),
            coi_final_Hz=float((coi[-1]-1)*FBASE),
            coi_signed_area_postevent_Hz_s=float(integral(t,FBASE*(coi-1))[-1]),
            grid_signed_area_min_Hz_s=float(f_area[-1].min()),
            grid_signed_area_max_Hz_s=float(f_area[-1].max()),
            window_area_boundary_identity_max_error_Hz_s=float(np.max(np.abs(f_area[-1]-boundary_area))),
            tail_coi_mean_Hz=float(np.mean((coi[tail]-1)*FBASE)),
            tail_coi_change_10s_Hz=float((coi[-1]-coi[np.flatnonzero(tail)[0]])*FBASE),
            tail_grid_absmax_Hz=float(np.max(np.abs(gf[tail]))),
            tail_all_rhs_mean_MW=float(np.mean(all_rhs[tail])),
            tail_load_delta_mean_MW=float(np.mean(load[tail])-load[0]),
            tail_net_loss_delta_mean_MW=float(np.mean(lossnet[tail])-lossnet[0]),
            tail_filter_loss_delta_mean_MW=float(np.mean(losses_filter[tail])-losses_filter[0]),
            tail_mechanical_delta_mean_MW=float(np.mean(p_mech[tail])-p_mech[0]),
            max_SG_dispersion_energy_MJ=float(dispersion.max()),
            max_incremental_kinetic_energy_MJ=float(inc_kin.max()),
            kinetic_coi_identity_max_error_MJ=float(np.max(np.abs(kinetic_identity_error))))
        summary.update(governor_endpoint_MJ=float(gov_boundary),
            governor_antiwindup_integral_MJ=float(gov_antiwindup),
            inferred_clipped_samples=inferred_clipped_samples,
            weighted_rotor_frequency_area_MJ=float(weighted_frequency_area),
            load_and_loss_increment_integral_MJ=integral_load_loss,
            governor_energy_area_identity_error_MJ=float(weighted_frequency_area-integrated_frequency_prediction),
            governor_reconstructed_ODE_max_rms_pu_s=gov_ode_max_rms,
            executed_energy_error_restarted_1p15s_to61s_MJ=restart_residual,
            tail_SG_frequency_spread_max_Hz=float(np.max(np.ptp(sf[tail][:,active],axis=1))))
        summaries.append(summary)
        for k,b in enumerate(BUSES):
            node_rows.append(dict(event=event,bus=int(b),rho=float(rho[k]),SG_rating_MVA=float(sn[k]),
                H_s=float(h[k]),M_MW_s=float(m[k]),
                filter_current_initial_squared=float(current2[0,k]),
                DC_initial_MJ=float(edc_node[0,k]),DC_final_change_MJ=float(edc_node[-1,k]-edc_node[0,k]),
                DC_sign_fit_slope_plus=slope(d_dc_node[:,k],dc_rhs_node[:,k],dcfit),
                grid_frequency_final_Hz=float(gf[-1,k]),
                grid_frequency_tail_mean_Hz=float(gf[tail,k].mean()),
                grid_frequency_signed_area_Hz_s=float(f_area[-1,k]),
                rotor_frequency_final_Hz=float(sf[-1,k]) if active[k] else None))
        for step in [1,2,4]:
            ix=np.arange(0,len(t),step)
            integ=integral(t[ix],all_rhs[ix],jumps)
            res=executed_energy[ix]-executed_energy[0]-integ
            quadrature_rows.append(dict(event=event,sample_step_s=dt*step,
                energy_error_final_MJ=float(res[-1]),energy_error_max_MJ=float(np.max(np.abs(res))),
                method="trapezoid; pre-jump endpoint replaced by one-sided linear extrapolation"))
        ix=np.unique(np.r_[np.arange(0,len(t),20),np.flatnonzero((t>=.98)&(t<=1.15)),len(t)-1])
        traces.append(pd.DataFrame(dict(event=event,time_s=t[ix],coi_Hz=(coi[ix]-1)*FBASE,
            kinetic_change_MJ=kin[ix]-kin[0],filter_change_MJ=ef[ix]-ef[0],
            DC_change_MJ=edc[ix]-edc[0],source_net_MW=all_rhs[ix],
            source_net_integral_MJ=int_all[ix],executed_energy_error_MJ=residual_exec[ix],
            positive_energy_error_MJ=residual_positive[ix],load_MW=load[ix],network_losses_MW=lossnet[ix],
            filter_losses_MW=losses_filter[ix],mechanical_MW=p_mech[ix],DC_MW=pdc[ix])))
    pd.DataFrame(summaries).to_csv(OUT/"IEEE39_BALANCE_SUMMARY.csv",index=False)
    pd.DataFrame(node_rows).to_csv(OUT/"IEEE39_BALANCE_NODES.csv",index=False)
    pd.DataFrame(quadrature_rows).to_csv(OUT/"IEEE39_BALANCE_QUADRATURE.csv",index=False)
    pd.DataFrame(gov_rows).to_csv(OUT/"IEEE39_GOVERNOR_BALANCE.csv",index=False)
    pd.concat(traces,ignore_index=True).to_csv(OUT/"IEEE39_BALANCE_TRACES.csv",index=False)
    software=tomllib.loads((ROOT/"reports/experiment_N/SOFTWARE_PROVENANCE.toml").read_text())
    library_root=Path(software["ieee39_source_path"]).parents[2]
    stock=library_root/"src/Library/Renewables/ComposableInverter.jl"
    stock_check=dict(path=str(stock),available=stock.exists())
    if stock.exists():
        stock_check.update(sha256=sha(stock),same_DC_sign="(P_ac - P_dc) / v_dc_state" in stock.read_text(),
            export_power_comment_present="power from converter to AC filter" in stock.read_text())
    parameter_inputs=[dict(path=str(p.relative_to(ROOT)).replace("\\","/"),sha256=sha(p)) for p in
        [candidate,ROOT/"reports/experiment_N/TABLE_N01_original_operating_point.csv",
         ROOT/"reports/experiment_D/inputs/machine.csv",ROOT/"reports/experiment_D/inputs/gov.csv",
         ROOT/"experiments/bnd_expQ2B/certified_search/validate_independent_pd.jl"]]
    audit=dict(scope="archived Q2B, no pure delay, no new simulation; numerical audit",
        source_hashes=source_hashes,input_files=input_files,candidate_sha256=expected,
        installed_stock_source=stock_check,parameter_inputs=parameter_inputs,
        python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,
        constants=dict(SBASE_MVA=SBASE,FBASE_Hz=FBASE,Rf=RF,Xf=XF,Cdc=CDC,window_s=WINDOW),
        physical_storage_gate="NOT_PASSED: executed DC sign gives K+Ef-Edc, not K+Ef+Edc",
        frozen_source_gate=all(x["pass_exact"] for x in source_hashes),
        outputs=["IEEE39_BALANCE_SUMMARY.csv","IEEE39_BALANCE_NODES.csv",
                 "IEEE39_BALANCE_QUADRATURE.csv","IEEE39_BALANCE_TRACES.csv","IEEE39_GOVERNOR_BALANCE.csv"])
    (OUT/"IEEE39_BALANCE_AUDIT.json").write_text(json.dumps(audit,indent=2),encoding="utf-8")
    print(pd.DataFrame(summaries)[["event","executed_energy_integral_error_final_MJ",
        "positive_energy_integral_error_max_MJ","dc_sign_fit_slope_plus",
        "coi_final_Hz","tail_coi_change_10s_Hz"]].to_string(index=False))


if __name__ == "__main__":
    main()
