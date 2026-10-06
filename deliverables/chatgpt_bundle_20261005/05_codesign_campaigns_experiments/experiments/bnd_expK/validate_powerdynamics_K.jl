using CSV, DataFrames, TOML, SHA, LinearAlgebra, Statistics
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_K")
const TABLES=joinpath(OUT,"tables")

function verified_candidate(path)
    isfile(path) && isfile(path*".sha256") || error("candidate is not frozen")
    hash=bytes2hex(sha256(read(path)))
    hash==strip(read(path*".sha256",String)) || error("candidate hash mismatch")
    c=TOML.parsefile(path)
    c["design_used_PowerDynamics"]==false || error("design purity failed")
    return c,hash
end
nom,nomhash=verified_candidate(joinpath(OUT,"Z_K_NOMINAL_FINAL.toml"))
robust_paths=filter(p->startswith(basename(p),"Z_K_ROBUST_beta_") &&
    endswith(p,".toml"),readdir(OUT;join=true))
for p in robust_paths; verified_candidate(p) end
include(joinpath(ROOT,"src","bnd_design_k","BNDDesignK.jl"))
const FROZEN_CONTEXT=BNDDesignK.design_context(ROOT)

# This import gate is intentionally below the hash verification above.
using NetworkDynamics, PowerDynamics, OrdinaryDiffEqRosenbrock, SciMLBase
include(joinpath(ROOT,"src","pd39","PD39.jl"))
using .PD39

function weighted_load_bus(nw,bus,rho,kp,ki)
    0<rho<=1 || error("load-bus construction requires a GFL share")
    pd=PD39.PD39Model
    data=pd.ieee39_data()
    row=data.bus[findfirst(==(bus),data.bus.bus),:]
    row.has_load || error("special builder requires colocated ZIP load")
    pd.set_Sbase!(pd.BASE_MVA);pd.set_fbase!(pd.BASE_FREQ)
    xf=0.03;vdc=2.5;cdc=1.25;ftau=300.0;fi=600.0
    gfl=pd.WeightedSimpleGFLDC(
        Xf=xf,Rf=0.01,PLL_Kp=kp,PLL_Ki=ki,PLL_τ_lpf=1/(ftau*2pi),
        CC1_KP=(xf/(2pi*pd.BASE_FREQ))*(fi*2pi),
        CC1_KI=(xf/(2pi*pd.BASE_FREQ))*(fi*2pi)^2/4,
        CC1_F=0,CC1_Fcoupl=0,C_dc=cdc,V_dc=vdc,
        kp_v_dc=vdc*cdc*(5*2pi),ki_v_dc=vdc*cdc*(5*2pi)^2/4,
        port_scale=rho,name=:gfl)
    load=deepcopy(pd.OfficialIEEE39.load)
    parts=if rho<1
        sg=deepcopy(bus==39 ? pd.OfficialIEEE39.uncontrolled_machine :
            pd.OfficialIEEE39.controlled_machine)
        pd.MTKBus(sg,gfl,load)
    else
        pd.MTKBus(gfl,load)
    end
    vertices,edges=pd.copy_network_components(nw)
    pf=pd.get_pfmodel(vertices[bus])
    replacement=pd.compile_bus(parts;name=Symbol("bus$(bus)"),vidx=bus,pf)
    pd.set_Sbase!();pd.set_fbase!()
    pd.set_default!(replacement,Symbol("busbar₊Vbase"),Float64(row.base_kv))
    pd.OfficialIEEE39.apply_csv_params!(replacement,data.load,bus)
    # The ZIP Vset is a free initialized parameter in the official example.
    # Recover the actual frozen value rather than accepting its default guess.
    lr=only(eachrow(FROZEN_CONTEXT.net.load_audit[
        FROZEN_CONTEXT.net.load_audit.bus .== bus,:]))
    vset=abs(FROZEN_CONTEXT.net.voltage[bus])*
        sqrt(lr.CSV_setpoint_load_MW/lr.initialized_load_MW)
    pd.set_default!(replacement,r"ZIPLoad₊Vset$",vset)
    if rho<1
        pd.OfficialIEEE39.apply_csv_params!(replacement,data.machine,bus)
        Bool(row.has_avr) && pd.OfficialIEEE39.apply_csv_params!(replacement,data.avr,bus)
        Bool(row.has_gov) && pd.OfficialIEEE39.apply_csv_params!(replacement,data.gov,bus)
        mrow=data.machine[findfirst(==(bus),data.machine.bus),:]
        pd.set_default!(replacement,r"Sn$",(1-rho)*Float64(mrow.Sn))
    end
    vertices[bus]=replacement
    out=pd.Network(vertices,edges)
    pd.set_jac_prototype!(out)
    out
end

function network_for(c)
    nw=PD39.baseline_network()
    for i in 1:10
        bus=29+i
        rho=Float64(c["rho"][i])
        if rho>0 && bus in (31,39)
            nw=weighted_load_bus(nw,bus,rho,Float64(c["Kp"][i]),Float64(c["Ki"][i]))
        elseif rho==1.0
            nw=PD39.replace_bus(nw,bus;template=PD39.simple_gfldc_template())
        elseif 0<rho<1
            nw=PD39.PD39Model.weighted_replacement_network(nw,bus,rho,1.0)
        end
    end
    vertices,edges=PD39.PD39Model.copy_network_components(nw)
    for i in 1:10
        rho=Float64(c["rho"][i]);rho>0 || continue
        bus=29+i
        set_default!(vertices[bus],r"PLL_Kp$",Float64(c["Kp"][i]))
        set_default!(vertices[bus],r"PLL_Ki$",Float64(c["Ki"][i]))
    end
    out=Network(vertices,edges)
    set_jac_prototype!(out)
    out
end

function frequency_state(state,bus)
    for sym in (:ctrld_gen₊machine₊ω,:machine₊ω,:gfl₊pll₊ω)
        try return Float64(state[VIndex(bus,sym)]) catch end
    end
    NaN
end
function voltage_pu(state,bus)
    try
        hypot(Float64(state[VIndex(bus,:busbar₊u_r)]),
            Float64(state[VIndex(bus,:busbar₊u_i)]))
    catch
        NaN
    end
end
function pulse_network(nw,bus,pulse)
    vertices,edges=PD39.PD39Model.copy_network_components(nw)
    defaults=get_defaults_dict(vertices[bus])
    ps=only([s for s in keys(defaults) if occursin("Pset",string(s))])
    qs=only([s for s in keys(defaults) if occursin("Qset",string(s))])
    p0=defaults[ps];q0=defaults[qs]
    affect=(u,p,ctx)->begin
        factor=ctx.t<1.1 ? 1+pulse : 1.0
        p[ps]=p0*factor;p[qs]=q0*factor
    end
    set_callback!(vertices[bus],PresetTimeComponentCallback([1.0,1.1],
        ComponentAffect(affect,(),(ps,qs))))
    out=Network(vertices,edges)
    set_jac_prototype!(out)
    out
end
function run_tds(nw,eq,c,event_bus,pulse)
    pn=pulse_network(nw,event_bus,pulse)
    prob=SciMLBase.ODEProblem(pn,eq.state,(0.0,12.0))
    sol=SciMLBase.solve(prob,OrdinaryDiffEqRosenbrock.Rodas5P();
        callback=get_callbacks(pn),initializealg=SciMLBase.NoInit(),
        saveat=0.02,abstol=1e-8,reltol=1e-8)
    times=Float64.(sol.t)
    retained=[i for i in 1:10 if Float64(c["epsilon"][i])>0]
    isempty(retained) && error("no retained SG COI")
    # Inertia values come from the frozen SG model, not an arbitrary weight.
    weights=[Float64(c["epsilon"][i])*FROZEN_CONTEXT.net.sg[29+i].op.parameters.inertia*
        FROZEN_CONTEXT.net.sg[29+i].op.parameters.rating_mva for i in retained]
    refs=[frequency_state(eq.state,29+i) for i in retained]
    all(isfinite,refs) || error("SG frequency state unavailable")
    coi=Float64[];vdev=Float64[]
    vref=[voltage_pu(eq.state,b) for b in 1:39]
    for t in times
        state=NetworkDynamics.NWState(sol,t)
        omega=[frequency_state(state,29+i) for i in retained]
        push!(coi,60dot(weights,omega.-refs)/sum(weights))
        push!(vdev,maximum(abs.([voltage_pu(state,b) for b in 1:39].-vref)))
    end
    dt=diff(times);dcoi=diff(coi)
    rates=[dcoi[i]/dt[i] for i in eachindex(dt) if dt[i]>0]
    (;times,coi,vdev,peak_frequency=maximum(abs.(coi)),
        peak_rocof=maximum(abs.(rates)),peak_voltage=maximum(vdev),
        finite=all(isfinite,coi)&&all(isfinite,rates)&&all(isfinite,vdev),
        success=SciMLBase.successful_retcode(sol.retcode))
end

front=CSV.read(joinpath(TABLES,"TABLE_K09_robust_frontier.csv"),DataFrame)
betag=1.6991206999182038e-6
rows=front[abs.(Float64.(front.beta_req).-betag).<1e-14,:]
selected=nrow(rows)>0 && rows.status[1]=="BEST_SAMPLED_CANDIDATE" ?
    first(filter(p->begin
        c,_=verified_candidate(p)
        abs(Float64(c["beta_req"])-betag)<1e-14
    end,robust_paths)) : nothing
candidates=[("nominal",nom,nomhash)]
if selected!==nothing
    c,h=verified_candidate(selected)
    push!(candidates,("robust_ExpG_beta",c,h))
end
worstbus=16
tpath=joinpath(TABLES,"TABLE_K14_transient_capacity.csv")
if isfile(tpath)
    tr=CSV.read(tpath,DataFrame)
    nrow(tr)>0 && (worstbus=Int(tr.bus[argmin(tr.deltaP_max_total_MW)]))
end
pdrows=NamedTuple[];tdsrows=NamedTuple[];traces=DataFrame[]
for (label,c,hash) in candidates
    println("K PowerDynamics validating ",label," hash=",hash);flush(stdout)
    try
        nw=network_for(c)
        eq=PD39.initialize_equilibrium(nw;sparse=false,check=:error)
        audit=PD39.stability_audit(eq.state)
        alpha=audit.max_real
        delta=alpha-Float64(c["spectral_abscissa_s_inv"])
        pass=eq.fixed_point && audit.finite && alpha<=-0.05 &&
            abs(delta)<=1e-3 && PD39._eq_residual(eq.state)<1e-8
        push!(pdrows,(candidate=label,status=pass ? "PASS" : "FAIL",
            alpha_analytic=Float64(c["spectral_abscissa_s_inv"]),
            alpha_PD=alpha,delta_alpha=delta,
            equilibrium_residual_inf=PD39._eq_residual(eq.state),
            finite_poles=length(audit.eigenvalues),
            numerical_gauges=length(audit.gauge_eigenvalues),
            fixed_point=eq.fixed_point,candidate_sha256=hash,error=""))
        # Run the requested nonlinear falsification pulses whenever a finite
        # initialized equilibrium exists, even if the spectral comparison
        # itself failed. Such TDS results cannot rescue a failed PD gate.
        if eq.fixed_point && audit.finite
            for bus in unique([16,worstbus]),pulse in (0.0005,0.001)
                try
                    t=run_tds(nw,eq,c,bus,pulse)
                    push!(tdsrows,(candidate=label,event_bus=bus,pulse_fraction=pulse,
                        status=t.success&&t.finite ? "EVALUATED" : "FAIL",
                        peak_frequency_Hz=t.peak_frequency,peak_RoCoF_Hz_s=t.peak_rocof,
                        peak_voltage_deviation_pu=t.peak_voltage,
                        finite=t.finite,candidate_sha256=hash,error=""))
                    push!(traces,DataFrame(candidate=fill(label,length(t.times)),
                        event_bus=fill(bus,length(t.times)),pulse_fraction=fill(pulse,length(t.times)),
                        time_s=t.times,COI_frequency_Hz=t.coi,
                        max_voltage_deviation_pu=t.vdev))
                catch err
                    push!(tdsrows,(candidate=label,event_bus=bus,pulse_fraction=pulse,
                        status="FAIL_EXCEPTION",peak_frequency_Hz=NaN,
                        peak_RoCoF_Hz_s=NaN,peak_voltage_deviation_pu=NaN,
                        finite=false,candidate_sha256=hash,error=sprint(showerror,err)))
                end
            end
        end
    catch err
        push!(pdrows,(candidate=label,status="FAIL_EXCEPTION",
            alpha_analytic=Float64(c["spectral_abscissa_s_inv"]),
            alpha_PD=NaN,delta_alpha=NaN,equilibrium_residual_inf=NaN,
            finite_poles=0,numerical_gauges=0,fixed_point=false,
            candidate_sha256=hash,error=sprint(showerror,err)))
    end
end
CSV.write(joinpath(TABLES,"TABLE_K12_powerdynamics_validation.csv"),DataFrame(pdrows))
if !isempty(tdsrows)
    tds=DataFrame(tdsrows)
    tds[!,:relative_frequency_scaling_error]=fill(NaN,nrow(tds))
    tds[!,:relative_voltage_scaling_error]=fill(NaN,nrow(tds))
    for label in unique(tds.candidate),bus in unique(tds.event_bus)
        ids=findall((tds.candidate.==label).&(tds.event_bus.==bus))
        length(ids)==2 || continue
        lo,hi=sort(ids;by=i->tds.pulse_fraction[i])
        if all(isfinite,tds.peak_frequency_Hz[[lo,hi]])
            ferr=abs(tds.peak_frequency_Hz[hi]/tds.peak_frequency_Hz[lo]-2)/2
            verr=abs(tds.peak_voltage_deviation_pu[hi]/tds.peak_voltage_deviation_pu[lo]-2)/2
            tds.relative_frequency_scaling_error[ids].=ferr
            tds.relative_voltage_scaling_error[ids].=verr
        end
    end
    CSV.write(joinpath(TABLES,"TABLE_K13_tds_validation.csv"),tds)
else
    CSV.write(joinpath(TABLES,"TABLE_K13_tds_validation.csv"),
        DataFrame(candidate=String[],event_bus=Int[],pulse_fraction=Float64[],status=String[]))
end
!isempty(traces) && CSV.write(joinpath(TABLES,"TABLE_K16_tds_trajectories.csv"),vcat(traces...))
println("K_PD_NOMINAL: ",first(pdrows).status)
