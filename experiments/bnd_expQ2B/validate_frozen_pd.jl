using CSV, DataFrames, LinearAlgebra, SHA, TOML, Pkg
using NetworkDynamics, PowerDynamics, OrdinaryDiffEqRosenbrock, SciMLBase

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_Q2B")
const CORE=joinpath(OUT,"CORE")
const CANDIDATE=joinpath(CORE,"Z_Q2B_SECURE_T05_FINAL.toml")
const SHA_EXPECTED="f53250383cdb906a2d0cb3130ca9c1e72c63b472717c5c82aea1d5838ac785ce"
bytes2hex(sha256(read(CANDIDATE)))==SHA_EXPECTED || error("frozen candidate SHA mismatch")
strip(read(CANDIDATE*".sha256",String))==SHA_EXPECTED || error("SHA sidecar mismatch")

include(joinpath(ROOT,"src","bnd_design_p","ExpP.jl"))
include(joinpath(ROOT,"src","bnd_expQ","LinearSecurity.jl"))
include(joinpath(ROOT,"src","bnd_expQ2B","FiniteWindow.jl"))
include(joinpath(ROOT,"src","bnd_model_expN","PDReferenceN.jl"))
const N=ExpP.PDExactDesignN
const PDRef=PDReferenceN

function event_network(nw,bus,deltaP_MW;start=1.0,base_mva=100.0)
    mdl=PDRef.PD39.PD39Model
    vertices,edges=mdl.copy_network_components(nw)
    defaults=get_defaults_dict(vertices[bus])
    pkeys=[s for s in keys(defaults) if occursin("Pset",string(s))]
    qkeys=[s for s in keys(defaults) if occursin("Qset",string(s))]
    length(pkeys)==1 && length(qkeys)==1 || error("load P/Q setpoints missing at bus $bus")
    pkey=only(pkeys);qkey=only(qkeys);p0=defaults[pkey];q0=defaults[qkey]
    dp=Float64(deltaP_MW)/base_mva
    affect=(u,p,ctx)->begin
        ctx.t>=start && (p[pkey]=p0-dp)
        p[qkey]=q0
    end
    cb=PresetTimeComponentCallback([start],ComponentAffect(affect,(),(pkey,qkey)))
    set_callback!(vertices[bus],cb)
    active=Network(vertices,edges);set_jac_prototype!(active)
    active,(;p0,q0,delta_p_pu=dp,deltaP_MW)
end

function wrapdiff(x)
    mod(x+pi,2pi)-pi
end

function add_vindex!(df,name,vals)
    df[!,Symbol(name)]=vals
end

function main()
    c=TOML.parsefile(CANDIDATE)
    rho=Float64.(c["rho"]);epsilon=Float64.(c["epsilon"])
    kp=Float64.(c["Kp"]);ki=Float64.(c["Ki"])
    ctx=N.design_context(ROOT)
    # Analytical design matrix is read-only after freeze; compare it to the
    # separately rebuilt PD network below.
    an=N.spectrum(ctx,rho,kp,ki)
    anwindows=FiniteWindow.design_metrics(ctx,an.model,rho;
        load_bus=16,disturbance_MW=100.0,windows=(0.2,0.5,1.0,2.0),
        dt_s=0.01,horizon_s=60.0,gauge_vector=N.gauge_vector)
    CSV.write(joinpath(CORE,"TABLE_Q2B_frozen_analytic_windows.csv"),anwindows.metrics)

    base=PDRef.frozen_baseline()
    nw=PDRef.build_architecture(base,rho,kp,ki)
    state=PDRef.trim_state(nw,base,rho,kp,ki)
    trim=PDRef.residual_audit(nw,state)
    pq=PDRef.direct_power_audit(state,base,rho)
    CSV.write(joinpath(CORE,"TABLE_Q2B_PD_component_trim.csv"),pq)
    maxP=maximum(Float64.(pq.max_P_error_pu));maxQ=maximum(Float64.(pq.max_Q_error_pu))
    load_err=0.0
    for bus in (31,39)
        row=only(eachrow(pq[pq.bus.==bus,:]));op=base.rows[bus]
        # The actual component loads stay as initialized in the frozen base.
        load_err=max(load_err,abs(row.ZIP_P_MW-100*Float64(state[VIndex(bus,:ZIPLoad₊P)]))/100,
            abs(row.ZIP_Q_Mvar-100*Float64(state[VIndex(bus,:ZIPLoad₊Q)]))/100)
        abs(row.ZIP_P_MW-100*Float64(base.state[VIndex(bus,:ZIPLoad₊P)]))<1e-9 ||
            error("bus $bus initialized ZIP P changed")
        abs(row.ZIP_Q_Mvar-100*Float64(base.state[VIndex(bus,:ZIPLoad₊Q)]))<1e-9 ||
            error("bus $bus initialized ZIP Q changed")
        abs(row.target_SG_P_MW+row.target_GFL_P_MW-op.P)<1e-9 || error("P share changed at $bus")
        abs(row.target_SG_Q_Mvar+row.target_GFL_Q_Mvar-op.Q)<1e-9 || error("Q share changed at $bus")
    end

    pdlin=linearize_network(state)
    lambda_pd=jacobian_eigenvals(pdlin)
    gauge_idx=argmin(abs.(lambda_pd));gauge=lambda_pd[gauge_idx]
    phys=lambda_pd[[j for j in eachindex(lambda_pd) if j!=gauge_idx]]
    length(phys)==length(an.lambda) || error("finite physical spectrum size mismatch")
    rr,cc=ExpP.linear_assignment(abs.(phys.-transpose(an.lambda)))
    pole_errors=abs.(phys[rr].-an.lambda[cc])
    alpha_pd=maximum(real.(phys));pole_err=maximum(pole_errors)

    # Independent PowerDynamics simulation of the frozen, sustained event.
    active,profile=event_network(nw,16,100.0;start=1.0,base_mva=100.0)
    dt=0.01;tend=61.0
    prob=SciMLBase.ODEProblem(active,state,(0.0,tend))
    sol=SciMLBase.solve(prob,OrdinaryDiffEqRosenbrock.Rodas5P();
        callback=get_callbacks(active),initializealg=SciMLBase.NoInit(),saveat=dt,
        abstol=1e-9,reltol=1e-9)
    simok=SciMLBase.successful_retcode(sol.retcode)
    simok || error("PowerDynamics sustained-step TDS failed: $(sol.retcode)")
    times=collect(0.0:dt:tend);nt=length(times)
    angles=zeros(10,nt);voltages=zeros(10,nt)
    rotor=zeros(10,nt);pll=zeros(10,nt)
    sgP=zeros(10,nt);sgQ=zeros(10,nt);gflP=zeros(10,nt);gflQ=zeros(10,nt)
    pdc=zeros(10,nt);vdc=zeros(10,nt);mech=zeros(10,nt);netP=zeros(10,nt)
    loadP=zeros(10,nt);loadQ=zeros(10,nt);kcl=zeros(10,nt)
    previous=zeros(10);theta=zeros(10);initial_angle=zeros(10);v0=zeros(10)
    bus_voltage_freq=zeros(10,nt)
    for (j,t) in enumerate(times)
        ss=NetworkDynamics.NWState(sol,t)
        for bus in 30:39
            k=bus-29
            ur=Float64(ss[VIndex(bus,:busbar₊u_r)]);ui=Float64(ss[VIndex(bus,:busbar₊u_i)])
            raw=atan(ui,ur)
            if j==1
                previous[k]=raw;theta[k]=raw;initial_angle[k]=raw;v0[k]=hypot(ur,ui)
            else
                theta[k]+=wrapdiff(raw-previous[k]);previous[k]=raw
            end
            angles[k,j]=theta[k]-initial_angle[k]
            voltages[k,j]=hypot(ur,ui)
            uprefix=bus==39 ? "machine₊" : "ctrld_gen₊machine₊"
            sn=Float64(ss[VIndex(bus,Symbol(uprefix*"Sn"))])
            psg=sn*Float64(ss[VIndex(bus,Symbol(uprefix*"P"))])
            qsg=sn*Float64(ss[VIndex(bus,Symbol(uprefix*"Q"))])
            omega=Float64(ss[VIndex(bus,Symbol(uprefix*"ω"))])
            pg=Float64(rho[k])*100*(ur*Float64(ss[VIndex(bus,:gfl₊filter₊i_f_r)])+
                ui*Float64(ss[VIndex(bus,:gfl₊filter₊i_f_i)]))
            qg=Float64(rho[k])*100*(ui*Float64(ss[VIndex(bus,:gfl₊filter₊i_f_r)])-
                ur*Float64(ss[VIndex(bus,:gfl₊filter₊i_f_i)]))
            sgP[k,j]=psg;sgQ[k,j]=qsg;gflP[k,j]=pg;gflQ[k,j]=qg
            rotor[k,j]=60*(omega-1.0)
            pll[k,j]=Float64(ss[VIndex(bus,:gfl₊pll₊Δω_rad_s)])/(2pi)
            pdc[k,j]=Float64(ss[VIndex(bus,:gfl₊P_dc)])
            vdc[k,j]=Float64(ss[VIndex(bus,:gfl₊v_dc_state)])
            mech[k,j]=Float64(ss[VIndex(bus,Symbol(uprefix*"τ_m"))])
            netP[k,j]=Float64(ss[VIndex(bus,:busbar₊P_MW)])
            lp=bus in (31,39) ? 100*Float64(ss[VIndex(bus,:ZIPLoad₊P)]) : 0.0
            lq=bus in (31,39) ? 100*Float64(ss[VIndex(bus,:ZIPLoad₊Q)]) : 0.0
            loadP[k,j]=lp;loadQ[k,j]=lq
            kcl[k,j]=psg+pg+lp-netP[k,j]
        end
    end
    for k in 1:10
        for j in 2:nt-1
            bus_voltage_freq[k,j]=(angles[k,j+1]-angles[k,j-1])/(4pi*dt)
        end
        bus_voltage_freq[k,1]=(angles[k,2]-angles[k,1])/(2pi*dt)
        bus_voltage_freq[k,end]=(angles[k,end]-angles[k,end-1])/(2pi*dt)
    end
    # Express time relative to the declared step and pad the pre-event window.
    trel=times.-1.0
    pre=collect(-4.0:dt:-1.01)
    tmetric=vcat(pre,trel)
    phmetric=hcat(zeros(10,length(pre)),angles)
    steady=zeros(10)
    n10=Int(round(10/dt))
    for k in 1:10
        steady[k]=(phmetric[k,end]-phmetric[k,end-n10])/(2pi*10.0)
    end
    signal=(;times_s=tmetric,phase_rad=phmetric,phase_jump_rad=zeros(10),
        F_inf_Hz=steady,poles=ComplexF64[],alpha=alpha_pd,A=zeros(1,1),B=zeros(1,1),
        C_bus=zeros(10,1),D_bus=zeros(10),condition_A=NaN,eigenvector_condition=NaN)
    win=FiniteWindow.window_metrics(signal;windows=(0.2,0.5,1.0,2.0),dt_s=dt,
        horizon_s=60.0,disturbance_MW=100.0)
    CSV.write(joinpath(CORE,"TABLE_Q2B_PD_TDS_metrics.csv"),DataFrame(win))

    df=DataFrame(time_s=times,event_time_s=trel)
    for k in 1:10
        bus=29+k
        add_vindex!(df,"phase_bus$(bus)_rad",angles[k,:])
        add_vindex!(df,"fSG_bus$(bus)_Hz",rotor[k,:])
        add_vindex!(df,"fPLL_bus$(bus)_Hz",pll[k,:])
        add_vindex!(df,"fBus_unfiltered_bus$(bus)_Hz",bus_voltage_freq[k,:])
        add_vindex!(df,"SG_P_bus$(bus)_MW",sgP[k,:]);add_vindex!(df,"SG_Q_bus$(bus)_Mvar",sgQ[k,:])
        add_vindex!(df,"GFL_P_bus$(bus)_MW",gflP[k,:]);add_vindex!(df,"GFL_Q_bus$(bus)_Mvar",gflQ[k,:])
        add_vindex!(df,"GFL_Pdc_bus$(bus)_pu",pdc[k,:]);add_vindex!(df,"GFL_Vdc_bus$(bus)_pu",vdc[k,:])
        add_vindex!(df,"SG_Pm_bus$(bus)_pu",mech[k,:]);add_vindex!(df,"KCL_P_error_bus$(bus)_MW",kcl[k,:])
        add_vindex!(df,"ZIP_P_bus$(bus)_MW",loadP[k,:]);add_vindex!(df,"ZIP_Q_bus$(bus)_Mvar",loadQ[k,:])
        add_vindex!(df,"V_bus$(bus)_pu",voltages[k,:])
    end
    df[!, :max_abs_fSG_Hz]=[maximum(abs.(rotor[:,j])) for j in 1:nt]
    df[!, :max_abs_fPLL_Hz]=[maximum(abs.(pll[:,j])) for j in 1:nt]
    df[!, :max_abs_fBus_unfiltered_Hz]=[maximum(abs.(bus_voltage_freq[:,j])) for j in 1:nt]
    CSV.write(joinpath(CORE,"TABLE_Q2B_PD_TDS_timeseries.csv"),df)

    # Put the PD trajectory onto the frozen causal T=0.5 s measurement used by
    # design. Compare its peak and settling diagnostics to the analytical one.
    pd05=only(filter(x->x.window_s==0.5,win))
    an05=only(filter(x->x.window_s==0.5,anwindows.metrics))
    fpeak=pd05.F_peak_Hz;ropeak=pd05.R_peak_Hz_s;finf=pd05.F_inf_Hz
    threshold=max(0.02*max(abs(finf),1e-9),1e-6)
    # Settled-to-asymptote means all bus-window frequencies remain within 2% of
    # their own last-10-s average for the remaining simulated horizon.
    f05=zeros(10,length(tmetric));lag=Int(round(.5/dt));z0=findfirst(==(0.0),tmetric)
    for j in z0:length(tmetric),k in 1:10
        jm=j-lag;θm=jm>=z0 ? phmetric[k,jm] : 0.0
        f05[k,j]=(phmetric[k,j]-θm)/(2pi*.5)
    end
    settle=NaN
    err_by_time=vec(maximum(abs.(f05.-reshape(steady,:,1)),dims=1))
    suffix_error=reverse(accumulate(max,reverse(err_by_time)))
    for j in z0:length(tmetric)
        if suffix_error[j]<=threshold
            settle=tmetric[j];break
        end
    end

    versions=Dict{String,Any}()
    deps=Pkg.dependencies()
    for name in ("PowerDynamics","NetworkDynamics")
        versions[name]=string(only(filter(x->x.name==name,collect(values(deps)))).version)
    end
    result=Dict{String,Any}(
        "candidate_sha256"=>SHA_EXPECTED,"post_freeze_pd_builds"=>1,
        "trim_residual_max"=>trim.maximum,"max_P_error_pu"=>maxP,
        "max_Q_error_pu"=>maxQ,"separate_load_error_pu"=>load_err,
        "component_bounds_pass"=>all(pq.bounds_status.=="PASS"),
        "gauge_eigenvalue"=>[real(gauge),imag(gauge)],
        "analytic_alpha"=>an.alpha,"pd_alpha"=>alpha_pd,
        "alpha_error"=>abs(an.alpha-alpha_pd),"maximum_relevant_pole_error"=>pole_err,
        "physical_pole_count"=>length(phys),"TDS_retcode"=>string(sol.retcode),
        "TDS_event_bus"=>16,"TDS_deltaP_MW"=>100.0,"TDS_start_s"=>1.0,
        "TDS_Fpeak_T05_Hz"=>fpeak,"TDS_Fsteady_T05_Hz"=>finf,
        "TDS_Rpeak_T05_Hz_s"=>ropeak,
        "analytic_Fpeak_T05_Hz"=>an05.F_peak_Hz,
        "analytic_Finf_Hz"=>an05.F_inf_Hz,
        "analytic_Rpeak_T05_Hz_s"=>an05.R_peak_Hz_s,
        "settling_to_steady_after_event_s"=>settle,
        "max_KCL_P_error_MW"=>maximum(abs.(kcl)),
        "voltage_min_pu"=>minimum(voltages),"voltage_max_pu"=>maximum(voltages),
        "frequency_limit_Hz"=>0.5,"rocof_limit_Hz_s"=>0.5,
        "frequency_pass"=>fpeak<=0.5,"rocof_pass"=>ropeak<=0.5,
        "PD_spectral_validation"=>(abs(an.alpha-alpha_pd)<=1e-4 && pole_err<=1e-4 &&
            sign(alpha_pd)==sign(an.alpha)),
        "PD_TDS_event_pass"=>(fpeak<=0.5 && ropeak<=0.5),
        "runtime_notes"=>"PD independent validation after candidate freeze; no gain retuning.",
        "versions"=>versions)
    open(joinpath(CORE,"PD_VALIDATION_RESULTS.toml"),"w") do io;TOML.print(io,result);end
    println("PD_TRIM max_residual=",trim.maximum," P_err=",maxP," Q_err=",maxQ)
    println("PD_SPECTRUM alpha_analytic=",an.alpha," alpha_PD=",alpha_pd,
        " alpha_err=",abs(an.alpha-alpha_pd)," pole_err=",pole_err,
        " poles=",length(phys)," gauge=",gauge)
    println("PD_TDS T05 Fpeak=",fpeak," Finf=",finf," Rpeak=",ropeak,
        " settle_s=",settle," freq_pass=",fpeak<=.5," rocof_pass=",ropeak<=.5)
end

main()
