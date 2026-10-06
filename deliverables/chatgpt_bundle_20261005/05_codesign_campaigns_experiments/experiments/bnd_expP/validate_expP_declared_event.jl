module ExpPDeclaredEventTDS

using SHA, TOML, CSV, DataFrames, LinearAlgebra
using NetworkDynamics, PowerDynamics, OrdinaryDiffEqRosenbrock, SciMLBase

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_P","P5")
const CANDIDATE=joinpath(OUT,"Z_P_NOMINAL_FINAL.toml")
bytes2hex(sha256(read(CANDIDATE)))==strip(read(CANDIDATE*".sha256",String)) ||
    error("P5 candidate hash mismatch before declared-event validation")
isfile(joinpath(OUT,"TABLE_P5_powerdynamics_validation.csv")) ||
    error("PowerDynamics spectrum gate must pass before declared-event TDS")
const candidate=TOML.parsefile(CANDIDATE)
include(joinpath(ROOT,"src","bnd_model_expN","PDReferenceN.jl"))
const PD39=PDReferenceN.PD39

function voltage(state,bus)
    hypot(Float64(state[VIndex(bus,:busbar₊u_r)]),
          Float64(state[VIndex(bus,:busbar₊u_i)]))
end

function injection38(state,rho)
    bus=38;r=rho[bus-29]
    sn=Float64(state[VIndex(bus,:ctrld_gen₊machine₊Sn)])
    sgp=sn*Float64(state[VIndex(bus,:ctrld_gen₊machine₊P)])
    sgq=sn*Float64(state[VIndex(bus,:ctrld_gen₊machine₊Q)])
    ur=Float64(state[VIndex(bus,:busbar₊u_r)]);ui=Float64(state[VIndex(bus,:busbar₊u_i)])
    ir=Float64(state[VIndex(bus,:gfl₊filter₊i_f_r)]);ii=Float64(state[VIndex(bus,:gfl₊filter₊i_f_i)])
    sgp,sgq,100r*(ur*ir+ui*ii),100r*(ui*ir-ur*ii)
end

function event_network(nw,bus,deltaP_MW;start=1.0,base_mva=100.0)
    vertices,edges=PD39.PD39Model.copy_network_components(nw)
    defaults=get_defaults_dict(vertices[bus])
    pkeys=[s for s in keys(defaults) if occursin("Pset",string(s))]
    qkeys=[s for s in keys(defaults) if occursin("Qset",string(s))]
    length(pkeys)==1 && length(qkeys)==1 || error("load P/Q setpoints absent at bus $bus")
    pkey=only(pkeys);qkey=only(qkeys);p0=defaults[pkey];q0=defaults[qkey]
    delta_p_pu=deltaP_MW/base_mva
    affect=(u,p,ctx)->begin
        ctx.t>=start && (p[pkey]=p0-delta_p_pu)
        p[qkey]=q0
    end
    callback=PresetTimeComponentCallback([start],
        ComponentAffect(affect,(),(pkey,qkey)))
    set_callback!(vertices[bus],callback)
    out=Network(vertices,edges);set_jac_prototype!(out)
    out,(;p0,q0,delta_p_pu,deltaP_MW)
end

function main()
    pd=CSV.read(joinpath(OUT,"TABLE_P5_powerdynamics_validation.csv"),DataFrame)
    only(pd.PD_validation)=="PASS" || error("P5 PD spectrum validation did not pass")
    scenario=TOML.parsefile(joinpath(ROOT,"experiments","bnd_expG","configs","DESIGN_SCENARIO_FROZEN.toml"))
    bus=Int(scenario["event_bus"]);deltaP=Float64(scenario["deltaP_MW"])
    rho=Float64.(candidate["rho"]);kp=Float64.(candidate["Kp"]);ki=Float64.(candidate["Ki"])
    base=PDReferenceN.frozen_baseline()
    nw=PDReferenceN.build_architecture(base,rho,kp,ki)
    state=PDReferenceN.trim_state(nw,base,rho,kp,ki)
    active,profile=event_network(nw,bus,deltaP;base_mva=Float64(scenario["base_MVA"]))
    dt=0.01;tend=60.0
    prob=SciMLBase.ODEProblem(active,state,(0.0,tend))
    sol=SciMLBase.solve(prob,OrdinaryDiffEqRosenbrock.Rodas5P();
        callback=get_callbacks(active),initializealg=SciMLBase.NoInit(),
        saveat=dt,abstol=1e-9,reltol=1e-9)
    SciMLBase.successful_retcode(sol.retcode) || error("declared event TDS failed: $(sol.retcode)")
    times=collect(0.0:dt:tend);omega0=Float64(state[VIndex(38,:ctrld_gen₊machine₊ω)])
    v0=[voltage(state,b) for b in 1:39]
    frames=NamedTuple[];freq=Float64[];rocof=Float64[];maxv=Float64[]
    for t in times
        s=NetworkDynamics.NWState(sol,t)
        f=60*(Float64(s[VIndex(38,:ctrld_gen₊machine₊ω)])-omega0)
        sgp,sgq,gfp,gfq=injection38(s,rho)
        volts=[voltage(s,b) for b in 1:39]
        push!(freq,f);push!(maxv,maximum(abs.(volts.-v0)))
        push!(frames,(;time_s=t,event_bus=bus,delta_P_MW=deltaP,
            COI_frequency_deviation_Hz=f,RoCoF_Hz_s=NaN,
            SG38_P_MW=sgp,SG38_Q_Mvar=sgq,GFL38_P_MW=gfp,GFL38_Q_Mvar=gfq,
            bus38_voltage_pu=volts[38],max_voltage_deviation_pu=maximum(abs.(volts.-v0)),
            voltage_min_pu=minimum(volts),voltage_max_pu=maximum(volts),
            GFL30_PLL_frequency_deviation_Hz=
              Float64(s[VIndex(30,:gfl₊pll₊Δω_rad_s)])/(2pi)))
    end
    rocof=vcat(0.0,diff(freq)./dt)
    frames=[merge(frames[i],(;RoCoF_Hz_s=rocof[i])) for i in eachindex(frames)]
    after=findall(times.>=1.0)
    freqpeak=maximum(abs.(freq[after]));rocofpeak=maximum(abs.(rocof[after]))
    threshold=max(0.02*freqpeak,1e-9);settle=NaN
    for j in after
        if maximum(abs.(freq[j:end]))<=threshold
            settle=times[j]-1.0;break
        end
    end
    rocof_limit=Float64(scenario["rocof_limit_Hz_s"])
    freq_limit=Float64(scenario["frequency_limit_Hz"])
    result=(;event_bus=bus,delta_P_MW=deltaP,pulse_profile="sustained Pset step at 1.0 s",
      reactive_load_change_MWAr=0.0,load_setpoint_delta_pu=profile.delta_p_pu,
      sample_dt_s=dt,horizon_s=tend,solver_status=string(sol.retcode),
      peak_COI_frequency_Hz=freqpeak,peak_RoCoF_Hz_s=rocofpeak,
      peak_voltage_deviation_pu=maximum(maxv),voltage_min_pu=minimum(r.voltage_min_pu for r in frames),
      voltage_max_pu=maximum(r.voltage_max_pu for r in frames),settling_time_s=settle,
      rocof_limit_Hz_s=rocof_limit,frequency_limit_Hz=freq_limit,
      rocof_pass=rocofpeak<=rocof_limit,frequency_pass=freqpeak<=freq_limit,
      overall_pass=rocofpeak<=rocof_limit && freqpeak<=freq_limit,
      tail_status=isfinite(settle) ? "2_PERCENT_SETTLED_WITHIN_HORIZON" : "NOT_SETTLED_WITHIN_60S",
      current_limit_claim="NOT_MADE; no validated hard GFL current limiter")
    CSV.write(joinpath(OUT,"TABLE_P5_declared_event_validation.csv"),DataFrame([result]))
    CSV.write(joinpath(OUT,"TABLE_P5_declared_event_trajectories.csv"),DataFrame(frames))
    println("DECLARED_EVENT_BUS=",bus," DELTAP_MW=",deltaP," ROCOF=",rocofpeak,
        " FREQ=",freqpeak," OVERALL=",result.overall_pass," TAIL=",result.tail_status)
end

main()

end
