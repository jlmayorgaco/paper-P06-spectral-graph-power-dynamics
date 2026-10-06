module ExpPTDS

using SHA, TOML
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_P","P5")
const CANDIDATE=joinpath(OUT,"Z_P_NOMINAL_FINAL.toml")
bytes2hex(sha256(read(CANDIDATE)))==strip(read(CANDIDATE*".sha256",String)) ||
    error("candidate SHA mismatch before TDS")
isfile(joinpath(OUT,"TABLE_P5_powerdynamics_validation.csv")) ||
    error("PowerDynamics spectral gate has not passed")
candidate=TOML.parsefile(CANDIDATE)

# Nonlinear PowerDynamics is loaded only after the frozen spectral gate exists.
using CSV, DataFrames, LinearAlgebra, NetworkDynamics, PowerDynamics
using OrdinaryDiffEqRosenbrock, SciMLBase
if !isdefined(@__MODULE__,:PDReferenceN)
    include(joinpath(ROOT,"src","bnd_model_expN","PDReferenceN.jl"))
end
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
    ur=Float64(state[VIndex(bus,:busbar₊u_r)])
    ui=Float64(state[VIndex(bus,:busbar₊u_i)])
    ir=Float64(state[VIndex(bus,:gfl₊filter₊i_f_r)])
    ii=Float64(state[VIndex(bus,:gfl₊filter₊i_f_i)])
    gfp=100r*(ur*ir+ui*ii)
    gfq=100r*(ui*ir-ur*ii)
    sgp,sgq,gfp,gfq
end

function pulse_network(nw,bus,pulse)
    vertices,edges=PD39.PD39Model.copy_network_components(nw)
    defaults=get_defaults_dict(vertices[bus])
    pkeys=[s for s in keys(defaults) if occursin("Pset",string(s))]
    qkeys=[s for s in keys(defaults) if occursin("Qset",string(s))]
    length(pkeys)==1 && length(qkeys)==1 || error("load setpoint absent at bus $bus")
    ps=only(pkeys);qs=only(qkeys)
    p0=defaults[ps];q0=defaults[qs]
    affect=(u,p,ctx)->begin
        factor=ctx.t<1.1 ? 1+pulse : 1.0
        p[ps]=p0*factor;p[qs]=q0*factor
    end
    callback=PresetTimeComponentCallback([1.0,1.1],
        ComponentAffect(affect,(),(ps,qs)))
    set_callback!(vertices[bus],callback)
    out=Network(vertices,edges)
    set_jac_prototype!(out)
    out,abs(p0)*100pulse
end

function one_run(nw,state,rho,bus,pulse)
    active,delta_p=pulse_network(nw,bus,pulse)
    prob=SciMLBase.ODEProblem(active,state,(0.0,60.0))
    sol=SciMLBase.solve(prob,OrdinaryDiffEqRosenbrock.Rodas5P();
        callback=get_callbacks(active),initializealg=SciMLBase.NoInit(),
        saveat=0.1,abstol=1e-9,reltol=1e-9)
    success=SciMLBase.successful_retcode(sol.retcode)
    success || error("TDS failed at bus $bus pulse $pulse: $(sol.retcode)")
    times=collect(0.0:0.1:60.0)
    ω0=Float64(state[VIndex(38,:ctrld_gen₊machine₊ω)])
    v0=[voltage(state,b) for b in 1:39]
    frames=NamedTuple[]
    freq=Float64[];vdev=Float64[]
    for t in times
        s=NetworkDynamics.NWState(sol,t)
        f=60*(Float64(s[VIndex(38,:ctrld_gen₊machine₊ω)])-ω0)
        vals=[voltage(s,b) for b in 1:39]
        dv=maximum(abs.(vals.-v0))
        sgp,sgq,gfp,gfq=injection38(s,rho)
        push!(freq,f);push!(vdev,dv)
        push!(frames,(;time_s=t,event_bus=bus,pulse_fraction=pulse,
            delta_P_MW=delta_p,COI_frequency_deviation_Hz=f,
            bus38_voltage_pu=vals[38],max_voltage_deviation_pu=dv,
            voltage_min_pu=minimum(vals),voltage_max_pu=maximum(vals),
            SG38_P_MW=sgp,SG38_Q_Mvar=sgq,
            GFL38_P_MW=gfp,GFL38_Q_Mvar=gfq,
            GFL38_vdc_pu=Float64(s[VIndex(38,:gfl₊v_dc_state)]),
            GFL30_PLL_frequency_deviation_Hz=
                Float64(s[VIndex(30,:gfl₊pll₊Δω_rad_s)])/(2pi)))
    end
    rocof=vcat(0.0,diff(freq)./0.1)
    frames=[merge(frames[i],(;RoCoF_Hz_s=rocof[i])) for i in eachindex(frames)]
    peak=maximum(abs.(freq))
    threshold=max(0.02*peak,1e-9)
    settling=NaN
    for j in eachindex(times)
        times[j]<1.1 && continue
        if maximum(abs.(freq[j:end]))<=threshold
            settling=times[j]-1.1;break
        end
    end
    metrics=(;event_bus=bus,pulse_fraction=pulse,delta_P_MW=delta_p,
        peak_COI_frequency_Hz=peak,peak_RoCoF_Hz_s=maximum(abs.(rocof)),
        peak_voltage_deviation_pu=maximum(vdev),settling_time_s=settling,
        voltage_min_pu=minimum(r.voltage_min_pu for r in frames),
        voltage_max_pu=maximum(r.voltage_max_pu for r in frames),
        peak_RoCoF_per_MW=maximum(abs.(rocof))/delta_p,
        peak_frequency_per_MW=peak/delta_p,
        solver_status=string(sol.retcode))
    metrics,frames,freq,vdev
end

function main()
    pd=CSV.read(joinpath(OUT,"TABLE_P5_powerdynamics_validation.csv"),DataFrame)
    pd.PD_validation[1]=="PASS" || error("PD spectral validation failed")
    rho=Float64.(candidate["rho"]);kp=Float64.(candidate["Kp"]);ki=Float64.(candidate["Ki"])
    base=PDReferenceN.frozen_baseline()
    nw=PDReferenceN.build_architecture(base,rho,kp,ki)
    state=PDReferenceN.trim_state(nw,base,rho,kp,ki)
    rows=NamedTuple[];traces=NamedTuple[]
    for bus in (8,16,29)
        responses=Any[]
        for pulse in (1e-5,2e-5)
            println("TDS bus=",bus," pulse=",pulse);flush(stdout)
            metrics,frames,freq,vdev=one_run(nw,state,rho,bus,pulse)
            push!(responses,(metrics,frames,freq,vdev))
        end
        freq_error=maximum(abs.(responses[2][3].-2 .* responses[1][3]))/
            max(maximum(abs.(responses[2][3])),1e-12)
        voltage_error=maximum(abs.(responses[2][4].-2 .* responses[1][4]))/
            max(maximum(abs.(responses[2][4])),1e-12)
        for resp in responses
            push!(rows,merge(resp[1],(;frequency_scaling_error=freq_error,
                voltage_scaling_error=voltage_error,
                TDS_validation=(freq_error<0.05 && voltage_error<0.05) ? "PASS" : "FAIL")))
            append!(traces,resp[2])
        end
        println("TDS_SCALING bus=",bus," frequency=",freq_error,
            " voltage=",voltage_error);flush(stdout)
    end
    CSV.write(joinpath(OUT,"TABLE_P5_tds_validation.csv"),DataFrame(rows))
    CSV.write(joinpath(OUT,"TABLE_P5_tds_trajectories.csv"),DataFrame(traces))
    println("TDS_DONE rows=",length(rows)," pass=",count(r->r.TDS_validation=="PASS",rows))
end

main()

end
