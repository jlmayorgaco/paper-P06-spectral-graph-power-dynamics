using CSV, DataFrames, LinearAlgebra, SHA, TOML
using NetworkDynamics, PowerDynamics, OrdinaryDiffEqRosenbrock, SciMLBase

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUT = joinpath(ROOT, "reports", "experiment_Q", "Q1A")
mkpath(OUT)
include(joinpath(ROOT, "src", "bnd_design_p", "ExpP.jl"))
include(joinpath(ROOT, "src", "bnd_model_expN", "PDReferenceN.jl"))
include(joinpath(ROOT, "src", "bnd_expQ", "FrequencyMetrics.jl"))
using .ExpP, .PDReferenceN, .FrequencyMetrics

function add_shared_event(nw, bus, p0, q0; start=1.0, base_mva=100.0)
    vertices, edges = PDReferenceN.PD39.PD39Model.copy_network_components(nw)
    defaults = get_defaults_dict(vertices[bus])
    pkey = only([s for s in keys(defaults) if occursin("Pset", string(s))])
    qkey = only([s for s in keys(defaults) if occursin("Qset", string(s))])
    amp = Ref(0.0)
    dur = Ref(Inf)
    affect = (u, p, ctx) -> begin
        active = ctx.t >= start && ctx.t < start + dur[]
        p[pkey] = active ? p0 - amp[] / base_mva : p0
        p[qkey] = q0
    end
    cb = PresetTimeComponentCallback([start, start+0.1],
        ComponentAffect(affect, (), (pkey,qkey)))
    set_callback!(vertices[bus], cb)
    active_net = Network(vertices,edges)
    set_jac_prototype!(active_net)
    active_net,amp,dur
end

function read_frequencies(sol, times, rho, dt; half_window=5, degree=3)
    ntime, nb = length(times), 10
    sg = fill(NaN,ntime,nb)
    pll = fill(NaN,ntime,nb)
    volt = Matrix{ComplexF64}(undef,ntime,nb)
    for (k,t) in enumerate(times)
        s = NetworkDynamics.NWState(sol,t)
        for bus in 30:39
            i=bus-29
            ur=Float64(s[VIndex(bus,:busbar₊u_r)])
            ui=Float64(s[VIndex(bus,:busbar₊u_i)])
            volt[k,i]=complex(ur,ui)
            if rho[i] < 1
                key=bus==39 ? :machine₊ω : :ctrld_gen₊machine₊ω
                sg[k,i]=sg_frequency_deviation(Float64(s[VIndex(bus,key)]))
            end
            if rho[i] > 0
                pll[k,i]=pll_frequency_deviation(
                    Float64(s[VIndex(bus,:gfl₊pll₊Δω_rad_s)]))
            end
        end
    end
    busfreq=fill(NaN,ntime,nb)
    busrocof=fill(NaN,ntime,nb)
    sg_rocof=fill(NaN,ntime,nb)
    pll_rocof=fill(NaN,ntime,nb)
    for i in 1:nb
        busfreq[:,i].=bus_frequency_deviation(volt[:,i],dt;half_window,degree)
        busrocof[:,i].=sg_polynomial_derivative(busfreq[:,i],dt;half_window,degree)
        any(isfinite,sg[:,i]) && (sg_rocof[:,i].=sg_polynomial_derivative(sg[:,i],dt;half_window,degree))
        any(isfinite,pll[:,i]) && (pll_rocof[:,i].=sg_polynomial_derivative(pll[:,i],dt;half_window,degree))
    end
    (;sg,pll,busfreq,busrocof,sg_rocof,pll_rocof,volt)
end

function peak_and_location(a::AbstractMatrix, buses; mask=trues(size(a,1)))
    best = (-Inf, NaN, NaN)
    tids=findall(mask)
    for i in axes(a,2)
        finite = filter(k -> isfinite(a[tids[k],i]), eachindex(tids))
        isempty(finite) && continue
        valid_tids=tids[finite]
        vals = abs.(a[valid_tids,i])
        j = argmax(vals)
        if vals[j] > best[1]
            best = (vals[j], buses[i], valid_tids[j])
        end
    end
    best
end

function main()
    candidate_path=joinpath(ROOT,"reports","experiment_P","P5","Z_P_NOMINAL_FINAL.toml")
    candidate_sha=bytes2hex(sha256(read(candidate_path)))
    candidate_sha==strip(read(candidate_path*".sha256",String)) ||
        error("frozen ExpP candidate hash mismatch")
    c=TOML.parsefile(candidate_path)
    rho,kp,ki=Float64.(c["rho"]),Float64.(c["Kp"]),Float64.(c["Ki"])
    base=PDReferenceN.frozen_baseline()
    nw=PDReferenceN.build_architecture(base,rho,kp,ki)
    state=PDReferenceN.trim_state(nw,base,rho,kp,ki)
    scenario=TOML.parsefile(joinpath(ROOT,"experiments","bnd_expG","configs","DESIGN_SCENARIO_FROZEN.toml"))
    event_bus=Int(scenario["event_bus"]);base_mva=Float64(scenario["base_MVA"])
    vertices,_=PDReferenceN.PD39.PD39Model.copy_network_components(nw)
    defaults=get_defaults_dict(vertices[event_bus])
    pkey=only([s for s in keys(defaults) if occursin("Pset",string(s))])
    qkey=only([s for s in keys(defaults) if occursin("Qset",string(s))])
    active,amp,dur=add_shared_event(nw,event_bus,defaults[pkey],defaults[qkey];base_mva)
    dt=0.01
    specs=[("step_1MW",1.0,Inf,15.0),
           ("step_10MW",10.0,Inf,15.0),
           ("step_100MW",100.0,Inf,60.0),
           ("pulse_100MW_0p1s",100.0,0.1,15.0)]
    rows=NamedTuple[]
    trajectories=DataFrame[]
    for (label,mw,duration,horizon) in specs
        amp[]=mw;dur[]=duration
        sol=SciMLBase.solve(SciMLBase.ODEProblem(active,deepcopy(state),(0.0,horizon)),
            OrdinaryDiffEqRosenbrock.Rodas5P();callback=get_callbacks(active),
            initializealg=SciMLBase.NoInit(),saveat=dt,abstol=1e-9,reltol=1e-9)
        SciMLBase.successful_retcode(sol.retcode) || error("$label TDS failed: $(sol.retcode)")
        times=collect(0.0:dt:horizon)
        freq=read_frequencies(sol,times,rho,dt)
        after=times .>= 1.0
        primary=hcat(freq.sg,freq.pll)
        primary_rocof=hcat(freq.sg_rocof,freq.pll_rocof)
        fp,fb,ft=peak_and_location(primary,repeat(collect(30:39),2);mask=after)
        rp,rb,rt=peak_and_location(primary_rocof,repeat(collect(30:39),2);mask=after)
        bp,bb,bt=peak_and_location(freq.busfreq,collect(30:39);mask=after)
        br,_,_=peak_and_location(freq.busrocof,collect(30:39);mask=after)
        sgpeak,sgbus,_=peak_and_location(freq.sg,collect(30:39);mask=after)
        pllpeak,pllbus,_=peak_and_location(freq.pll,collect(30:39);mask=after)
        sgr,_,_=peak_and_location(freq.sg_rocof,collect(30:39);mask=after)
        pllr,_,_=peak_and_location(freq.pll_rocof,collect(30:39);mask=after)
        f_limit=Float64(scenario["frequency_limit_Hz"])
        r_limit=Float64(scenario["rocof_limit_Hz_s"])
        primary_f_pass=fp<=f_limit
        primary_r_pass=rp<=r_limit
        bus_f_pass=bp<=f_limit
        bus_r_pass=br<=r_limit
        metric_split=(primary_f_pass != bus_f_pass) || (primary_r_pass != bus_r_pass) ||
                     ((sgpeak<=f_limit) != (pllpeak<=f_limit)) ||
                     ((sgr<=r_limit) != (pllr<=r_limit))
        tailthreshold=max(0.02*fp,1e-9)
        settle=NaN
        tids=findall(after)
        for j in tids
            vals=primary[j:end,:]
            finitevals=abs.(vals[isfinite.(vals)])
            !isempty(finitevals) && maximum(finitevals)<=tailthreshold &&
                (settle=times[j]-1.0;break)
        end
        push!(rows,(;scenario=label,disturbance_bus=event_bus,delta_P_MW=mw,
            profile=duration==Inf ? "sustained" : "pulse_0.1s",
            horizon_s=horizon,dt_s=dt,solver_status=string(sol.retcode),
            primary_output="max local SG rotor / GFL PLL frequency deviation at buses 30:39",
            primary_peak_Hz=fp,primary_peak_bus=Int(fb),primary_peak_time_s=times[Int(ft)],
            primary_RoCoF_Hz_s=rp,primary_RoCoF_bus=Int(rb),primary_RoCoF_time_s=times[Int(rt)],
            SG_rotor_peak_Hz=sgpeak,SG_rotor_peak_bus=Int(sgbus),
            SG_rotor_RoCoF_Hz_s=sgr,PLL_peak_Hz=pllpeak,PLL_peak_bus=Int(pllbus),
            PLL_RoCoF_Hz_s=pllr,bus_frequency_peak_Hz=bp,bus_frequency_peak_bus=Int(bb),
            bus_frequency_peak_time_s=times[Int(bt)],bus_frequency_RoCoF_Hz_s=br,
            Fmax_pass=primary_f_pass,Rmax_pass=primary_r_pass,
            SG_decision=(sgpeak<=f_limit && sgr<=r_limit),
            PLL_decision=(pllpeak<=f_limit && pllr<=r_limit),
            bus_frequency_decision=(bus_f_pass && bus_r_pass),
            metric_decision_split=metric_split,primary_settling_s=settle))
        traj=DataFrame(time_s=times)
        for i in 1:10
            b=29+i
            traj[!,Symbol("bus$(b)_SG_df_Hz")]=freq.sg[:,i]
            traj[!,Symbol("bus$(b)_PLL_df_Hz")]=freq.pll[:,i]
            traj[!,Symbol("bus$(b)_bus_df_Hz")]=freq.busfreq[:,i]
            traj[!,Symbol("bus$(b)_SG_RoCoF_Hz_s")]=freq.sg_rocof[:,i]
            traj[!,Symbol("bus$(b)_PLL_RoCoF_Hz_s")]=freq.pll_rocof[:,i]
            traj[!,Symbol("bus$(b)_bus_RoCoF_Hz_s")]=freq.busrocof[:,i]
        end
        traj[!,:scenario]=fill(label,nrow(traj))
        push!(trajectories,traj)
        println("Q1A_CASE=",label," primary_peak=",fp," primary_RoCoF=",rp,
            " PLL_peak=",pllpeak," bus_peak=",bp," bus_RoCoF=",br," split=",metric_split)
    end
    table=DataFrame(rows)
    CSV.write(joinpath(OUT,"TABLE_Q1A_frequency_metric_cases.csv"),table)
    CSV.write(joinpath(OUT,"TABLE_Q1A_frequency_trajectories.csv"),vcat(trajectories...))
    cutoff=savgol_cutoff_hz(dt;half_window=5,degree=3)
    decision_split=any(table.metric_decision_split)
    status=decision_split ? "METRIC_DEPENDENT" : "CONSISTENT"
    result=Dict("stage"=>"Q1A","status"=>status,"candidate_sha256"=>candidate_sha,
        "frequency_output_definition"=>"max over all available local SG rotor and GFL PLL frequency deviations at generator buses 30:39",
        "rotor_definition"=>"60*(omega_pu-1) Hz",
        "pll_definition"=>"Delta_omega_rad_s/(2*pi) Hz; state is relative to nominal rotating frame",
        "bus_definition"=>"d(unwrapped arg(V_bus))/dt/(2*pi), Savitzky-Golay derivative, 11 samples at 100 Hz, cubic local polynomial",
        "bus_filter_support_s"=>0.1,"bus_filter_minus3dB_derivative_bandwidth_Hz"=>cutoff,
        "decision_metric_split"=>decision_split,"case_count"=>nrow(table),
        "amplitude_MW"=>[1.0,10.0,100.0,100.0],"pulse_duration_s"=>0.1,
        "PowerDynamics_calls_in_design"=>0,"independent_PD_TDS_runs"=>length(specs))
    ExpP.write_json(joinpath(OUT,"Q1A_RESULTS.json"),result)
    write(joinpath(OUT,"STAGE_SUMMARY.md"),"""
    # Q1A — Frequency-output definitions and metric audit

    Status: **$(status)**. Frozen ExpP candidate only; no design variable was changed.

    - Primary metric frozen for Q: maximum absolute deviation among every present SG rotor-speed output and every present local GFL PLL-frequency output on generator buses 30–39. Rotor speed uses `60(ω−1)` Hz; the installed PLL state is rad/s, converted by `Δω/(2π)`.
    - Independent bus diagnostic: derivative of unwrapped bus-voltage phase, using a cubic 11-sample Savitzky–Golay local polynomial at 100 Hz. Window support 0.1 s; measured -3 dB derivative-response bandwidth `$(cutoff)` Hz. It injects no plant signal.
    - PD traces cover sustained 1, 10, and 100 MW steps and a 100 MW, 0.1 s bus-16 pulse. Per-family peaks, RoCoF, limit decisions, and bus/time locations are tabulated.
    - Frequency metric classification is `$(status)` because frequency families $(decision_split ? "produce at least one different inherited-limit decision" : "have the same pass/fail decision in every tested scenario"). The primary rotor/PLL metric is fixed by the specification; bus frequency remains an independent diagnostic.
    """)
end

main()
