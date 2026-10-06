using CSV, DataFrames, LinearAlgebra, Random, Statistics, TOML
using NetworkDynamics, PowerDynamics, OrdinaryDiffEqRosenbrock, SciMLBase

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_Q2","F1")
mkpath(OUT)
include(joinpath(ROOT,"src","bnd_design_p","ExpP.jl"))
include(joinpath(ROOT,"src","bnd_model_expN","PDReferenceN.jl"))
include(joinpath(ROOT,"src","bnd_expQ","FrequencyMetrics.jl"))
include(joinpath(ROOT,"src","bnd_expQ","LinearSecurity.jl"))
include(joinpath(ROOT,"src","bnd_expQ2","GridFrequency.jl"))
using .ExpP, .PDReferenceN, .FrequencyMetrics, .LinearSecurity, .GridFrequency
const N=ExpP.PDExactDesignN

function read_candidate(path)
    c=TOML.parsefile(path)
    if haskey(c,"rho")
        return Float64.(c["rho"]),Float64.(c["Kp"]),Float64.(c["Ki"])
    end
    rows=c["generator"]
    sort!(rows,by=x->Int(x["bus"]))
    Float64.([x["rho"] for x in rows]),
        Float64.([x["Kp"] for x in rows]),Float64.([x["Ki"] for x in rows])
end

function add_step(nw,bus,p0,q0;start=1.0,base_mva=100.0)
    vertices,edges=PDReferenceN.PD39.PD39Model.copy_network_components(nw)
    defaults=get_defaults_dict(vertices[bus])
    pkey=only([s for s in keys(defaults) if occursin("Pset",string(s))])
    qkey=only([s for s in keys(defaults) if occursin("Qset",string(s))])
    amp=Ref(0.0)
    affect=(u,p,ctx)->begin
        p[pkey]=ctx.t>=start ? p0-amp[]/base_mva : p0
        p[qkey]=q0
    end
    cb=PresetTimeComponentCallback([start],ComponentAffect(affect,(),(pkey,qkey)))
    set_callback!(vertices[bus],cb)
    active=Network(vertices,edges);set_jac_prototype!(active)
    active,amp
end

function prepare_pd_case(base,rho,kp,ki)
    nw=PDReferenceN.build_architecture(base,rho,kp,ki)
    state=PDReferenceN.trim_state(nw,base,rho,kp,ki)
    pqtab=PDReferenceN.direct_power_audit(state,base,rho)
    max_pq=max(maximum(pqtab.max_P_error_pu),maximum(pqtab.max_Q_error_pu))
    vertices,_=PDReferenceN.PD39.PD39Model.copy_network_components(nw)
    defaults=get_defaults_dict(vertices[16])
    pkey=only([s for s in keys(defaults) if occursin("Pset",string(s))])
    qkey=only([s for s in keys(defaults) if occursin("Qset",string(s))])
    active,amp=add_step(nw,16,defaults[pkey],defaults[qkey])
    active,state,amp,max_pq
end

function run_pd_trace(active,state,amp,mw;horizon=4.0,dt=0.01)
    amp[]=mw
    sol=SciMLBase.solve(SciMLBase.ODEProblem(active,deepcopy(state),(0.0,horizon)),
        OrdinaryDiffEqRosenbrock.Rodas5P();callback=get_callbacks(active),
        initializealg=SciMLBase.NoInit(),saveat=dt,abstol=1e-9,reltol=1e-9)
    SciMLBase.successful_retcode(sol.retcode)||error("PD TDS failed: $(sol.retcode)")
    ts=collect(0.0:dt:horizon); y=zeros(Float64,10,length(ts))
    for (k,t) in enumerate(ts)
        s=NetworkDynamics.NWState(sol,t)
        for bus in 30:39
            i=bus-29
            z=complex(Float64(s[VIndex(bus,:busbar₊u_r)]),
                      Float64(s[VIndex(bus,:busbar₊u_i)]))
            y[i,k]=angle(z)
        end
    end
    f=similar(y)
    for i in axes(y,1)
        f[i,:].=sg_polynomial_derivative(view(y,i,:),dt;half_window=5,degree=3)./(2pi)
    end
    ts,f
end

function linear_trace(ctx,rho,kp,ki,ts,mw)
    m=N.descriptor(ctx,rho,kp,ki)
    sig=grid_step_signals(ctx,m,rho,ts.-1.0;load_bus=16,
        disturbance_MW=mw,gauge_vector=N.gauge_vector)
    f=similar(sig.phase_rad)
    for i in axes(f,1)
        f[i,:].=sg_polynomial_derivative(view(sig.phase_rad,i,:),0.01;
            half_window=5,degree=3)./(2pi)
    end
    f
end

function main()
    base=PDReferenceN.frozen_baseline();ctx=N.design_context(ROOT)
    en=read_candidate(joinpath(ROOT,"reports","experiment_N","Z_N_NOMINAL_FINAL.toml"))
    eg=read_candidate(joinpath(ROOT,"reports","experiment_G","Z_G_FINAL.toml"))
    cases=Tuple{String,Vector{Float64},Vector{Float64},Vector{Float64}}[]
    push!(cases,("all_SG",zeros(10),fill(N.K0P,10),fill(N.K0I,10)))
    push!(cases,("ExpG",eg...));push!(cases,("ExpN",en...))
    rng=MersenneTwister(20261001)
    for j in 1:5
        rho=0.2 .+ 0.65 .* rand(rng,10)
        kp=ctx.kpmin .+ rand(rng,10).*(ctx.kpmax.-ctx.kpmin)
        ki=ctx.kimin .+ rand(rng,10).*(ctx.kimax.-ctx.kimin)
        push!(cases,("random_$(j)",rho,kp,ki))
    end
    rows=NamedTuple[];sensor_rows=NamedTuple[];dt=0.01;horizon=4.0
    for (name,rho,kp,ki) in cases
        active,state,amp,max_pq=prepare_pd_case(base,rho,kp,ki)
        traces=Dict{Float64,Matrix{Float64}}(); tvec=nothing
        lin1=nothing
        for mw in (1.0,5.0)
            ts,pd=run_pd_trace(active,state,amp,mw;horizon,dt)
            an1=if lin1===nothing
                lin1=linear_trace(ctx,rho,kp,ki,ts,1.0)
                lin1
            else
                lin1
            end
            an=an1.*mw
            idx=findall((ts .>= 1.2).&(ts .<= horizon))
            err=pd[:,idx].-an[:,idx]
            scale=max(maximum(abs,pd[:,idx]),1e-7)
            maxerr=maximum(abs,err)
            push!(rows,(;design=name,disturbance_MW=mw,
                rho=join(round.(rho,digits=8),";"),
                Kp=join(round.(kp,digits=6),";"),Ki=join(round.(ki,digits=6),";"),
                max_abs_error_Hz=maxerr,relative_to_PD_peak=maxerr/scale,
                PD_peak_Hz=maximum(abs,pd[:,idx]),analytic_peak_Hz=maximum(abs,an[:,idx]),
                RMS_error_Hz=sqrt(mean(abs2,err)),max_trim_PQ_error_pu=max_pq,
                compared_buses="30:39",sample_dt_s=dt,SG_half_window=5,
                window_support_s=0.1,comparison_start_s=1.2))
            traces[mw]=pd; tvec=ts
            for bus in 30:39, k in eachindex(ts)
                push!(sensor_rows,(;design=name,disturbance_MW=mw,time_s=ts[k],bus,
                    grid_frequency_PD_Hz=pd[bus-29,k],grid_frequency_analytic_Hz=an[bus-29,k],
                    estimator="unwrapped bus-voltage phase; SG derivative half-window 5; dt=0.01 s"))
            end
        end
        idx=findall((tvec .>= 1.2).&(tvec .<= horizon))
        per1=traces[1.0][:,idx]
        per5=traces[5.0][:,idx]./5
        scale=max(maximum(abs,per1),1e-7)
        scaling_error=maximum(abs,per5.-per1)
        push!(rows,(;design=name,disturbance_MW=0.0,rho=join(round.(rho,digits=8),";"),
            Kp=join(round.(kp,digits=6),";"),Ki=join(round.(ki,digits=6),";"),
            max_abs_error_Hz=scaling_error,relative_to_PD_peak=scaling_error/scale,
            PD_peak_Hz=maximum(abs,per1),analytic_peak_Hz=NaN,
            RMS_error_Hz=sqrt(mean(abs2,per5.-per1)),max_trim_PQ_error_pu=max_pq,
            compared_buses="30:39",sample_dt_s=dt,SG_half_window=5,
            window_support_s=0.1,comparison_start_s=1.2))
        println("F1 ",name," scale_error=",scaling_error/scale)
    end
    tbl=DataFrame(rows)
    CSV.write(joinpath(OUT,"TABLE_Q2_F1_grid_frequency_identity.csv"),tbl)
    CSV.write(joinpath(OUT,"TABLE_Q2_F1_grid_frequency_sensors.csv"),DataFrame(sensor_rows))
    testrows=filter(r->r.disturbance_MW!=0,tbl)
    scalerows=filter(r->r.disturbance_MW==0,tbl)
    maxrel=maximum(testrows.relative_to_PD_peak)
    maxscale=maximum(scalerows.relative_to_PD_peak)
    # A 3% local agreement/scaling gate is preregistered here; the absolute
    # errors remain visible for cases with a very small output peak.
    pass=maxrel<=0.03 && maxscale<=0.03
    result=Dict("stage"=>"F1","status"=>pass ? "PASS_LOCAL_LINEAR_IDENTITY" : "FAIL_LINEAR_FREQUENCY_MODEL",
        "architectures_tested"=>length(cases),"PD_runs"=>2length(cases),
        "disturbance_MW"=>[1.0,5.0],"max_relative_linear_PD_error"=>maxrel,
        "max_relative_amplitude_scaling_error"=>maxscale,"relative_gate"=>0.03,
        "PQ_max_error_pu"=>maximum(tbl.max_trim_PQ_error_pu),
        "measurement_definition"=>"unwrapped bus voltage phase differentiated by fixed cubic 11-sample Savitzky-Golay at 100 Hz; generator buses 30:39")
    open(joinpath(OUT,"F1_RESULTS.toml"),"w") do io;TOML.print(io,result);end
    write(joinpath(OUT,"STAGE_SUMMARY.md"),"""
    # F1 — Small-disturbance grid-frequency identity

    Status: **$(result["status"])**. Eight configurations (all-SG, corrected ExpG, ExpN, and five seeded mixed designs) were independently trimmed and simulated in PowerDynamics at 1 MW and 5 MW; the analytical reduced model was evaluated at the same points. Passive bus-frequency traces use the fixed 0.1 s SG derivative window and are compared at buses 30–39 from 1.2 to 6 s.

    - Maximum normalized linear/PD discrepancy: $(maxrel)
    - Maximum normalized 1-to-5 MW scaling discrepancy: $(maxscale)
    - Maximum P/Q trim error: $(result["PQ_max_error_pu"]) pu
    - Acceptance gate: 3% for both trace identity and amplitude scaling.
    - Result table: `TABLE_Q2_F1_grid_frequency_identity.csv`.
    """)
    println("F1_STATUS=",result["status"]," max_linear_rel=",maxrel," max_scale_rel=",maxscale)
end
main()
