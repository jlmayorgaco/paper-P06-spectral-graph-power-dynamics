using CSV, DataFrames, TOML, LinearAlgebra
using NetworkDynamics, PowerDynamics, OrdinaryDiffEqRosenbrock, SciMLBase

module Single
include(joinpath(@__DIR__,"run_event_screen.jl"))
end

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
BLAS.set_num_threads(1)

function actuator_metric(sol,ts,rho)
    slack=Inf;winning_bus=0;winning_time=NaN;winning_state=""
    for t in ts
        t<1.0 && continue
        te=abs(t-1.0)<1e-12 ? t+1e-9 : t
        ss=abs(t-1.0)<1e-12 ? Single.consistent_event_observation(sol,te) : NWState(sol,te)
        for b in 30:38
            i=b-29
            rho[i]<1 || continue
            xg=Float64(ss[VIndex(b,:ctrld_gen₊gov₊xg1)])
            lo=Float64(ss[VIndex(b,:ctrld_gen₊gov₊V_min)])
            hi=Float64(ss[VIndex(b,:ctrld_gen₊gov₊V_max)])
            vr=Float64(ss[VIndex(b,:ctrld_gen₊avr₊vr)])
            vlo=Float64(ss[VIndex(b,:ctrld_gen₊avr₊vr_min)])
            vhi=Float64(ss[VIndex(b,:ctrld_gen₊avr₊vr_max)])
            vals=((xg-lo)/(hi-lo),(hi-xg)/(hi-lo),(vr-vlo)/(vhi-vlo),(vhi-vr)/(vhi-vlo))
            v,ix=findmin(vals)
            if v<slack
                slack=v;winning_bus=b;winning_time=t
                winning_state=("governor_lower","governor_upper","AVR_lower","AVR_upper")[ix]
            end
        end
    end
    (;slack,winning_bus,winning_time,winning_state)
end

function main()
    probe_name=length(ARGS)==0 ? "L3_FROZEN_EVENT_GRADIENT_PROBES.csv" : ARGS[1]
    probes=CSV.read(joinpath(@__DIR__,probe_name),DataFrame)
    base=Single.P.frozen_baseline()
    outpath=joinpath(@__DIR__,"L3_EVENT_ACTUATOR_SENSITIVITY.csv")
    rows=isfile(outpath) ? [NamedTuple(r) for r in eachrow(CSV.read(outpath,DataFrame))] : NamedTuple[]
    completed=Set(String(r.candidate_id) for r in rows)
    dt=0.01;tspan=(0.0,61.0);ts=collect(0.0:dt:61.0)
    for p in eachrow(probes)
        id=String(p.candidate_id)
        id in completed && continue
        path=joinpath(@__DIR__,String(p.path))
        d=TOML.parsefile(path)
        rho,kp,ki=Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"])
        println("EVENT_GRADIENT_START ",id);flush(stdout)
        nw=Single.P.build_architecture(base,rho,kp,ki)
        state=Single.P.trim_state(nw,base,rho,kp,ki)
        trim=Single.P.residual_audit(nw,state)
        trim.maximum<1e-8 || error("equilibrium residual failed for $id")
        active=Single.event_network(nw,16,100.0)
        started=time();lastlog=Ref(started)
        watch=DiscreteCallback((u,t,int)->time()-lastlog[]>20,int->begin
            println("EVENT_GRADIENT_PROGRESS ",id," t=",int.t);flush(stdout)
            lastlog[]=time()
            time()-started>240 ? terminate!(int) : SciMLBase.u_modified!(int,false)
        end;save_positions=(false,false))
        sol=solve(ODEProblem(active,state,tspan),Rodas5P();
            callback=CallbackSet(get_callbacks(active),watch),
            initializealg=SciMLBase.NoInit(),saveat=dt,
            abstol=1e-9,reltol=1e-9,maxiters=1_000_000)
        complete=SciMLBase.successful_retcode(sol.retcode) && sol.t[end]>=61-1e-8
        metric=complete ? actuator_metric(sol,ts,rho) :
            (;slack=NaN,winning_bus=0,winning_time=NaN,winning_state="INCOMPLETE")
        push!(rows,(;candidate_id=id,gain_kind=p.gain_kind,bus=p.bus,
            delta_log_gain=p.delta_log_gain,complete,retcode=string(sol.retcode),
            actuator_slack=metric.slack,winning_bus=metric.winning_bus,
            winning_time_s=metric.winning_time,winning_state=metric.winning_state,
            residual_inf=trim.maximum,runtime_s=time()-started,
            status=complete ? "ONE_EVENT_ACTUATOR_GRADIENT_PROBE_ONLY" : "INCOMPLETE"))
        CSV.write(outpath,DataFrame(rows))
        println("EVENT_GRADIENT_DONE ",last(rows));flush(stdout)
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
