using CSV,DataFrames,TOML,SHA,LinearAlgebra,NetworkDynamics,PowerDynamics,OrdinaryDiffEqRosenbrock,SciMLBase
include(joinpath(@__DIR__,"..","nonlinear_codesign_20261001","PDPhysicalReference.jl"))
include(joinpath(@__DIR__,"..","bnd_expQ2B","certified_search","ConsistentEventObservation.jl"))
const P=PDPhysicalReference
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(@__DIR__,"baseline_reproduction")
const DEST=joinpath(OUT,"independent_pd")
mkpath(DEST);BLAS.set_num_threads(1)
say(args...)=(println(args...);flush(stdout))

function event_network(nw,bus,delta)
    vertices,edges=P.PD39.PD39Model.copy_network_components(nw)
    defaults=get_defaults_dict(vertices[bus]);keys0=collect(keys(defaults))
    pk=only(filter(s->occursin("Pset",string(s)),keys0))
    qk=only(filter(s->occursin("Qset",string(s)),keys0))
    p0=defaults[pk];q0=defaults[qk]
    affect=(u,p,ctx)->begin;p[pk]=p0-delta/100;p[qk]=q0;end
    set_callback!(vertices[bus],PresetTimeComponentCallback([1.],ComponentAffect(affect,(),(pk,qk))))
    active=Network(vertices,edges);set_jac_prototype!(active);active
end

function main()
    say("PD_BASE_BUILD")
    base=P.frozen_baseline();template=nothing;rows=NamedTuple[]
    labels=isempty(ARGS) ? ["baseline","joint_final","fixed_gains_final"] : ARGS
    for label in labels
        # Checkpoint files are updated during optimization. Consume them only
        # after the main run freezes both arms and completes its source audit.
        path=joinpath(OUT,label*".toml");d=TOML.parsefile(path)
        d["dc_convention"]=="physical_supply" || error("Wrong DC contract")
        rho=Float64.(d["rho"]);kp=Float64.(d["Kp"]);ki=Float64.(d["Ki"])
        all(0 .<rho.<1) || error("This validation reuses one interior architecture")
        say("PD_CANDIDATE_BUILD ",label)
        template===nothing && (template=P.build_architecture(base,rho,kp,ki))
        nw=template;state=P.trim_state(nw,base,rho,kp,ki)
        audit=P.residual_audit(nw,state);audit.maximum<1e-8 || error("PD trim failed")
        pq=P.direct_power_audit(state,base,rho);CSV.write(joinpath(DEST,label*"_trim.csv"),pq)
        all(pq.bounds_status.=="PASS") || error("PD trim limits failed")
        maximum(pq.max_P_error_pu)<1e-9 || error("PD sharing mismatch")
        linear=linearize_network(state);poles=jacobian_eigenvals(linear)
        iz=argmin(abs.(poles));abs(poles[iz])<1e-5 || error("Gauge not found")
        physical=poles[setdiff(eachindex(poles),[iz])];alpha=maximum(real.(physical))
        CSV.write(joinpath(DEST,label*"_poles.csv"),DataFrame(real=real.(poles),imag=imag.(poles)))
        active=event_network(nw,8,100.);dt=.005;started=time();lastlog=Ref(started)
        watch=DiscreteCallback((u,t,i)->time()-lastlog[]>20,i->begin
            say("PD_PROGRESS ",label," t=",i.t," wall=",time()-started);lastlog[]=time()
            time()-started>240 ? terminate!(i) : SciMLBase.u_modified!(i,false)
        end;save_positions=(false,false))
        say("PD_EVENT_START ",label," alpha=",alpha)
        sol=solve(ODEProblem(active,state,(0.,61.)),Rodas5P();
            callback=CallbackSet(get_callbacks(active),watch),initializealg=SciMLBase.NoInit(),
            saveat=dt,abstol=1e-10,reltol=1e-10,maxiters=1000000)
        ok=SciMLBase.successful_retcode(sol.retcode)&&sol.t[end]>=61-1e-8
        ok || error("PD incomplete $label: $(sol.retcode) t=$(sol.t[end])")
        ts=collect(0.:dt:61.);nt=length(ts);phase=zeros(39,nt);voltage=similar(phase)
        previous=zeros(39);initial=zeros(39);angle=zeros(39)
        for (j,t) in enumerate(ts)
            ss=abs(t-1)<1e-12 ? consistent_event_observation(sol,t+1e-9) : NWState(sol,t)
            for b in 1:39
                ur=Float64(ss[VIndex(b,:busbar₊u_r)]);ui=Float64(ss[VIndex(b,:busbar₊u_i)])
                raw=atan(ui,ur)
                if j==1
                    previous[b]=raw;initial[b]=raw;angle[b]=raw
                else
                    angle[b]+=mod(raw-previous[b]+pi,2pi)-pi;previous[b]=raw
                end
                phase[b,j]=angle[b]-initial[b];voltage[b,j]=hypot(ur,ui)
            end
        end
        lag=round(Int,.5/dt);freq=zeros(size(phase));roc=similar(freq)
        for j in 1:nt,b in 1:39
            p1=j>lag ? phase[b,j-lag] : 0.;p2=j>2lag ? phase[b,j-2lag] : 0.
            freq[b,j]=(phase[b,j]-p1)/pi
            roc[b,j]=(phase[b,j]-2p1+p2)/(pi/2)
        end
        fi=argmax(abs.(freq));ri=argmax(abs.(roc))
        row=(;label,sha256=bytes2hex(sha256(read(path))),complete=ok,alpha,
            trim_residual=audit.maximum,F=maximum(abs,freq),R=maximum(abs,roc),
            Fbus=fi[1],Ftime=ts[fi[2]]-1,Rbus=ri[1],Rtime=ts[ri[2]]-1,
            Vmin=minimum(voltage),Vmax=maximum(voltage),elapsed_s=time()-started)
        push!(rows,row);CSV.write(joinpath(DEST,"events.csv"),DataFrame(rows))
        CSV.write(joinpath(DEST,label*"_trajectory.csv"),DataFrame(time_after_event_s=ts.-1,
            Fmax_Hz=vec(maximum(abs.(freq),dims=1)),Rmax_Hz_s=vec(maximum(abs.(roc),dims=1)),
            Vmin_pu=vec(minimum(voltage,dims=1)),Vmax_pu=vec(maximum(voltage,dims=1))))
        say("PD_DONE ",row)
    end
end

try
    main()
catch err
    open(joinpath(DEST,"failure.toml"),"w") do io
        TOML.print(io,Dict("error"=>sprint(showerror,err)))
    end
    rethrow()
end


