# Independent acceptance uses the compiled PowerDynamics DAE, not ReducedDAE.
include(joinpath(@__DIR__,"..","codesign_validation_20261001","validate_pd.jl"))
include(joinpath(@__DIR__,"PDPhysicalReference.jl"))
function validate_final()
    candidate=abspath(ARGS[1]);d=TOML.parsefile(candidate)
    rho=Float64.(d["rho"]);kp=Float64.(d["Kp"]);ki=Float64.(d["Ki"])
    label=length(ARGS)>1 ? ARGS[2] : "PD_final"
    dest=joinpath(ROOT,"reports","nonlinear_codesign_20261001",label);mkpath(dest)
    reference=get(d,"dc_convention","legacy")=="physical_supply" ? PDPhysicalReference : PDRef
    base=reference.frozen_baseline();nw=reference.build_architecture(base,rho,kp,ki)
    state=reference.trim_state(nw,base,rho,kp,ki);audit=reference.residual_audit(nw,state)
    println("FINAL_PD_TRIM ",audit);flush(stdout)
    rows=NamedTuple[];dt=.005;T=.5;lag=round(Int,T/dt)
    cases=length(ARGS)>2 && ARGS[3]=="smoke" ? [(8,100.)] : [(b,d) for b in (8,16,29) for d in (100.,-100.)]
    for (bus,delta) in cases
        println("FINAL_PD_START ",bus," ",delta);flush(stdout)
        active=event_network(nw,bus,delta);started=time()
        watch=DiscreteCallback((u,t,int)->time()-started>180.,int->terminate!(int);save_positions=(false,false))
        elapsed=@elapsed sol=SciMLBase.solve(SciMLBase.ODEProblem(active,state,(0.,61.)),Rodas5P();
            callback=CallbackSet(get_callbacks(active),watch),initializealg=SciMLBase.NoInit(),
            saveat=dt,abstol=1e-9,reltol=1e-9,maxiters=250000)
        ok=SciMLBase.successful_retcode(sol.retcode) && sol.t[end]>=61.0-1e-8
        if !ok
            open(joinpath(dest,"incomplete_$(bus)_$(Int(delta)).toml"),"w") do io
                TOML.print(io,Dict("complete"=>false,"last_time"=>sol.t[end],"retcode"=>string(sol.retcode)))
            end
            println("FINAL_PD_INCOMPLETE ",bus," ",delta);flush(stdout);continue
        end
        times=collect(0.:dt:61.);phase=zeros(39,length(times));volt=similar(phase)
        previous=zeros(39);initial=zeros(39);unwrap=zeros(39)
        for (j,t) in enumerate(times)
            ss=NWState(sol,abs(t-1.)<1e-12 ? t+1e-9 : t)
            for b in 1:39
                ur=Float64(ss[VIndex(b,:busbar₊u_r)]);ui=Float64(ss[VIndex(b,:busbar₊u_i)])
                raw=atan(ui,ur)
                if j==1;previous[b]=raw;initial[b]=raw;unwrap[b]=raw
                else;unwrap[b]+=mod(raw-previous[b]+pi,2pi)-pi;previous[b]=raw;end
                phase[b,j]=unwrap[b]-initial[b];volt[b,j]=hypot(ur,ui)
            end
        end
        f=zeros(size(phase));r=zeros(size(phase))
        for j in eachindex(times),b in 1:39
            a=phase[b,j];p1=j>lag ? phase[b,j-lag] : 0.;p2=j>2lag ? phase[b,j-2lag] : 0.
            f[b,j]=(a-p1)/(2pi*T);r[b,j]=(a-2p1+p2)/(2pi*T^2)
        end
        fi=argmax(abs.(f));ri=argmax(abs.(r))
        row=(;bus,delta,complete=ok,Fpeak_Hz=maximum(abs,f),Rpeak_Hz_s=maximum(abs,r),
            F_bus=fi[1],F_time=times[fi[2]]-1,R_bus=ri[1],R_time=times[ri[2]]-1,
            Vmin=minimum(volt),Vmax=maximum(volt),runtime_s=elapsed,
            F_pass=maximum(abs,f)<=.5,R_pass=maximum(abs,r)<=.5,
            V_pass=minimum(volt)>=.9 && maximum(volt)<=1.1)
        push!(rows,row);CSV.write(joinpath(dest,"events.csv"),DataFrame(rows))
        trace=DataFrame(time_after_event_s=times.-1,Fmax_Hz=vec(maximum(abs.(f),dims=1)),
            Rmax_Hz_s=vec(maximum(abs.(r),dims=1)),Vmin=vec(minimum(volt,dims=1)),Vmax=vec(maximum(volt,dims=1)))
        CSV.write(joinpath(dest,"trajectory_$(bus)_$(Int(delta)).csv"),trace)
        println("FINAL_PD_RESULT ",row);flush(stdout)
    end
    open(joinpath(dest,"provenance.toml"),"w") do io
        TOML.print(io,Dict("candidate_sha256"=>bytes2hex(sha256(read(candidate))),
            "trim_residual"=>audit.maximum,"completed_cases"=>length(rows),"monitored_buses"=>collect(1:39),
            "scope"=>"Independent compiled PowerDynamics DAE; finite horizon and numerical resolution only. No current/energy or infinite-horizon stability certificate."))
    end
end
validate_final()
