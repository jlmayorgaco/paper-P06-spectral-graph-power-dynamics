include("setup.jl")
using NetworkDynamics,PowerDynamics,OrdinaryDiffEqRosenbrock,SciMLBase
include(joinpath(ROOT,"src","bnd_model_expN","PDReferenceN.jl"))
const PDRef=PDReferenceN

function event_network(nw,bus,deltaP;duration=Inf)
    mdl=PDRef.PD39.PD39Model
    vertices,edges=mdl.copy_network_components(nw)
    defaults=get_defaults_dict(vertices[bus])
    pk=only([s for s in keys(defaults) if occursin("Pset",string(s))])
    qk=only([s for s in keys(defaults) if occursin("Qset",string(s))])
    p0=defaults[pk];q0=defaults[qk]
    affect=(u,p,ctx)->begin
      p[pk]=p0-(ctx.t<1+duration ? deltaP/100 : 0.0);p[qk]=q0
    end
    breaks=isfinite(duration) ? [1.0,1.0+duration] : [1.0]
    set_callback!(vertices[bus],PresetTimeComponentCallback(breaks,ComponentAffect(affect,(),(pk,qk))))
    active=Network(vertices,edges);set_jac_prototype!(active);active
end

"""Locate the installed hard anti-windup boundaries without changing limits.
The installed vector field uses strict >/<. A tiny outward numerical nudge
selects the saturated side at an outward crossing; inward crossings are free.
"""
function limiter_callbacks(nw,state,rho,events;guard=1e-10)
    symbols=string.(NetworkDynamics.SII.variable_symbols(nw));callbacks=Any[]
    for bus in 30:38
      rho[bus-29]<1 || continue
      for (kind,xname,lname,hname) in (("governor",:ctrld_gen₊gov₊xg1,:ctrld_gen₊gov₊V_min,:ctrld_gen₊gov₊V_max),
          ("AVR",:ctrld_gen₊avr₊vr,:ctrld_gen₊avr₊vr_min,:ctrld_gen₊avr₊vr_max))
        idx=findfirst(==(string(VIndex(bus,xname))),symbols)
        idx===nothing && error("limiter state missing")
        for (side,limit) in (("lower",Float64(state[VIndex(bus,lname)])),("upper",Float64(state[VIndex(bus,hname)])))
          let idx=idx,side=side,limit=limit,bus=bus,kind=kind
            condition=(u,t,int)->u[idx]-limit
            affect=int->begin
              int.u[idx]=limit+(side=="upper" ? 1 : -1)*guard*max(1,abs(limit))
              push!(events,(;time_s=Float64(int.t),bus,kind,side,limit,guard))
              println("LIMIT_EVENT ",last(events));flush(stdout)
              SciMLBase.u_modified!(int,true)
            end
            cb=side=="upper" ? ContinuousCallback(condition,affect,nothing;abstol=guard/100,reltol=0,save_positions=(false,false)) :
              ContinuousCallback(condition,nothing,affect;abstol=guard/100,reltol=0,save_positions=(false,false))
            push!(callbacks,cb)
          end
        end
      end
    end
    callbacks
end

function main()
    path=isempty(ARGS) ? joinpath(@__DIR__,"frozen","original_candidate.toml") : abspath(ARGS[1])
    label=length(ARGS)>1 ? ARGS[2] : "original"
    dest=joinpath(OUT,"PD_"*label);mkpath(dest)
    d=TOML.parsefile(path);sha=bytes2hex(sha256(read(path)))
    rho=Float64.(d["rho"]);kp=Float64.(d["Kp"]);ki=Float64.(d["Ki"])
    an=N.spectrum(CTX,rho,kp,ki)
    println("PD_BUILD ",label," sha=",sha);flush(stdout)
    base=PDRef.frozen_baseline();nw=PDRef.build_architecture(base,rho,kp,ki)
    state=PDRef.trim_state(nw,base,rho,kp,ki)
    trim=PDRef.residual_audit(nw,state);pq=PDRef.direct_power_audit(state,base,rho)
    CSV.write(joinpath(dest,"component_trim.csv"),pq)
    pl=linearize_network(state);lp=jacobian_eigenvals(pl)
    # Exactly one rotational symmetry is removed; all other finite poles remain.
    ig=argmin(abs.(lp));phys=lp[setdiff(eachindex(lp),[ig])]
    length(phys)==length(an.lambda) || error("physical spectrum dimension mismatch")
    rr,cc=ExpP.linear_assignment(abs.(phys.-transpose(an.lambda)))
    pe=maximum(abs.(phys[rr].-an.lambda[cc]));ap=maximum(real.(phys))
    spectrum=Dict("candidate_sha256"=>sha,"trim_residual"=>trim.maximum,
      "max_P_error_pu"=>maximum(pq.max_P_error_pu),"max_Q_error_pu"=>maximum(pq.max_Q_error_pu),
      "trim_bounds_pass"=>all(pq.bounds_status.=="PASS"),"analytic_alpha"=>an.alpha,
      "pd_alpha"=>ap,"max_matched_pole_error"=>pe,"removed_gauge_abs"=>abs(lp[ig]),
      "physical_poles"=>length(phys),"nominal_margin_pass"=>ap<=-0.05)
    open(joinpath(dest,"spectrum.toml"),"w") do io;TOML.print(io,spectrum);end
    println("PD_SPECTRUM ",spectrum);flush(stdout)
    cases=[(name="small_16",bus=16,d=1.0,duration=Inf,dt=.01,tol=1e-9),
      (name="step_16",bus=16,d=100.0,duration=Inf,dt=.01,tol=1e-9),
      (name="step_29",bus=29,d=100.0,duration=Inf,dt=.01,tol=1e-9)]
    if label!="original"
      append!(cases,[(name="step_8",bus=8,d=100.0,duration=Inf,dt=.01,tol=1e-9),
       (name="negative_16",bus=16,d=-100.0,duration=Inf,dt=.01,tol=1e-9),
       (name="stress_29",bus=29,d=150.0,duration=Inf,dt=.01,tol=1e-9),
       (name="pulse_16",bus=16,d=100.0,duration=.1,dt=.01,tol=1e-9),
       (name="refined_16",bus=16,d=100.0,duration=Inf,dt=.005,tol=1e-10),
       (name="step_8_relaxed",bus=8,d=100.0,duration=Inf,dt=.01,tol=1e-7)])
    end
    event_limits=length(ARGS)>2 && ARGS[3]=="event_limits"
    if event_limits
      cases=[(name="step_8_event_limits",bus=8,d=100.0,duration=Inf,dt=.01,tol=1e-9),
        (name="step_8_event_limits_refined",bus=8,d=100.0,duration=Inf,dt=.005,tol=1e-10)]
    end
    if length(ARGS)>2 && ARGS[3]=="refined_bus8"
      cases=[(name="refined_8",bus=8,d=100.0,duration=Inf,dt=.005,tol=1e-10)]
    end
    if length(ARGS)>2 && ARGS[3]=="small_inputs"
      cases=[(name="small_8",bus=8,d=1.0,duration=Inf,dt=.01,tol=1e-9),
             (name="small_29",bus=29,d=1.0,duration=Inf,dt=.01,tol=1e-9)]
    end
    if label=="repaired"
      sort!(cases;by=c->c.name=="step_8" ? 0 : 1,alg=MergeSort)
    end
    eventpath=joinpath(dest,"events.csv")
    rows=isfile(eventpath) ? NamedTuple.(eachrow(CSV.read(eventpath,DataFrame))) : NamedTuple[]
    done=Set(String(r.event) for r in rows)
    for case in cases
      case.name in done && continue
      case.name=="step_8_relaxed" && any(r->r.event=="step_8",rows) && continue
      println("PD_EVENT_START ",label," ",case.name);flush(stdout)
      active=event_network(nw,case.bus,case.d;duration=case.duration)
      dt=case.dt;tend=61.0
      started=time();lastlog=Ref(started)
      limitevents=NamedTuple[]
      limits=event_limits ? limiter_callbacks(active,state,rho,limitevents;guard=case.dt==.005 ? 1e-11 : 1e-10) : Any[]
      watch=DiscreteCallback((u,t,integ)->time()-lastlog[]>20,
        integ->begin
          println("PD_PROGRESS ",case.name," t=",integ.t," dt=",integ.dt," wall_s=",time()-started);flush(stdout)
          lastlog[]=time()
          if time()-started>180
            terminate!(integ)
          else
            SciMLBase.u_modified!(integ,false)
          end
        end;save_positions=(false,false))
      elapsed=@elapsed sol=SciMLBase.solve(SciMLBase.ODEProblem(active,state,(0.0,tend)),Rodas5P();
        callback=CallbackSet(get_callbacks(active),watch,limits...),initializealg=SciMLBase.NoInit(),saveat=dt,
        abstol=case.tol,reltol=case.tol,maxiters=250000)
      ok=SciMLBase.successful_retcode(sol.retcode) && sol.t[end]>=tend-1e-8
      if !isempty(limitevents);CSV.write(joinpath(dest,case.name*"_limit_events.csv"),DataFrame(limitevents));end
      if !ok
        failure=Dict("event"=>case.name,"retcode"=>string(sol.retcode),"last_time_s"=>sol.t[end],
          "runtime_s"=>elapsed,"candidate_sha256"=>sha,"status"=>"INCOMPLETE_TDS",
          "physical_instability_proven"=>false,"tolerance"=>case.tol)
        open(joinpath(dest,case.name*"_incomplete.toml"),"w") do io;TOML.print(io,failure);end
        println("PD_EVENT_INCOMPLETE ",failure);flush(stdout);continue
      end
      times=collect(0.0:dt:tend);nt=length(times)
      angles=zeros(10,nt);volt=zeros(10,nt);curr=zeros(10,nt);dc=zeros(10,nt)
      rotor=zeros(10,nt);pll=zeros(10,nt);kcl=zeros(10,nt)
      previous=zeros(10);theta=zeros(10);initial=zeros(10);i0=zeros(10)
      for (j,t) in enumerate(times)
        # Use the post-event side at exact parameter discontinuities.
        te=abs(t-1.0)<1e-12 || (isfinite(case.duration)&&abs(t-1-case.duration)<1e-12) ? t+1e-9 : t
        ss=NWState(sol,te)
        for bus in 30:39
          k=bus-29;ur=Float64(ss[VIndex(bus,:busbar₊u_r)]);ui=Float64(ss[VIndex(bus,:busbar₊u_i)])
          raw=atan(ui,ur)
          if j==1;previous[k]=raw;theta[k]=raw;initial[k]=raw
          else;theta[k]+=mod(raw-previous[k]+pi,2pi)-pi;previous[k]=raw;end
          angles[k,j]=theta[k]-initial[k];volt[k,j]=hypot(ur,ui)
          ir=Float64(ss[VIndex(bus,:gfl₊filter₊i_f_r)]);ii=Float64(ss[VIndex(bus,:gfl₊filter₊i_f_i)])
          j==1 && (i0[k]=hypot(ir,ii));curr[k,j]=hypot(ir,ii)/i0[k]
          dc[k,j]=Float64(ss[VIndex(bus,:gfl₊v_dc_state)])
          prefix=bus==39 ? "machine₊" : "ctrld_gen₊machine₊"
          rotor[k,j]=rho[k]<1 ? 60*(Float64(ss[VIndex(bus,Symbol(prefix*"ω"))])-1) : NaN
          pll[k,j]=Float64(ss[VIndex(bus,:gfl₊pll₊Δω_rad_s)])/(2pi)
          sg=rho[k]<1 ? Float64(ss[VIndex(bus,Symbol(prefix*"Sn"))])*Float64(ss[VIndex(bus,Symbol(prefix*"P"))]) : 0.0
          gf=100rho[k]*(ur*ir+ui*ii)
          load=bus in (31,39) ? 100Float64(ss[VIndex(bus,:ZIPLoad₊P)]) : 0.0
          kcl[k,j]=sg+gf+load-Float64(ss[VIndex(bus,:busbar₊P_MW)])
        end
      end
      nlag=Int(round(4/dt));nevent=Int(round(1/dt))
      tm=collect(-nlag:Int(round(60/dt))).*dt
      phase=hcat(zeros(10,nlag-nevent),angles)
      steady=(phase[:,end]-phase[:,end-Int(round(10/dt))])/(20pi)
      sig=(times_s=tm,phase_rad=phase,F_inf_Hz=steady)
      metrics=FiniteWindow.window_metrics(sig;windows=(.2,.5,1.,2.),dt_s=dt,horizon_s=60.0)
      CSV.write(joinpath(dest,case.name*"_windows.csv"),DataFrame(metrics))
      m=only(filter(x->x.window_s==.5,metrics))
      lag=Int(round(.5/dt));freq=zeros(10,length(tm));roc=similar(freq);fill!(roc,0)
      z0=findfirst(==(0.0),tm)
      for j in z0:length(tm)
        p1=j-lag>=z0 ? phase[:,j-lag] : zeros(10)
        p2=j-2lag>=z0 ? phase[:,j-2lag] : zeros(10)
        freq[:,j]=(phase[:,j]-p1)/pi;roc[:,j]=(phase[:,j]-2p1+p2)/(pi/2)
      end
      # Compare small-signal waveforms using the corrected analytic phase.
      am=FiniteWindow.design_metrics(CTX,an.model,rho;load_bus=case.bus,disturbance_MW=case.d,
        windows=(.5,),dt_s=dt,horizon_s=60.0,gauge_vector=N.gauge_vector)
      af=zeros(10,length(tm));apad=hcat(zeros(10,nlag-Int(round(1/dt))),am.signals.phase_rad)
      for j in z0:length(tm)
        p1=j-lag>=z0 ? apad[:,j-lag] : zeros(10)
        af[:,j]=(apad[:,j]-p1)/pi
      end
      small_error=isfinite(case.duration) ? NaN : maximum(abs.(freq-af))/max(maximum(abs,af),1e-12)
      lastfreq=maximum(abs,steady)
      taildrift=maximum(abs.((phase[:,end]-phase[:,end-Int(round(5/dt))])/(10pi)-steady))
      row=(;event=case.name,bus=case.bus,commanded_Pset_MW=case.d,duration_s=case.duration,
        Fpeak_Hz=m.F_peak_Hz,Rpeak_Hz_s=m.R_peak_Hz_s,last10s_Hz=lastfreq,
        last5_vs10_Hz=taildrift,F_pass=m.F_peak_Hz<=.5,R_pass=m.R_peak_Hz_s<=.5,
        Vmin=minimum(volt),Vmax=maximum(volt),current_ratio_max=maximum(curr),
        Vdc_min=minimum(dc),Vdc_max=maximum(dc),KCL_error_MW=maximum(abs,kcl),
        linear_relative_waveform_error=small_error,dt,solver_tolerance=case.tol,runtime_s=elapsed,
        retcode=string(sol.retcode))
      push!(rows,row);CSV.write(joinpath(dest,"events.csv"),DataFrame(rows))
      df=DataFrame(time_after_event_s=tm,Fmax_Hz=vec(maximum(abs.(freq),dims=1)),
        Rmax_Hz_s=vec(maximum(abs.(roc),dims=1)),Flinear_max_Hz=vec(maximum(abs.(af),dims=1)))
      for k in 1:10;df[!,Symbol("f_bus$(29+k)_Hz")]=freq[k,:];end
      CSV.write(joinpath(dest,case.name*"_trajectory.csv"),df)
      println("PD_EVENT_DONE ",label," ",row);flush(stdout)
    end
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    main()
end
