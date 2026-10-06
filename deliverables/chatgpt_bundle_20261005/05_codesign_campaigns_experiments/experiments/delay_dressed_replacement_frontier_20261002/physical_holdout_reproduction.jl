using CSV,DataFrames,TOML,SHA,LinearAlgebra,NetworkDynamics,PowerDynamics,OrdinaryDiffEqRosenbrock,SciMLBase
include(joinpath(@__DIR__,"..","nonlinear_codesign_20261001","PDPhysicalReference.jl"))
include(joinpath(@__DIR__,"..","bnd_expQ2B","certified_search","ConsistentEventObservation.jl"))
const P=PDPhysicalReference
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const DEST=joinpath(@__DIR__,"baseline_reproduction","independent_pd")
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

candidate=TOML.parsefile(joinpath(@__DIR__,"inputs","prior_joint_input.toml"))
rho=Float64.(candidate["rho"]);kp=Float64.(candidate["Kp"]);ki=Float64.(candidate["Ki"])
base=P.frozen_baseline();nw=P.build_architecture(base,rho,kp,ki);state=P.trim_state(nw,base,rho,kp,ki)
P.residual_audit(nw,state).maximum<1e-8 || error("physical equilibrium residual failed")
events=[("joint_bus8_minus100",8,-100.),("joint_bus16_plus100",16,100.),("joint_bus16_minus100",16,-100.),
        ("joint_bus29_plus100",29,100.),("joint_bus29_minus100",29,-100.)]
rows=NamedTuple[];dt=.005;ts=collect(0.:dt:61.);nt=length(ts)
for (label,bus,delta) in events
    active=event_network(nw,bus,delta);started=time();lastlog=Ref(started)
    watch=DiscreteCallback((u,t,i)->time()-lastlog[]>20,i->begin
        say("PD_HOLDOUT_PROGRESS ",label," t=",i.t," wall=",time()-started);lastlog[]=time()
        time()-started>240 ? terminate!(i) : SciMLBase.u_modified!(i,false)
    end;save_positions=(false,false))
    say("PD_HOLDOUT_START ",label)
    sol=solve(ODEProblem(active,state,(0.,61.)),Rodas5P();
        callback=CallbackSet(get_callbacks(active),watch),initializealg=SciMLBase.NoInit(),
        saveat=dt,abstol=1e-10,reltol=1e-10,maxiters=1000000)
    ok=SciMLBase.successful_retcode(sol.retcode)&&sol.t[end]>=61-1e-8
    ok || error("PD holdout incomplete $label: $(sol.retcode) t=$(sol.t[end])")
    phase=zeros(39,nt);voltage=similar(phase);previous=zeros(39);initial=zeros(39);angle=zeros(39)
    for (j,t) in enumerate(ts)
        ss=abs(t-1)<1e-12 ? consistent_event_observation(sol,t+1e-9) : NWState(sol,t)
        for b in 1:39
            ur=Float64(ss[VIndex(b,:busbar₊u_r)]);ui=Float64(ss[VIndex(b,:busbar₊u_i)])
            raw=atan(ui,ur)
            if j==1;previous[b]=raw;initial[b]=raw;angle[b]=raw
            else;angle[b]+=mod(raw-previous[b]+pi,2pi)-pi;previous[b]=raw;end
            phase[b,j]=angle[b]-initial[b];voltage[b,j]=hypot(ur,ui)
        end
    end
    lag=round(Int,.5/dt);freq=zeros(size(phase));roc=similar(freq)
    for j in 1:nt,b in 1:39
        p1=j>lag ? phase[b,j-lag] : 0.;p2=j>2lag ? phase[b,j-2lag] : 0.
        freq[b,j]=(phase[b,j]-p1)/pi;roc[b,j]=(phase[b,j]-2p1+p2)/(pi/2)
    end
    keep=findall(ts .>= 1.0-1e-12);fi=argmax(abs.(freq[:,keep]));ri=argmax(abs.(roc[:,keep]))
    row=(;label,bus,delta_load_MW=delta,complete=true,retcode=string(sol.retcode),last_time_s=sol.t[end],
        F_post_event_Hz=maximum(abs.(freq[:,keep])),Fbus=mod1(fi[1],39),Ftime_s=ts[keep[fi[2]]]-1,
        RoCoF_post_event_Hz_s=maximum(abs.(roc[:,keep])),Rbus=mod1(ri[1],39),Rtime_s=ts[keep[ri[2]]]-1,
        Vmin_post_event_pu=minimum(voltage[:,keep]),Vmax_post_event_pu=maximum(voltage[:,keep]),elapsed_s=time()-started)
    push!(rows,row);CSV.write(joinpath(DEST,"TABLE_D00_PD_HOLDOUTS.csv"),DataFrame(rows))
    CSV.write(joinpath(DEST,label*"_trajectory.csv"),DataFrame(time_after_event_s=ts.-1,
        Fmax_Hz=vec(maximum(abs.(freq),dims=1)),RoCoFmax_Hz_s=vec(maximum(abs.(roc),dims=1)),
        Vmin_pu=vec(minimum(voltage,dims=1)),Vmax_pu=vec(maximum(voltage,dims=1))))
    say("PD_HOLDOUT_DONE ",row)
end
