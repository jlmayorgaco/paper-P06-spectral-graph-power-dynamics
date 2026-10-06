using CSV, DataFrames, LinearAlgebra, SHA, TOML
using NetworkDynamics, PowerDynamics, OrdinaryDiffEqRosenbrock, SciMLBase
const ROOT=normpath(joinpath(@__DIR__,"..",".."));const CORE=joinpath(ROOT,"reports","experiment_Q2B","CORE")
const CANDIDATE=joinpath(CORE,"Z_Q2B_SECURE_T05_FINAL.toml")
const EXPECTED="f53250383cdb906a2d0cb3130ca9c1e72c63b472717c5c82aea1d5838ac785ce"
bytes2hex(sha256(read(CANDIDATE)))==EXPECTED || error("frozen candidate SHA mismatch")
include(joinpath(ROOT,"src","bnd_design_p","ExpP.jl"))
include(joinpath(ROOT,"src","bnd_expQ","LinearSecurity.jl"))
include(joinpath(ROOT,"src","bnd_expQ2B","FiniteWindow.jl"))
include(joinpath(ROOT,"src","bnd_model_expN","PDReferenceN.jl"))
const Ref=PDReferenceN;const PD39=PDReferenceN.PD39

function event_network(nw,bus,delta;start=1.0)
    vertices,edges=PD39.PD39Model.copy_network_components(nw)
    defaults=get_defaults_dict(vertices[bus]);pk=only([x for x in keys(defaults) if occursin("Pset",string(x))]);qk=only([x for x in keys(defaults) if occursin("Qset",string(x))])
    p0=defaults[pk];q0=defaults[qk];dp=Float64(delta)/100
    affect=(u,p,ctx)->begin;ctx.t>=start && (p[pk]=p0-dp);p[qk]=q0;end
    cb=PresetTimeComponentCallback([start],ComponentAffect(affect,(),(pk,qk)))
    set_callback!(vertices[bus],cb);active=Network(vertices,edges);set_jac_prototype!(active)
    active
end
wrapdiff(x)=mod(x+pi,2pi)-pi

function simulate_case(base,nw,rho,kp,ki,bus,delta)
    active=event_network(nw,bus,delta)
    state=Ref.trim_state(active,base,rho,kp,ki)
    dt=.01;times=collect(0.0:dt:61.0)
    prob=SciMLBase.ODEProblem(active,state,(0.0,61.0))
    sol=SciMLBase.solve(prob,OrdinaryDiffEqRosenbrock.Rodas5P();callback=get_callbacks(active),
        initializealg=SciMLBase.NoInit(),saveat=dt,abstol=1e-9,reltol=1e-9)
    SciMLBase.successful_retcode(sol.retcode)||error("OOS TDS failed bus=$bus MW=$delta retcode=$(sol.retcode)")
    angles=zeros(10,length(times));prev=zeros(10);theta=zeros(10);initial=zeros(10)
    sgf=zeros(10,length(times));pllf=similar(sgf);vdc=similar(sgf)
    for (j,t) in enumerate(times)
        ss=NetworkDynamics.NWState(sol,t)
        for b in 30:39
            k=b-29;ur=Float64(ss[VIndex(b,:busbar₊u_r)]);ui=Float64(ss[VIndex(b,:busbar₊u_i)])
            a=atan(ui,ur)
            if j==1;prev[k]=a;theta[k]=a;initial[k]=a
            else;theta[k]+=wrapdiff(a-prev[k]);prev[k]=a;end
            angles[k,j]=theta[k]-initial[k]
            prefix=b==39 ? "machine₊" : "ctrld_gen₊machine₊"
            sgf[k,j]=60*(Float64(ss[VIndex(b,Symbol(prefix*"ω"))])-1)
            pllf[k,j]=Float64(ss[VIndex(b,:gfl₊pll₊Δω_rad_s)])/(2pi)
            vdc[k,j]=Float64(ss[VIndex(b,:gfl₊v_dc_state)])
        end
    end
    tmetric=vcat(collect(-4.0:dt:-1.01),times.-1.0)
    phase=hcat(zeros(10,length(tmetric)-length(times)),angles)
    n10=Int(round(10/dt));steady=[(phase[k,end]-phase[k,end-n10])/(2pi*10) for k in 1:10]
    sig=(;times_s=tmetric,phase_rad=phase,phase_jump_rad=zeros(10),F_inf_Hz=steady,
        poles=ComplexF64[],alpha=NaN,A=zeros(1,1),B=zeros(1,1),C_bus=zeros(10,1),
        D_bus=zeros(10),condition_A=NaN,eigenvector_condition=NaN)
    m=only(FiniteWindow.window_metrics(sig;windows=(.5,),dt_s=dt,horizon_s=60,
        disturbance_MW=delta))
    (;bus,delta,retcode=string(sol.retcode),F_peak_Hz=m.F_peak_Hz,
      F_inf_Hz=m.F_inf_Hz,R_peak_Hz_s=m.R_peak_Hz_s,
      F_peak_bus=m.F_peak_bus,F_peak_time_s=m.F_peak_time_s,
      R_peak_bus=m.R_peak_bus,R_peak_time_s=m.R_peak_time_s,
      max_abs_SG_rotor_Hz=maximum(abs.(sgf)),max_abs_PLL_Hz=maximum(abs.(pllf)),
      vdc_min_pu=minimum(vdc),vdc_max_pu=maximum(vdc),
      pass=(m.F_peak_Hz<=.5 && m.R_peak_Hz_s<=.5),times,phase,sgf,pllf,vdc)
end

function main()
    c=TOML.parsefile(CANDIDATE);rho=Float64.(c["rho"]);kp=Float64.(c["Kp"]);ki=Float64.(c["Ki"])
    base=Ref.frozen_baseline();nw=Ref.build_architecture(base,rho,kp,ki)
    cases=((8,100.0),(29,100.0),(16,25.0),(16,50.0))
    results=NamedTuple[];traces=DataFrame[]
    for (bus,delta) in cases
        println("OOS_START bus=",bus," MW=",delta)
        r=simulate_case(base,nw,rho,kp,ki,bus,delta)
        push!(results,(;bus=r.bus,delta_MW=r.delta,retcode=r.retcode,
            F_peak_Hz=r.F_peak_Hz,F_inf_Hz=r.F_inf_Hz,R_peak_Hz_s=r.R_peak_Hz_s,
            F_peak_bus=r.F_peak_bus,F_peak_time_s=r.F_peak_time_s,
            R_peak_bus=r.R_peak_bus,R_peak_time_s=r.R_peak_time_s,
            max_abs_SG_rotor_Hz=r.max_abs_SG_rotor_Hz,max_abs_PLL_Hz=r.max_abs_PLL_Hz,
            vdc_min_pu=r.vdc_min_pu,vdc_max_pu=r.vdc_max_pu,pass=r.pass))
        df=DataFrame(time_s=r.times,event_time_s=r.times.-1)
        phase_start=size(r.phase,2)-length(r.times)+1
        for k in 1:10
            b=29+k;df[!,Symbol("phase_bus$(b)_rad")]=collect(view(r.phase,k,phase_start:size(r.phase,2)))
            df[!,Symbol("fSG_bus$(b)_Hz")]=collect(view(r.sgf,k,:))
            df[!,Symbol("fPLL_bus$(b)_Hz")]=collect(view(r.pllf,k,:))
            df[!,Symbol("Vdc_bus$(b)_pu")]=collect(view(r.vdc,k,:))
        end
        df[!,:event_bus]=fill(bus,nrow(df));df[!,:event_MW]=fill(delta,nrow(df));push!(traces,df)
        CSV.write(joinpath(CORE,"TABLE_Q2B_PD_out_of_sample.csv"),DataFrame(results))
        CSV.write(joinpath(CORE,"TABLE_Q2B_PD_out_of_sample_timeseries.csv"),vcat(traces...))
        println("OOS_DONE bus=",bus," MW=",delta," Fpeak=",r.F_peak_Hz,
            " Finf=",r.F_inf_Hz," Rpeak=",r.R_peak_Hz_s," pass=",r.pass)
    end
    CSV.write(joinpath(CORE,"TABLE_Q2B_PD_out_of_sample.csv"),DataFrame(results))
    CSV.write(joinpath(CORE,"TABLE_Q2B_PD_out_of_sample_timeseries.csv"),vcat(traces...))
end
main()
