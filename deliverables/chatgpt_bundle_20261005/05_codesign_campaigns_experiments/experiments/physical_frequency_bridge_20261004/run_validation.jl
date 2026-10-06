using CSV, DataFrames, LinearAlgebra, TOML, SHA, NetworkDynamics, PowerDynamics,
    OrdinaryDiffEqRosenbrock, SciMLBase, Dates
const OUT=@__DIR__
const ROOT=normpath(joinpath(OUT,"..",".."))
const CAND=joinpath(ROOT,"reports","experiment_Q2B","CERTIFIED_SEARCH","Z_LOCAL_SECURE_FINAL.toml")
const EXPECTED="66513a8d3a1b4cc37c9d6f2b1a5824371fed2d0c5621d2ccdb7540490fe59124"
bytes2hex(sha256(read(CAND)))==EXPECTED || error("Candidate changed")
bytes2hex(sha256(read(joinpath(OUT,"PREREGISTRATION.json"))))==strip(read(joinpath(OUT,"PREREGISTRATION.sha256"),String)) || error("Protocol changed")
include(joinpath(ROOT,"src","bnd_model_expN","PDReferenceN.jl"))
const P=PDReferenceN
include("PhysicalBridge.jl")
include(joinpath(ROOT,"experiments","bnd_expQ2B","certified_search","ConsistentEventObservation.jl"))
BLAS.set_num_threads(1)
write_toml(name,x)=open(io->TOML.print(io,x),joinpath(OUT,name),"w")
val(s,b,n)=Float64(s[VIndex(b,Symbol(n))])
function event_network(nw,bus,delta)
    vertices,edges=P.PD39.PD39Model.copy_network_components(nw)
    defaults=get_defaults_dict(vertices[bus])
    pk=only(filter(s->occursin("Pset",string(s)),collect(keys(defaults))))
    qk=only(filter(s->occursin("Qset",string(s)),collect(keys(defaults))))
    p0=defaults[pk];q0=defaults[qk]
    affect=(u,p,ctx)->begin;p[pk]=p0-delta/100;p[qk]=q0;end
    set_callback!(vertices[bus],PresetTimeComponentCallback([1.],ComponentAffect(affect,(),(pk,qk))))
    active=Network(vertices,edges);set_jac_prototype!(active);active
end
function matrix_audit(nwo,so,nwp,sp)
    lo=linearize_network(so);lp=linearize_network(sp)
    ko=string.(lo.sym);kp=string.(lp.sym)
    perm=[only(findall(==(k),kp)) for k in ko]
    jo=Matrix(lo.A);jp=Matrix(lp.A)[perm,perm]
    n=length(ko);sgn=[occursin("gfl₊v_dc_state",k) ? -1. : 1. for k in ko]
    S=Diagonal(sgn);rel=norm(jp-S*jo*S)/norm(jo)
    mo=lo.M isa UniformScaling ? Matrix{Float64}(I,n,n) : Matrix(lo.M)
    mp=lp.M isa UniformScaling ? Matrix{Float64}(I,n,n) : Matrix(lp.M)[perm,perm]
    massrel=norm(mp-S*mo*S)/max(1,norm(mo))
    eo=jacobian_eigenvals(lo);ep=jacobian_eigenvals(lp)
    all(isfinite,eo)&&all(isfinite,ep)||error("Unexpected nonfinite reduced poles")
    ig=argmin(abs.(ep));alpha=maximum(real.(ep[setdiff(eachindex(ep),[ig])]))
    eo_sorted=sort(eo;by=z->(real(z),imag(z)));ep_sorted=sort(ep;by=z->(real(z),imag(z)))
    E=abs.(eo.-transpose(ep));haus=max(maximum(minimum(E;dims=1)),maximum(minimum(E;dims=2)))
    CSV.write(joinpath(OUT,"LEGACY_POLES.csv"),DataFrame(real=real.(eo),imag=imag.(eo)))
    CSV.write(joinpath(OUT,"PHYSICAL_POLES.csv"),DataFrame(real=real.(ep),imag=imag.(ep)))
    CSV.write(joinpath(OUT,"JACOBIAN_STATE_MAP.csv"),DataFrame(state=ko,physical_index=perm,similarity_sign=sgn))
    reso=P.residual_audit(nwo,so).maximum;resp=P.residual_audit(nwp,sp).maximum
    pq=P.direct_power_audit(sp,BASE[],RHO[]);CSV.write(joinpath(OUT,"PHYSICAL_TRIM_PQ.csv"),pq)
    x=Dict("legacy_residual"=>reso,"physical_residual"=>resp,"jacobian_similarity_relative"=>rel,
      "mass_similarity_relative"=>massrel,"pole_hausdorff_absolute"=>haus,"finite_pole_count"=>length(ep),
      "gauge_abs"=>abs(ep[ig]),"physical_alpha"=>alpha,"max_P_error_pu"=>maximum(pq.max_P_error_pu),
      "max_Q_error_pu"=>maximum(pq.max_Q_error_pu),"candidate_sha256"=>EXPECTED)
    write_toml("PARITY.toml",x);println("PARITY ",x);flush(stdout)
    reso<1e-10&&resp<1e-10&&rel<1e-10&&massrel<1e-12&&haus<1e-5&&maximum(pq.max_P_error_pu)<1e-10&&maximum(pq.max_Q_error_pu)<1e-10 || error("BLOCKED_BASELINE_PARITY")
end
const BASE=Ref{Any}();const RHO=Ref{Any}()
function run_case(model,nw,state,rho,case)
    name=model*"_"*case.name
    println("CASE_START ",name);flush(stdout)
    active=event_network(nw,case.bus,case.d);dt=.005;started=time();last=Ref(time())
    watch=DiscreteCallback((u,t,i)->time()-last[]>30,i->begin
      println("PROGRESS $name t=$(i.t) wall=$(time()-started)");flush(stdout);last[]=time()
      time()-started>600 ? SciMLBase.terminate!(i) : SciMLBase.u_modified!(i,false)
    end;save_positions=(false,false))
    sol=solve(ODEProblem(active,state,(0.,61.)),Rodas5P();callback=CallbackSet(get_callbacks(active),watch),
      initializealg=SciMLBase.NoInit(),saveat=dt,abstol=1e-10,reltol=1e-10,maxiters=2000000)
    if !(SciMLBase.successful_retcode(sol.retcode)&&sol.t[end]>60.999)
      write_toml(name*"_INCOMPLETE.toml",Dict("retcode"=>string(sol.retcode),"last_time"=>sol.t[end],"physical_instability_proven"=>false))
      return (;model,event=case.name,complete=false,Fpeak=NaN,Rpeak=NaN,frequency_pass=false,rocof_pass=false,voltage_min=NaN,voltage_max=NaN,vdc_min=NaN,vdc_max=NaN,energy_RHS_error_MW=NaN,power_balance_error_MW=NaN,governor_upper_margin=NaN,governor_lower_margin=NaN,elapsed_s=time()-started)
    end
    ts=collect(0.:dt:61.);nt=length(ts);theta=zeros(10,nt);voltage=similar(theta);sgf=fill(NaN,10,nt)
    pllf=zeros(10,nt);sgP=zeros(10,nt);gfP=zeros(10,nt);dc=zeros(10,nt);pdc=zeros(10,nt)
    sgQ=zeros(10,nt);gfQ=zeros(10,nt);mech=zeros(10,nt);xg1=fill(NaN,10,nt);xg2=copy(xg1)
    govlow=copy(xg1);govhigh=copy(xg1)
    loads=zeros(nt);losses=zeros(nt);balance=zeros(nt);K=zeros(nt);Ef=zeros(nt);Edc=zeros(nt)
    lossf=zeros(nt);dotK=zeros(nt);dotEf=zeros(nt);dotEdc=zeros(nt);input=zeros(nt);eres=zeros(nt)
    rawlast=zeros(10);initial=zeros(10);angle=zeros(10)
    loadbuses=[b for b in 1:39 if try;val(state,b,"ZIPLoad₊P");true;catch;false;end]
    syms=NetworkDynamics.SII.variable_symbols(active);lookup=Dict(string(s)=>j for (j,s) in enumerate(syms))
    mass=active.mass_matrix;md=mass isa UniformScaling ? ones(length(syms)) : diag(mass)
    derivative(rhs,b,n)=rhs[lookup[string(VIndex(b,Symbol(n)))]]/md[lookup[string(VIndex(b,Symbol(n)))]]
    rhs=similar(uflat(state));sign_dc=model=="physical" ? 1. : -1.
    for (j,t) in enumerate(ts)
      te=abs(t-1)<1e-12 ? t+1e-9 : t
      ss=te==t ? NWState(sol,te) : consistent_event_observation(sol,te)
      active(rhs,uflat(ss),pflat(ss),te)
      for b in 30:39
        k=b-29;ur=val(ss,b,"busbar₊u_r");ui=val(ss,b,"busbar₊u_i");raw=atan(ui,ur)
        if j==1;rawlast[k]=raw;angle[k]=raw;initial[k]=raw;else;angle[k]+=mod(raw-rawlast[k]+pi,2pi)-pi;rawlast[k]=raw;end
        theta[k,j]=angle[k]-initial[k];voltage[k,j]=hypot(ur,ui)
        ir=val(ss,b,"gfl₊filter₊i_f_r");ii=val(ss,b,"gfl₊filter₊i_f_i");scale=100rho[k]
        gfP[k,j]=scale*(ur*ir+ui*ii);gfQ[k,j]=scale*(ui*ir-ur*ii)
        pllf[k,j]=val(ss,b,"gfl₊pll₊Δω_rad_s")/(2pi)
        dc[k,j]=val(ss,b,"gfl₊v_dc_state");pdc[k,j]=val(ss,b,"gfl₊P_dc")
        rf=val(ss,b,"gfl₊filter₊Rf");xf=val(ss,b,"gfl₊filter₊Xf");cdc=val(ss,b,"gfl₊C_dc")
        Ef[j]+=.5scale*xf/(2pi*60)*(ir^2+ii^2);Edc[j]+=.5scale*cdc*dc[k,j]^2
        lossf[j]+=scale*rf*(ir^2+ii^2)
        dotEf[j]+=scale*xf/(2pi*60)*(ir*derivative(rhs,b,"gfl₊filter₊i_f_r")+ii*derivative(rhs,b,"gfl₊filter₊i_f_i"))
        dotEdc[j]+=scale*cdc*dc[k,j]*derivative(rhs,b,"gfl₊v_dc_state")
        if rho[k]<1
          pre=b==39 ? "machine₊" : "ctrld_gen₊machine₊"
          sn=val(ss,b,pre*"Sn");h=val(ss,b,pre*"H");w=val(ss,b,pre*"ω")
          # This frozen campaign has stator resistance and explicit machine damping zero.
          val(ss,b,pre*"R_s")==0&&val(ss,b,pre*"D")==0 || error("Unregistered SG dissipation")
          sgf[k,j]=60*(w-1);sgP[k,j]=sn*val(ss,b,pre*"P");sgQ[k,j]=sn*val(ss,b,pre*"Q")
          mech[k,j]=sn*val(ss,b,pre*"τ_m")*w;K[j]+=h*sn*w^2
          dotK[j]+=2h*sn*w*derivative(rhs,b,pre*"ω")
          if b!=39
            xg1[k,j]=val(ss,b,"ctrld_gen₊gov₊xg1");xg2[k,j]=val(ss,b,"ctrld_gen₊gov₊xg2")
            govlow[k,j]=xg1[k,j]-val(ss,b,"ctrld_gen₊gov₊V_min")
            govhigh[k,j]=val(ss,b,"ctrld_gen₊gov₊V_max")-xg1[k,j]
          end
        end
      end
      loads[j]=-100sum(val(ss,b,"ZIPLoad₊P") for b in loadbuses)
      losses[j]=sum(val(ss,b,"busbar₊P_MW") for b in 1:39)
      balance[j]=sum(sgP[:,j])+sum(gfP[:,j])-loads[j]-losses[j]
      input[j]=sum(mech[:,j])+sum(100 .*rho.*pdc[:,j])-loads[j]-losses[j]-lossf[j]
      eres[j]=dotK[j]+dotEf[j]+sign_dc*dotEdc[j]-input[j]
    end
    lag=Int(round(.5/dt));f=zeros(size(theta));r=similar(f)
    for j in 1:nt
      p1=j>lag ? theta[:,j-lag] : zeros(10);p2=j>2lag ? theta[:,j-2lag] : zeros(10)
      f[:,j]=(theta[:,j]-p1)/pi;r[:,j]=(theta[:,j]-2p1+p2)/(pi/2)
    end
    df=DataFrame(time_s=ts,event_time_s=ts.-1.,load_MW=loads,network_losses_MW=losses,
      filter_losses_MW=lossf,balance_error_MW=balance,rotor_energy_MJ=K,filter_energy_MJ=Ef,
      dc_energy_MJ=Edc,aggregate_energy_MJ=K.+Ef.+sign_dc.*Edc,net_input_MW=input,
      dotK_MW=dotK,dotEf_MW=dotEf,dotEdc_MW=dotEdc,energy_RHS_error_MW=eres)
    for k in 1:10
      b=k+29
      for (key,arr) in (("grid_Hz",f),("rocof_Hz_s",r),("SG_Hz",sgf),("PLL_Hz",pllf),("V_pu",voltage),("SG_MW",sgP),("GFL_MW",gfP),("SG_Mvar",sgQ),("GFL_Mvar",gfQ),("SG_mechanical_MW",mech),("Vdc_pu",dc),("Pdc_pu",pdc),("gov_xg1",xg1),("gov_xg2",xg2),("gov_upper_margin",govhigh),("gov_lower_margin",govlow))
        df[!,Symbol("$(key)_bus$b")]=arr[k,:]
      end
    end
    CSV.write(joinpath(OUT,name*"_SENSORS.csv"),df)
    peak=maximum(abs,f);roc=maximum(abs,r)
    row=(;model,event=case.name,complete=true,Fpeak=peak,Rpeak=roc,frequency_pass=peak<=.5,rocof_pass=roc<=.5,
      voltage_min=minimum(voltage),voltage_max=maximum(voltage),vdc_min=minimum(dc),vdc_max=maximum(dc),
      energy_RHS_error_MW=maximum(abs,eres),power_balance_error_MW=maximum(abs,balance),
      governor_upper_margin=minimum(filter(isfinite,vec(govhigh))),governor_lower_margin=minimum(filter(isfinite,vec(govlow))),elapsed_s=time()-started)
    println("CASE_DONE ",row);flush(stdout);row
end
function main()
    c=TOML.parsefile(CAND);rho=Float64.(c["rho"]);kp=Float64.(c["Kp"]);ki=Float64.(c["Ki"])
    write_toml("RUN_ENVIRONMENT.toml",Dict("started_UTC"=>string(now(UTC)),"julia_version"=>string(VERSION),"PowerDynamics_version"=>string(pkgversion(PowerDynamics)),"NetworkDynamics_version"=>string(pkgversion(NetworkDynamics)),"runner_sha256"=>bytes2hex(sha256(read(@__FILE__))),"model_sha256"=>bytes2hex(sha256(read(joinpath(OUT,"PhysicalBridge.jl"))))))
    println("BUILD_BASELINE");flush(stdout);base=P.frozen_baseline();BASE[]=base;RHO[]=rho
    println("BUILD_LEGACY");flush(stdout);old=P.build_architecture(base,rho,kp,ki);so=P.trim_state(old,base,rho,kp,ki)
    println("BUILD_PHYSICAL");flush(stdout);new=PhysicalBridge.build_architecture(base,rho,kp,ki);sp=P.trim_state(new,base,rho,kp,ki)
    matrix_audit(old,so,new,sp)
    cases=[(name="bus16_1MW",bus=16,d=1.),(name="bus16_100MW",bus=16,d=100.),(name="bus8_100MW",bus=8,d=100.),(name="bus29_100MW",bus=29,d=100.)]
    rows=NamedTuple[]
    push!(rows,run_case("legacy",old,so,rho,cases[1]));CSV.write(joinpath(OUT,"EVENT_RESULTS.csv"),DataFrame(rows))
    for case in cases
      push!(rows,run_case("physical",new,sp,rho,case));CSV.write(joinpath(OUT,"EVENT_RESULTS.csv"),DataFrame(rows))
    end
    write_toml("EXECUTION_COMPLETE.toml",Dict("completed_UTC"=>string(now(UTC)),"new_runs"=>length(rows),"all_complete"=>all(r.complete for r in rows),"no_tuning"=>true))
end
main()
