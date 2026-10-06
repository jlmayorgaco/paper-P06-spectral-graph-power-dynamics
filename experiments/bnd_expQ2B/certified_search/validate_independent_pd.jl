using CSV, DataFrames, LinearAlgebra, TOML, SHA, NetworkDynamics, PowerDynamics,
    OrdinaryDiffEqRosenbrock, SciMLBase
const ROOT=normpath(joinpath(@__DIR__,"..","..",".."))
const OUT=joinpath(get(ENV,"EXPQ2B_CERT_OUT",joinpath(ROOT,"reports","experiment_Q2B","CERTIFIED_SEARCH")),"PD")
mkpath(OUT)
const CAND=joinpath(dirname(OUT),"Z_LOCAL_SECURE_FINAL.toml")
const HASH=strip(read(CAND*".sha256",String))
bytes2hex(sha256(read(CAND)))==HASH || error("Candidate hash mismatch")
include(joinpath(ROOT,"src","bnd_model_expN","PDReferenceN.jl"))
const P=PDReferenceN
const BUILT=Ref{Any}(nothing)
const SOLUTIONS=Dict{String,Any}()
include("RiccatiCertificate.jl")
include("ConsistentEventObservation.jl")
BLAS.set_num_threads(1)
function pdkey(sym)
    m=match(r"VIndex\((\d+), :(.*)\)",string(sym));m===nothing&&error("Unknown PD coordinate $sym")
    bus=parse(Int,m[1]);name=m[2];tail=last(split(name,"₊"))
    occursin("₊gov₊",name)&&return (bus,"SG","gov_"*tail)
    occursin("₊avr₊",name)&&return (bus,"SG","avr_"*tail)
    if occursin("machine₊",name)
        names=Dict("ψ″_q"=>"psi2q","ψ″_d"=>"psi2d","E′_d"=>"Epd","E′_q"=>"Epq","ω"=>"omega","δ"=>"delta")
        return (bus,"SG","machine_"*names[tail])
    end
    names=Dict("γ_q"=>"gamma_q","γ_d"=>"gamma_d","θ"=>"theta","Δω_rad_s"=>"delta_omega_rad_s","Δω_i_rad_s"=>"delta_omega_i_rad_s")
    (bus,"GFL",get(names,tail,tail))
end
function independent_matrix(sys,state)
    mass=sys.M isa UniformScaling ? ones(size(sys.A,1)) : diag(sys.M)
    dyn=findall(mass.==1);alg=findall(mass.==0);J=Matrix(sys.A)
    Ar=J[dyn,dyn]-J[dyn,alg]*(J[alg,alg]\J[alg,dyn])
    keys=pdkey.(sys.sym[dyn]);lookup=Dict(k=>i for (i,k) in enumerate(keys))
    inventory=CSV.read(joinpath(dirname(OUT),"FROZEN_STATE_MAP.csv"),DataFrame)
    wanted=[(Int(r.bus),String(r.kind),String(r.state_name)) for r in eachrow(inventory)]
    perm=[lookup[k] for k in wanted];Ar=Ar[perm,perm]
    g=zeros(length(perm))
    for (i,(bus,kind,name)) in enumerate(wanted)
        if name in ("machine_delta","theta");g[i]=1.;end
        name=="i_f_i" && (g[i]=Float64(state[VIndex(bus,:gfl₊filter₊i_f_r)]))
        name=="i_f_r" && (g[i]=-Float64(state[VIndex(bus,:gfl₊filter₊i_f_i)]))
    end
    Q=nullspace(reshape(g/norm(g),1,:));A=Q'*Ar*Q
    CSV.write(joinpath(OUT,"INDEPENDENT_PHYSICAL_A.csv"),DataFrame(A,:auto))
    A
end
function event_network(nw,bus,delta;duration=Inf)
    vertices,edges=P.PD39.PD39Model.copy_network_components(nw)
    defaults=get_defaults_dict(vertices[bus])
    pk=only(filter(s->occursin("Pset",string(s)),collect(keys(defaults))))
    qk=only(filter(s->occursin("Qset",string(s)),collect(keys(defaults))))
    p0=defaults[pk];q0=defaults[qk]
    affect=(u,p,ctx)->begin;p[pk]=p0-(ctx.t<1+duration ? delta/100 : 0.);p[qk]=q0;end
    set_callback!(vertices[bus],PresetTimeComponentCallback(isfinite(duration) ? [1.,1+duration] : [1.],ComponentAffect(affect,(),(pk,qk))))
    active=Network(vertices,edges);set_jac_prototype!(active);active
end
function run()
    c=TOML.parsefile(CAND);rho=c["rho"];kp=c["Kp"];ki=c["Ki"]
    println("INDEPENDENT_PD_BUILD hash=$HASH");flush(stdout)
    if BUILT[]===nothing
        base=P.frozen_baseline();nw=P.build_architecture(base,rho,kp,ki)
        state=P.trim_state(nw,base,rho,kp,ki);BUILT[]=(base,nw,state)
    else
        base,nw,state=BUILT[]
    end
    trim=P.residual_audit(nw,state);pq=P.direct_power_audit(state,base,rho)
    CSV.write(joinpath(OUT,"COMPONENT_PQ.csv"),pq)
    linear=linearize_network(state);lp=jacobian_eigenvals(linear)
    ig=argmin(abs.(lp));physical=lp[setdiff(eachindex(lp),[ig])]
    poles=CSV.read(joinpath(dirname(OUT),"FROZEN_ANALYTICAL_POLES.csv"),DataFrame)
    analytic=complex.(poles.real,poles.imag);length(physical)==length(analytic)||error("Pole count mismatch")
    # Symmetric nearest spectral distance; all poles and multiplicities are
    # also retained for independent matching/review.
    E=abs.(physical.-transpose(analytic));pole_error=max(maximum(minimum(E;dims=1)),maximum(minimum(E;dims=2)))
    spectrum=Dict("candidate_sha256"=>HASH,"alpha_PD"=>maximum(real.(physical)),
        "alpha_analytic"=>maximum(real.(analytic)),"alpha_error"=>abs(maximum(real.(physical))-maximum(real.(analytic))),
        "complete_spectrum_Hausdorff_error"=>pole_error,"physical_pole_count"=>length(physical),
        "gauge_abs"=>abs(lp[ig]),"trim_residual"=>trim.maximum,"max_P_error_pu"=>maximum(pq.max_P_error_pu),
        "max_Q_error_pu"=>maximum(pq.max_Q_error_pu),"trim_bounds_pass"=>all(pq.bounds_status.=="PASS"),
        "margin_pass"=>maximum(real.(physical))<=-.05)
    open(joinpath(OUT,"SPECTRUM.toml"),"w") do io;TOML.print(io,spectrum);end
    CSV.write(joinpath(OUT,"ALL_PD_POLES.csv"),DataFrame(real=real.(lp),imag=imag.(lp)))
    println("PD_SPECTRUM ",spectrum);flush(stdout)
    spectrum["margin_pass"]&&pole_error<1e-5 || error("Independent spectrum gate failed")
    Apd=independent_matrix(linear,state)
    Aan=Matrix(CSV.read(joinpath(dirname(OUT),"FROZEN_PHYSICAL_A.csv"),DataFrame))
    matrix_relative_error=norm(Apd-Aan)/norm(Aan)
    println("PD_MATRIX_RELATIVE_ERROR ",matrix_relative_error);flush(stdout)
    br=riccati_certificate(Apd;beta=Float64(c["beta_requirement"]),outdir=OUT)
    br["strict_negative_LMI"] || error("PD robustness certificate failed")
    spectrum["matrix_relative_error"]=matrix_relative_error;spectrum["beta_certified"]=true
    open(joinpath(OUT,"SPECTRUM.toml"),"w") do io;TOML.print(io,spectrum);end
    loadbuses=[b for b in 1:39 if try;Float64(state[VIndex(b,:ZIPLoad₊P)]);true;catch;false;end]
    for b in (31,39)
        abs(Float64(state[VIndex(b,:ZIPLoad₊P)])-Float64(base.state[VIndex(b,:ZIPLoad₊P)]))<1e-10||error("ZIP load changed")
        abs(Float64(state[VIndex(b,:ZIPLoad₊Q)])-Float64(base.state[VIndex(b,:ZIPLoad₊Q)]))<1e-10||error("ZIP Q changed")
    end
    cases=[(name="bus16_100MW",bus=16,d=100.,duration=Inf),
        (name="bus16_1MW",bus=16,d=1.,duration=Inf),
        (name="bus16_25MW",bus=16,d=25.,duration=Inf),
        (name="bus16_50MW",bus=16,d=50.,duration=Inf),
        (name="bus8_100MW",bus=8,d=100.,duration=Inf),
        (name="bus29_100MW",bus=29,d=100.,duration=Inf),
        (name="bus16_pulse",bus=16,d=100.,duration=.1)]
    rows=NamedTuple[]
    for case in cases
        println("PD_EVENT_START $(case.name)");flush(stdout)
        active=event_network(nw,case.bus,case.d;duration=case.duration);dt=.005;started=time();last=Ref(time())
        watch=DiscreteCallback((u,t,i)->time()-last[]>30,i->begin
            println("PD_PROGRESS $(case.name) t=$(i.t) wall=$(time()-started)");flush(stdout);last[]=time()
                time()-started>300 ? SciMLBase.terminate!(i) : SciMLBase.u_modified!(i,false)
        end;save_positions=(false,false))
        sol=get!(SOLUTIONS,case.name) do
            solve(ODEProblem(active,state,(0.,61.)),Rodas5P();callback=CallbackSet(get_callbacks(active),watch),
                initializealg=SciMLBase.NoInit(),saveat=dt,abstol=1e-10,reltol=1e-10,maxiters=2000000)
        end
        ok=SciMLBase.successful_retcode(sol.retcode)&&sol.t[end]>60.999
        if !ok
            open(joinpath(OUT,"$(case.name)_INCOMPLETE.toml"),"w") do io
                TOML.print(io,Dict("retcode"=>string(sol.retcode),"last_time"=>sol.t[end],"physical_instability_proven"=>false))
            end
            case.name=="bus16_100MW" && error("Primary nonlinear validation incomplete")
            continue
        end
        ts=collect(0.:dt:61.);nt=length(ts);theta=zeros(10,nt);voltage=similar(theta);sgf=fill(NaN,10,nt)
        pllf=zeros(10,nt);sgP=zeros(10,nt);gfP=zeros(10,nt);dc=zeros(10,nt);pdc=zeros(10,nt)
        sgQ=zeros(10,nt);gfQ=zeros(10,nt);mech=zeros(10,nt)
        loads=zeros(nt);losses=zeros(nt);balance=zeros(nt);rawlast=zeros(10);initial=zeros(10);angle=zeros(10)
        for (j,t) in enumerate(ts)
            te=abs(t-1)<1e-12||(isfinite(case.duration)&&abs(t-1-case.duration)<1e-12) ? t+1e-9 : t
            ss=te==t ? NWState(sol,te) : consistent_event_observation(sol,te)
            for b in 30:39
                k=b-29;ur=Float64(ss[VIndex(b,:busbar₊u_r)]);ui=Float64(ss[VIndex(b,:busbar₊u_i)]);raw=atan(ui,ur)
                if j==1;rawlast[k]=raw;angle[k]=raw;initial[k]=raw;else;angle[k]+=mod(raw-rawlast[k]+pi,2pi)-pi;rawlast[k]=raw;end
                theta[k,j]=angle[k]-initial[k];voltage[k,j]=hypot(ur,ui)
                ir=Float64(ss[VIndex(b,:gfl₊filter₊i_f_r)]);ii=Float64(ss[VIndex(b,:gfl₊filter₊i_f_i)])
                gfP[k,j]=100rho[k]*(ur*ir+ui*ii);gfQ[k,j]=100rho[k]*(ui*ir-ur*ii)
                pllf[k,j]=Float64(ss[VIndex(b,:gfl₊pll₊Δω_rad_s)])/(2pi)
                dc[k,j]=Float64(ss[VIndex(b,:gfl₊v_dc_state)]);pdc[k,j]=Float64(ss[VIndex(b,:gfl₊P_dc)])
                if rho[k]<1
                    prefix=b==39 ? "machine₊" : "ctrld_gen₊machine₊"
                    sgf[k,j]=60*(Float64(ss[VIndex(b,Symbol(prefix*"ω"))])-1)
                    sgP[k,j]=Float64(ss[VIndex(b,Symbol(prefix*"Sn"))])*Float64(ss[VIndex(b,Symbol(prefix*"P"))])
                    sgQ[k,j]=Float64(ss[VIndex(b,Symbol(prefix*"Sn"))])*Float64(ss[VIndex(b,Symbol(prefix*"Q"))])
                    mech[k,j]=Float64(ss[VIndex(b,Symbol(prefix*"Sn"))])*Float64(ss[VIndex(b,Symbol(prefix*"τ_m"))])*Float64(ss[VIndex(b,Symbol(prefix*"ω"))])
                end
            end
            loads[j]=-100sum(Float64(ss[VIndex(b,:ZIPLoad₊P)]) for b in loadbuses)
            losses[j]=sum(Float64(ss[VIndex(b,:busbar₊P_MW)]) for b in 1:39)
            balance[j]=sum(sgP[:,j])+sum(gfP[:,j])-loads[j]-losses[j]
        end
        lag=Int(round(.5/dt));f=zeros(size(theta));r=similar(f)
        for j in 1:nt
            p1=j>lag ? theta[:,j-lag] : zeros(10);p2=j>2lag ? theta[:,j-2lag] : zeros(10)
            f[:,j]=(theta[:,j]-p1)/pi;r[:,j]=(theta[:,j]-2p1+p2)/(pi/2)
        end
        steady=(theta[:,end]-theta[:,end-Int(round(10/dt))])/(20pi)
        peak=maximum(abs,f);roc=maximum(abs,r)
        pdc_error=maximum(abs.(pdc.-pdc[:,1]))
        row=(;event=case.name,bus=case.bus,MW=case.d,duration=case.duration,Fpeak=peak,Rpeak=roc,
            Fsteady=maximum(abs,steady),frequency_pass=peak<=.5,rocof_pass=roc<=.5,
            voltage_min=minimum(voltage),voltage_max=maximum(voltage),vdc_min=minimum(dc),vdc_max=maximum(dc),
            Pdc_change=pdc_error,power_balance_error_MW=maximum(abs,balance),elapsed_s=time()-started,retcode=string(sol.retcode))
        push!(rows,row);CSV.write(joinpath(OUT,"EVENTS.csv"),DataFrame(rows))
        df=DataFrame(time_s=ts,event_time_s=ts.-1.,load_MW=loads,network_losses_MW=losses,balance_error_MW=balance)
        for k in 1:10
            b=k+29
            for (name,arr) in (("grid_Hz",f),("rocof_Hz_s",r),("SG_Hz",sgf),("PLL_Hz",pllf),("V_pu",voltage),("SG_MW",sgP),("GFL_MW",gfP),("SG_Mvar",sgQ),("GFL_Mvar",gfQ),("SG_mechanical_MW",mech),("Vdc_pu",dc),("Pdc_pu",pdc))
                df[!,Symbol("$(name)_bus$b")]=arr[k,:]
            end
        end
        CSV.write(joinpath(OUT,"$(case.name)_SENSORS.csv"),df)
        println("PD_EVENT_DONE ",row);flush(stdout)
        if case.name=="bus16_100MW" && !(peak<=.5&&roc<=.5)
            error("Frozen candidate failed the primary nonlinear event; no retuning")
        end
    end
end
try;run();catch err;showerror(stdout,err,catch_backtrace());println();end
for line in eachline(stdin)
    strip(line)=="EXIT" && break
    try;include_string(Main,line,"validation_console");catch err;showerror(stdout,err,catch_backtrace());println();end
    println("COMMAND_DONE");flush(stdout)
end
