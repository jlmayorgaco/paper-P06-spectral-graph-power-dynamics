include("ReducedDAE.jl")
using .ReducedDAE, LinearAlgebra, CSV, DataFrames, TOML, Sockets
const R=ReducedDAE
BLAS.set_num_threads(1);const CTX=R.N.design_context(R.ROOT)
const CFG=TOML.parsefile(joinpath(@__DIR__,"problem.toml"))
const CASES=[(;bus=Int(s["bus"]),delta=Float64(s["delta_MW"])) for s in CFG["scenarios"]]
const DEST=joinpath(R.OUT,get(ENV,"BND_SEARCH_SUBDIRECTORY","search"));mkpath(DEST)
const COUNT=Ref(0)
const GUARD=CFG["minimum_actuator_fraction_slack"]
const NBUNDLE=parse(Int,get(ENV,"BND_MODAL_BUNDLE_SIZE","1"))

function modal_bundle(rho,kp,ki)
    NBUNDLE==1 && return [R.N.simple_mode_sensitivities(CTX,rho,kp,ki)]
    s=R.N.spectrum(CTX,rho,kp,ki)
    vals=eigen(s.quotient'*s.model.Ared*s.quotient).values
    ids=sort([j for j in eachindex(vals) if imag(vals[j])>=-1e-8];by=j->real(vals[j]),rev=true)
    [R.N.simple_mode_sensitivities(CTX,rho,kp,ki;mode=j) for j in ids[1:min(NBUNDLE,length(ids))]]
end

function evaluate(p;grad=false,gradient_horizon=20.,gradient_dt=.05)
    rho=p[1:10];kp=clamp.(exp.(p[11:20]),CTX.kpmin,CTX.kpmax);ki=clamp.(exp.(p[21:30]),CTX.kimin,CTX.kimax)
    COUNT[]+=1;tag=lpad(COUNT[],4,'0')*(grad ? "_grad" : "_eval")
    open(joinpath(DEST,tag*"_parameters.toml"),"w") do io
        TOML.print(io,Dict("parameters"=>p,"modal_bundle_size"=>NBUNDLE,
            "problem_instance"=>CFG,"gradient_horizon_s"=>gradient_horizon,"gradient_step_s"=>gradient_dt))
    end
    c=Float64[];jac=Vector{Float64}[];rows=NamedTuple[]
    fm=CFG["frequency_limit_Hz"];rm=CFG["rocof_limit_Hz_s"]
    vlo=CFG["voltage_min_pu"];vhi=CFG["voltage_max_pu"];window=CFG["window_s"]
    decay=CFG["nominal_decay_margin_s_inv"]
    modes=modal_bundle(rho,kp,ki)
    condition=maximum(m.condition for m in modes)
    condition<1e6 || error("Ill-conditioned simple poles require a spectral-subspace audit")
    for modal in modes
        mg=real.(vcat(modal.rho,modal.Kp.*kp,modal.Ki.*ki))/decay
        push!(c,(real(modal.lambda)+decay)/decay);push!(jac,mg)
    end
    for case in CASES
        println("CASE ",tag," bus=",case.bus," delta=",case.delta);flush(stdout)
        m=R.model(CTX,rho,kp,ki;case.bus,case.delta,dc_convention=Symbol(get(CFG,"dc_convention","legacy")))
        if grad
            # Early transients guide the proposal; every accepted trial is independently
            # checked to 60 s. A late active constraint requires horizon enrichment.
            trj=R.tangent_simulate(m;horizon=gradient_horizon,dt=gradient_dt,monitor_buses=collect(1:39))
            trj.limiter_fraction<0 && error("smooth sensitivity invalid: limiter crossed")
            peaks=R.tangent_peaks(trj;window)
            append!(c,[peaks[1].peak/fm-1,peaks[2].peak/rm-1,(vlo-trj.vmin)/.1,(trj.vmax-vhi)/.1,(GUARD-trj.limiter_fraction)/.1])
            append!(jac,[peaks[1].gradient/fm,peaks[2].gradient/rm,-trj.gvmin/.1,trj.gvmax/.1,-trj.limiter_gradient/.1])
            push!(rows,(;case...,F=peaks[1].peak,R=peaks[2].peak,Ftime=peaks[1].time,Rtime=peaks[2].time,
                Vmin=trj.vmin,Vmax=trj.vmax,limiter_fraction=trj.limiter_fraction,
                limiter_bus=trj.limiter_bus,limiter_kind=trj.limiter_kind,limiter_time=trj.limiter_time,limiter_side=trj.limiter_side))
        else
            sol=R.simulate(m;horizon=CFG["horizon_s"],dt=.01,tol=2e-8,wall_limit=70.,rotating_frame=true)
            met=R.metrics(m,sol;dt=.01,window,monitor_buses=collect(1:39))
            push!(rows,(;case...,met...))
            if sol.ok
                append!(c,[met.Fpeak_Hz/fm-1,met.Rpeak_Hz_s/rm-1,(vlo-met.Vmin)/.1,(met.Vmax-vhi)/.1,(GUARD-met.limiter_fraction)/.1])
            else
                append!(c,fill(10.,5))
            end
        end
        CSV.write(joinpath(DEST,tag*".csv"),DataFrame(rows))
    end
    if grad
        A=reduce(vcat,transpose.(jac))
        CSV.write(joinpath(DEST,tag*"_jacobian.csv"),DataFrame(A,:auto))
        open(joinpath(DEST,tag*"_constraints.toml"),"w") do io
            TOML.print(io,Dict("constraints"=>c,"parameters"=>p,"modal_condition"=>condition,"modal_bundle_size"=>NBUNDLE,
                "gradient_horizon_s"=>gradient_horizon,"gradient_step_s"=>gradient_dt))
        end
        return vcat(c,vec(A),condition)
    end
    c
end

if abspath(PROGRAM_FILE)==(@__FILE__)
server=listen(ip"127.0.0.1",0);println("READY ",getsockname(server)[2]);flush(stdout)
sock=accept(server)
while isopen(sock)
    line=readline(sock);line=="QUIT" && break
    try
        cmd,vals=split(line,';');p=parse.(Float64,split(vals,','))
        cmd in ("EVAL","GRAD","GRADFULL") || error("unknown oracle command")
        out=evaluate(p;grad=cmd!="EVAL",gradient_horizon=cmd=="GRADFULL" ? CFG["horizon_s"] : 20.)
        println(sock,join(out,','));flush(sock)
    catch ex
        println(sock,"ERR ",replace(sprint(showerror,ex),'\n'=>' '));flush(sock)
    end
end
close(sock);close(server)
end
