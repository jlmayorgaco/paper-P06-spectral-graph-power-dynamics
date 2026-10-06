using LinearAlgebra, CSV, DataFrames, TOML, SHA, Statistics
BLAS.set_num_threads(1)
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_design_k","BNDDesignK.jl"))
using .BNDDesignK
const OUT=joinpath(ROOT,"reports","experiment_K")
const TABLES=joinpath(OUT,"tables")
function frozen(name)
    path=joinpath(OUT,name)
    strip(read(path*".sha256",String))==bytes2hex(sha256(read(path))) || error("freeze mismatch")
    TOML.parsefile(path)
end
ctx=design_context(ROOT)
candidates=Dict("nominal"=>frozen("Z_K_NOMINAL_FINAL.toml"),
    "robust_ExpG_beta"=>frozen("Z_K_ROBUST_beta_1p6991206999182038em6.toml"))
traces=CSV.read(joinpath(TABLES,"TABLE_K16_tds_trajectories.csv"),DataFrame)
rows=NamedTuple[];out=NamedTuple[]
for (label,c) in candidates
    ev=evaluate(ctx,Float64.(c["epsilon"]),Float64.(c["Kp"]),Float64.(c["Ki"]);vectors=true)
    Cf=zeros(ev.model.n_dynamic);wt=0.0
    for (i,b) in enumerate(ctx.buses)
        ev.epsilon[i]>0 || continue
        row=only(eachrow(ev.model.state_map[(ev.model.state_map.bus.==b).&
            (ev.model.state_map.kind.=="SG"),:]))
        h=ctx.net.sg[b].op.parameters.inertia*ctx.net.sg[b].op.parameters.rating_mva
        w=ev.epsilon[i]*h;wt+=w
        Cf[Int(row.last)-1]=60w
    end
    Cf./=wt
    Cq=transpose(ev.Q)*Cf
    for bus in unique(Int.(traces.event_bus[traces.candidate.==label]))
        lr=only(eachrow(ctx.net.load_audit[ctx.net.load_audit.bus.==bus,:]))
        S=complex(lr.initialized_load_MW,lr.initialized_load_Mvar)/100
        Ii=conj(S/ctx.net.voltage[bus])
        iv=zeros(size(ev.model.Gy,1));iv[2bus-1]=real(Ii);iv[2bus]=imag(Ii)
        Bd=transpose(ev.Q)*(ev.model.B*(ev.model.Gy\iv))
        q=[dot(Cq,ev.right[:,j])*dot(ev.left[:,j],Bd)/
            dot(ev.left[:,j],ev.right[:,j]) for j in eachindex(ev.lambda)]
        step(t)=t<=0 ? 0.0 : real(sum(abs(ev.lambda[j])<1e-10 ? q[j]*t :
            q[j]*expm1(ev.lambda[j]*t)/ev.lambda[j] for j in eachindex(q)))
        for pulse in unique(Float64.(traces.pulse_fraction[
                (traces.candidate.==label).&(traces.event_bus.==bus)]))
            sub=traces[(traces.candidate.==label).&(traces.event_bus.==bus).&
                (traces.pulse_fraction.==pulse),:]
            linear=[pulse*(step(t-1.0)-step(t-1.1)) for t in sub.time_s]
            nonlinear=Float64.(sub.COI_frequency_Hz)
            maxerr=maximum(abs.(linear.-nonlinear))
            maxopp=maximum(abs.(-linear.-nonlinear))
            ref=max(maximum(abs.(nonlinear)),1e-12)
            push!(rows,(candidate=label,event_bus=bus,pulse_fraction=pulse,
                frozen_load_P_MW=lr.initialized_load_MW,
                frozen_load_Q_Mvar=lr.initialized_load_Mvar,
                max_linear_frequency_Hz=maximum(abs.(linear)),
                max_nonlinear_frequency_Hz=maximum(abs.(nonlinear)),
                max_abs_mismatch_Hz=maxerr,relative_mismatch=maxerr/ref,
                opposite_sign_relative_mismatch=maxopp/ref,
                status="POSTFREEZE_MODEL_MISMATCH_DIAGNOSTIC"))
            append!(out,[(candidate=label,event_bus=bus,pulse_fraction=pulse,
                time_s=Float64(sub.time_s[i]),linear_frequency_Hz=linear[i],
                nonlinear_frequency_Hz=nonlinear[i]) for i in 1:nrow(sub)])
        end
    end
end
CSV.write(joinpath(TABLES,"TABLE_K17_linear_nonlinear_mismatch.csv"),DataFrame(rows))
CSV.write(joinpath(TABLES,"TABLE_K18_linear_nonlinear_traces.csv"),DataFrame(out))
println("K_LINEAR_NONLINEAR_MAX_RELATIVE_MISMATCH: ",maximum(x->x.relative_mismatch,rows))
