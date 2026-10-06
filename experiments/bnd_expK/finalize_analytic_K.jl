using LinearAlgebra, CSV, DataFrames, SHA, TOML, Statistics
BLAS.set_num_threads(1)
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_design_k","BNDDesignK.jl"))
using .BNDDesignK
const OUT=joinpath(ROOT,"reports","experiment_K")
const TABLES=joinpath(OUT,"tables")
ctx=design_context(ROOT)

function frozen(path)
    isfile(path) && isfile(path*".sha256") || error("missing freeze: $path")
    hash=bytes2hex(sha256(read(path)))
    hash==strip(read(path*".sha256",String)) || error("hash mismatch: $path")
    return TOML.parsefile(path),hash
end
nom,nomhash=frozen(joinpath(OUT,"Z_K_NOMINAL_FINAL.toml"))
ev=evaluate(ctx,Float64.(nom["epsilon"]),Float64.(nom["Kp"]),Float64.(nom["Ki"]);vectors=true)

# The frozen comparison artifacts are read only after the ExpK candidate hashes
# exist. Prior inspection of summary text is disclosed in the candidate files.
erows=CSV.read(joinpath(ROOT,"reports","experiment_E","tables",
    "TABLE_E14_provisional_per_generator_design.csv"),DataFrame)
grows=CSV.read(joinpath(ROOT,"reports","experiment_G","tables",
    "TABLE_G13_final_per_generator_design.csv"),DataFrame)
sort!(erows,:bus);sort!(grows,:bus)
ee=evaluate(ctx,1 .- Float64.(erows.rho),Float64.(erows.Kp),Float64.(erows.Ki))
gg=evaluate(ctx,Float64.(grows.epsilon_retained),Float64.(grows.Kp),Float64.(grows.Ki))
comparison=DataFrame(case=["ExpK_nominal","ExpE_provisional","ExpG_final"],
    support_mask=[ev.support,ee.support,gg.support],
    retained_SG_MW=[ev.retained_mw,ee.retained_mw,gg.retained_mw],
    GFL_MW=[ev.gfl_mw,ee.gfl_mw,gg.gfl_mw],
    exact_ExpK_alpha_s_inv=[ev.alpha,ee.alpha,gg.alpha],
    meets_margin=[ev.alpha<=-0.05,ee.alpha<=-0.05,gg.alpha<=-0.05],
    KKT_certified=[false,false,false],
    provenance=["frozen_candidate","read_after_ExpK_hash","read_after_ExpK_hash"])
CSV.write(joinpath(TABLES,"TABLE_K11_blinded_comparison.csv"),comparison)

# Local dual values are defined only on continuous coordinates of this fixed
# support. Missing buses are architecture choices and receive no derivative.
dual=NamedTuple[]
active=ev.active
grads=[fixed_support_gradient(ctx,ev,j) for j in active]
interior=[i for i in eachindex(ctx.buses) if 0<ev.epsilon[i]<1]
y=zeros(length(active))
if length(interior)>0 && length(active)>0
    G=[real(grads[j].epsilon[i]) for i in interior,j in eachindex(active)]
    # Nonnegative least-squares is not invoked: this is the unconstrained
    # pseudoinverse diagnostic, clipped only for reporting feasibility of y.
    trial=pinv(G)*(-ctx.power[interior])
    if all(trial .>= 0);y=trial end
end
for (i,b) in enumerate(ctx.buses)
    if ev.epsilon[i]==0
        push!(dual,(bus=b,epsilon=0.0,P_MW=ctx.power[i],
            spectral_shadow_price=NaN,rocof_shadow_price=0.0,
            Psi_MW=NaN,multiplier_sum=sum(y),
            status="DISCRETE_REMOVAL_NO_CONTINUOUS_DERIVATIVE"))
    else
        spectral=sum(y[j]*real(grads[j].epsilon[i]) for j in eachindex(active))
        psi=-spectral-ctx.power[i]
        push!(dual,(bus=b,epsilon=ev.epsilon[i],P_MW=ctx.power[i],
            spectral_shadow_price=spectral,rocof_shadow_price=0.0,
            Psi_MW=psi,multiplier_sum=sum(y),
            status="DIAGNOSTIC_MULTIPLIER_NOT_KKT_CERTIFIED"))
    end
end
CSV.write(joinpath(TABLES,"TABLE_K06_dual_survival_values.csv"),DataFrame(dual))

"Frequency and RoCoF step transfer from a one MW active-load injection."
function transient(ev,event_bus)
    weights=Float64[];states=Int[]
    for (i,b) in enumerate(ctx.buses)
        ev.epsilon[i]>0 || continue
        row=only(eachrow(ev.model.state_map[(ev.model.state_map.bus .== b) .&
            (ev.model.state_map.kind .== "SG"),:]))
        h=ctx.net.sg[b].op.parameters.inertia*ctx.net.sg[b].op.parameters.rating_mva
        push!(weights,ev.epsilon[i]*h)
        push!(states,Int(row.last)-1)
    end
    total=sum(weights)
    total>0 || return nothing
    Cf=zeros(ev.model.n_dynamic)
    for (state,w) in zip(states,weights);Cf[state]=60w/total end
    V=ctx.net.voltage[event_bus]
    Ii=conj((1/100)/V)
    iv=zeros(size(ev.model.Gy,1));iv[2event_bus-1]=real(Ii);iv[2event_bus]=imag(Ii)
    Bd=ev.model.B*(ev.model.Gy\iv)
    Cq=transpose(ev.Q)*Cf;Bq=transpose(ev.Q)*Bd
    R=ev.right;L=ev.left;λ=ev.lambda
    q=[dot(Cq,R[:,j])*dot(L[:,j],Bq)/dot(L[:,j],R[:,j]) for j in eachindex(λ)]
    times=collect(0.0:0.025:80.0)
    rocof=[real(sum(q[j]*exp(λ[j]*t) for j in eachindex(λ))) for t in times]
    frequency=[real(sum(abs(λ[j])<1e-10 ? q[j]*t : q[j]*expm1(λ[j]*t)/λ[j]
        for j in eachindex(λ))) for t in times]
    peakR=maximum(abs.(rocof));peakF=maximum(abs.(frequency))
    settling=NaN
    threshold=0.02peakF
    for k in eachindex(times)
        if maximum(abs.(frequency[k:end]))<=threshold
            settling=times[k];break
        end
    end
    (;event_bus,times,rocof,frequency,peakR,peakF,nadir=minimum(frequency),
        settling,coi_inertia=total,residues=q,eigenvalues=λ,
        deltaP_max_rocof=0.5/max(peakR,eps()),
        deltaP_max_frequency=0.5/max(peakF,eps()))
end

loadbuses=sort(unique(Int.(ctx.net.load_audit.bus)))
transientrows=NamedTuple[]
responses=Dict{Int,Any}()
for b in loadbuses
    tr=transient(ev,b)
    tr===nothing && continue
    responses[b]=tr
    push!(transientrows,(candidate="nominal",bus=b,
        peak_RoCoF_Hz_s_per_MW=tr.peakR,
        peak_frequency_Hz_per_MW=tr.peakF,
        nadir_Hz_per_MW=tr.nadir,settling_s=tr.settling,
        deltaP_max_RoCoF_MW=tr.deltaP_max_rocof,
        deltaP_max_frequency_MW=tr.deltaP_max_frequency,
        deltaP_max_total_MW=min(tr.deltaP_max_rocof,tr.deltaP_max_frequency)))
end
robust_file=joinpath(OUT,"Z_K_ROBUST_beta_1p6991206999182038em6.toml")
if isfile(robust_file)
    robust,_=frozen(robust_file)
    rev=evaluate(ctx,Float64.(robust["epsilon"]),Float64.(robust["Kp"]),
        Float64.(robust["Ki"]);vectors=true)
    for b in loadbuses
        rt=transient(rev,b)
        rt===nothing && continue
        push!(transientrows,(candidate="robust_ExpG_beta",bus=b,
            peak_RoCoF_Hz_s_per_MW=rt.peakR,
            peak_frequency_Hz_per_MW=rt.peakF,
            nadir_Hz_per_MW=rt.nadir,settling_s=rt.settling,
            deltaP_max_RoCoF_MW=rt.deltaP_max_rocof,
            deltaP_max_frequency_MW=rt.deltaP_max_frequency,
            deltaP_max_total_MW=min(rt.deltaP_max_rocof,rt.deltaP_max_frequency)))
    end
end
CSV.write(joinpath(TABLES,"TABLE_K14_transient_capacity.csv"),DataFrame(transientrows))
worst=sort(filter(x->x.candidate=="nominal",transientrows);
    by=x->x.deltaP_max_total_MW)[1]
tr=responses[worst.bus]
CSV.write(joinpath(TABLES,"TABLE_K15_analytic_transient_trace.csv"),
    DataFrame(time_s=tr.times,frequency_Hz_per_MW=tr.frequency,
        RoCoF_Hz_s_per_MW=tr.rocof,event_bus=fill(worst.bus,length(tr.times))))

# Port Schur decomposition at the limiting pole. The direct and self-energy
# derivatives are computed from the same exact 20-coordinate port closure.
explain=NamedTuple[]
for j in ev.active
    pole=ev.lambda[j]
    pc=BNDDesignK.BNDDesignG.ClosureSpectrum.port_closure(ctx.net,pole,
        1 .- ev.epsilon,ev.kp,ev.ki;derivatives=true)
    C=pc.C;Cs=pc.Cs
    F=svd(C);v=F.V[:,end]
    busindex=argmax([norm(v[2k-1:2k]) for k in eachindex(ctx.buses)])
    kk=collect(2busindex-1:2busindex)
    rr=setdiff(collect(1:size(C,1)),kk)
    R=C[rr,rr]\Matrix{ComplexF64}(I,length(rr),length(rr))
    Ckr=C[kk,rr];Crk=C[rr,kk]
    Γ=-Ckr*R*Crk
    Teff=C[kk,kk]+Γ
    Fs=svd(Teff);u=Fs.U[:,end];w=Fs.V[:,end]
    function split(dC)
        dΓ=-dC[kk,rr]*R*Crk-Ckr*R*dC[rr,kk]+Ckr*R*dC[rr,rr]*R*Crk
        return dC[kk,kk],dΓ
    end
    direct_s,self_s=split(Cs)
    den=dot(u,(direct_s+self_s)*w)
    cancellation=sum(norm(C[kk,2k-1:2k]*R[findall(in(2k-1:2k),rr),:])
        for k in eachindex(ctx.buses) if k!=busindex) # descriptive only
    for (i,b) in enumerate(ctx.buses)
        ev.epsilon[i]>0 || continue
        direct,self=split(-pc.d_rho[i])
        ds_direct=-dot(u,direct*w)/den
        ds_self=-dot(u,self*w)/den
        push!(explain,(pole_real=real(pole),pole_imag=imag(pole),
            pivot_bus=ctx.buses[busindex],parameter_bus=b,
            direct_dlambda_real=real(ds_direct),self_energy_dlambda_real=real(ds_self),
            total_dlambda_real=real(ds_direct+ds_self),
            self_energy_share=abs(ds_self)/max(abs(ds_direct)+abs(ds_self),eps()),
            port_schur_sigma_min=Fs.S[end],
            cancellation_ratio=cancellation/max(norm(Γ),eps()),
            status="EXACT_LOCAL_SCHUR_DERIVATIVE"))
    end
end
CSV.write(joinpath(TABLES,"TABLE_K10_BND_explanation.csv"),DataFrame(explain))

println("EXP_K_POSTFREEZE_COMPARISON: E_alpha=",ee.alpha," G_alpha=",gg.alpha)
println("WORST_DISTURBANCE_BUS: ",worst.bus)
println("ROCOF_UNIT_GAIN: ",worst.peak_RoCoF_Hz_s_per_MW)
println("FREQUENCY_UNIT_GAIN: ",worst.peak_frequency_Hz_per_MW)
println("DELTA_P_MAX_TOTAL: ",worst.deltaP_max_total_MW)
