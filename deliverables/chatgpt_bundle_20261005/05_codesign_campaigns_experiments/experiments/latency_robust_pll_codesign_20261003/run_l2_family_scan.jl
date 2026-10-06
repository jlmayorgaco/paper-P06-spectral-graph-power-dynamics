using CSV, DataFrames, LinearAlgebra, TOML

module Roots
include(joinpath(@__DIR__, "..", "optimal_latency_margin_20261002", "run_fast_root_locus.jl"))
end

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const PREV=joinpath(ROOT,"experiments","optimal_latency_margin_20261002")
BLAS.set_num_threads(1)

function modevector(L,s,tau_ms)
    D=Roots.MegaOracle.DC.delta_matrix(L,s,fill(tau_ms/1000,10))
    v=ComplexF64.(L.Q*svd(D).V[:,end]);v/norm(v)
end

function root_crossing(m,s_seed;lo=34.0,hi=45.0,initial=40.0)
    tau=initial;s=ComplexF64(s_seed)
    for j in 1:50
        r=Roots.refine(m,s,fill(tau/1000,10);tol=1e-10)
        r.converged || return nothing
        alpha_err=real(r.s)+0.05
        if abs(alpha_err)<2e-6
            return (;tau_ms=tau,s=r.s,residual=r.sv,steps=j)
        end
        slope=real(r.root_tau_sensitivity)/1000
        (isfinite(slope) && slope>0.05) || return nothing
        change=clamp(-alpha_err/slope,-0.25,0.25)
        next=clamp(tau+change,lo,hi)
        abs(next-tau)<1e-8 && return nothing
        tau=next;s=r.s+change/1000*r.root_tau_sensitivity
    end
    nothing
end

function main()
    ctx=Roots.MegaOracle.DC.R.N.design_context(ROOT)
    z=TOML.parsefile(joinpath(@__DIR__,"designs","Z_zero_delay_tuned.toml"))
    n=TOML.parsefile(joinpath(@__DIR__,"designs","N_nominal.toml"))
    rho=Float64.(z["rho"]);z["rho"]==n["rho"] || error("rho mismatch")
    L0=CSV.read(joinpath(@__DIR__,"L0_REPRODUCTION.csv"),DataFrame)
    templates=Dict{String,Vector{ComplexF64}}()
    for id in ("Z","N")
        d=id=="Z" ? z : n
        L=Roots.MegaOracle.DC.linearization(ctx,rho,Float64.(d["Kp"]),Float64.(d["Ki"]))
        row=only(eachrow(L0[L0.design_id.==id,:]))
        templates[id]=modevector(L,row.critical_root_real+im*row.critical_root_imag,row.tau_crit_local_ms)
    end
    prev=CSV.read(joinpath(PREV,"T08_FAST_MODE_PROVENANCE.csv"),DataFrame)
    seedrows=prev[(abs.(prev.tau_ms.-40.0).<1e-8),:]
    seeds=ComplexF64.(seedrows.root_real.+im.*seedrows.root_imag)
    rows=NamedTuple[]; active=NamedTuple[]
    for eta in 0.0:0.1:1.0
        kp=exp.((1-eta).*log.(Float64.(z["Kp"])).+eta.*log.(Float64.(n["Kp"])))
        ki=exp.((1-eta).*log.(Float64.(z["Ki"])).+eta.*log.(Float64.(n["Ki"])))
        L=Roots.MegaOracle.DC.linearization(ctx,rho,kp,ki)
        m=Roots.reduced_model(L)
        roots=NamedTuple[]
        for (seed_id,seed) in enumerate(seeds)
            r40=Roots.refine(m,seed,fill(0.040,10);tol=1e-10)
            r40.converged || continue
            any(abs(r40.s-x.s40)<1e-3 for x in roots) && continue
            cross=root_crossing(m,r40.s)
            cross===nothing && continue
            v=modevector(L,cross.s,cross.tau_ms)
            mz=abs(dot(templates["Z"],v))^2
            mn=abs(dot(templates["N"],v))^2
            label=mz>mn ? "Z_template_like" : "N_template_like"
            push!(roots,(;s40=r40.s,seed_id,tau_ms=cross.tau_ms,real=real(cross.s),
                imag=imag(cross.s),frequency_hz=imag(cross.s)/(2pi),residual=cross.residual,
                MAC_Z=mz,MAC_N=mn,label))
        end
        sort!(roots,by=x->x.tau_ms)
        for (j,x) in enumerate(roots)
            push!(rows,(;eta,rank_by_local_crossing=j,root_seed_id=x.seed_id,
                local_crossing_ms=x.tau_ms,critical_real=x.real,frequency_hz=x.frequency_hz,
                small_factor_residual=x.residual,MAC_Z=x.MAC_Z,MAC_N=x.MAC_N,
                family_label=x.label,status="LOCAL_EXACT_CHARACTERISTIC_NOT_FULL_CONTOUR"))
        end
        if !isempty(roots)
            x=first(roots)
            push!(active,(;eta,local_earliest_crossing_ms=x.tau_ms,
                active_family_label=x.label,active_frequency_hz=x.frequency_hz,
                competing_family_gap_ms=length(roots)>1 ? roots[2].tau_ms-x.tau_ms : NaN,
                discovered_distinct_fast_roots=length(roots),
                status="LOCAL_ROOT_CATALOGUE_NOT_COMPLETE_CONTOUR_CERTIFICATE"))
        else
            push!(active,(;eta,local_earliest_crossing_ms=NaN,
                active_family_label="INDETERMINATE",active_frequency_hz=NaN,
                competing_family_gap_ms=NaN,discovered_distinct_fast_roots=0,
                status="NO_ROOTS_REFINED"))
        end
        CSV.write(joinpath(@__DIR__,"L2_FAMILY_SCAN.csv"),DataFrame(rows))
        CSV.write(joinpath(@__DIR__,"L2_ACTIVE_FAMILY_SCAN.csv"),DataFrame(active))
        println("L2_SCAN_ETA ",eta," roots=",length(roots)," active=",last(active));flush(stdout)
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
