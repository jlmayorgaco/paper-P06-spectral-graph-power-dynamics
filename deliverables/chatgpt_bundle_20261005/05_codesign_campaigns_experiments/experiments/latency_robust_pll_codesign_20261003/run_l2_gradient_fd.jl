using CSV, DataFrames, LinearAlgebra, TOML

module Roots
include(joinpath(@__DIR__, "..", "optimal_latency_margin_20261002", "run_fast_root_locus.jl"))
end

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
BLAS.set_num_threads(1)

function root_crossing(ctx,rho,kp,ki,s0,tau0)
    L=Roots.MegaOracle.DC.linearization(ctx,rho,kp,ki)
    m=Roots.reduced_model(L)
    t=tau0;s=ComplexF64(s0)
    for j in 1:15
        r=Roots.refine(m,s,fill(t/1000,10);tol=1e-12)
        r.converged || return nothing
        err=real(r.s)+0.05
        abs(err)<2e-8 && return (;tau_ms=t,s=r.s,residual=r.sv)
        slope=real(r.root_tau_sensitivity)/1000
        (isfinite(slope) && slope>0.05) || return nothing
        change=clamp(-err/slope,-0.2,0.2)
        t+=change;s=r.s+change/1000*r.root_tau_sensitivity
        30<t<45 || return nothing
    end
    nothing
end

function main()
    ctx=Roots.MegaOracle.DC.R.N.design_context(ROOT)
    baseline=CSV.read(joinpath(@__DIR__,"L0_REPRODUCTION.csv"),DataFrame)
    grad=CSV.read(joinpath(@__DIR__,"L2_BASELINE_ANALYTIC_GRADIENTS.csv"),DataFrame)
    tests=Dict("Z"=>[(:Kp,35),(:Kp,36),(:Ki,35),(:Ki,36)],
               "N"=>[(:Kp,30),(:Ki,30),(:Kp,37),(:Ki,37),(:Kp,39)])
    out=NamedTuple[]
    for id in ("Z","N")
        d=TOML.parsefile(joinpath(@__DIR__,"designs",id=="Z" ?
            "Z_zero_delay_tuned.toml" : "N_nominal.toml"))
        rho,kp,ki=Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"])
        row=only(eachrow(baseline[baseline.design_id.==id,:]))
        tau0=Float64(row.tau_crit_local_ms)
        s0=ComplexF64(row.critical_root_real+im*row.critical_root_imag)
        for (kind,bus) in tests[id]
            i=bus-29;k=kind==:Kp ? kp[i] : ki[i];h=0.001*k
            results=Dict{String,Any}()
            for (side,sign) in (("plus",1),("minus",-1))
                kp1=copy(kp);ki1=copy(ki)
                if kind==:Kp;kp1[i]+=sign*h;else;ki1[i]+=sign*h;end
                results[side]=root_crossing(ctx,rho,kp1,ki1,s0,tau0)
            end
            gp=only(eachrow(grad[(grad.design_id.==id).&(grad.bus.==bus),:]))
            analytic=kind==:Kp ? gp.analytic_d_tau_margin_ms_d_Kp : gp.analytic_d_tau_margin_ms_d_Ki
            ok=results["plus"]!==nothing && results["minus"]!==nothing
            fd=ok ? (results["plus"].tau_ms-results["minus"].tau_ms)/(2h) : NaN
            relative=ok ? abs(fd-analytic)/max(abs(analytic),abs(fd),1e-12) : NaN
            signagree=ok && (sign(fd)==sign(analytic))
            push!(out,(;design_id=id,bus,gain_kind=String(kind),gain_value=k,perturbation=h,
                analytic_margin_gradient_ms_per_gain=analytic,
                FD_margin_gradient_ms_per_gain=fd,
                plus_local_crossing_ms=ok ? results["plus"].tau_ms : NaN,
                minus_local_crossing_ms=ok ? results["minus"].tau_ms : NaN,
                relative_error=relative,sign_agreement=signagree,
                status=ok && relative<=0.05 && signagree ? "LOCAL_GRADIENT_VALIDATED" :
                    "FAIL_OR_SWITCH_OR_NONCONVERGENCE",
                scope="LOCAL_ROOT_BRANCH;BASELINE_BOUNDARY_FULL_CONTOUR_VERIFIED;PERTURBATION_COUNTS_PENDING"))
            CSV.write(joinpath(@__DIR__,"L2_TAU_MARGIN_GRADIENT_VALIDATION.csv"),DataFrame(out))
            println("L2_FD ",last(out));flush(stdout)
        end
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
