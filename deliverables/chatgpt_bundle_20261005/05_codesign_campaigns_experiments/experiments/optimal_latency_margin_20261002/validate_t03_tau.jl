using CSV, DataFrames, LinearAlgebra, TOML

module FastRoots
include(joinpath(@__DIR__,"run_fast_root_locus.jl"))
end

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const MEGA=joinpath(ROOT,"experiments","analytical_delay_codesign_mega_20261002")
BLAS.set_num_threads(1)

function main()
    ctx=FastRoots.MegaOracle.DC.R.N.design_context(ROOT)
    crossings=CSV.read(joinpath(@__DIR__,"T02_PRECISE_CROSSINGS.csv"),DataFrame)
    derivatives=CSV.read(joinpath(@__DIR__,"T03_CRITICAL_MODE_FAMILIES.csv"),DataFrame)
    rows=NamedTuple[]
    for point in eachrow(crossings)
        id=String(point.design_id)
        path=id=="seed_875" ? joinpath(MEGA,"seed_uniform_875.toml") :
             joinpath(MEGA,"M1_ZERO_DELAY_DESIGN.toml")
        d=TOML.parsefile(path)
        L=FastRoots.MegaOracle.DC.linearization(ctx,Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"]))
        model=FastRoots.reduced_model(L)
        s=point.critical_root_real+im*point.critical_root_imag
        tau=fill(point.local_root_crossing_ms/1000,length(L.Ai))
        for i in (1,6,7,8,10)
            h=1e-6
            tp=copy(tau);tm=copy(tau);tp[i]+=h;tm[i]-=h
            rp=FastRoots.refine(model,s,tp;tol=1e-11)
            rm=FastRoots.refine(model,s,tm;tol=1e-11)
            fd=(rp.s-rm.s)/(2h)
            analytical=only(derivatives.d_real_lambda_d_tau_i_per_s[
                (derivatives.design_id.==id).&(derivatives.bus.==29+i)])
            rel=abs(real(fd)-analytical)/max(abs(real(fd)),abs(analytical),1e-12)
            push!(rows,(;design_id=id,bus=29+i,analytic_real=analytical,
                centered_FD_real=real(fd),absolute_error=abs(real(fd)-analytical),
                relative_error=rel,plus_converged=rp.converged,minus_converged=rm.converged,
                perturbation_s=h,status=(rp.converged&&rm.converged&&rel<=0.01) ?
                    "NUMERICALLY_VALIDATED" : "INDETERMINATE"))
            CSV.write(joinpath(@__DIR__,"T03_TAU_DERIVATIVE_VALIDATION.csv"),DataFrame(rows))
            println("T03_TAU_FD ",last(rows));flush(stdout)
        end
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
