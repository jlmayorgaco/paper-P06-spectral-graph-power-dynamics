using CSV, DataFrames, TOML, LinearAlgebra
module Scan
include(joinpath(@__DIR__, "run_l2_family_scan.jl"))
end
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
BLAS.set_num_threads(1)

function main()
    length(ARGS)==1 || error("usage: export_linear_dde.jl design.toml")
    design=abspath(ARGS[1]); id=splitext(basename(design))[1]
    d=TOML.parsefile(design)
    ctx=Scan.Roots.MegaOracle.DC.R.N.design_context(ROOT)
    L=Scan.Roots.MegaOracle.DC.linearization(ctx,Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"]))
    out=joinpath(@__DIR__,"linear_dde",id);mkpath(out)
    CSV.write(joinpath(out,"A0.csv"),DataFrame(L.A0,:auto))
    CSV.write(joinpath(out,"B.csv"),DataFrame(L.B,:auto))
    CSV.write(joinpath(out,"C.csv"),DataFrame(L.C,:auto))
    roots=CSV.read(joinpath(@__DIR__,"evaluations",id,"ROOTS.csv"),DataFrame)
    firstroot=roots[1,:]
    critical=Float64(firstroot.local_crossing_ms)
    seed=ComplexF64(firstroot.critical_real+im*firstroot.critical_imag)
    rows=NamedTuple[]
    model=Scan.Roots.reduced_model(L)
    for factor in (0.90,0.98,1.02)
        tau_ms=factor*critical
        r=Scan.Roots.refine(model,seed,fill(tau_ms/1000,10);tol=1e-11)
        r.converged || error("root continuation failed at $factor")
        D=Scan.Roots.MegaOracle.DC.delta_matrix(L,r.s,fill(tau_ms/1000,10))
        v=svd(D).V[:,end];v/=norm(v)
        CSV.write(joinpath(out,"mode_$(replace(string(factor),"."=>"p"))_real.csv"),DataFrame(x=real.(v)))
        CSV.write(joinpath(out,"mode_$(replace(string(factor),"."=>"p"))_imag.csv"),DataFrame(x=imag.(v)))
        push!(rows,(;factor,tau_ms,root_real=real(r.s),root_imag=imag(r.s),
            normalized_characteristic_residual=norm(D*v)/max(1,norm(D)*norm(v)),
            status="EXACT_CHARACTERISTIC_ROOT_WITH_NUMERICAL_RESIDUAL"))
    end
    CSV.write(joinpath(out,"roots.csv"),DataFrame(rows))
    println("EXPORTED_LINEAR_DDE ",id," dimension=",size(L.A0,1))
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
