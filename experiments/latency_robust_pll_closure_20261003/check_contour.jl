using CSV, DataFrames, TOML, LinearAlgebra
include(joinpath(@__DIR__,"evaluate_design.jl"))

function main_contour()
    length(ARGS)>=3 || error("usage: check_contour.jl design.toml tau_ms [tau_ms...]")
    d=TOML.parsefile(abspath(ARGS[1]))
    ctx=Scan.Roots.MegaOracle.DC.R.N.design_context(ROOT)
    L=Scan.Roots.MegaOracle.DC.linearization(ctx,Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"]))
    id=splitext(basename(ARGS[1]))[1]
    rows=NamedTuple[]
    for raw in ARGS[2:end]
        tau_ms=parse(Float64,raw)
        c=classify(L,tau_ms)
        push!(rows,merge((;design_id=id,tau_ms),c))
        println("CONTOUR_CHECK ",last(rows));flush(stdout)
    end
    CSV.write(joinpath(@__DIR__,"evaluations",id,"ADDITIONAL_COUNTS.csv"),DataFrame(rows))
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main_contour()
