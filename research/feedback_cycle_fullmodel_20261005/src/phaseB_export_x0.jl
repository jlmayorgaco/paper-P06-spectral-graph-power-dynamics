using LinearAlgebra, CSV, DataFrames, TOML
const EXP = joinpath(@__DIR__, "..", "..", "..", "experiments")
include(joinpath(EXP, "graph_gsp_codesign_20261003", "DelayedEvents.jl"))
const R = DelayedEvents.R
ctx = R.N.design_context(R.ROOT)
d = TOML.parsefile(joinpath(EXP, "graph_gsp_codesign_20261003", "baseline.toml"))
kp, ki = Float64.(d["Kp"]), Float64.(d["Ki"])
for r in (0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875)
    m = R.model(ctx, fill(r, 10), kp, ki; dc_convention=:physical_supply)
    CSV.write(joinpath(@__DIR__, "..", "raw", "jac", "x0_rho$(replace(string(r), "." => "p")).csv"), DataFrame(x0=m.x0))
end
