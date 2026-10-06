# Phase A2/B/C: export exact Schur-eliminated Jacobians (Fx) of the frozen model for rho = 0 (all-SG), a rho grid, and the baseline.
# Delay loop (PLL error injection) is NOT in Fx; it is added in Python with exact exponentials.
using LinearAlgebra, CSV, DataFrames, TOML
const EXP = joinpath(@__DIR__, "..", "..", "..", "experiments")
include(joinpath(EXP, "graph_gsp_codesign_20261003", "DelayedEvents.jl"))
const R = DelayedEvents.R
ctx = R.N.design_context(R.ROOT)
d = TOML.parsefile(joinpath(EXP, "graph_gsp_codesign_20261003", "baseline.toml"))
kp, ki = Float64.(d["Kp"]), Float64.(d["Ki"])
out = joinpath(@__DIR__, "..", "raw", "jac"); mkpath(out)
for r in [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875]
    rho = fill(r, 10)
    m = R.model(ctx, rho, kp, ki; bus=8, delta=0.0, dc_convention=:physical_supply)
    de = R.derivatives(m.x0, m)
    tag = replace(string(r), "." => "p")
    CSV.write(joinpath(out, "Fx_rho$(tag).csv"), DataFrame(de.Fx, :auto))
    # state map: device kind and site per state
    rows = NamedTuple[]
    for i in 1:10
        for (k, ix) in enumerate(m.sgidx[i]); push!(rows, (; index=ix, site=29 + i, kind="SG", local_pos=k)); end
        for (k, ix) in enumerate(m.gfidx[i]); push!(rows, (; index=ix, site=29 + i, kind="GFL", local_pos=k)); end
    end
    CSV.write(joinpath(out, "statemap_rho$(tag).csv"), DataFrame(rows))
    println("rho=", r, " n=", length(m.x0), " rhs_inf=", begin dx = zeros(length(m.x0)); R.rhs!(dx, m.x0, m, 0.0); norm(dx, Inf) end, " eigmax=", maximum(real.(eigvals(de.Fx))))
end
