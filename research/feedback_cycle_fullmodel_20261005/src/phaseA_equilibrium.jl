# Phase A1: equilibrium parity of the frozen model (read-only use of ReducedDAE)
using LinearAlgebra, CSV, DataFrames, TOML
const EXP = joinpath(@__DIR__, "..", "..", "..", "experiments")
include(joinpath(EXP, "graph_gsp_codesign_20261003", "DelayedEvents.jl"))
const R = DelayedEvents.R
d = TOML.parsefile(joinpath(EXP, "graph_gsp_codesign_20261003", "baseline.toml"))
rho, kp, ki = Float64.(d["rho"]), Float64.(d["Kp"]), Float64.(d["Ki"])
ctx = R.N.design_context(R.ROOT)
m = R.model(ctx, rho, kp, ki; bus=8, delta=0.0, dc_convention=:physical_supply)
x0 = copy(m.x0); dx = zeros(length(x0)); R.rhs!(dx, x0, m, 0.0)
v = R.voltage(x0, m; allbus=true); vm = hypot.(v[1:2:end], v[2:2:end])
rows = [("n_states", length(x0)), ("rhs_inf_norm", norm(dx, Inf)), ("vmin_pu", minimum(vm)), ("vmax_pu", maximum(vm)), ("n_buses", length(vm))]
out = joinpath(@__DIR__, "..", "derived", "TABLE02a_equilibrium_julia.csv")
CSV.write(out, DataFrame(item=first.(rows), value=last.(rows))); foreach(println, rows)
