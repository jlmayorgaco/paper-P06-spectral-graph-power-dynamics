# Verify source/data closure in an extracted bundle, without rerunning events.
using LinearAlgebra, CSV, DataFrames, TOML, ForwardDiff
root = isempty(ARGS) ? normpath(joinpath(@__DIR__, "..", "..")) : abspath(ARGS[1])
include(joinpath(root, "experiments", "graph_gsp_codesign_20261003", "DelayedEvents.jl"))
const R = DelayedEvents.R
BLAS.set_num_threads(1)
ctx = R.N.design_context(root)
out = joinpath(root, "experiments", "interaction_decision_20261004")
d = TOML.parsefile(joinpath(out, "designs", "corrected.toml"))
rho, kp, ki = Float64.(d["rho"]), Float64.(d["Kp"]), Float64.(d["Ki"])
m = R.model(ctx, rho, kp, ki; dc_convention=:physical_supply)
f = zeros(length(m.x0)); R.rhs!(f, m.x0, m, 0.)
A = ForwardDiff.jacobian(m.x0) do x
    dx = similar(x); R.rhs!(dx, x, m, 0.); dx
end
C = ForwardDiff.jacobian(x -> DelayedEvents.detector(x, m), m.x0)
names = ("Adev", "Bv", "Cs", "Cf", "Ds", "Y", "Etheta", "Hv", "Bp", "Bi")
p = Dict(n => Matrix(CSV.read(joinpath(root, "experiments", "graph_gsp_codesign_20261003", "model", n*".csv"), DataFrame)) for n in names)
w = repeat(1 .-rho, inner=2)
Vk = -(p["Y"] + Diagonal(w)*p["Ds"]) \ (Diagonal(w)*p["Cs"] + Diagonal(1 .-w)*p["Cf"])
Cp = p["Etheta"] + p["Hv"]*Vk
Ap = p["Adev"] + p["Bv"]*Vk + (p["Bp"]*Diagonal(kp) + p["Bi"]*Diagonal(ki))*Cp
result = Dict("equilibrium_residual" => norm(f, Inf), "jacobian_relative_error" => norm(A-Ap)/norm(A),
              "detector_relative_error" => norm(C-Cp)/norm(C), "GFL_MW" => dot(ctx.power, rho),
              "SG_MW" => dot(ctx.power, 1 .-rho), "julia_version" => string(VERSION))
@assert result["equilibrium_residual"] < 1e-7
@assert result["jacobian_relative_error"] < 1e-8
@assert result["detector_relative_error"] < 1e-8
TOML.print(stdout, result; sorted=true)
