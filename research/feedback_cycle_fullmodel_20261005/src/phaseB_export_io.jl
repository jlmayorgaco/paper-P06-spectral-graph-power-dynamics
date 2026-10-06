# Linear I/O of the frozen model: load-admittance inputs at buses 8,16,29 (per MW), outputs = all-bus complex voltage (re,im).
using LinearAlgebra, CSV, DataFrames, TOML
const EXP = joinpath(@__DIR__, "..", "..", "..", "experiments")
include(joinpath(EXP, "graph_gsp_codesign_20261003", "DelayedEvents.jl"))
const R = DelayedEvents.R
ctx = R.N.design_context(R.ROOT)
d = TOML.parsefile(joinpath(EXP, "graph_gsp_codesign_20261003", "baseline.toml"))
kp, ki = Float64.(d["Kp"]), Float64.(d["Ki"])
out = joinpath(@__DIR__, "..", "raw", "jac")
for r in (0.0, 0.875)
    rho = fill(r, 10); tag = replace(string(r), "." => "p")
    m0 = R.model(ctx, rho, kp, ki; bus=8, delta=0.0, dc_convention=:physical_supply)
    de = R.derivatives(m0.x0, m0); lift = vcat(m0.net.lift, Matrix{Float64}(I, 20, 20))
    CSV.write(joinpath(out, "Vall_x_rho$(tag).csv"), DataFrame(lift * de.vx, :auto))
    CSV.write(joinpath(out, "Vall_0_rho$(tag).csv"), DataFrame(v=lift * de.v))
    dx0 = zeros(length(m0.x0)); R.rhs!(dx0, m0.x0, m0, 0.0); v0 = R.voltage(m0.x0, m0; allbus=true)
    cols = Dict{String,Vector{Float64}}(); lin = Float64[]
    for bus in (8, 16, 29), mw in (1.0, 0.1)
        mb = R.model(ctx, rho, kp, ki; bus, delta=mw, dc_convention=:physical_supply)
        dxb = zeros(length(m0.x0)); R.rhs!(dxb, m0.x0, mb, 0.0)
        cols["dx_bus$(bus)_$(mw)"] = (dxb - dx0) / mw; cols["dv_bus$(bus)_$(mw)"] = (R.voltage(m0.x0, mb; allbus=true) - v0) / mw
    end
    for bus in (8, 16, 29)
        push!(lin, norm(cols["dx_bus$(bus)_1.0"] - cols["dx_bus$(bus)_0.1"]) / max(norm(cols["dx_bus$(bus)_1.0"]), 1e-30))
        CSV.write(joinpath(out, "inp_bus$(bus)_rho$(tag).csv"), DataFrame(dx=cols["dx_bus$(bus)_0.1"], dummy=0.0)[:, [:dx]])
        CSV.write(joinpath(out, "inpv_bus$(bus)_rho$(tag).csv"), DataFrame(dv=cols["dv_bus$(bus)_0.1"]))
    end
    println("rho=", r, " linearity (rel diff 1MW vs 0.1MW)=", lin)
end
