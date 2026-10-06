"""Freeze ExpE generator discovery and independent gain bounds; no stability calls."""

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "src", "bnd_design_e", "PhysicalData.jl"))
using .PhysicalData

result = freeze_design_domain(ROOT)
replaceable = result.generators[result.generators.replaceable .== true, :]
println("EXP_E_DOMAIN_STATUS: ", result.status)
println("REPLACEABLE_GENERATOR_COUNT: ", size(replaceable, 1))
println("REPLACEABLE_BUSES: ", join(replaceable.bus, ","))
println("TOTAL_INITIAL_SG_MW: ", sum(replaceable.SG_dispatch_initial_MW))
println("INDEPENDENT_PLL_GAINS: YES")
println("DESIGN_DOMAIN_SHA256: ", result.sha256)
println("POWERDYNAMICS_IMPORTED: NO")
