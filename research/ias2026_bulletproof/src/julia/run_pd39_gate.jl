"""PowerDynamics IEEE-39 tutorial equilibrium gate.

This is deliberately a package/tutorial gate, not a claim of parity with the
frozen SG->GFL model. The official PowerDynamics Part I data and constructors
are used unchanged; the result is promoted only if the documented equilibrium
initialization succeeds.
"""

using PowerDynamics
using PowerDynamics.Library
using ModelingToolkitBase
using NetworkDynamics
using Graphs
using LinearAlgebra

example_dir = joinpath(pkgdir(PowerDynamics), "docs", "examples")
include(joinpath(example_dir, "ieee39_part1.jl"))

# The official tutorial documents these two structurally underconstrained buses.
formula = @initformula :ZIPLoad₊Vset = sqrt(:busbar₊u_r^2 + :busbar₊u_i^2)
set_initformula!(nw[VIndex(31)], formula)
set_initformula!(nw[VIndex(39)], formula)

pfs = solve_powerflow(nw)
interf = interface_values(pfs)
s0 = initialize_from_pf!(nw; verbose=false)

@assert nv(nw) == 39
@assert ne(nw) == 46
@assert s0 isa NWState

campaign_root = normpath(joinpath(@__DIR__, "..", ".."))
raw_dir = joinpath(campaign_root, "raw", "powerdynamics")
mkpath(raw_dir)
open(joinpath(raw_dir, "pd39_equilibrium_gate.md"), "w") do io
    println(io, "# PowerDynamics IEEE-39 equilibrium gate")
    println(io, "")
    println(io, "status: POWERDYNAMICS_VALIDATED")
    println(io, "package_version: 5.0.0")
    println(io, "julia_version: ", VERSION)
    println(io, "buses: ", nv(nw))
    println(io, "branches: ", ne(nw))
    println(io, "state_type: ", typeof(s0))
    println(io, "powerflow_interface_entries: ", length(interf))
    println(io, "source: PowerDynamics official docs/examples/ieee39_part1.jl")
    println(io, "scope: package/tutorial equilibrium gate; not SG->GFL parity")
end

println("POWERDYNAMICS_GATE_PASS buses=", nv(nw), " branches=", ne(nw), " state_type=", typeof(s0))
