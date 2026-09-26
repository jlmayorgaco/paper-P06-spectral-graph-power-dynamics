using CSV
using DataFrames
using Dates
using PowerDynamics
using TOML

include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
using .PD39

const OUT_DIR = joinpath(@__DIR__, "..", "..", "results", "pd39", "single_replacement")
mkpath(OUT_DIR)

base_nw = baseline_network()
base_eq = initialize_equilibrium(base_nw; sparse = false)
base_stab = stability_audit(base_eq.state)
rows = NamedTuple[]

for bus in CANDIDATE_SG_BUSES
    nw = replace_buses(base_nw, [bus])
    eq = initialize_equilibrium(nw; pfs = base_eq.pfs, sparse = false)
    if !eq.powerflow_finite || !eq.state_finite || !eq.fixed_point
        push!(rows, (
            bus = bus,
            equilibrium_status = "equilibrium_failed",
            dynamic_margin = NaN,
            max_real = NaN,
            stable = false,
            margin_loss = NaN,
        ))
        continue
    end
    stab = stability_audit(eq.state)
    push!(rows, (
        bus = bus,
        equilibrium_status = "ok",
        dynamic_margin = stab.dynamic_margin,
        max_real = stab.max_real,
        stable = stab.stable,
        margin_loss = base_stab.dynamic_margin - stab.dynamic_margin,
    ))
end

result = DataFrame(rows)
CSV.write(joinpath(OUT_DIR, "single_replacement.csv"), result)
CSV.write(joinpath(OUT_DIR, "static_nodes.csv"), static_weak_nodes())
open(joinpath(OUT_DIR, "metadata.toml"), "w") do io
    TOML.print(io, Dict(
        "timestamp_utc" => string(now(UTC)),
        "baseline_margin" => base_stab.dynamic_margin,
        "gauge_tol" => 1e-8,
        "margin_tol" => 1e-8,
    ))
end

println("PD39 single-replacement campaign complete")
