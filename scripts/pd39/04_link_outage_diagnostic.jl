using CSV
using Dates
using DataFrames
using PowerDynamics
using TOML

include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
using .PD39

const OUT_DIR = joinpath(@__DIR__, "..", "..", "results", "pd39", "link_outage")
mkpath(OUT_DIR)

base_nw = baseline_network()
base_eq = initialize_equilibrium(base_nw; sparse = false)
base_stab = stability_audit(base_eq.state)
result = dynamic_weak_links(base_nw; base_margin = base_stab.dynamic_margin)
result = leftjoin(result, static_weak_links(), on = :link)
CSV.write(joinpath(OUT_DIR, "link_outage_results.csv"), result)
open(joinpath(OUT_DIR, "metadata.toml"), "w") do io
    TOML.print(io, Dict(
        "timestamp_utc" => string(now(UTC)),
        "branch_count" => nrow(result),
        "base_margin" => base_stab.dynamic_margin,
        "outage_rule" => "one_at_a_time_branch_removal_no_redispatch",
    ))
end

println("PD39 link-outage diagnostic complete")
