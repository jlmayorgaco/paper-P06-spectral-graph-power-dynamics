using CSV
using DataFrames
using Dates
using PowerDynamics
using TOML

include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
using .PD39

const OUT_DIR = joinpath(@__DIR__, "..", "..", "results", "pd39", "portfolio_campaign")
mkpath(OUT_DIR)

base_nw = baseline_network()
base_eq = initialize_equilibrium(base_nw; sparse = false)
base_eq.powerflow_finite && base_eq.state_finite && base_eq.fixed_point ||
    error("baseline qualification must pass before campaign")

results = evaluate_campaign(
    base_nw;
    pfs = base_eq.pfs,
    portfolios = candidate_portfolios(),
    scenarios = uncertainty_scenarios(),
    margin_target = ROBUST_MARGIN_TARGET,
)
CSV.write(joinpath(OUT_DIR, "portfolio_scenario_results.csv"), results)
CSV.write(joinpath(OUT_DIR, "static_nodes.csv"), static_weak_nodes())
CSV.write(joinpath(OUT_DIR, "static_links.csv"), static_weak_links())

design = minimum_intervention_design(results; margin_target = ROBUST_MARGIN_TARGET)
open(joinpath(OUT_DIR, "design_summary.toml"), "w") do io
    TOML.print(io, Dict(
        "timestamp_utc" => string(now(UTC)),
        "portfolio_count" => length(candidate_portfolios()),
        "scenario_count" => length(uncertainty_scenarios()),
        "margin_target" => ROBUST_MARGIN_TARGET,
        "robust_design_found" => !isnothing(design),
        "selected_portfolio" => isnothing(design) ? "" : design.portfolio,
        "selected_intervention_count" => isnothing(design) ? -1 : design.intervention_count,
        "selected_robust_margin" => isnothing(design) ? NaN : design.robust_margin,
    ))
end

println("PD39 portfolio campaign complete: robust_design_found=$( !isnothing(design) )")
