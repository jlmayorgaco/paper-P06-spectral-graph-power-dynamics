using CSV
using DataFrames
using Dates
using Statistics
using TOML

const ROOT = joinpath(@__DIR__, "..", "..")
const RESULT_DIR = joinpath(ROOT, "results", "pd39")
const OUT_DIR = joinpath(RESULT_DIR, "summary")
mkpath(OUT_DIR)

function average_ranks(values)
    order = sortperm(values)
    ranks = zeros(Float64, length(values))
    i = 1
    while i <= length(order)
        j = i
        while j < length(order) && values[order[j + 1]] == values[order[i]]
            j += 1
        end
        ranks[order[i:j]] .= (i + j) / 2
        i = j + 1
    end
    return ranks
end

function spearman(x, y)
    rx = average_ranks(Float64.(x))
    ry = average_ranks(Float64.(y))
    return cor(rx, ry)
end

portfolio = CSV.read(joinpath(RESULT_DIR, "portfolio_campaign", "portfolio_scenario_results.csv"), DataFrame)
portfolio_summary = combine(
    groupby(portfolio, :portfolio),
    :intervention_count => first => :intervention_count,
    :dynamic_margin => minimum => :robust_margin,
    :robust_feasible_for_scenario => all => :robust_feasible,
    :equilibrium_status => (x -> all(x .== "ok")) => :all_equilibria_ok,
)
sort!(portfolio_summary, [:intervention_count, :portfolio])
CSV.write(joinpath(OUT_DIR, "portfolio_summary.csv"), portfolio_summary)

by_cardinality = combine(
    groupby(portfolio_summary, :intervention_count),
    nrow => :portfolio_count,
    :robust_feasible => sum => :robust_feasible_count,
    :robust_margin => maximum => :best_robust_margin,
)
sort!(by_cardinality, :intervention_count)
CSV.write(joinpath(OUT_DIR, "portfolio_by_cardinality.csv"), by_cardinality)

nonempty = filter(row -> row.intervention_count > 0 && row.robust_feasible, portfolio_summary)
sort!(nonempty, [:intervention_count, :robust_margin, :portfolio], rev = [false, true, false])
CSV.write(joinpath(OUT_DIR, "best_nonempty_portfolios.csv"), nonempty[1:min(10, nrow(nonempty)), :])

single = CSV.read(joinpath(RESULT_DIR, "single_replacement", "single_replacement.csv"), DataFrame)
static_nodes = CSV.read(joinpath(RESULT_DIR, "single_replacement", "static_nodes.csv"), DataFrame)
nodes = leftjoin(single, select(static_nodes, :bus, :static_weakness, :electrical_strength), on = :bus)
baseline_summary = TOML.parsefile(joinpath(RESULT_DIR, "baseline", "baseline.toml"))
baseline_margin = Float64(baseline_summary["dynamic_margin"])
nodes[!, :dynamic_margin_loss] = baseline_margin .- nodes.dynamic_margin
CSV.write(joinpath(OUT_DIR, "dynamic_node_scores.csv"), nodes)
node_rho = spearman(nodes.static_weakness, nodes.dynamic_margin_loss)

links = CSV.read(joinpath(RESULT_DIR, "link_outage", "link_outage_results.csv"), DataFrame)
link_ok = filter(row -> row.equilibrium_status == "ok", links)
CSV.write(joinpath(OUT_DIR, "dynamic_link_scores_qualified.csv"), link_ok)
link_rho = nrow(link_ok) >= 3 ? spearman(link_ok.static_weakness, link_ok.margin_loss) : NaN

scenario_summary = combine(
    groupby(portfolio, :scenario),
    :dynamic_margin => minimum => :minimum_margin,
    :dynamic_margin => maximum => :maximum_margin,
    :robust_feasible_for_scenario => sum => :feasible_rows,
)
CSV.write(joinpath(OUT_DIR, "scenario_summary.csv"), scenario_summary)

best_nonempty = nonempty[1, :]
summary = Dict(
    "timestamp_utc" => string(now(UTC)),
    "portfolio_rows" => nrow(portfolio),
    "portfolio_count" => nrow(portfolio_summary),
    "scenario_count" => length(unique(portfolio.scenario)),
    "portfolio_failed_rows" => sum(portfolio.equilibrium_status .!= "ok"),
    "portfolio_feasible_rows" => sum(portfolio.robust_feasible_for_scenario),
    "null_design_portfolio" => portfolio_summary.portfolio[1],
    "null_design_margin" => portfolio_summary.robust_margin[1],
    "best_nonempty_portfolio" => best_nonempty.portfolio,
    "best_nonempty_intervention_count" => best_nonempty.intervention_count,
    "best_nonempty_robust_margin" => best_nonempty.robust_margin,
    "single_replacement_node_spearman" => node_rho,
    "link_outage_rows" => nrow(links),
    "link_outage_qualified_rows" => nrow(link_ok),
    "link_outage_failed_rows" => nrow(links) - nrow(link_ok),
    "link_outage_spearman" => link_rho,
)
open(joinpath(OUT_DIR, "summary.toml"), "w") do io
    TOML.print(io, summary)
end

println("PD39 summary complete")
println("node_spearman=$(node_rho)")
println("link_spearman=$(link_rho)")
println("best_nonempty=$(best_nonempty.portfolio) margin=$(best_nonempty.robust_margin)")
