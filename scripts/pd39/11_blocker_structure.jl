using CSV
using DataFrames

const ROOT = joinpath(@__DIR__, "..", "..")
const DISC = CSV.read(joinpath(ROOT, "results", "pd39", "portfolio_campaign", "portfolio_scenario_results.csv"), DataFrame)
const CAND = CSV.read(joinpath(ROOT, "results", "pd39", "model_qualification", "candidate_table.csv"), DataFrame)
const MW = Dict(Int(r.bus) => Float64(r.dispatch_p_mw) for r in eachrow(CAND))
const OUT = joinpath(ROOT, "results", "PD39_PORTFOLIO_BLOCKER_STRUCTURE.csv")
mkpath(dirname(OUT))

portfolios = unique(DISC.portfolio)
portfolios = sort(collect(portfolios), by = p -> (p == "none" ? 0 : count(==(';'), p) + 1, p == "none" ? 0 : parse(Int, replace(p, ";" => ""))))

function buses_of(p)
    p == "none" ? Int[] : parse.(Int, split(p, ";"))
end

rows = NamedTuple[]
for p in portfolios
    s = filter(r -> r.portfolio == p, DISC)
    nominal = only(filter(r -> r.scenario == "nominal", s))
    m9 = minimum(s.dynamic_margin)
    alpha_nom = nominal.max_real
    pass_stability = all(s.stable)
    pass_robust = all(s.robust_feasible_for_scenario)
    buses = buses_of(p)
    proper = filter(q -> q != p && (q == "none" || length(buses_of(q)) < length(buses)), portfolios)
    proper_stable = all(all(filter(r -> r.portfolio == q, DISC).stable) for q in proper)
    proper_robust = all(all(filter(r -> r.portfolio == q, DISC).robust_feasible_for_scenario) for q in proper)
    push!(rows, (portfolio = p, cardinality = length(buses), converted_mw = sum((get(MW, b, 0.0) for b in buses); init = 0.0),
        nominal_alpha = alpha_nom, m_9 = m9,
        pass_stability = pass_stability, pass_0_05_robustness = pass_robust,
        minimal_true_blocker = !pass_stability && proper_stable,
        minimal_0_05_blocker = !pass_robust && proper_robust))
end
CSV.write(OUT, DataFrame(rows))
println("wrote ", OUT, " rows=", length(rows))
