using CSV
using DataFrames

const ROOT = joinpath(@__DIR__, "..", "..")
const RESULTS = joinpath(ROOT, "results")

function best_row(q, strategy)
    isempty(q) && return nothing
    if strategy == "P0_penetration_only"
        return q[argmax(Float64.(q.converted_mw)), :]
    elseif strategy == "P1_nominal_small_signal"
        return q[argmax(Float64.(q.converted_mw)), :]
    elseif strategy == "P2_static_strength"
        # Frozen rule: maximize the best available candidate-bus strength;
        # break ties by converted MW. SCR/gSCR are unavailable and excluded.
        vals = Float64.(q.min_candidate_static_strength)
        vals[.!isfinite.(vals)] .= -Inf
        k = sortperm(1:nrow(q), by = i -> (-vals[i], -Float64(q.converted_mw[i])))[1]
        return q[k, :]
    elseif strategy == "P3_robust_dynamic"
        return q[argmax(Float64.(q.m_9)), :]
    end
    error("unknown strategy")
end

base = CSV.read(joinpath(RESULTS, "PD39_CLASSICAL_BASELINES.csv"), DataFrame)
disc = CSV.read(joinpath(RESULTS, "pd39", "summary", "portfolio_summary.csv"), DataFrame)
disc = select(disc, :portfolio, :robust_margin => :m_9, :robust_feasible, :all_equilibria_ok)
base = leftjoin(base, disc; on = :portfolio, makeunique = true)
base.m_9 = coalesce.(base.m_9, base.m_9_1)
base.robust_feasible = coalesce.(base.robust_feasible, false)
base.all_equilibria_ok = coalesce.(base.all_equilibria_ok, false)

targets = [("25pct", 1155.0), ("50pct", 2310.0), ("75pct", 3465.0), ("100pct", 4620.0)]
strategies = ["P0_penetration_only", "P1_nominal_small_signal", "P2_static_strength", "P3_robust_dynamic"]
rows = NamedTuple[]
for (target_id, target_mw) in targets, strategy in strategies
    feasible = filter(r -> Float64(r.converted_mw) >= target_mw - 1e-9 && r.all_equilibria_ok, base)
    if strategy == "P1_nominal_small_signal"
        feasible = filter(r -> Float64(r.nominal_alpha) <= -0.05, feasible)
    elseif strategy == "P3_robust_dynamic"
        feasible = filter(r -> r.robust_feasible, feasible)
    end
    r = best_row(feasible, strategy)
    if r === nothing
        push!(rows, (target = target_id, target_mw = target_mw, strategy = strategy,
            portfolio = "", converted_mw = NaN, cardinality = NaN, nominal_alpha = NaN,
            m_9 = NaN, static_strength = NaN, status = "no_feasible_discovery_design"))
    else
        push!(rows, (target = target_id, target_mw = target_mw, strategy = strategy,
            portfolio = r.portfolio, converted_mw = r.converted_mw, cardinality = r.cardinality,
            nominal_alpha = r.nominal_alpha, m_9 = r.m_9,
            static_strength = r.min_candidate_static_strength, status = "selected"))
    end
end
CSV.write(joinpath(RESULTS, "PD39_PLANNER_SELECTIONS.csv"), DataFrame(rows))
println("planner selections complete")
