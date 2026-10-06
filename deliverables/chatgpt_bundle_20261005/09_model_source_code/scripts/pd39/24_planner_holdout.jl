using CSV
using DataFrames
using Statistics

include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
using .PD39

const ROOT = joinpath(@__DIR__, "..", "..")
const RESULTS = joinpath(ROOT, "results")
portfolio_from_string(s) = s == "" || s == "none" ? Int[] : parse.(Int, split(s, ";"))
function line_delta(r)
    [Float64(r[Symbol("delta_line_$(lpad(i, 2, '0'))")]) for i in 1:46]
end

plans = CSV.read(joinpath(RESULTS, "PD39_PLANNER_SELECTIONS.csv"), DataFrame)
plans = filter(r -> r.status == "selected", plans)
conditions = CSV.read(joinpath(RESULTS, "PD39_HOLDOUT_CONDITIONS.csv"), DataFrame)
unique_portfolios = unique(plans.portfolio)
detail = NamedTuple[]
for portfolio_text in unique_portfolios
    portfolio = portfolio_from_string(portfolio_text)
    for cr in eachrow(conditions)
        nw = build_confirmatory_network(portfolio;
            controller_delta = (Float64(cr.delta_pll), Float64(cr.delta_xf), Float64(cr.delta_cc)),
            load_delta = (Float64(cr.delta_load_p), Float64(cr.delta_load_q)),
            ibr_delta = Float64(cr.delta_ibr_p), branch_delta = line_delta(cr), bounds = :primary)
        r = try
            run_margin_case(nw)
        catch err
            (status = "failed", error_type = "exception", error_message = sprint(showerror, err),
             alpha = NaN, margin = NaN, stable = false, equilibrium_status = "equilibrium_failed",
             equilibrium_residual = NaN)
        end
        push!(detail, (portfolio = portfolio_text, condition = cr.condition,
            status = r.status, error_type = r.error_type, error_message = r.error_message,
            alpha = r.alpha, margin = r.margin, stable = r.stable,
            robust = isfinite(r.margin) && r.margin >= 0.05,
            equilibrium_status = r.equilibrium_status, equilibrium_residual = r.equilibrium_residual))
    end
end
dd = DataFrame(detail)
CSV.write(joinpath(RESULTS, "PD39_PLANNER_HOLDOUT.csv"), dd)
sumrows = NamedTuple[]
for p in unique(dd.portfolio)
    q = filter(r -> r.portfolio == p, dd)
    a = filter(isfinite, Float64.(q.alpha))
    push!(sumrows, (portfolio = p, conditions = nrow(q), ok_cases = sum(q.status .== "ok"),
        robust_cases = sum(q.robust), robust_coverage = mean(q.robust),
        worst_alpha = isempty(a) ? NaN : maximum(a), median_alpha = isempty(a) ? NaN : median(a),
        worst_margin = isempty(a) ? NaN : -maximum(a),
        failed_cases = sum(q.status .!= "ok")))
end
summary = DataFrame(sumrows)
CSV.write(joinpath(RESULTS, "PD39_PLANNER_HOLDOUT_SUMMARY.csv"), summary)
CSV.write(joinpath(RESULTS, "PD39_PLANNER_HOLDOUT_SELECTIONS.csv"), plans)
println("planner holdout complete")
