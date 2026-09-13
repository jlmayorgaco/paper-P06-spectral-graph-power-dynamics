using DataFrames
using NetworkDynamics

export ROBUST_MARGIN_TARGET, candidate_portfolios, evaluate_portfolio,
    evaluate_campaign, minimum_intervention_design

const ROBUST_MARGIN_TARGET = 0.05

"Enumerate all candidate subsets in increasing cardinality and bit-mask order."
function candidate_portfolios()
    portfolios = Vector{Vector{Int}}()
    for mask in 0:(2 ^ length(CANDIDATE_SG_BUSES) - 1)
        subset = Int[]
        for (bit, bus) in enumerate(CANDIDATE_SG_BUSES)
            ((mask >> (bit - 1)) & 1 == 1) && push!(subset, bus)
        end
        push!(portfolios, subset)
    end
    return portfolios
end

"Evaluate one SG→IBR portfolio under one fixed controller scenario."
function evaluate_portfolio(base_nw, portfolio; pfs = nothing, scenario = first(uncertainty_scenarios()),
                            margin_target = ROBUST_MARGIN_TARGET)
    template = simple_gfldc_template(
        pll_scale = scenario.pll_scale,
        filter_scale = scenario.filter_scale,
        current_control_scale = scenario.current_control_scale,
    )
    nw = replace_buses(base_nw, portfolio; template = template)
    eq = initialize_equilibrium(nw; pfs = pfs, sparse = false)
    stab = stability_audit(eq.state)
    return (
        portfolio = join(portfolio, ";"),
        replaced_buses = copy(portfolio),
        intervention_count = length(portfolio),
        scenario = scenario.id,
        powerflow_finite = eq.powerflow_finite,
        equilibrium_finite = eq.state_finite,
        fixed_point = eq.fixed_point,
        dynamic_margin = stab.dynamic_margin,
        max_real = stab.max_real,
        stable = stab.stable,
        robust_feasible_for_scenario = eq.state_finite && eq.fixed_point &&
            stab.stable && stab.dynamic_margin >= margin_target,
    )
end

"Evaluate the full preregistered portfolio × uncertainty grid."
function evaluate_campaign(base_nw; pfs = nothing, portfolios = candidate_portfolios(),
                           scenarios = uncertainty_scenarios(), margin_target = ROBUST_MARGIN_TARGET)
    rows = NamedTuple[]
    for portfolio in portfolios
        for scenario in scenarios
            push!(rows, evaluate_portfolio(
                base_nw,
                portfolio;
                pfs = pfs,
                scenario = scenario,
                margin_target = margin_target,
            ))
        end
    end
    return DataFrame(rows)
end

"Select the minimum-cardinality robustly feasible design."
function minimum_intervention_design(results; margin_target = ROBUST_MARGIN_TARGET)
    grouped = combine(
        groupby(results, :portfolio),
        :intervention_count => first => :intervention_count,
        :dynamic_margin => minimum => :robust_margin,
        :robust_feasible_for_scenario => all => :robust_feasible,
    )
    feasible = filter(row -> row.robust_feasible && row.robust_margin >= margin_target, grouped)
    isempty(feasible) && return nothing
    sort!(feasible, [:intervention_count, :robust_margin], rev = [false, true])
    return feasible[1, :]
end
