using CSV
using DataFrames
using Statistics

const ROOT = joinpath(@__DIR__, "..", "..")
const RESULTS = joinpath(ROOT, "results")
disc = CSV.read(joinpath(RESULTS, "pd39", "portfolio_campaign", "portfolio_scenario_results.csv"), DataFrame)
cands = CSV.read(joinpath(RESULTS, "pd39", "model_qualification", "candidate_table.csv"), DataFrame)

function group_worst(q)
    isempty(q) ? (alpha = NaN, margin = NaN, scenario = "") : begin
        r = q[argmax(Float64.(q.max_real)), :]
        (alpha = Float64(r.max_real), margin = Float64(r.dynamic_margin), scenario = r.scenario)
    end
end

v8 = collect(Int.(cands.bus))
candidate_rows = NamedTuple[]
for c in eachrow(cands)
    bus = Int(c.bus)
    s = filter(r -> r.portfolio == string(bus), disc)
    pred = filter(r -> count(==(';'), r.portfolio) == 6 && occursin(";$(bus);", ";" * r.portfolio * ";"), disc)
    g = group_worst(pred)
    push!(candidate_rows, (bus = bus, dispatch_p_mw = c.dispatch_p_mw, machine_rating_mva = c.machine_rating_mva,
        voltage_setpoint_pu = c.voltage_setpoint_pu,
        singleton_nominal_alpha = isempty(s) ? NaN : Float64(first(s).max_real),
        singleton_nominal_margin = isempty(s) ? NaN : Float64(first(s).dynamic_margin),
        singleton_m9 = isempty(s) ? NaN : minimum(Float64.(s.dynamic_margin)),
        worst_7of8_margin_when_retained = g.margin,
        worst_7of8_alpha_when_retained = g.alpha,
        worst_7of8_scenario_when_retained = g.scenario))
end
CSV.write(joinpath(RESULTS, "PD39_TABLE_T2_CANDIDATES.csv"), DataFrame(candidate_rows))

predrows = NamedTuple[]
for missing in v8
    retained = setdiff(v8, [missing])
    p = join(retained, ";")
    q = filter(r -> r.portfolio == p, disc)
    g = group_worst(q)
    mw = sum(Float64.(filter(r -> Int(r.bus) in retained, cands).dispatch_p_mw))
    push!(predrows, (missing_sg = missing, portfolio = p, cardinality = 7,
        converted_mw = mw, worst_alpha = g.alpha, worst_margin = g.margin,
        critical_scenario = g.scenario))
end
CSV.write(joinpath(RESULTS, "PD39_TABLE_T3_7OF8_V8.csv"), DataFrame(predrows))

if isfile(joinpath(RESULTS, "PD39_PLANNER_SELECTIONS.csv"))
    CSV.write(joinpath(RESULTS, "PD39_TABLE_T6_TRANSITION_PARETO.csv"),
        CSV.read(joinpath(RESULTS, "PD39_PLANNER_SELECTIONS.csv"), DataFrame))
end
println("tables complete")
