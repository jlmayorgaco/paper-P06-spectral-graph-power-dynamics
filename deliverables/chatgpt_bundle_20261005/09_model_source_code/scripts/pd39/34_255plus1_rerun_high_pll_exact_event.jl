using CSV
using DataFrames

include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
using .PD39

const ROOT = joinpath(@__DIR__, "..", "..")
const RESULTS = joinpath(ROOT, "results")
const DISCOVERY = CSV.read(joinpath(RESULTS, "pd39", "portfolio_campaign", "portfolio_scenario_results.csv"), DataFrame)
const STATIC = CSV.read(joinpath(RESULTS, "PD39_CLASSICAL_BASELINES.csv"), DataFrame)
const AUDIT = CSV.read(joinpath(RESULTS, "PD39_255PLUS1_NUMERICAL_TRUTH_AUDIT.csv"), DataFrame)
const V8 = portfolio_string(CANDIDATE_SG_BUSES)
const SCENARIOS = uncertainty_scenarios()
const HIGH_PLL = filter(s -> s.pll_scale == 1.2, SCENARIOS)

parse_portfolio(s::AbstractString) = isempty(s) || s == "none" ? Int[] : parse.(Int, split(s, ";"))
scenario_controller(s) = (s.pll_scale - 1, s.filter_scale - 1, s.current_control_scale - 1)

function best_worst_7of8()
    d = combine(groupby(DISCOVERY, :portfolio), :dynamic_margin => minimum => :m9)
    d = innerjoin(d, STATIC[:, [:portfolio, :cardinality]], on = :portfolio)
    d = filter(r -> r.cardinality == 7, d)
    best = first(sort(collect(eachrow(d)), by = r -> (-r.m9, r.portfolio)))
    worst = first(sort(collect(eachrow(d)), by = r -> (r.m9, r.portfolio)))
    return (best = String(best.portfolio), worst = String(worst.portfolio))
end

function linear_alpha(p, sid)
    r = filter(r -> String(r.portfolio) == p && String(r.scenario) == sid && String(r.audit_layer) == "independent", eachrow(AUDIT))
    isempty(r) ? NaN : Float64(first(r).reference_alpha)
end

function trace_path(label, scenario)
    joinpath(RESULTS, "PD39_255PLUS1_TDS_TRACE_$(label)_$(scenario.id).csv")
end

bw = best_worst_7of8()
cases = [("V8", V8), ("best_7of8", bw.best), ("worst_7of8", bw.worst)]
rows = NamedTuple[]
println("Exact high-PLL TDS rerun: load_bus=$(largest_load_bus()), pulse_window=1.0-1.1 s")
for (label, p) in cases, s in HIGH_PLL
    println("B exact TDS ", label, " / ", s.id)
    flush(stdout)
    la = linear_alpha(p, s.id)
    try
        td = simulate_tds_case(parse_portfolio(p); pulse = 0.01,
            controller_delta = scenario_controller(s), bounds = :discovery,
            tspan = (0.0, 20.0), saveat = 0.01)
        tp = trace_path(label, s)
        CSV.write(tp, DataFrame(time = td.t, frequency_spread_hz = td.frequency_spread_hz,
            voltage_pu = td.voltage, voltage_deviation_pu = td.voltage_deviation_pu))
        push!(rows, (case_id = label, portfolio = p, scenario = s.id, pulse = 0.01,
            pulse_start_s = 1.0, pulse_end_s = 1.1, load_bus = largest_load_bus(),
            status = "ok", error = "", linear_alpha = la,
            estimated_rate_s_inv = td.estimated_rate_s_inv,
            sign_consistent = isfinite(la) && isfinite(td.estimated_rate_s_inv) && sign(la) == sign(td.estimated_rate_s_inv),
            max_frequency_spread_hz = td.max_frequency_spread_hz,
            max_voltage_deviation_pu = td.max_voltage_deviation_pu,
            frequency_settling_s = td.frequency_settling_s,
            voltage_settling_s = td.voltage_settling_s,
            voltage_min_pu = td.voltage_min_pu, voltage_max_pu = td.voltage_max_pu,
            trace = tp))
    catch err
        push!(rows, (case_id = label, portfolio = p, scenario = s.id, pulse = 0.01,
            pulse_start_s = 1.0, pulse_end_s = 1.1, load_bus = largest_load_bus(),
            status = "failed", error = sprint(showerror, err), linear_alpha = la,
            estimated_rate_s_inv = NaN, sign_consistent = false,
            max_frequency_spread_hz = NaN, max_voltage_deviation_pu = NaN,
            frequency_settling_s = NaN, voltage_settling_s = NaN,
            voltage_min_pu = NaN, voltage_max_pu = NaN, trace = ""))
    end
end
out = DataFrame(rows)
nrow(out) == 12 || error("expected exactly 12 TDS rows")
CSV.write(joinpath(RESULTS, "PD39_255PLUS1_HIGH_PLL_TDS.csv"), out)
println("exact high-PLL TDS complete rows=", nrow(out), " ok=", count(out.status .== "ok"), " sign_consistent=", count(out.sign_consistent))
