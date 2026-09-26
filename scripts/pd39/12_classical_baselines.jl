using CSV
using DataFrames
using Graphs
using LinearAlgebra
using Statistics

include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
using .PD39

const ROOT = joinpath(@__DIR__, "..", "..")
const OUT = joinpath(ROOT, "results", "PD39_CLASSICAL_BASELINES.csv")
mkpath(dirname(OUT))

data = ieee39_data()
candidate = candidate_table()
machine = data.machine
bus = data.bus
portfolio_table = CSV.read(joinpath(ROOT, "results", "PD39_PORTFOLIO_BLOCKER_STRUCTURE.csv"), DataFrame)
disc = CSV.read(joinpath(ROOT, "results", "pd39", "portfolio_campaign", "portfolio_scenario_results.csv"), DataFrame)
static_nodes = static_weak_nodes(data)
static_by_bus = Dict(Int(r.bus) => r for r in eachrow(static_nodes))

function buses_of(p)
    p == "none" ? Int[] : parse.(Int, split(p, ";"))
end

# Effective resistance is calculated on the undirected reactance-strength graph.
L = zeros(Float64, nrow(bus), nrow(bus))
for r in eachrow(data.branch)
    i, j = Int(r.src_bus), Int(r.dst_bus)
    w = 1 / max(abs(Float64(r.X)), eps(Float64))
    L[i, i] += w; L[j, j] += w; L[i, j] -= w; L[j, i] -= w
end
Lplus = pinv(L)
effective_to_slack = [Lplus[i, i] + Lplus[SLACK_BUS, SLACK_BUS] - 2Lplus[i, SLACK_BUS] for i in 1:nrow(bus)]

all_machine_mw = Dict(Int(r.bus) => 100.0 * Float64(bus[findfirst(bus.bus .== r.bus), :P]) for r in eachrow(machine))
all_machine_mva = Dict(Int(r.bus) => Float64(r.Sn) for r in eachrow(machine))
all_machine_inertia = Dict(Int(r.bus) => Float64(r.H) * Float64(r.Sn) for r in eachrow(machine))
total_machine_mw = sum(values(all_machine_mw))
total_machine_mva = sum(values(all_machine_mva))
total_inertia = sum(values(all_machine_inertia))

rows = NamedTuple[]
for r in eachrow(portfolio_table)
    p = r.portfolio
    replaced = buses_of(p)
    ds = [candidate[candidate.bus .== b, :][1, :] for b in replaced]
    static_weak = [Float64(static_by_bus[b].static_weakness) for b in replaced]
    strengths = [Float64(static_by_bus[b].electrical_strength) for b in replaced]
    vset = [Float64(bus[findfirst(bus.bus .== b), :V]) for b in replaced]
    p_mw = sum((100.0 * Float64(bus[findfirst(bus.bus .== b), :P]) for b in replaced); init = 0.0)
    mva = sum((all_machine_mva[b] for b in replaced); init = 0.0)
    remain = setdiff(keys(all_machine_mw), replaced)
    nominal = only(filter(x -> x.portfolio == p && x.scenario == "nominal", disc))
    push!(rows, (portfolio = p, cardinality = length(replaced), converted_mw = p_mw,
        ibr_mva = mva, ibr_penetration_mw = p_mw / total_machine_mw,
        ibr_penetration_mva = mva / total_machine_mva,
        remaining_sg_mw = sum(all_machine_mw[b] for b in remain; init = 0.0),
        remaining_sg_mva = sum(all_machine_mva[b] for b in remain; init = 0.0),
        remaining_inertia_mva_s = sum(all_machine_inertia[b] for b in remain; init = 0.0),
        total_inertia_mva_s = total_inertia,
        min_candidate_static_strength = isempty(strengths) ? NaN : minimum(strengths),
        mean_candidate_static_weakness = isempty(static_weak) ? NaN : mean(static_weak),
        sum_candidate_static_weakness = sum(static_weak; init = 0.0),
        min_candidate_voltage_setpoint = isempty(vset) ? NaN : minimum(vset),
        mean_effective_resistance_to_slack = isempty(replaced) ? NaN : mean(effective_to_slack[replaced]),
        max_effective_resistance_to_slack = isempty(replaced) ? NaN : maximum(effective_to_slack[replaced]),
        nominal_alpha = nominal.max_real, m_9 = r.m_9,
        scr = NaN, gscr = NaN,
        scr_status = "unavailable: no frozen short-circuit/converter fault-current model",
        gscr_status = "unavailable: generalized-SCR assumptions/data not established"))
end
CSV.write(OUT, DataFrame(rows))
println("wrote ", OUT, " rows=", length(rows), " total_machine_mw=", total_machine_mw,
        " total_machine_mva=", total_machine_mva, " total_inertia_mva_s=", total_inertia)
