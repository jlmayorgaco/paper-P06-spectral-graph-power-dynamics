using CSV
using DataFrames
using Dates
using Distributed
using PowerDynamics
using TOML

include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
using .PD39

const OUT_DIR = get(ENV, "PD39_OUT_DIR", joinpath(@__DIR__, "..", "..", "results", "pd39", "portfolio_campaign"))
mkpath(OUT_DIR)

desired_workers = try
    parse(Int, get(ENV, "PD39_WORKERS", "4"))
catch
    4
end
desired_workers = max(1, desired_workers)
current_workers = max(0, nprocs() - 1)
if current_workers < desired_workers
    addprocs(desired_workers - current_workers; exeflags = "--project=$(Base.active_project())")
end

@everywhere begin
    if !isdefined(Main, :PD39)
        include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
    end
    using .PD39
    using CSV
    using DataFrames
end

@everywhere const PD39_WORKER_CACHE = Ref{Any}(nothing)

@everywhere function run_pd39_chunk(portfolio_chunk, scenarios, chunk_id, out_dir)
    if isnothing(PD39_WORKER_CACHE[])
        base_nw = baseline_network()
        base_eq = initialize_equilibrium(base_nw; sparse = false)
        base_eq.powerflow_finite && base_eq.state_finite && base_eq.fixed_point ||
            error("baseline qualification failed on worker")
        PD39_WORKER_CACHE[] = (base_nw, base_eq.pfs)
    end
    base_nw, pfs = PD39_WORKER_CACHE[]
    rows = NamedTuple[]
    for portfolio in portfolio_chunk
        for scenario in scenarios
            push!(rows, evaluate_portfolio(
                base_nw,
                portfolio;
                pfs = pfs,
                scenario = scenario,
                margin_target = ROBUST_MARGIN_TARGET,
            ))
        end
    end
    path = joinpath(out_dir, "chunk_$(chunk_id).csv")
    CSV.write(path, DataFrame(rows))
    return path
end

portfolios = candidate_portfolios()
scenarios = uncertainty_scenarios()
scenario_limit = try parse(Int, get(ENV, "PD39_SCENARIO_COUNT", "0")) catch; 0 end
scenario_limit > 0 && (scenarios = scenarios[1:min(scenario_limit, length(scenarios))])
chunk_size = try parse(Int, get(ENV, "PD39_CHUNK_SIZE", "16")) catch; 16 end
chunk_size = max(1, chunk_size)
chunks = [portfolios[i:min(i + chunk_size - 1, end)] for i in 1:chunk_size:length(portfolios)]
chunk_limit = try parse(Int, get(ENV, "PD39_MAX_CHUNKS", "0")) catch; 0 end
chunk_limit > 0 && (chunks = chunks[1:min(chunk_limit, length(chunks))])

println("PD39 parallel campaign: workers=$(max(0, nprocs() - 1)), chunks=$(length(chunks)), scenarios=$(length(scenarios))")
expected_rows = length(scenarios)
completed = Dict{Int, String}()
pending = Tuple{Int, Vector{Vector{Int}}}[]
for (chunk_id, chunk) in enumerate(chunks)
    path = joinpath(OUT_DIR, "chunk_$(chunk_id).csv")
    if isfile(path)
        try
            if nrow(CSV.read(path, DataFrame)) == length(chunk) * expected_rows
                completed[chunk_id] = path
                continue
            end
        catch
            # Treat an unreadable or partial checkpoint as pending.
        end
    end
    push!(pending, (chunk_id, chunk))
end
println("PD39 checkpoints: completed=$(length(completed)), pending=$(length(pending))")
new_paths = pmap(pending) do item
    chunk_id, chunk = item
    run_pd39_chunk(chunk, scenarios, chunk_id, OUT_DIR)
end
for ((chunk_id, _), path) in zip(pending, new_paths)
    completed[chunk_id] = path
end
paths = [completed[chunk_id] for chunk_id in 1:length(chunks)]

frames = [CSV.read(path, DataFrame) for path in paths]
results = vcat(frames...; cols = :union)
CSV.write(joinpath(OUT_DIR, "portfolio_scenario_results.csv"), results)
CSV.write(joinpath(OUT_DIR, "static_nodes.csv"), static_weak_nodes())
CSV.write(joinpath(OUT_DIR, "static_links.csv"), static_weak_links())

design = minimum_intervention_design(results; margin_target = ROBUST_MARGIN_TARGET)
open(joinpath(OUT_DIR, "design_summary.toml"), "w") do io
    TOML.print(io, Dict(
        "timestamp_utc" => string(now(UTC)),
        "portfolio_count" => sum(length, chunks),
        "scenario_count" => length(scenarios),
        "worker_count" => max(0, nprocs() - 1),
        "chunk_count" => length(chunks),
        "margin_target" => ROBUST_MARGIN_TARGET,
        "attempted_rows" => nrow(results),
        "failed_rows" => count(results.equilibrium_status .== "failed"),
        "robust_design_found" => !isnothing(design),
        "selected_portfolio" => isnothing(design) ? "" : design.portfolio,
        "selected_intervention_count" => isnothing(design) ? -1 : design.intervention_count,
        "selected_robust_margin" => isnothing(design) ? NaN : design.robust_margin,
    ))
end

println("PD39 portfolio campaign complete: rows=$(nrow(results)), robust_design_found=$( !isnothing(design) )")
