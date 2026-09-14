using CSV
using DataFrames
using Distributed

include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
using .PD39

const ROOT = joinpath(@__DIR__, "..", "..")
const RESULTS = joinpath(ROOT, "results")
const HOLDOUT = CSV.read(joinpath(RESULTS, "PD39_HOLDOUT_CONDITIONS.csv"), DataFrame)
const OUT_DIR = joinpath(RESULTS, "PD39_255PLUS1_HOLDOUT_CHUNKS")
mkpath(OUT_DIR)

desired_workers = try parse(Int, get(ENV, "PD39_WORKERS", "4")) catch; 4 end
desired_workers = max(1, desired_workers)
if nprocs() - 1 < desired_workers
    addprocs(desired_workers - (nprocs() - 1); exeflags = "--project=$(Base.active_project())")
end

const CONDITIONS = [NamedTuple(r) for r in eachrow(HOLDOUT)]
const PORTFOLIOS = candidate_portfolios()
const CHUNK_SIZE = try parse(Int, get(ENV, "PD39_HOLDOUT_CHUNK_SIZE", "8")) catch; 8 end
const CHUNKS = [PORTFOLIOS[i:min(i + CHUNK_SIZE - 1, end)] for i in 1:CHUNK_SIZE:length(PORTFOLIOS)]

@everywhere begin
    if !isdefined(Main, :PD39)
        include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
    end
    using .PD39
    using CSV
    using DataFrames
end

@everywhere function holdout_row(portfolio, c)
    b = [Float64(c[Symbol("delta_line_$(lpad(i, 2, '0'))")]) for i in 1:46]
    result = try
        nw = build_confirmatory_network(portfolio;
            controller_delta = (Float64(c.delta_pll), Float64(c.delta_xf), Float64(c.delta_cc)),
            load_delta = (Float64(c.delta_load_p), Float64(c.delta_load_q)),
            ibr_delta = Float64(c.delta_ibr_p), branch_delta = b, bounds = :primary)
        run_margin_case(nw)
    catch err
        (status = "failed", error_type = "build_exception", error_message = sprint(showerror, err),
         equilibrium_status = "equilibrium_failed", alpha = NaN, margin = NaN,
         stable = false, equilibrium_residual = NaN)
    end
    p = portfolio_string(portfolio)
    meta = candidate_table()
    return (portfolio = p, cardinality = length(portfolio), converted_mw = sum(meta[findall(in(portfolio), meta.bus), :dispatch_p_mw]),
        condition = String(c.condition), seed = Int(c.seed), status = result.status,
        error_type = result.error_type, error_message = result.error_message,
        equilibrium_status = result.equilibrium_status, alpha = Float64(result.alpha),
        margin = Float64(result.margin), stable = Bool(result.stable),
        robust = isfinite(Float64(result.margin)) && Float64(result.margin) >= 0.05,
        equilibrium_residual = Float64(result.equilibrium_residual))
end

@everywhere function run_holdout_chunk(portfolio_chunk, conditions, chunk_id, out_dir)
    rows = NamedTuple[]
    for p in portfolio_chunk, c in conditions
        push!(rows, holdout_row(p, c))
    end
    path = joinpath(out_dir, "chunk_$(chunk_id).csv")
    CSV.write(path, DataFrame(rows))
    return path
end

println("PD39 255+1 complete holdout census: workers=$(nprocs()-1), portfolios=$(length(PORTFOLIOS)), conditions=$(length(CONDITIONS)), chunks=$(length(CHUNKS))")
completed = Dict{Int,String}()
pending = Tuple{Int,Vector{Vector{Int}}}[]
for (id, chunk) in enumerate(CHUNKS)
    path = joinpath(OUT_DIR, "chunk_$(id).csv")
    if isfile(path)
        try
            if nrow(CSV.read(path, DataFrame)) == length(chunk) * length(CONDITIONS)
                completed[id] = path
                continue
            end
        catch
        end
    end
    push!(pending, (id, chunk))
end
println("checkpoints completed=$(length(completed)) pending=$(length(pending))")
paths = pmap(pending) do item
    id, chunk = item
    println("C chunk ", id, " portfolios=", length(chunk))
    run_holdout_chunk(chunk, CONDITIONS, id, OUT_DIR)
end
for ((id, _), path) in zip(pending, paths)
    completed[id] = path
end
all_paths = [completed[i] for i in 1:length(CHUNKS)]
census = vcat([CSV.read(p, DataFrame) for p in all_paths]...; cols = :union)
expected = length(PORTFOLIOS) * length(CONDITIONS)
nrow(census) == expected || error("expected $(expected) rows, got $(nrow(census))")
length(unique(string.(census.portfolio) .* "|" .* string.(census.condition))) == expected || error("portfolio-condition keys are not unique")
length(unique(census.portfolio)) == 256 || error("portfolio count mismatch")
length(unique(census.condition)) == 24 || error("condition count mismatch")
sort!(census, [:condition, :cardinality, :portfolio])
CSV.write(joinpath(RESULTS, "PD39_255PLUS1_HOLDOUT_CENSUS.csv"), census)

function buses(p)
    isempty(p) || p == "none" ? Int[] : parse.(Int, split(p, ";"))
end

function proper_subsets(p)
    bs = buses(p)
    out = String[]
    for i in eachindex(bs)
        sub = [bs[j] for j in eachindex(bs) if j != i]
        push!(out, portfolio_string(sub))
    end
    return out
end

function condition_summary(df)
    rows = NamedTuple[]
    for cond in sort(unique(String.(df.condition)))
        d = filter(r -> String(r.condition) == cond, eachrow(df))
        byp = Dict(String(r.portfolio) => r for r in d)
        h0 = String[]
        h005 = String[]
        for p in keys(byp)
            r = byp[p]
            proper = proper_subsets(p)
            all_stable_proper = all(haskey(byp, q) && String(byp[q].status) == "ok" && Float64(byp[q].alpha) < 0 for q in proper)
            all_robust_proper = all(haskey(byp, q) && String(byp[q].status) == "ok" && Float64(byp[q].alpha) <= -0.05 for q in proper)
            if String(r.status) == "ok" && Float64(r.alpha) >= 0 && all_stable_proper
                push!(h0, p)
            end
            if String(r.status) == "ok" && Float64(r.alpha) > -0.05 && all_robust_proper
                push!(h005, p)
            end
        end
        sort!(h0); sort!(h005)
        push!(rows, (condition = cond, attempts = nrow(d), ok = count(String.(d.status) .== "ok"),
            failed = count(String.(d.status) .!= "ok"), h0_count = length(h0),
            h0_portfolios = join(h0, "|"), h005_count = length(h005),
            h005_portfolios = join(h005, "|")))
    end
    out = DataFrame(rows)
    CSV.write(joinpath(RESULTS, "PD39_255PLUS1_HOLDOUT_CONDITION_SUMMARY.csv"), out)
    return out
end

summary = condition_summary(census)
println("C census complete rows=", nrow(census), " H0_total=", sum(summary.h0_count), " H0.05_total=", sum(summary.h005_count))
