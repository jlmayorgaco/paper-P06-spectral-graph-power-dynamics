using CSV
using DataFrames
using LinearAlgebra
using Statistics

include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
using .PD39

const ROOT = joinpath(@__DIR__, "..", "..")
const RESULTS = joinpath(ROOT, "results")
const SCALE = vcat(fill(0.10, 3), fill(0.05, 3), fill(0.10, 46))
const TARGETS = [("rho_0", 0.0), ("rho_0.05", -0.05)]
const N_DIRECTIONS = parse(Int, get(ENV, "PD39_RADIUS_DIRECTIONS", "128"))
const N_BISECT = parse(Int, get(ENV, "PD39_RADIUS_BISECT", "10"))
const SEED = parse(Int, get(ENV, "PD39_RADIUS_SEED", "39025"))

function selected_portfolios()
    d = CSV.read(joinpath(RESULTS, "PD39_SELECTION_SET.csv"), DataFrame)
    ids = split(get(ENV, "PD39_RADIUS_IDS", "V8,missing_30,missing_34,missing_35,missing_36,missing_37,missing_38,representative_25pct,representative_50pct,representative_75pct,bus37_singleton,all_SG"), ',')
    rows = NamedTuple[]
    for id in ids
        q = filter(r -> r.selection_id == id, d)
        isempty(q) && continue
        r = first(q)
        portfolio = r.portfolio == "none" ? Int[] : parse.(Int, split(r.portfolio, ";"))
        push!(rows, (selection_id = id, portfolio = portfolio, cardinality = r.cardinality,
            converted_mw = r.converted_mw))
    end
    return rows
end

function perturbation_parts(x, t)
    z = t .* Float64.(x) .* SCALE
    return z[1:3], z[4:5], z[6], z[7:end]
end

function evaluate(portfolio, x, t)
    c, load, ibr, lines = perturbation_parts(x, t)
    nw = build_confirmatory_network(portfolio; controller_delta = c,
        load_delta = load, ibr_delta = ibr, branch_delta = lines, bounds = :primary)
    return run_margin_case(nw)
end

function one_direction(portfolio, x, target_name, target; direction)
    base = evaluate(portfolio, x, 0.0)
    if base.status != "ok" || !isfinite(base.alpha)
        return (status = "baseline_not_estimable", radius_inf = NaN, radius_l2 = NaN,
            t_boundary = NaN, feasible_alpha = NaN, violating_alpha = NaN,
            boundary_alpha = NaN, equilibrium_residual = NaN, direction = direction,
            delta_star = "", physical_delta_star = "", active_coordinates = "", target = target_name)
    end
    base.alpha >= target && return (status = "already_violating_at_baseline", radius_inf = 0.0,
        radius_l2 = 0.0, t_boundary = 0.0, feasible_alpha = base.alpha,
        violating_alpha = base.alpha, boundary_alpha = base.alpha,
        equilibrium_residual = base.equilibrium_residual, direction = direction,
        delta_star = join(string.(zeros(52)), ";"), physical_delta_star = join(string.(zeros(52)), ";"),
        active_coordinates = "baseline",
        target = target_name)
    endpoint = evaluate(portfolio, x, 1.0)
    if endpoint.status != "ok" || !isfinite(endpoint.alpha) || endpoint.alpha < target
        return (status = endpoint.status == "ok" ? "no_boundary_in_primary_box" : "endpoint_not_estimable",
            radius_inf = NaN, radius_l2 = NaN, t_boundary = NaN,
            feasible_alpha = base.alpha, violating_alpha = endpoint.alpha,
            boundary_alpha = NaN, equilibrium_residual = endpoint.equilibrium_residual,
            direction = direction, delta_star = "", physical_delta_star = "",
            active_coordinates = "", target = target_name)
    end
    lo = 0.0
    hi_t = 1.0
    feasible = base
    violating = endpoint
    for _ in 1:N_BISECT
        mid = (lo + hi_t) / 2
        trial = evaluate(portfolio, x, mid)
        if trial.status == "ok" && isfinite(trial.alpha) && trial.alpha < target
            lo = mid
            feasible = trial
        else
            hi_t = mid
            violating = trial
        end
    end
    z = lo .* Float64.(x)
    physical_z = z .* SCALE
    act = sortperm(abs.(z); rev = true)[1:min(5, length(z))]
    return (status = "boundary_found", radius_inf = lo * maximum(abs.(x)),
        radius_l2 = lo * norm(x), t_boundary = lo,
        feasible_alpha = feasible.alpha, violating_alpha = violating.alpha,
        boundary_alpha = (feasible.alpha + violating.alpha) / 2,
        equilibrium_residual = feasible.equilibrium_residual, direction = direction,
        delta_star = join(string.(z), ";"), physical_delta_star = join(string.(physical_z), ";"),
        active_coordinates = join(string.(act), ";"),
        target = target_name)
end

function run_portfolio(row, dirs)
    out = NamedTuple[]
    for (j, x) in enumerate(eachrow(dirs))
        direction = Float64.(collect(x)) .* 2 .- 1
        for (target_name, target) in TARGETS
            r = one_direction(row.portfolio, direction, target_name, target; direction = j)
            push!(out, merge((selection_id = row.selection_id, portfolio = portfolio_string(row.portfolio),
                cardinality = row.cardinality, converted_mw = row.converted_mw), r))
        end
        if j % 8 == 0
            println(row.selection_id, " direction ", j, "/", size(dirs, 1))
        end
    end
    return out
end

dirs = maximin_lhs(N_DIRECTIONS, 52; seed = SEED, candidates = 128)
out = NamedTuple[]
for row in selected_portfolios()
    append!(out, run_portfolio(row, dirs))
    CSV.write(joinpath(RESULTS, "PD39_STRUCTURED_RADIUS_SEARCH.csv"), DataFrame(out))
end

if !isempty(out)
    d = DataFrame(out)
    rows = NamedTuple[]
    for id in unique(d.selection_id), target in unique(d.target)
        q = filter(r -> r.selection_id == id && r.target == target && r.status == "boundary_found", d)
        if isempty(q)
            push!(rows, (selection_id = id, target = target, directions = N_DIRECTIONS,
                boundaries_found = 0, smallest_boundary_found = NaN, l2_at_smallest = NaN,
                best_direction = NaN, feasible_alpha = NaN, violating_alpha = NaN))
        else
            k = argmin(Float64.(q.radius_inf))
            r = q[k]
            push!(rows, (selection_id = id, target = target, directions = N_DIRECTIONS,
                boundaries_found = nrow(q), smallest_boundary_found = r.radius_inf,
                l2_at_smallest = r.radius_l2, best_direction = r.direction,
                feasible_alpha = r.feasible_alpha, violating_alpha = r.violating_alpha))
        end
    end
    CSV.write(joinpath(RESULTS, "PD39_STRUCTURED_RADIUS_SUMMARY.csv"), DataFrame(rows))
end
println("structured-radius search complete")
