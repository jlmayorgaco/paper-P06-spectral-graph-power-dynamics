using CSV
using DataFrames
using Random

include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
using .PD39

const ROOT = joinpath(@__DIR__, "..", "..")
const DISC = CSV.read(joinpath(ROOT, "results", "pd39", "portfolio_campaign", "portfolio_scenario_results.csv"), DataFrame)
const BLOCK = CSV.read(joinpath(ROOT, "results", "PD39_PORTFOLIO_BLOCKER_STRUCTURE.csv"), DataFrame)
const CAND = CSV.read(joinpath(ROOT, "results", "pd39", "model_qualification", "candidate_table.csv"), DataFrame)
const MW = Dict(Int(r.bus) => Float64(r.dispatch_p_mw) for r in eachrow(CAND))
const PCAND = sum(values(MW))

const HOLDOUT = joinpath(ROOT, "results", "PD39_HOLDOUT_CONDITIONS.csv")
const SELECTION = joinpath(ROOT, "results", "PD39_SELECTION_SET.csv")
mkpath(dirname(HOLDOUT))

# Exact frozen 52-coordinate box: 3 control + 3 operating + 46 branch coordinates.
x = maximin_lhs(24, 52; seed = 39024, candidates = 256)
rows = NamedTuple[]
for i in 1:size(x, 1)
    c = 0.20 .* x[i, 1:3] .- 0.10
    o = 0.10 .* x[i, 4:6] .- 0.05
    b = 0.20 .* x[i, 7:end] .- 0.10
    nt = (condition = "C$(lpad(i, 2, '0'))", seed = 39024,
          delta_pll = c[1], delta_xf = c[2], delta_cc = c[3],
          delta_load_p = o[1], delta_load_q = o[2], delta_ibr_p = o[3])
    names = Symbol.("delta_line_" .* lpad.(string.(1:46), 2, '0'))
    push!(rows, merge(nt, NamedTuple{Tuple(names)}(Tuple(b))))
end
CSV.write(HOLDOUT, DataFrame(rows))

function buses_of(p)
    p == "none" ? Int[] : parse.(Int, split(p, ";"))
end

function bitkey(p)
    p == "none" ? 0 : parse(Int, replace(p, ";" => ""))
end

function choose_rep(target_pct)
    target = target_pct * PCAND
    d = filter(r -> r.cardinality < 8 && r.pass_stability, BLOCK)
    sort(d, [:converted_mw, :m_9, :portfolio])
    score(r) = (abs(r.converted_mw - target), -r.m_9, bitkey(r.portfolio))
    return first(sort(collect(eachrow(d)), by = score))
end

selected = Dict{String,NamedTuple}()
function add_selection!(id, p, reason)
    r = only(filter(x -> x.portfolio == p, BLOCK))
    selected[id] = (selection_id = id, portfolio = p, cardinality = r.cardinality,
        converted_mw = r.converted_mw, nominal_alpha = r.nominal_alpha, m_9 = r.m_9,
        selection_reason = reason)
end

for bus in CAND.bus
    p = join([b for b in CAND.bus if b != bus], ";")
    add_selection!("missing_$(bus)", p, "all eight 7-of-8 predecessors")
end
add_selection!("V8", "30;32;33;34;35;36;37;38", "full transition fixed placement")
add_selection!("all_SG", "none", "all-SG reference")
add_selection!("bus37_singleton", "37", "discovery best nonempty singleton")
for pct in (0.25, 0.50, 0.75)
    r = choose_rep(pct)
    add_selection!("representative_$(Int(round(100pct)))pct", r.portfolio,
                   "discovery-only nearest active-MW target, then max m_9")
end

# At most four closest nominal-alpha pairs per cardinality, selected only from
# discovery values; this bounds the confirmatory representative set without
# selecting on any new radius outcome.
pair_rows = NamedTuple[]
for k in 1:7
    d = filter(r -> r.cardinality == k && r.pass_stability, BLOCK)
    pairs = NamedTuple[]
    for i in 1:nrow(d)-1, j in i+1:nrow(d)
        a, b = d[i, :], d[j, :]
        mwgap = abs(a.converted_mw - b.converted_mw)
        adiff = abs(a.nominal_alpha - b.nominal_alpha)
        if mwgap <= 0.10 * PCAND && adiff <= 0.02
            push!(pairs, (cardinality = k, portfolio_a = a.portfolio, portfolio_b = b.portfolio,
                          mw_difference = mwgap, alpha_difference = adiff))
        end
    end
    sort!(pairs, by = r -> (r.alpha_difference, r.mw_difference, bitkey(r.portfolio_a), bitkey(r.portfolio_b)))
    for (j, pair) in enumerate(pairs[1:min(4, length(pairs))])
        pair_id = "K$(k)_P$(j)"
        push!(pair_rows, merge(pair, (pair_id = pair_id,)))
        add_selection!("pair_K$(k)_P$(j)_A", pair.portfolio_a, "matched-alpha pair $(pair_id) A")
        add_selection!("pair_K$(k)_P$(j)_B", pair.portfolio_b, "matched-alpha pair $(pair_id) B")
    end
end

sel = DataFrame(collect(values(selected)))
CSV.write(SELECTION, sel)
CSV.write(joinpath(ROOT, "results", "PD39_MATCHED_ALPHA_PAIRS.csv"), DataFrame(pair_rows))
println("wrote ", HOLDOUT, " rows=", nrow(DataFrame(rows)))
println("wrote ", SELECTION, " rows=", nrow(sel), " pairs=", length(pair_rows))
println("candidate_MW=", PCAND)
