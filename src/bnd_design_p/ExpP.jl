module ExpP

using LinearAlgebra, SHA, TOML, CSV, DataFrames
include(joinpath(@__DIR__, "..", "bnd_model_expN", "PDExactDesignN.jl"))
using .PDExactDesignN

export PDExactDesignN, ROOT, read_candidate, verify_expN_freeze,
       reset_counters!, counter_snapshot, design_descriptor, design_spectrum,
       port_transfer, nodal_operator, linear_assignment, write_json

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const _COUNTERS = Dict{Symbol,Int}(:descriptor => 0, :spectrum => 0,
                                   :mode_sensitivity => 0)

function reset_counters!()
    for key in keys(_COUNTERS)
        _COUNTERS[key] = 0
    end
    nothing
end

counter_snapshot() = copy(_COUNTERS)

function design_descriptor(ctx, rho, kp, ki)
    _COUNTERS[:descriptor] += 1
    PDExactDesignN.descriptor(ctx, rho, kp, ki)
end

function design_spectrum(ctx, rho, kp, ki)
    _COUNTERS[:spectrum] += 1
    PDExactDesignN.spectrum(ctx, rho, kp, ki)
end

function read_candidate(root::AbstractString=ROOT)
    p = joinpath(root, "reports", "experiment_N", "Z_N_NOMINAL_FINAL.toml")
    isfile(p) && isfile(p * ".sha256") || error("ExpN frozen candidate missing")
    digest = bytes2hex(sha256(read(p)))
    digest == strip(read(p * ".sha256", String)) || error("ExpN candidate SHA mismatch")
    candidate = TOML.parsefile(p)
    candidate["model_sha"] == verify_expN_freeze(root) || error("candidate/model freeze mismatch")
    candidate, digest
end

"""Check every source/input hash recorded in the immutable ExpN freeze."""
function verify_expN_freeze(root::AbstractString=ROOT)
    freeze_path = joinpath(root, "reports", "experiment_N", "MODEL_FREEZE.json")
    text = read(freeze_path, String)
    m = match(r"\"MODEL_SHA\"\s*:\s*\"([0-9a-f]{64})\"", text)
    m === nothing && error("MODEL_SHA absent from ExpN freeze")
    pairs = collect(eachmatch(r"\"((?:src|reports)/[^\"]+)\"\s*:\s*\"([0-9a-f]{64})\"", text))
    isempty(pairs) && error("no inputs found in ExpN freeze")
    checked = 0
    for pair in pairs
        rel, expected = pair.captures
        path = joinpath(root, split(rel, '/')...)
        isfile(path) || error("ExpN frozen input missing: $rel")
        bytes2hex(sha256(read(path))) == expected || error("ExpN frozen input changed: $rel")
        checked += 1
    end
    checked == 17 || error("ExpN input count differs from freeze: $checked")
    m.captures[1]
end

port_transfer(J, s) = J.D + J.C * ((s * I(size(J.A, 1)) - J.A) \ J.B)

"""Return the ExpN raw-port closure T(s)=Ystatic+Σ share_i Π_i Yraw_i(s) Π_iᵀ.

The component output is current entering the device. In the paper convention
current injection into the network is Yinj=-Yraw, hence this is exactly
Ystatic-ΣΠYinjΠᵀ.
"""
function nodal_operator(ctx, rho, kp, ki, s)
    model = design_descriptor(ctx, rho, kp, ki)
    T = ComplexF64.(ctx.net.y_static)
    for b in model.blocks
        ix = (2b.bus - 1):(2b.bus)
        T[ix, ix] .+= b.share .* port_transfer(b.J, s)
    end
    T
end

"""Square minimum-cost assignment using the O(n³) shortest augmenting path method."""
function linear_assignment(cost::AbstractMatrix{<:Real})
    n, m = size(cost)
    n == m || throw(DimensionMismatch("assignment matrix must be square"))
    u = zeros(Float64, n + 1)
    v = zeros(Float64, n + 1)
    p = zeros(Int, n + 1)
    way = zeros(Int, n + 1)
    for i in 1:n
        p[1] = i
        j0 = 1
        minv = fill(Inf, n + 1)
        used = falses(n + 1)
        while true
            used[j0] = true
            i0 = p[j0]
            delta = Inf
            j1 = 0
            for j in 2:(n + 1)
                if !used[j]
                    cur = Float64(cost[i0, j - 1]) - u[i0 + 1] - v[j]
                    if cur < minv[j]
                        minv[j] = cur
                        way[j] = j0
                    end
                    if minv[j] < delta
                        delta = minv[j]
                        j1 = j
                    end
                end
            end
            for j in 1:(n + 1)
                if used[j]
                    p[j] != 0 && (u[p[j] + 1] += delta)
                    v[j] -= delta
                else
                    minv[j] -= delta
                end
            end
            j0 = j1
            p[j0] == 0 && break
        end
        while true
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
            j0 == 1 && break
        end
    end
    rows = zeros(Int, n)
    cols = zeros(Int, n)
    for j in 2:(n + 1)
        rows[p[j]] = p[j]
        cols[p[j]] = j - 1
    end
    rows, cols
end

json_escape(s::AbstractString) = "\"" * replace(String(s), "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n", "\r" => "\\r", "\t" => "\\t") * "\""
function json_string(x)
    if x === nothing
        "null"
    elseif x isa Bool
        x ? "true" : "false"
    elseif x isa AbstractString
        json_escape(x)
    elseif x isa Integer
        string(x)
    elseif x isa AbstractFloat
        isfinite(x) ? repr(x) : "null"
    elseif x isa AbstractVector || x isa Tuple
        "[" * join(json_string.(collect(x)), ",") * "]"
    elseif x isa AbstractDict
        ks = sort!(collect(keys(x)); by=string)
        "{" * join((json_escape(string(k)) * ":" * json_string(x[k]) for k in ks), ",") * "}"
    else
        json_escape(string(x))
    end
end

function write_json(path::AbstractString, data)
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, json_string(data), "\n")
    end
    path
end

end
