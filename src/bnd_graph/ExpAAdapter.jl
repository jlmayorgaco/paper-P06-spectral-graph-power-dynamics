module ExpAAdapter

using LinearAlgebra
using TOML
using ..DynamicSelfEnergy: StateSpaceSelfEnergy, sigma, sigma_derivative

export BNDOperatorBundle, bundle_from_realization, save_expA_bundle,
       load_expA_bundle, required_metadata_keys, missing_metadata

const required_metadata_keys = ["source", "retained_state_names", "units",
    "sign_convention", "operating_point", "model_version", "exact_or_approximate",
    "valid_frequency_band"]

struct BNDOperatorBundle
    M::Matrix{Float64}
    L::Matrix{Float64}
    sigma::Union{Nothing,Function}
    sigma_derivative::Union{Nothing,Function}
    metadata::Dict{String,Any}
    sigma_available::Bool
    realization::Union{Nothing,StateSpaceSelfEnergy}
end

function BNDOperatorBundle(M::AbstractMatrix, L::AbstractMatrix,
                           sigma_fun, derivative_fun, metadata::AbstractDict;
                           sigma_available::Bool=(sigma_fun !== nothing),
                           realization=nothing)
    size(M, 1) == size(M, 2) || throw(ArgumentError("M must be square"))
    size(L) == size(M) || throw(ArgumentError("L and M must have equal dimensions"))
    sigma_available && sigma_fun === nothing &&
        throw(ArgumentError("sigma_available=true requires a callable Σ(s)"))
    md = Dict{String,Any}(string(k) => v for (k, v) in metadata)
    return BNDOperatorBundle(Matrix{Float64}(M), Matrix{Float64}(L),
        sigma_available ? sigma_fun : nothing,
        sigma_available ? derivative_fun : nothing,
        md, sigma_available, realization)
end

function bundle_from_realization(M, L, sys::StateSpaceSelfEnergy, metadata)
    return BNDOperatorBundle(M, L,
        s -> sigma(sys, s), s -> sigma_derivative(sys, s), metadata;
        sigma_available=true, realization=sys)
end

missing_metadata(bundle::BNDOperatorBundle) =
    filter(k -> !haskey(bundle.metadata, k), required_metadata_keys)

matrix_rows(A::AbstractMatrix) = [collect(A[i, :]) for i in axes(A, 1)]

function _toml_value(x)
    if x isa AbstractDict
        return Dict(string(k) => _toml_value(v) for (k, v) in x)
    elseif x isa Tuple
        return [_toml_value(v) for v in x]
    elseif x isa AbstractVector
        return [_toml_value(v) for v in x]
    elseif x isa Symbol
        return string(x)
    elseif x === nothing
        return ""
    else
        return x
    end
end

"""Write the versioned ExpA interchange format. Function closures are never serialized."""
function save_expA_bundle(path::AbstractString, bundle::BNDOperatorBundle)
    data = Dict{String,Any}(
        "schema_version" => "1.0",
        "sigma_available" => bundle.sigma_available,
        "M" => matrix_rows(bundle.M),
        "L" => matrix_rows(bundle.L),
        "metadata" => _toml_value(bundle.metadata))
    if bundle.sigma_available
        bundle.realization isa StateSpaceSelfEnergy ||
            throw(ArgumentError("the portable v1 format requires a state-space realization"))
        sys = bundle.realization
        data["realization"] = Dict{String,Any}(
            "D0" => matrix_rows(sys.D0),
            "Ac" => matrix_rows(sys.Ac),
            "B" => matrix_rows(sys.B),
            "C" => matrix_rows(sys.C),
            "channel_pairs" => [collect(p) for p in sys.channel_pairs],
            "offdiag_state_indices" => sys.offdiag_state_indices,
            "offdiag_base_scale" => sys.offdiag_base_scale)
    end
    open(path, "w") do io
        TOML.print(io, data)
    end
    return path
end

function _matrix_or_empty(raw, nrows::Int, ncols::Int)
    isempty(raw) && return zeros(Float64, nrows, ncols)
    A = reduce(vcat, (reshape(Float64.(row), 1, :) for row in raw))
    size(A) == (nrows, ncols) || throw(ArgumentError("serialized matrix has wrong dimensions"))
    return A
end

"""Load an ExpA bundle. If only T_exact(s) was exported, Σ remains unavailable."""
function load_expA_bundle(path::AbstractString)
    data = TOML.parsefile(path)
    get(data, "schema_version", "") == "1.0" ||
        throw(ArgumentError("unsupported ExpA bundle schema version"))
    Mraw = data["M"]
    n = length(Mraw)
    M = _matrix_or_empty(Mraw, n, n)
    L = _matrix_or_empty(data["L"], n, n)
    metadata = Dict{String,Any}(string(k) => v for (k, v) in get(data, "metadata", Dict()))
    available = Bool(get(data, "sigma_available", false))
    if !available
        return BNDOperatorBundle(M, L, nothing, nothing, metadata;
                                 sigma_available=false, realization=nothing)
    end
    raw = data["realization"]
    D0 = _matrix_or_empty(raw["D0"], n, n)
    Acraw = raw["Ac"]
    nc = length(Acraw)
    Ac = _matrix_or_empty(Acraw, nc, nc)
    Braw = raw["B"]
    B = isempty(Braw) ? zeros(Float64, nc, n) : _matrix_or_empty(Braw, nc, n)
    Craw = raw["C"]
    C = isempty(Craw) ? zeros(Float64, n, nc) : _matrix_or_empty(Craw, n, nc)
    pairs = [Tuple(Int.(p)) for p in get(raw, "channel_pairs", Vector{Vector{Int}}())]
    sys = StateSpaceSelfEnergy(D0, Ac, B, C;
        channel_pairs=pairs,
        offdiag_state_indices=Int.(get(raw, "offdiag_state_indices", Int[])),
        offdiag_base_scale=Float64(get(raw, "offdiag_base_scale", 1.0)))
    return bundle_from_realization(M, L, sys, metadata)
end

end
