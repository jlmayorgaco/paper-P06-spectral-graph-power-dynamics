module AnalyticDeviceModel

using LinearAlgebra

export rank_audit, port_transfer, json_string, write_json

"Scale-aware numerical rank audit retaining the entire singular spectrum."
function rank_audit(A::AbstractMatrix; rtol=nothing)
    s = svdvals(Matrix{Float64}(A))
    scale = isempty(s) ? 0.0 : first(s)
    reltol = rtol === nothing ? max(size(A)...) * eps(Float64) : Float64(rtol)
    tol = reltol * scale
    return (rank=count(>(tol), s), tolerance=tol, relative_tolerance=reltol,
            singular_values=s, residual=isempty(s) ? 0.0 : norm(s[(count(>(tol),s)+1):end]))
end

"Return the small-signal terminal admittance/current map for a state-space device."
function port_transfer(A, B, C, D, s::Number)
    return D + C * ((s * I(size(A,1)) - A) \ B)
end

json_escape(s::AbstractString) = replace(String(s), "\\" => "\\\\", "\"" => "\\\"",
                                         "\n" => "\\n", "\r" => "\\r", "\t" => "\\t")
function json_value(x)
    x === nothing && return "null"
    x === missing && return "null"
    x isa Bool && return x ? "true" : "false"
    x isa AbstractString && return "\"$(json_escape(x))\""
    x isa Symbol && return json_value(String(x))
    x isa Integer && return string(x)
    x isa AbstractFloat && return isfinite(x) ? repr(Float64(x)) : "null"
    x isa Complex && return json_value(Dict("real"=>real(x),"imag"=>imag(x)))
    x isa NamedTuple && return json_value(Dict(string(k)=>v for (k,v) in pairs(x)))
    x isa AbstractDict && return "{" * join(("\"$(json_escape(string(k)))\":" * json_value(v)
        for (k,v) in sort!(collect(pairs(x)); by=p->string(first(p)))), ",") * "}"
    x isa Tuple && return json_value(collect(x))
    x isa AbstractArray && return "[" * join(json_value.(collect(x)), ",") * "]"
    return json_value(string(x))
end

function json_string(x; indent=2)
    # Compact JSON is deterministic and sufficient for machine-readable results.
    return json_value(x) * "\n"
end

write_json(path, x) = open(path,"w") do io; write(io,json_string(x)); end

end
