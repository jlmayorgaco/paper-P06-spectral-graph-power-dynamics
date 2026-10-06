module PhysicalGraph

using LinearAlgebra
using CSV
using DataFrames

export assemble_ybus, kron_reduce, kron_identity_audit, passive_port_graph,
       spectral_summary, generator_ports

"""Assemble the conventional nodal admittance matrix from the frozen pi-line data.

The archived PowerDynamics pi-line reports currents with the opposite sign to
the conventional passive Ybus convention. This routine returns the conventional
matrix, so a passive branch satisfies I_drawn = Ybus * V. `r_src` is the source
side ideal-transformer ratio and the destination ratio is one in the frozen
case. Line shunts are attached to their named ends.
"""
function assemble_ybus(branches::AbstractDataFrame, nbus::Integer)
    nbus > 0 || throw(ArgumentError("nbus must be positive"))
    Y = zeros(ComplexF64, nbus, nbus)
    for row in eachrow(branches)
        src, dst = Int(row.src_bus), Int(row.dst_bus)
        1 <= src <= nbus && 1 <= dst <= nbus ||
            throw(ArgumentError("branch bus index outside 1:nbus"))
        src != dst || throw(ArgumentError("self-loop branch is unsupported"))
        z = complex(Float64(row.R), Float64(row.X))
        abs(z) > 0 || throw(ArgumentError("zero series impedance on $src-$dst"))
        tap = Float64(row.r_src)
        tap > 0 || throw(ArgumentError("nonpositive source tap on $src-$dst"))
        y = inv(z)
        ysrc = complex(Float64(row.G_src), Float64(row.B_src))
        ydst = complex(Float64(row.G_dst), Float64(row.B_dst))
        Y[src,src] += tap^2 * (y + ysrc)
        Y[src,dst] -= tap * y
        Y[dst,src] -= tap * y
        Y[dst,dst] += y + ydst
    end
    return Y
end

"""Kron-reduce a complex nodal operator onto the retained ports."""
function kron_reduce(Y::AbstractMatrix{<:Complex}, retained::AbstractVector{<:Integer})
    n, m = size(Y)
    n == m || throw(DimensionMismatch("Y must be square"))
    p = Int.(retained)
    length(unique(p)) == length(p) || throw(ArgumentError("duplicate retained port"))
    all(i -> 1 <= i <= n, p) || throw(ArgumentError("retained port out of range"))
    eliminated = setdiff(collect(1:n), p)
    Ypp = Y[p,p]
    isempty(eliminated) && return (Yport=Matrix(Ypp), retained=p,
        eliminated=eliminated, interior_condition=1.0)
    Yii = Y[eliminated,eliminated]
    κ = cond(Yii)
    isfinite(κ) || throw(LinearAlgebra.SingularException(0))
    Yport = Ypp - Y[p,eliminated] * (Yii \ Y[eliminated,p])
    return (Yport=Matrix(Yport), retained=p, eliminated=eliminated,
        interior_condition=κ)
end

"""Check port-current equivalence against a directly solved zero-interior-injection network."""
function kron_identity_audit(Y, reduction; trials=5)
    p, i = reduction.retained, reduction.eliminated
    isempty(i) && return (relative_residual=0.0, max_abs_residual=0.0,
        trials=trials)
    Vp = Matrix{ComplexF64}(undef, length(p), trials)
    for c in 1:trials, r in eachindex(p)
        Vp[r,c] = complex(sin(0.37*r + 0.61*c), cos(0.23*r - 0.47*c))
    end
    Vi = -(Y[i,i] \ (Y[i,p] * Vp))
    I_direct = Y[p,p] * Vp + Y[p,i] * Vi
    I_reduced = reduction.Yport * Vp
    err = norm(I_direct - I_reduced)
    scale = max(norm(I_direct), norm(I_reduced), eps(Float64))
    return (relative_residual=err/scale,
        max_abs_residual=maximum(abs, I_direct-I_reduced), trials=trials)
end

"""Build the preregistered passive generator-port conductance graph.

The Kron reduction is performed on the full complex passive branch Ybus first.
Taking its Hermitian part afterwards preserves the port dissipation operator;
it is not equivalent in general to Kron-reducing only the original real part.
"""
function passive_port_graph(Y, retained)
    red = kron_reduce(Y, retained)
    Yp = red.Yport
    Lc = (Yp + Yp') / 2
    residual = norm(Lc - Lc', Inf) / max(norm(Lc, Inf), eps(Float64))
    # The reduced passive-port Hermitian part is real symmetric up to roundoff.
    Lreal = Matrix{Float64}(real.((Lc + Lc')/2))
    return (Lc=Lreal, Yport=Yp, reduction=red,
        hermitian_residual=residual,
        imaginary_leakage=norm(imag.(Lc), Inf) / max(norm(Lc, Inf), eps(Float64)))
end

function spectral_summary(A; relative_tolerance=1e-10)
    H = Hermitian((A + A')/2)
    λ = eigvals(H)
    scale = max(maximum(abs, λ; init=0.0), eps(Float64))
    tol = relative_tolerance * scale
    λmin, λmax = minimum(λ), maximum(λ)
    nzero = count(x -> abs(x) <= tol, λ)
    positive = minimum(λ) > tol
    psd = minimum(λ) >= -tol
    κ = positive ? λmax/λmin : Inf
    return (lambda_min=λmin, lambda_max=λmax, condition_number=κ,
        zero_mode_count=nzero, psd=psd, spd=positive,
        symmetry_residual=norm(A-A',Inf)/max(norm(A,Inf),eps(Float64)),
        eigenvalues=λ, tolerance=tol)
end

generator_ports(buses::AbstractDataFrame) =
    Int.(buses.bus[Bool.(buses.has_gen)]) |> sort

end
