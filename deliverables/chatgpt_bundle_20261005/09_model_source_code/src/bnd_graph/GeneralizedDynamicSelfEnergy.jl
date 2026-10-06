module GeneralizedDynamicSelfEnergy

using LinearAlgebra

export AngleOperator, pi_q, pi_q_derivative, psi, psi_derivative,
       retained_operator, generalized_reconstruction_error,
       derivative_check

"""Exact angle-level operator exported by the Experiment-A partition.

The fields are the actual blocks of the algebraically reduced dynamic matrix.
`M` is the second-order metric obtained from `q̇ = Ckin v`; it is not presumed
to be a physical inertia matrix.
"""
struct AngleOperator
    M::Matrix{Float64}
    D0::Matrix{Float64}
    L0::Matrix{Float64}
    Avc::Matrix{Float64}
    Acc::Matrix{Float64}
    Acq::Matrix{Float64}
    Acv::Matrix{Float64}
end

function AngleOperator(M, D0, L0, Avc, Acc, Acq, Acv)
    n = size(M, 1)
    size(M) == (n, n) || throw(ArgumentError("M must be square"))
    size(D0) == size(L0) == (n, n) || throw(ArgumentError("D0/L0 dimensions disagree"))
    size(Avc, 1) == n || throw(ArgumentError("Avc row dimension must match retained q"))
    size(Acc, 1) == size(Acc, 2) == size(Avc, 2) ||
        throw(ArgumentError("Acc/Avc condensed dimensions disagree"))
    size(Acq) == (size(Acc, 1), n) || throw(ArgumentError("Acq dimensions disagree"))
    size(Acv) == (size(Acc, 1), n) || throw(ArgumentError("Acv dimensions disagree"))
    return AngleOperator(Matrix{Float64}(M), Matrix{Float64}(D0), Matrix{Float64}(L0),
                         Matrix{Float64}(Avc), Matrix{Float64}(Acc),
                         Matrix{Float64}(Acq), Matrix{Float64}(Acv))
end

function pi_q(op::AngleOperator, s::Number)
    Kcv = op.Acv * op.M
    Q = ComplexF64(s) * I - op.Acc
    return -op.Avc * (Q \ (op.Acq + s * Kcv))
end

function pi_q_derivative(op::AngleOperator, s::Number)
    Kcv = op.Acv * op.M
    Q = ComplexF64(s) * I - op.Acc
    X = Q \ (op.Acq + s * Kcv)
    return op.Avc * ((Q \ X) - (Q \ Kcv))
end

psi(op::AngleOperator, LG::AbstractMatrix, s::Number) =
    s * op.D0 + (op.L0 - LG) + pi_q(op, s)

psi_derivative(op::AngleOperator, s::Number) = op.D0 + pi_q_derivative(op, s)

retained_operator(op::AngleOperator, s::Number) =
    s^2 * op.M + s * op.D0 + op.L0 + pi_q(op, s)

function generalized_reconstruction_error(op::AngleOperator, LG::AbstractMatrix, s::Number)
    T = retained_operator(op, s)
    reconstructed = s^2 * op.M + LG + psi(op, LG, s)
    return norm(T - reconstructed) / max(norm(T), eps(Float64))
end

"""Compare the analytic derivative with a 5-point symmetric finite difference."""
function derivative_check(op::AngleOperator, s::Number)
    z = ComplexF64(s)
    h = 2e-5 * max(1.0, abs(z))
    numeric = (-pi_q(op, z + 2h) + 8pi_q(op, z + h) -
               8pi_q(op, z - h) + pi_q(op, z - 2h)) / (12h)
    analytic = pi_q_derivative(op, z)
    return norm(analytic - numeric) / max(norm(analytic), norm(numeric), eps(Float64))
end

end
