module PolePredictor

using LinearAlgebra
using ..DynamicSelfEnergy: StateSpaceSelfEnergy, sigma, sigma_derivative,
                          diagonal_only, augmented_matrix
using ..IntermodalSelfEnergy: schur_self_energy

export uncoupled_modal_pole, pole_derivative, predict_pole_shift,
       track_augmented_pole, complex_newton

function pole_derivative(sys::StateSpaceSelfEnergy, s::Number, k::Integer)
    S = sigma(sys, s)
    Sp = sigma_derivative(sys, s)
    return 2s + S[k, k] + s * Sp[k, k]
end

function uncoupled_modal_pole(sys::StateSpaceSelfEnergy, Lambda::AbstractVector,
                              k::Integer; sign::Int=1, mass=I)
    dsys = diagonal_only(sys)
    A = augmented_matrix(mass, Diagonal(Lambda), dsys)
    E = eigen(A)
    vals = E.values
    target_imag = sqrt(max(Lambda[k], 0.0))
    candidates = filter(x -> sign * imag(x) > 1e-8, vals)
    isempty(candidates) && throw(ArgumentError("no oscillatory pole found for mode $k"))
    target = sign * target_imag
    n = length(Lambda)
    controller_states = [2n + r for r in eachindex(dsys.channel_pairs)
                         if dsys.channel_pairs[r] == (k, k)]
    support_indices = vcat([k, n + k], controller_states)
    scores = Float64[]
    for value in candidates
        j = argmin(abs.(E.values .- value))
        vector = E.vectors[:, j]
        participation = sum(abs2, vector[support_indices]) / sum(abs2, vector)
        distance = abs(imag(value) - target) + 0.02 * abs(real(value))
        push!(scores, distance + 100.0 * (1 - participation))
    end
    return candidates[argmin(scores)]
end

"""One Newton solve with analytic derivative; used for independent NEP checks."""
function complex_newton(f, df, s0::Number; tol::Real=1e-12, maxiter::Integer=80)
    s = ComplexF64(s0)
    for _ in 1:maxiter
        fv = f(s)
        dv = df(s)
        abs(dv) > eps(Float64) || throw(ArgumentError("Newton derivative is singular"))
        step = fv / dv
        s -= step
        abs(step) <= tol * max(1.0, abs(s)) && return s
    end
    throw(ErrorException("complex Newton iteration did not converge"))
end

function predict_pole_shift(sys::StateSpaceSelfEnergy, Lambda::AbstractVector,
                            s0::Number, k::Integer)
    T = s0^2 .* Matrix{ComplexF64}(I, length(Lambda), length(Lambda)) .+
        s0 .* sigma(sys, s0) + Diagonal(Lambda)
    gamma = schur_self_energy(T, k).gamma
    derivative = pole_derivative(sys, s0, k)
    shift = -gamma / derivative
    return (gamma=gamma, derivative=derivative, shift=shift, predicted=s0 + shift)
end

"""Continue an eigenvalue branch using proximity and normalized eigenvector overlap."""
function track_augmented_pole(A::AbstractMatrix, previous_value::Number,
                              previous_vector::AbstractVector; positive_imag::Bool=true)
    E = eigen(A)
    candidates = findall(i -> positive_imag ? imag(E.values[i]) > 1e-8 : true,
                         eachindex(E.values))
    isempty(candidates) && throw(ArgumentError("no eligible eigenvalue for continuation"))
    vp = previous_vector / norm(previous_vector)
    scores = map(candidates) do i
        vi = E.vectors[:, i] / norm(E.vectors[:, i])
        overlap = abs(dot(vp, vi))
        distance = abs(E.values[i] - previous_value) / max(1.0, abs(previous_value))
        distance + 0.35 * (1 - overlap)
    end
    j = candidates[argmin(scores)]
    return (value=E.values[j], vector=E.vectors[:, j], overlap=abs(dot(vp, E.vectors[:, j] / norm(E.vectors[:, j]))))
end

end
