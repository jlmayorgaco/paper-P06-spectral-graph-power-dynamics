module SyntheticSystems

using LinearAlgebra
using Random
using ..GraphBasis: generalized_graph_basis
using ..DynamicSelfEnergy: StateSpaceSelfEnergy, with_offdiagonal_scale

export SyntheticCase, build_synthetic_cases, build_optional_s5,
       build_resonance_case, nodal_realization, modal_realization, case_table_row

struct SyntheticCase
    name::String
    seed::Int
    description::String
    M::Matrix{Float64}
    L::Matrix{Float64}
    Lambda::Vector{Float64}
    basis
    modal_sys::StateSpaceSelfEnergy
    nodal_sys::StateSpaceSelfEnergy
    modal_frequencies_hz::Vector{Float64}
end

function metric_factor(n::Int)
    R = Matrix{Float64}(I, n, n)
    for i in 1:n
        R[i, i] = 0.9 + 0.12 * i
        for j in i+1:n
            R[i, j] = 0.025 * sin(0.7i + 1.3j)
        end
    end
    return R
end

function transfer_realization(D0::AbstractMatrix, modal_frequencies_hz::Vector{Float64};
                              seed::Int,
                              diagonal_gain::Real=0.08,
                              offdiagonal_gain::Real=0.0,
                              offdiagonal_edges=nothing,
                              pole_frequencies_hz=[0.2, 0.7, 2.0, 5.0, 20.0])
    n = length(modal_frequencies_hz)
    rng = MersenneTwister(seed)
    edges = if offdiagonal_edges === nothing
        [(i, j) for i in 1:n for j in 1:n if i != j && offdiagonal_gain != 0]
    else
        collect(offdiagonal_edges)
    end
    channels = [(i, i) for i in 1:n]
    append!(channels, edges)
    nc = length(channels)
    rates = zeros(Float64, nc)
    B = zeros(Float64, nc, n)
    C = zeros(Float64, n, nc)
    offdiag_indices = Int[]
    for (r, (i, j)) in enumerate(channels)
        fp = pole_frequencies_hz[mod1(r, length(pole_frequencies_hz))]
        rate = 2π * fp
        rates[r] = rate
        target_f = i == j ? modal_frequencies_hz[i] :
                   sqrt(max(modal_frequencies_hz[i] * modal_frequencies_hz[j], 0.0))
        target_w = 2π * target_f
        if i == j
            magnitude = diagonal_gain * (0.96 + 0.08 * rand(rng))
            residue = magnitude * hypot(rate, target_w)
        else
            push!(offdiag_indices, r)
            directed_scale = 0.82 + 0.36 * rand(rng)
            residue = offdiagonal_gain * directed_scale * hypot(rate, target_w)
            # A small, deterministic asymmetry avoids assuming reciprocal coupling.
            isodd(i + j) && (residue = -residue)
        end
        b = sqrt(abs(residue))
        c = sign(residue) * b
        B[r, j] = b
        C[i, r] = c
    end
    sys = StateSpaceSelfEnergy(D0, Matrix(Diagonal(-rates)), B, C;
        channel_pairs=channels, offdiag_state_indices=offdiag_indices,
        offdiag_base_scale=1.0)
    return sys
end

function nodal_realization(modal_sys::StateSpaceSelfEnergy, R::AbstractMatrix)
    return StateSpaceSelfEnergy(R' * modal_sys.D0 * R, modal_sys.Ac,
        modal_sys.B * R, R' * modal_sys.C;
        channel_pairs=modal_sys.channel_pairs,
        offdiag_state_indices=modal_sys.offdiag_state_indices,
        offdiag_base_scale=modal_sys.offdiag_base_scale)
end

modal_realization(sys::StateSpaceSelfEnergy) = sys

function make_case(name, seed, description, frequencies_hz, D0, sys)
    n = length(frequencies_hz)
    lambda = (2π .* frequencies_hz).^2
    R = metric_factor(n)
    M = R' * R
    L = R' * Diagonal(lambda) * R
    basis = generalized_graph_basis(M, L; zero_tol=1e-12, tol_deg=1e-7)
    return SyntheticCase(String(name), Int(seed), String(description), M, Matrix(L),
        basis.Lambda, basis, sys, nodal_realization(sys, R), copy(frequencies_hz))
end

function build_synthetic_cases(config=nothing)
    seedcfg = config === nothing ? Dict{String,Any}() : get(config, "seeds", Dict{String,Any}())
    freqcfg = config === nothing ? Dict{String,Any}() : get(config, "graph_modal_frequencies_hz", Dict{String,Any}())
    gaincfg = config === nothing ? Dict{String,Any}() : get(config, "controller_gains", Dict{String,Any}())
    directcfg = config === nothing ? Dict{String,Any}() : get(config, "direct_diagonal", Dict{String,Any}())
    controller_poles = config === nothing ? [0.2,0.7,2.0,5.0,20.0] :
        Float64.(get(config, "controller_pole_frequencies_hz", [0.2,0.7,2.0,5.0,20.0]))
    resonance_poles = config === nothing ? [3.0,8.0,0.2,20.0,5.0] :
        Float64.(get(config, "S4_controller_pole_frequencies_hz", [3.0,8.0,0.2,20.0,5.0]))
    default_f = [0.0, 0.8, 2.1, 4.5, 8.0]
    default_d = [0.15, 0.20, 0.23, 0.27, 0.31]
    fbase = Float64.(get(freqcfg, "S0_S1_S2", default_f))
    f3 = Float64.(get(freqcfg, "S3", [0.0, 0.8, 0.81, 4.5, 8.0]))
    f4 = Float64.(get(freqcfg, "S4", [0.0, 3.0, 3.005, 6.0, 9.0]))
    seedof(key, default) = Int(get(seedcfg, key, default))
    gainof(key, default) = Float64(get(gaincfg, key, default))
    diagonal(key, default) = Float64.(get(directcfg, key, default))
    d0 = Matrix(Diagonal(diagonal("S0", default_d)))

    # S0: all controller paths are diagonal in the generalized graph basis.
    s0sys = transfer_realization(d0, fbase; seed=seedof("S0",301), diagonal_gain=gainof("S0_diagonal",0.10),
                                 pole_frequencies_hz=controller_poles)
    S0 = make_case("S0", seedof("S0",301), "Commuting dynamic baseline: Σ̂(s) diagonal for all s.", fbase, d0, s0sys)

    # S1: static symmetric, non-proportional modal damping, with no controller states.
    d1 = Matrix(Diagonal(diagonal("S1", default_d)))
    static_edges = config === nothing ? [(1,2,0.025),(2,3,0.045),(2,4,-0.030),(3,5,0.035)] :
        get(config, "static_nonproportional_edges", [[1,2,0.025],[2,3,0.045],[2,4,-0.030],[3,5,0.035]])
    for edge in static_edges
        i,j,a = Int(edge[1]), Int(edge[2]), Float64(edge[3])
        d1[i, j] = a
        d1[j, i] = a
    end
    s1sys = StateSpaceSelfEnergy(d1, zeros(0, 0), zeros(0, 5), zeros(5, 0))
    S1 = make_case("S1", seedof("S1",302), "Static symmetric non-proportional modal self-energy.", fbase, d1, s1sys)

    # S2: stable first-order controller paths with poles spanning 0.2–20 Hz.
    d2 = Matrix(Diagonal(diagonal("S2", default_d)))
    s2sys = transfer_realization(d2, fbase; seed=seedof("S2",303),
        diagonal_gain=gainof("S2_diagonal",0.10), offdiagonal_gain=gainof("S2_offdiagonal",0.045),
        pole_frequencies_hz=controller_poles)
    S2 = make_case("S2", seedof("S2",303), "Dynamic state-space self-energy with frequency-dependent intermodal entries.", fbase, d2, s2sys)

    # S3: the off-diagonal state-space paths are scaled by ε; diagonal paths stay fixed.
    d3 = Matrix(Diagonal(diagonal("S3", [0.18, 0.21, 0.24, 0.29, 0.34])))
    s3sys = transfer_realization(d3, f3; seed=seedof("S3",304),
        diagonal_gain=gainof("S3_diagonal",0.075), offdiagonal_gain=gainof("S3_offdiagonal",0.25),
        pole_frequencies_hz=controller_poles)
    S3 = make_case("S3", seedof("S3",304), "Weak-to-strong physically realized off-diagonal coupling sweep.", f3, d3, s3sys)

    # S4: two distinct but closely spaced graph frequencies with weak reciprocal paths.
    d4 = Matrix(Diagonal(diagonal("S4", fill(0.18, 5))))
    s4sys = transfer_realization(d4, f4; seed=seedof("S4",305),
                                 diagonal_gain=gainof("S4_diagonal",0.015),
                                 offdiagonal_gain=gainof("S4_offdiagonal",0.005),
                                 offdiagonal_edges=[(2, 3), (3, 2)],
                                 pole_frequencies_hz=resonance_poles)
    S4 = make_case("S4", seedof("S4",305), "Near-resonant nondegenerate graph-mode pair with small dynamic coupling.", f4, d4, s4sys)
    return [S0, S1, S2, S3, S4]
end

"""Optional static counterexample: positive diagonal Hermitian entries,
but an indefinite total Hermitian operator.
"""
function build_optional_s5(; seed::Int=306, cross_gain::Real=1.2)
    f = [0.0, 1.0, 2.0, 4.0, 7.0]
    D = Matrix{Float64}(I, 5, 5)
    D[1, 2] = D[2, 1] = cross_gain
    sys = StateSpaceSelfEnergy(D, zeros(0, 0), zeros(0, 5), zeros(5, 0))
    return make_case("S5", seed, "Optional positive diagonal but indefinite Hermitian counterexample.", f, D, sys)
end

function build_resonance_case(f3::Real=3.005; base_frequency_hz::Real=3.0,
                              seed::Int=305,
                              diagonal_gain::Real=0.015,
                              offdiagonal_gain::Real=0.005,
                              pole_frequencies_hz=[3.0, 8.0, 0.2, 20.0, 5.0])
    f0 = Float64(base_frequency_hz)
    f = [0.0, f0, Float64(f3), 2 * f0, 3 * f0]
    D = Matrix(Diagonal(fill(0.18, 5)))
    sys = transfer_realization(D, f; seed=seed, diagonal_gain=diagonal_gain,
        offdiagonal_gain=offdiagonal_gain, offdiagonal_edges=[(2, 3), (3, 2)],
        pole_frequencies_hz=pole_frequencies_hz)
    return make_case("S4", seed,
        "Near-resonant nondegenerate graph-mode pair with small dynamic coupling.",
        f, D, sys)
end

function case_table_row(case::SyntheticCase)
    return (case=case.name, m=size(case.M, 1), nc=size(case.modal_sys.Ac, 1),
            seed=case.seed, min_eig_M=minimum(eigvals(Symmetric(case.M))),
            zero_modes=length(case.basis.zero_modes),
            Ac_spectral_abscissa=size(case.modal_sys.Ac, 1) == 0 ? nothing :
                maximum(real.(eigvals(case.modal_sys.Ac))), description=case.description)
end

end
