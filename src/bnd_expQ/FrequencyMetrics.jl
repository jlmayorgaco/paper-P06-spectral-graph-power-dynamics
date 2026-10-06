module FrequencyMetrics

using LinearAlgebra

export sg_frequency_deviation, pll_frequency_deviation, unwrap_phase,
       sg_polynomial_derivative, savgol_cutoff_hz, bus_frequency_deviation,
       primary_frequency_peak

"""Synchronous-machine rotor speed `omega` is per-unit on the 60 Hz base."""
sg_frequency_deviation(omega_pu; f0_hz=60.0) = Float64(f0_hz) * (Float64(omega_pu) - 1.0)

"""Installed SimpleGFLDC Δω state is rad/s relative to its frame."""
pll_frequency_deviation(delta_omega_rad_s) = Float64(delta_omega_rad_s) / (2pi)

function unwrap_phase(theta::AbstractVector{<:Real})
    n = length(theta)
    out = Vector{Float64}(undef, n)
    n == 0 && return out
    out[1] = Float64(theta[1])
    for k in 2:n
        d = Float64(theta[k]) - Float64(theta[k-1])
        d = mod(d + pi, 2pi) - pi
        out[k] = out[k-1] + d
    end
    out
end

"""
Savitzky–Golay local-polynomial differentiator.

The polynomial is fit in sample-index coordinates and its first derivative is
converted to per-second units. At endpoints, the same least-squares fit uses
the nearest complete one-sided window. The default 11-sample window at 100 Hz
has 0.1 s support; its -3 dB derivative-response bandwidth is computed by
`savgol_cutoff_hz`, rather than inferred from the window length.
"""
function sg_polynomial_derivative(y::AbstractVector{<:Real}, dt::Real;
                                  half_window::Integer=5, degree::Integer=3)
    dt > 0 || throw(ArgumentError("dt must be positive"))
    half_window >= 1 || throw(ArgumentError("half_window must be positive"))
    degree >= 1 || throw(ArgumentError("degree must be positive"))
    n = length(y)
    n >= degree + 1 || throw(ArgumentError("not enough samples for polynomial degree"))
    out = Vector{Float64}(undef, n)
    window = 2half_window + 1
    for k in 1:n
        lo = clamp(k-half_window, 1, max(1, n-window+1))
        hi = min(n, lo+window-1)
        idx = lo:hi
        τ = Float64.(idx .- k)
        X = hcat((τ .^ j for j in 0:degree)...)
        coef = X \ Float64.(y[idx])
        out[k] = coef[2] / Float64(dt)
    end
    out
end

function savgol_cutoff_hz(dt::Real; half_window::Integer=5, degree::Integer=3,
                          points::Integer=20001)
    dt > 0 || throw(ArgumentError("dt must be positive"))
    points >= 100 || throw(ArgumentError("frequency grid too small"))
    τ = collect(-half_window:half_window)
    X = hcat((Float64.(τ) .^ j for j in 0:degree)...)
    # Derivative at sample zero is the second row of the least-squares map.
    w = (X'X) \ X'
    deriv_weights = vec(w[2, :]) ./ Float64(dt)
    nyquist = 0.5 / Float64(dt)
    fs = range(max(nyquist/(points-1), 1e-6), nyquist; length=points)
    response = [abs(sum(deriv_weights[k] * cis(-2pi*f*τ[k]*dt)
                        for k in eachindex(τ))) / (2pi*f) for f in fs]
    target = inv(sqrt(2.0))
    k = findfirst(<=(target), response)
    k === nothing ? nyquist : Float64(fs[k])
end

function bus_frequency_deviation(voltage::AbstractVector{<:Complex}, dt::Real;
                                 half_window::Integer=5, degree::Integer=3)
    phase = unwrap_phase(angle.(voltage))
    sg_polynomial_derivative(phase, dt; half_window, degree) ./ (2pi)
end

function primary_frequency_peak(sg::AbstractMatrix{<:Real},
                                pll::AbstractMatrix{<:Real})
    vals = vcat(vec(sg[isfinite.(sg)]), vec(pll[isfinite.(pll)]))
    isempty(vals) && return NaN
    maximum(abs, vals)
end

end
