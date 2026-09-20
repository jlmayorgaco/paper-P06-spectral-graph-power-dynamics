"""Nonlinear fixed-step DAE traces for the frozen custom same-model cases."""

using LinearAlgebra
using Printf
using Statistics

include(joinpath(@__DIR__, "run_true_same_model_ieee39.jl"))

const TDS_OUT = joinpath(CAMPAIGN, "raw", "tds_same_model")
mkpath(TDS_OUT)
const DT = 0.002
const T_END = 8.0
const PULSE_END = 0.10
const PULSE_FRACTION = 0.02
const PULSE_BUS = 3

function g_residual_disturbance(d::FrozenDAE, x, z, t)
    v = ComplexF64[voltage(d, z, bus) for bus in BUS_IDS]
    injection = zeros(ComplexF64, length(BUS_IDS))
    for slot in d.slots
        vv = v[slot.bus]
        contribution = slot.kind == :machine ? machine_injection(slot.machine, x[slot.start:slot.stop], vv) : gfl_injection(slot.gfl, x[slot.start:slot.stop], vv)
        injection[slot.bus] += contribution
    end
    for (bus, load) in LOAD_BY_BUS
        effective = bus == PULSE_BUS && t <= PULSE_END ? load + PULSE_FRACTION * real(load) : load
        injection[bus] -= conj(effective) / conj(v[bus])
    end
    residual = d.ybus * v - injection
    out = zeros(Float64, 2 * length(BUS_IDS))
    out[1:2:end] = real.(residual)
    out[2:2:end] = imag.(residual)
    out
end

function solve_algebraic(d::FrozenDAE, x, z_initial, t)
    z = copy(z_initial)
    for iteration in 1:16
        r = g_residual_disturbance(d, x, z, t)
        norm(r, Inf) < 1e-9 && return z, true, iteration, norm(r, Inf)
        nz = length(z)
        gz = zeros(Float64, nz, nz)
        for j in 1:nz
            h = cbrt(eps(Float64)) * max(abs(z[j]), 1.0)
            zp, zm = copy(z), copy(z)
            zp[j] += h
            zm[j] -= h
            gz[:, j] = (g_residual_disturbance(d, x, zp, t) - g_residual_disturbance(d, x, zm, t)) / (2h)
        end
        step = -(gz \ r)
        accepted = false
        current = norm(r, Inf)
        for scale in (1.0, 0.5, 0.25, 0.125, 0.0625)
            trial = z .+ scale .* step
            trial_norm = norm(g_residual_disturbance(d, x, trial, t), Inf)
            if isfinite(trial_norm) && trial_norm < current
                z = trial
                accepted = true
                break
            end
        end
        accepted || return z, false, iteration, current
    end
    z, false, 16, norm(g_residual_disturbance(d, x, z, t), Inf)
end

function rhs_with_algebraic(d::FrozenDAE, x, z_guess, t)
    z, ok, iterations, residual = solve_algebraic(d, x, z_guess, t)
    ok || error("TDS algebraic solve failed at t=$(t), residual=$(residual)")
    f_residual(d, x, z), z, iterations, residual
end

function signal_values(d::FrozenDAE, x, z)
    speeds = Float64[]
    for slot in d.slots
        slot.kind == :machine && push!(speeds, x[slot.start + 1] - 1.0)
    end
    mean_speed = isempty(speeds) ? NaN : mean(speeds)
    v3 = abs(voltage(d, z, 3))
    v30 = abs(voltage(d, z, 30))
    (mean_speed=mean_speed, tracked_state_1=x[1], v3=v3, v30=v30, state_norm=norm(x))
end

function seed_critical_mode(d::FrozenDAE, x, z)
    fx, fz, gx, gz = jacobians(d, x, z)
    a = fx - fz * (gz \ gx)
    spectrum = eigen(a)
    keep = findall(abs.(spectrum.values) .>= GAUGE_TOL)
    index = keep[argmax(real.(spectrum.values[keep]))]
    vector = real.(spectrum.vectors[:, index])
    norm(vector) < 1e-12 && (vector = imag.(spectrum.vectors[:, index]))
    x .+ 5e-2 .* vector ./ max(norm(vector), 1e-12)
end

function simulate_case(label, replaced::Set{Int}; mode_seed=false)
    v, _, _ = solve_powerflow()
    d, x = build_dae(v, replaced)
    z = zeros(Float64, 2 * length(BUS_IDS))
    z[1:2:end] = real.(v)
    z[2:2:end] = imag.(v)
    mode_seed && (x = seed_critical_mode(d, x, z))
    n_steps = Int(round(T_END / DT))
    times = collect(0.0:DT:T_END)
    speed = zeros(Float64, length(times))
    v3 = zeros(Float64, length(times))
    v30 = zeros(Float64, length(times))
    tracked = zeros(Float64, length(times))
    state_norm = zeros(Float64, length(times))
    algebraic_residual = zeros(Float64, length(times))
    max_newton = 0
    for k in 1:length(times)
        t = times[k]
        z, ok, nit, gres = solve_algebraic(d, x, z, t)
        ok || error("initial/output algebraic solve failed for $(label) at t=$(t)")
        max_newton = max(max_newton, nit)
        sig = signal_values(d, x, z)
        speed[k], tracked[k], v3[k], v30[k], state_norm[k] = sig.mean_speed, sig.tracked_state_1, sig.v3, sig.v30, sig.state_norm
        algebraic_residual[k] = gres
        k == length(times) && break
        k1, z1, _, _ = rhs_with_algebraic(d, x, z, t)
        k2, z2, _, _ = rhs_with_algebraic(d, x .+ (DT / 2) .* k1, z1, t + DT / 2)
        k3, z3, _, _ = rhs_with_algebraic(d, x .+ (DT / 2) .* k2, z2, t + DT / 2)
        k4, z4, _, _ = rhs_with_algebraic(d, x .+ DT .* k3, z3, t + DT)
        x = x .+ (DT / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
        z = z4
    end
    path = joinpath(TDS_OUT, "$(label)_trace.csv")
    open(path, "w") do io
        println(io, "time_s,mean_machine_speed_deviation_pu,tracked_state_1,bus3_voltage_pu,bus30_voltage_pu,state_norm,algebraic_residual")
        for k in eachindex(times)
            println(io, join((times[k], speed[k], tracked[k], v3[k], v30[k], state_norm[k], algebraic_residual[k]), ','))
        end
    end
    (times=times, speed=speed, tracked=tracked, v3=v3, v30=v30, state_norm=state_norm, algebraic_residual=algebraic_residual, max_newton=max_newton)
end

function detrended(values, times)
    tail = times .>= 0.2
    values[tail] .- mean(values[tail])
end

function prony_ar2(values, dt)
    n = length(values)
    n < 4 && return (alpha=NaN, frequency_hz=NaN)
    x1 = values[2:n-1]
    x2 = values[1:n-2]
    y = values[3:n]
    coeff = hcat(x1, x2) \ y
    roots = eigvals(ComplexF64[coeff[1] coeff[2]; 1.0 0.0])
    roots = roots[abs.(roots) .> 1e-10]
    isempty(roots) && return (alpha=NaN, frequency_hz=NaN)
    root = roots[argmax(abs.(imag.(log.(roots))))]
    lambda = log(root) / dt
    (alpha=real(lambda), frequency_hz=abs(imag(lambda)) / (2pi))
end

function band_dft(values, times)
    centered = values .- mean(values)
    frequencies = collect(0.3:0.005:1.5)
    best_f, best_amp = NaN, -Inf
    for f in frequencies
        projection = sum(centered .* exp.(-2pi * im * f .* times))
        amp = abs(projection)
        if amp > best_amp
            best_amp = amp
            best_f = f
        end
    end
    (frequency_hz=best_f, amplitude=best_amp)
end

function run_tds()
    cases = [("base", Set{Int}(), false), ("proper_30+33+35", Set([30, 33, 35]), false), ("repaired_30+33+35", Set([30, 33, 35]), false), ("v4_blocker", TARGET_BUSES, false), ("v4_blocker_mode_seed", TARGET_BUSES, true)]
    open(joinpath(TDS_OUT, "tds_same_model_status.csv"), "w") do io
        println(io, "case,declared_stable,prony_alpha_s_inv,prony_frequency_hz,dft_frequency_hz,dft_amplitude,tail_growth_ratio,max_state_norm,max_algebraic_residual,max_newton,tds_label")
        for (label, replaced, mode_seed) in cases
            trace = simulate_case(label, replaced; mode_seed=mode_seed)
            y = detrended(trace.tracked, trace.times)
            prony = prony_ar2(y, DT)
            tail_times = trace.times[trace.times .>= 0.2]
            dft = band_dft(y, tail_times)
            tail = trace.times .>= 3.0
            early = (trace.times .>= 0.2) .& (trace.times .<= 1.0)
            growth = maximum(abs.(trace.tracked[tail])) / max(maximum(abs.(trace.tracked[early])), 1e-12)
            declared = startswith(label, "v4_blocker") ? false : true
            tds_label = prony.alpha < 0.0 ? "STABLE" : (prony.alpha > 0.0 ? "UNSTABLE" : "INDETERMINATE")
            println(io, join((label, declared, prony.alpha, prony.frequency_hz, dft.frequency_hz, dft.amplitude, growth, maximum(trace.state_norm), maximum(abs.(trace.algebraic_residual)), trace.max_newton, tds_label), ','))
            @printf("TDS_SAME_MODEL case=%s declared_stable=%s prony_alpha=%.6e prony_f=%.6f dft_f=%.6f growth=%.6e label=%s\n", label, declared, prony.alpha, prony.frequency_hz, dft.frequency_hz, growth, tds_label)
        end
    end
end

run_tds()
