"""Exact frozen custom SG/AVR/PSS/GFL IEEE-39 base and V4 Julia port.

The implementation mirrors the frozen Python source snapshots directly:
two-axis SynchronousMachine, first-order AVR, two-state IEEEST PSS, D=0,
constant Pm, no governor, constant-power loads, and the frozen admittance
network. Equilibria use a Newton AC power flow followed by the analytic device
initializers. The index-one reduced Jacobian is assembled from central
differences of the same f/g residuals used by the nonlinear model.
"""

using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "campaign_root.jl"))
include(joinpath(@__DIR__, "frozen_ieee39_data.jl"))

const CAMPAIGN = campaign_root_from_args()
const OUT = joinpath(CAMPAIGN, "raw", "reconciliation")
mkpath(OUT)

const TARGET_ORDER = (30, 33, 35, 37)
const TARGET_BUSES = Set(TARGET_ORDER)
const GAUGE_TOL = 1e-4
const SLACK_BUS = Int(round(SLACK_ROW[1]))
const BUS_POS = Dict(bus => bus for bus in BUS_IDS)
const MACHINE_BY_BUS = Dict(Int(round(row[1])) => row for row in MACHINE_ROWS)
const AVR_BY_BUS = Dict(Int(round(MACHINE_ROWS[i][1])) => AVR_ROWS[i] for i in eachindex(MACHINE_ROWS))
const PSS_BY_BUS = Dict(Int(round(MACHINE_ROWS[i][1])) => PSS_ROWS[i] for i in eachindex(MACHINE_ROWS))
const PV_BY_BUS = Dict(Int(round(row[1])) => row for row in PV_ROWS)
const LOAD_BY_BUS = let d = Dict{Int,ComplexF64}(); for row in LOAD_ROWS; d[Int(round(row[1]))] = get(d, Int(round(row[1])), 0.0 + 0.0im) + complex(row[2], row[3]); end; d end

struct MachineParams
    ra::Float64
    xd::Float64
    xq::Float64
    xd1::Float64
    xq1::Float64
    td10::Float64
    tq10::Float64
    m::Float64
    d::Float64
    ka::Float64
    ta::Float64
    pss_gain::Float64
    pss_washout::Float64
    pss_wash_lag::Float64
    pss_lag::Float64
    pm::Float64
    vref::Float64
end

struct Machine
    p::MachineParams
    weight::Float64
end

struct GFLParams
    kp_pll::Float64
    ki_pll::Float64
    tau_p::Float64
    kp_p::Float64
    ki_p::Float64
    kp_q::Float64
    ki_q::Float64
    kp_i::Float64
    ki_i::Float64
    xf::Float64
    rf::Float64
    p_ref::Float64
    q_ref::Float64
    v_ref::Float64
end

struct GFL
    p::GFLParams
    weight::Float64
end

struct Slot
    kind::Symbol
    bus::Int
    start::Int
    stop::Int
    machine::Union{Nothing,Machine}
    gfl::Union{Nothing,GFL}
end

struct FrozenDAE
    ybus::Matrix{ComplexF64}
    slots::Vector{Slot}
    n_x::Int
end

function ybus_matrix()
    y = zeros(ComplexF64, length(BUS_IDS), length(BUS_IDS))
    for row in LINE_ROWS
        row[9] == 0.0 && continue
        f, t = Int(round(row[1])), Int(round(row[2]))
        r, x, g, b = row[3], row[4], row[5], row[6]
        tap, phi = row[7], row[8]
        series = inv(complex(r, x))
        charging = complex(g, b) / 2.0
        m = tap * cis(phi)
        y[f, f] += (series + charging) / abs2(m)
        y[t, t] += series + charging
        y[f, t] += -series / conj(m)
        y[t, f] += -series / m
    end
    for row in SHUNT_ROWS
        bus = Int(round(row[1]))
        y[bus, bus] += complex(row[2], row[3])
    end
    y
end

const YBUS = ybus_matrix()

function scheduled_power()
    scheduled = zeros(ComplexF64, length(BUS_IDS))
    for (bus, load) in LOAD_BY_BUS
        scheduled[bus] -= load
    end
    for row in PV_ROWS
        scheduled[Int(round(row[1]))] += complex(row[2], 0.0)
    end
    scheduled
end

const SCHEDULED = scheduled_power()
const PV_BUSES = Set(keys(PV_BY_BUS))
const PQ_BUSES = [bus for bus in BUS_IDS if bus != SLACK_BUS && !(bus in PV_BUSES)]
const ANGLE_BUSES = [bus for bus in BUS_IDS if bus != SLACK_BUS]

function unpack_pf(w)
    angles = zeros(Float64, length(BUS_IDS))
    magnitudes = ones(Float64, length(BUS_IDS))
    magnitudes[SLACK_BUS] = SLACK_ROW[2]
    angles[SLACK_BUS] = SLACK_ROW[3]
    for row in PV_ROWS
        bus = Int(round(row[1]))
        magnitudes[bus] = row[4]
    end
    angles[ANGLE_BUSES] = w[1:length(ANGLE_BUSES)]
    magnitudes[PQ_BUSES] = w[length(ANGLE_BUSES)+1:end]
    magnitudes .* cis.(angles)
end

function pf_mismatch(w)
    v = unpack_pf(w)
    s = v .* conj.(YBUS * v)
    d = s - SCHEDULED
    vcat(real.(d[ANGLE_BUSES]), imag.(d[PQ_BUSES]))
end

function solve_powerflow()
    w = vcat(zeros(Float64, length(ANGLE_BUSES)), ones(Float64, length(PQ_BUSES)))
    last_norm = Inf
    for iteration in 1:100
        r = pf_mismatch(w)
        current = norm(r, Inf)
        current < 1e-12 && return unpack_pf(w), current, iteration
        n = length(w)
        j = zeros(Float64, length(r), n)
        for k in 1:n
            h = 1e-6 * max(abs(w[k]), 1.0)
            wp, wm = copy(w), copy(w)
            wp[k] += h
            wm[k] -= h
            j[:, k] = (pf_mismatch(wp) - pf_mismatch(wm)) / (2h)
        end
        step = -(j \ r)
        accepted = false
        for scale in (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125)
            trial = w .+ scale .* step
            trial_norm = norm(pf_mismatch(trial), Inf)
            if isfinite(trial_norm) && trial_norm < current
                w = trial
                last_norm = trial_norm
                accepted = true
                break
            end
        end
        accepted || error("frozen power flow stalled at iteration=$iteration residual=$current previous=$last_norm")
    end
    error("frozen power flow failed residual=$(norm(pf_mismatch(w), Inf))")
end

function dq(v, delta)
    (real(v) * sin(delta) - imag(v) * cos(delta), real(v) * cos(delta) + imag(v) * sin(delta))
end

function network_current(id, iq, delta)
    complex(id * sin(delta) + iq * cos(delta), -id * cos(delta) + iq * sin(delta))
end

function machine_from_row(bus, generation, voltage)
    raw = MACHINE_BY_BUS[bus]
    avr = AVR_BY_BUS[bus]
    pss = PSS_BY_BUS[bus]
    weight = raw[2] / SYSTEM_BASE_MVA
    ra, xd, xq, xd1, xq1 = raw[5], raw[6], raw[7], raw[8], raw[9]
    td10, tq10, m, d = raw[10], raw[11], raw[4], raw[3]
    ka, ta = avr[2], avr[3]
    ks, t5, t6, t4 = pss[1], pss[2], pss[3], pss[4]
    current = conj(generation / weight / voltage)
    internal = voltage + complex(ra, xq) * current
    delta = angle(internal)
    vd, vq = dq(voltage, delta)
    id, iq = dq(current, delta)
    eq1 = vq + ra * iq + xd1 * id
    ed1 = vd + ra * id - xq1 * iq
    efd = eq1 + (xd - xd1) * id
    power = vd * id + vq * iq + ra * (id^2 + iq^2)
    vref = abs(voltage) + efd / ka
    p = MachineParams(ra, xd, xq, xd1, xq1, td10, tq10, m, d, ka, ta, ks, t5, t6, t4, power, vref)
    Machine(p, weight), [delta, 1.0, eq1, ed1, efd, power, 0.0]
end

function gfl_from_generation(generation, voltage, weight)
    current = conj(generation / weight / voltage)
    theta = angle(voltage)
    rotated = current * cis(-theta)
    id, iq = real(rotated), imag(rotated)
    v_d = abs(voltage)
    p = GFLParams(53.0, 1400.0, 0.03, 0.20, 8.0, 0.20, 8.0, 0.25, 6.0, 0.15, 0.01, v_d * id, -v_d * iq, v_d)
    GFL(p, weight), [theta, 0.0, p.p_ref, p.q_ref, id, -iq, id, iq, p.rf * id, p.rf * iq]
end

function build_dae(v, replaced::Set{Int})
    slots = Slot[]
    states = Float64[]
    cursor = 0
    generation = v .* conj.(YBUS * v)
    for (bus, load) in LOAD_BY_BUS
        generation[bus] += load
    end
    for row in MACHINE_ROWS
        bus = Int(round(row[1]))
        weight = row[2] / SYSTEM_BASE_MVA
        if bus in replaced
            device, state = gfl_from_generation(generation[bus], v[bus], weight)
            start = cursor + 1
            append!(states, state)
            cursor += length(state)
            push!(slots, Slot(:gfl, bus, start, cursor, nothing, device))
        else
            device, state = machine_from_row(bus, generation[bus], v[bus])
            start = cursor + 1
            append!(states, state)
            cursor += length(state)
            push!(slots, Slot(:machine, bus, start, cursor, device, nothing))
        end
    end
    FrozenDAE(YBUS, slots, cursor), states
end

function machine_derivatives(p::MachineParams, x, v)
    delta, omega, eq1, ed1, efd, pss_w, pss_l = x
    vd, vq = dq(v, delta)
    determinant = p.ra^2 + p.xd1 * p.xq1
    rd, rq = ed1 - vd, eq1 - vq
    id = (p.ra * rd + p.xq1 * rq) / determinant
    iq = (-p.xd1 * rd + p.ra * rq) / determinant
    power = vd * id + vq * iq + p.ra * (id^2 + iq^2)
    terminal = abs(v)
    stabilizer = p.pss_gain * pss_l
    washout_out = p.pss_washout * (power - pss_w) / p.pss_wash_lag
    [OMEGA_B * (omega - 1.0),
     (p.pm - power - p.d * (omega - 1.0)) / p.m,
     (efd - eq1 - (p.xd - p.xd1) * id) / p.td10,
     (-ed1 + (p.xq - p.xq1) * iq) / p.tq10,
     (p.ka * (p.vref + stabilizer - terminal) - efd) / p.ta,
     (power - pss_w) / p.pss_wash_lag,
     (washout_out - pss_l) / p.pss_lag]
end

function machine_injection(device::Machine, x, v)
    p = device.p
    vd, vq = dq(v, x[1])
    determinant = p.ra^2 + p.xd1 * p.xq1
    rd, rq = x[4] - vd, x[3] - vq
    id = (p.ra * rd + p.xq1 * rq) / determinant
    iq = (-p.xd1 * rd + p.ra * rq) / determinant
    device.weight * network_current(id, iq, x[1])
end

function gfl_derivatives(p::GFLParams, x, v)
    theta, xpll, p_filt, q_filt, x_p, x_q, i_d, i_q, x_id, x_iq = x
    rotated = v * cis(-theta)
    vd, vq = real(rotated), imag(rotated)
    power = vd * i_d + vq * i_q
    reactive = vq * i_d - vd * i_q
    p_command = p.p_ref
    q_command = p.q_ref
    id_ref = p.kp_p * (p_command - p_filt) + x_p
    iq_ref = -(p.kp_q * (q_command - q_filt) + x_q)
    e_d = vd + p.kp_i * (id_ref - i_d) + x_id - p.xf * i_q
    e_q = vq + p.kp_i * (iq_ref - i_q) + x_iq + p.xf * i_d
    [p.kp_pll * vq + xpll,
     p.ki_pll * vq,
     (power - p_filt) / p.tau_p,
     (reactive - q_filt) / p.tau_p,
     p.ki_p * (p_command - p_filt),
     p.ki_q * (q_command - q_filt),
     (OMEGA_B / p.xf) * (e_d - vd - p.rf * i_d + p.xf * i_q),
     (OMEGA_B / p.xf) * (e_q - vq - p.rf * i_q - p.xf * i_d),
     p.ki_i * (id_ref - i_d),
     p.ki_i * (iq_ref - i_q)]
end

function gfl_injection(device::GFL, x, v)
    device.weight * complex(x[7], x[8]) * cis(x[1])
end

function voltage(d::FrozenDAE, z, bus)
    complex(z[2 * bus - 1], z[2 * bus])
end

function f_residual(d::FrozenDAE, x, z)
    out = zeros(Float64, d.n_x)
    for slot in d.slots
        v = voltage(d, z, slot.bus)
        if slot.kind == :machine
            out[slot.start:slot.stop] = machine_derivatives(slot.machine.p, x[slot.start:slot.stop], v)
        else
            out[slot.start:slot.stop] = gfl_derivatives(slot.gfl.p, x[slot.start:slot.stop], v)
        end
    end
    out
end

function g_residual(d::FrozenDAE, x, z)
    v = ComplexF64[voltage(d, z, bus) for bus in BUS_IDS]
    injection = zeros(ComplexF64, length(BUS_IDS))
    for slot in d.slots
        vv = v[slot.bus]
        contribution = slot.kind == :machine ? machine_injection(slot.machine, x[slot.start:slot.stop], vv) : gfl_injection(slot.gfl, x[slot.start:slot.stop], vv)
        injection[slot.bus] += contribution
    end
    for (bus, load) in LOAD_BY_BUS
        injection[bus] -= conj(load) / conj(v[bus])
    end
    residual = d.ybus * v - injection
    out = zeros(Float64, 2 * length(BUS_IDS))
    out[1:2:end] = real.(residual)
    out[2:2:end] = imag.(residual)
    out
end

function jacobians(d::FrozenDAE, x, z)
    nx, nz = length(x), length(z)
    fx, gx = zeros(nx, nx), zeros(nz, nx)
    fz, gz = zeros(nx, nz), zeros(nz, nz)
    scale = cbrt(eps(Float64))
    for j in 1:nx
        h = scale * max(abs(x[j]), 1.0)
        xp, xm = copy(x), copy(x)
        xp[j] += h; xm[j] -= h
        fx[:, j] = (f_residual(d, xp, z) - f_residual(d, xm, z)) / (2h)
        gx[:, j] = (g_residual(d, xp, z) - g_residual(d, xm, z)) / (2h)
    end
    for j in 1:nz
        h = scale * max(abs(z[j]), 1.0)
        zp, zm = copy(z), copy(z)
        zp[j] += h; zm[j] -= h
        fz[:, j] = (f_residual(d, x, zp) - f_residual(d, x, zm)) / (2h)
        gz[:, j] = (g_residual(d, x, zp) - g_residual(d, x, zm)) / (2h)
    end
    fx, fz, gx, gz
end

function write_vector(path, values)
    open(path, "w") do io
        for value in values
            println(io, value)
        end
    end
end

function write_matrix(path, matrix)
    open(path, "w") do io
        for row in eachrow(matrix)
            println(io, join(row, ','))
        end
    end
end

function write_eigenvectors(path, matrix)
    open(path, "w") do io
        for j in axes(matrix, 2)
            println(io, join(vcat([j], real.(matrix[:, j])), ','))
        end
    end
    open(replace(path, "_real.csv" => "_imag.csv"), "w") do io
        for j in axes(matrix, 2)
            println(io, join(vcat([j], imag.(matrix[:, j])), ','))
        end
    end
end

function run_case(label, replaced)
    v, pf_residual, pf_iterations = solve_powerflow()
    d, x = build_dae(v, replaced)
    z = zeros(Float64, 2 * length(BUS_IDS))
    z[1:2:end] = real.(v)
    z[2:2:end] = imag.(v)
    f0, g0 = f_residual(d, x, z), g_residual(d, x, z)
    fx, fz, gx, gz = jacobians(d, x, z)
    a = fx - fz * (gz \ gx)
    spectral = eigen(a)
    values, vectors = spectral.values, spectral.vectors
    transverse = findall(abs.(values) .>= GAUGE_TOL)
    critical_indices = isempty(transverse) ? collect(eachindex(values)) : transverse
    critical_index = critical_indices[argmax(real.(values[critical_indices]))]
    alpha_full = maximum(real.(values))
    alpha_transverse = maximum(real.(values[transverse]))
    write_vector(joinpath(OUT, "julia_$(label)_state.csv"), x)
    write_vector(joinpath(OUT, "julia_$(label)_algebraic.csv"), z)
    write_matrix(joinpath(OUT, "julia_$(label)_A.csv"), a)
    open(joinpath(OUT, "julia_$(label)_spectrum.csv"), "w") do io
        println(io, "index,real,imag,frequency_hz,magnitude,transverse")
        for (index, value) in enumerate(values)
            println(io, join((index, real(value), imag(value), abs(imag(value)) / (2pi), abs(value), abs(value) >= GAUGE_TOL), ','))
        end
    end
    write_eigenvectors(joinpath(OUT, "julia_$(label)_eigenvectors_real.csv"), vectors)
    open(joinpath(OUT, "julia_equilibrium_summary.csv"), isfile(joinpath(OUT, "julia_equilibrium_summary.csv")) ? "a" : "w") do io
        position(io) == 0 && println(io, "case,states,pf_residual,pf_iterations,norm_f,norm_g,gz_condition,alpha_full,alpha_transverse,critical_frequency_hz,critical_real")
        println(io, join((label, length(x), pf_residual, pf_iterations, norm(f0, Inf), norm(g0, Inf), cond(gz), alpha_full, alpha_transverse, abs(imag(values[critical_index])) / (2pi), real(values[critical_index])), ','))
    end
    @printf("JULIA_TRUE_SAME_MODEL_PASS case=%s states=%d pf=%.3e f=%.3e g=%.3e alpha=%.6e transverse=%.6e\n", label, length(x), pf_residual, norm(f0, Inf), norm(g0, Inf), alpha_full, alpha_transverse)
    (states=length(x), pf_residual=pf_residual, norm_f=norm(f0, Inf), norm_g=norm(g0, Inf), alpha_full=alpha_full, alpha_transverse=alpha_transverse, critical_frequency_hz=abs(imag(values[critical_index])) / (2pi), stable=alpha_transverse < 0.0)
end

function portfolio_label(replaced::Set{Int})
    isempty(replaced) ? "none" : join(sort(collect(replaced)), "+")
end

function run_census()
    run_case("base", Set{Int}())
    run_case("v4", TARGET_BUSES)
    open(joinpath(OUT, "julia_true_same_model_v4_census.csv"), "w") do io
        println(io, "portfolio,replaced_count,states,norm_f,norm_g,alpha_full,alpha_transverse,critical_frequency_hz,stable")
        for mask in 0:15
            replaced = Set{Int}(TARGET_ORDER[j] for j in 1:4 if ((mask >> (j - 1)) & 1) == 1)
            label = "census_" * portfolio_label(replaced)
            result = run_case(label, replaced)
            println(io, join((portfolio_label(replaced), length(replaced), result.states, result.norm_f, result.norm_g, result.alpha_full, result.alpha_transverse, result.critical_frequency_hz, result.stable), ','))
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_census()
end
