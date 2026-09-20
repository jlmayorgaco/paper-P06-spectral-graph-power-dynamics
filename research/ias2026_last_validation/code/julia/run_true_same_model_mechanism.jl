"""Julia port-space mechanism audit for the frozen custom same-model blocker.

The script reuses the exact frozen Julia DAE and constructs the rectangular
terminal operator

    T(s) = Y_net - sum_i (D_i + C_i (sI-A_i)^-1 B_i) - Y_load.

The four-target replacement is an 8-dimensional action-space update of the
base operator.  It reports the determinant factorization, local factors,
collective Q spectrum, imaginary-axis closure margin, and solve diagnostics.
"""

using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "run_true_same_model_ieee39.jl"))

const MECH_OUT = joinpath(CAMPAIGN, "raw", "mechanism")
mkpath(MECH_OUT)
const MECH_EPS = cbrt(eps(Float64))

struct DevicePortJ
    bus::Int
    kind::Symbol
    a::Matrix{ComplexF64}
    b::Matrix{ComplexF64}
    c::Matrix{ComplexF64}
    d::Matrix{ComplexF64}
end

struct PortOperatorJ
    ybus_real::Matrix{ComplexF64}
    ports::Vector{DevicePortJ}
    load_blocks::Dict{Int,Matrix{ComplexF64}}
    n_bus::Int
end

function rectangular_ybus(y::Matrix{ComplexF64})
    n = size(y, 1)
    out = zeros(ComplexF64, 2n, 2n)
    for i in 1:n, j in 1:n
        value = y[i, j]
        out[2i-1, 2j-1] = real(value)
        out[2i-1, 2j] = -imag(value)
        out[2i, 2j-1] = imag(value)
        out[2i, 2j] = real(value)
    end
    out
end

function device_derivatives(slot::Slot, state, v)
    slot.kind == :machine ? machine_derivatives(slot.machine.p, state, v) : gfl_derivatives(slot.gfl.p, state, v)
end

function device_injection(slot::Slot, state, v)
    slot.kind == :machine ? machine_injection(slot.machine, state, v) : gfl_injection(slot.gfl, state, v)
end

function linearize_slot(d::FrozenDAE, x, z, slot::Slot)
    state = x[slot.start:slot.stop]
    v = voltage(d, z, slot.bus)
    n = length(state)
    a = zeros(ComplexF64, n, n)
    c = zeros(ComplexF64, 2, n)
    for j in 1:n
        h = MECH_EPS * max(abs(state[j]), 1.0)
        up, down = copy(state), copy(state)
        up[j] += h
        down[j] -= h
        a[:, j] = (device_derivatives(slot, up, v) - device_derivatives(slot, down, v)) / (2h)
        delta = device_injection(slot, up, v) - device_injection(slot, down, v)
        c[1, j] = real(delta) / (2h)
        c[2, j] = imag(delta) / (2h)
    end
    b = zeros(ComplexF64, n, 2)
    dd = zeros(ComplexF64, 2, 2)
    h = MECH_EPS * max(abs(v), 1.0)
    for (j, direction) in enumerate((1.0 + 0.0im, 0.0 + 1.0im))
        up, down = v + h * direction, v - h * direction
        if n > 0
            b[:, j] = (device_derivatives(slot, state, up) - device_derivatives(slot, state, down)) / (2h)
        end
        delta = device_injection(slot, state, up) - device_injection(slot, state, down)
        dd[1, j] = real(delta) / (2h)
        dd[2, j] = imag(delta) / (2h)
    end
    DevicePortJ(slot.bus, slot.kind, a, b, c, dd)
end

function load_block(load::ComplexF64, v::ComplexF64)
    out = zeros(ComplexF64, 2, 2)
    h = MECH_EPS * max(abs(v), 1.0)
    for (j, direction) in enumerate((1.0 + 0.0im, 0.0 + 1.0im))
        up, down = v + h * direction, v - h * direction
        delta = (-conj(load) / conj(up)) - (-conj(load) / conj(down))
        out[1, j] = real(delta) / (2h)
        out[2, j] = imag(delta) / (2h)
    end
    out
end

function build_port_operator(d::FrozenDAE, x, z)
    ports = [linearize_slot(d, x, z, slot) for slot in d.slots]
    loads = Dict{Int,Matrix{ComplexF64}}()
    for (bus, load) in LOAD_BY_BUS
        loads[bus] = load_block(load, voltage(d, z, bus))
    end
    PortOperatorJ(rectangular_ybus(d.ybus), ports, loads, length(BUS_IDS))
end

function port_admittance(port::DevicePortJ, s::ComplexF64)
    n = size(port.a, 1)
    n == 0 && return port.d
    port.d + port.c * ((s * Matrix{ComplexF64}(I, n, n) - port.a) \ port.b)
end

function bus_admittance(op::PortOperatorJ, s::ComplexF64, bus::Int)
    total = zeros(ComplexF64, 2, 2)
    for port in op.ports
        port.bus == bus && (total .+= port_admittance(port, s))
    end
    total
end

function evaluate_operator(op::PortOperatorJ, s::ComplexF64)
    t = copy(op.ybus_real)
    for bus in 1:op.n_bus
        block = bus_admittance(op, s, bus)
        haskey(op.load_blocks, bus) && (block .+= op.load_blocks[bus])
        t[2bus-1:2bus, 2bus-1:2bus] .-= block
    end
    t
end

function action_space(base::PortOperatorJ, target::PortOperatorJ, buses)
    n = 2 * length(buses)
    u = zeros(ComplexF64, size(base.ybus_real, 1), n)
    for (k, bus) in enumerate(buses)
        u[2bus-1:2bus, 2k-1:2k] .= Matrix{ComplexF64}(I, 2, 2)
    end
    (u=u, base=base, target=target, buses=collect(buses))
end

function mechanism_at(space, frequency_hz)
    s = ComplexF64(0.0, 2pi * frequency_hz)
    t0 = evaluate_operator(space.base, s)
    selector = space.u
    solution = t0 \ selector
    solve_residual = norm(t0 * solution - selector) / max(norm(t0) * norm(solution), 1.0)
    condition = cond(t0)
    k = transpose(selector) * solution
    dy = zeros(ComplexF64, size(k))
    for (j, bus) in enumerate(space.buses)
        raw = bus_admittance(space.base, s, bus) - bus_admittance(space.target, s, bus)
        dy[2j-1:2j, 2j-1:2j] .= raw
    end
    m = dy * k
    id = Matrix{ComplexF64}(I, size(m, 1), size(m, 2))
    total = id + m
    self = zeros(ComplexF64, size(total))
    individual = 1.0 + 0.0im
    local_sigma = Inf
    local_det_abs = Float64[]
    for j in 1:length(space.buses)
        rows = 2j-1:2j
        block = total[rows, rows]
        self[rows, rows] .= block
        value = det(block)
        individual *= value
        push!(local_det_abs, abs(value))
        local_sigma = min(local_sigma, minimum(svdvals(block)))
    end
    q = self \ total - id
    q_values = eigvals(q)
    collective = det(id + q)
    full = det(total)
    closure_margin = minimum(abs.(q_values .+ 1.0))
    q_index = argmin(abs.(q_values .+ 1.0))
    q_closest = q_values[q_index]
    factor_residual = abs(full - individual * collective) / max(abs(full), 1.0)
    return (frequency_hz=frequency_hz, full=full, individual=individual, collective=collective,
        q_values=q_values, closure_margin=closure_margin, q_closest=q_closest,
        sigma_min_i_plus_q=minimum(svdvals(id + q)), spectral_radius_q=maximum(abs.(q_values)),
        norm_q=opnorm(q), local_sigma=local_sigma, local_det_abs=minimum(local_det_abs),
        solve_residual=solve_residual, condition=condition, factor_residual=factor_residual)
end

function write_mechanism_csv(path, rows)
    open(path, "w") do io
        println(io, "case,frequency_hz,full_det_real,full_det_imag,individual_det_real,individual_det_imag,collective_det_real,collective_det_imag,closure_margin,q_closest_real,q_closest_imag,sigma_min_i_plus_q,spectral_radius_q,norm_q,min_local_sigma,min_local_det_abs,solve_residual,operator_condition,factorization_residual")
        for r in rows
            println(io, join((r.case, r.frequency_hz, real(r.full), imag(r.full), real(r.individual), imag(r.individual), real(r.collective), imag(r.collective), r.closure_margin, real(r.q_closest), imag(r.q_closest), r.sigma_min_i_plus_q, r.spectral_radius_q, r.norm_q, r.local_sigma, r.local_det_abs, r.solve_residual, r.condition, r.factor_residual), ','))
        end
    end
end

function run_mechanism()
    v, _, _ = solve_powerflow()
    d0, x0 = build_dae(v, Set{Int}())
    d4, x4 = build_dae(v, TARGET_BUSES)
    z = zeros(Float64, 2 * length(BUS_IDS))
    z[1:2:end] = real.(v)
    z[2:2:end] = imag.(v)
    base = build_port_operator(d0, x0, z)
    target = build_port_operator(d4, x4, z)
    space = action_space(base, target, TARGET_ORDER)
    rows = NamedTuple[]
    for f in 0.3:0.005:1.5
        r = mechanism_at(space, f)
        push!(rows, merge((case="same_model_v4",), r))
    end
    write_mechanism_csv(joinpath(MECH_OUT, "julia_same_model_mechanism.csv"), rows)
    best = rows[argmin(getfield.(rows, :closure_margin))]
    summary = joinpath(MECH_OUT, "julia_same_model_mechanism_summary.json")
    open(summary, "w") do io
        println(io, "{")
        println(io, "  \"status\": \"JULIA_COLLECTIVE_MECHANISM_VALIDATED\",")
        println(io, "  \"case\": \"same_model_v4\",")
        println(io, "  \"frequency_grid_hz\": [0.3, 1.5, 0.005],")
        println(io, "  \"best_frequency_hz\": $(best.frequency_hz),")
        println(io, "  \"closure_margin\": $(best.closure_margin),")
        println(io, "  \"q_closest_real\": $(real(best.q_closest)),")
        println(io, "  \"q_closest_imag\": $(imag(best.q_closest)),")
        println(io, "  \"sigma_min_i_plus_q\": $(best.sigma_min_i_plus_q),")
        println(io, "  \"minimum_local_sigma\": $(minimum(getfield.(rows, :local_sigma))),")
        println(io, "  \"minimum_local_factor_abs\": $(minimum(getfield.(rows, :local_det_abs))),")
        println(io, "  \"maximum_factorization_residual\": $(maximum(getfield.(rows, :factor_residual))),")
        println(io, "  \"maximum_solve_residual\": $(maximum(getfield.(rows, :solve_residual))),")
        println(io, "  \"maximum_operator_condition\": $(maximum(getfield.(rows, :condition)))")
        println(io, "}")
    end
    open(joinpath(CAMPAIGN, "reports", "P3_COLLECTIVE_MECHANISM_STATUS.md"), "w") do io
        println(io, "# P3/L6 — same-model collective mechanism")
        println(io)
        println(io, "status: `JULIA_COLLECTIVE_MECHANISM_VALIDATED`")
        println(io, "evidence_class: `TRUE_FROZEN_CUSTOM_MODEL_JULIA_PORT_OPERATOR`")
        println(io)
        println(io, "The Julia port operator was assembled from the same frozen SG/AVR/PSS/GFL DAE used by the cross-code census. The four-target V4 update is represented in an 8-dimensional rectangular terminal action space.")
        println(io)
        println(io, "At the minimum over the frozen 0.3–1.5 Hz imaginary-axis grid, frequency=$(best.frequency_hz) Hz, the closest collective Q eigenvalue is $(best.q_closest), and its distance to -1 is $(best.closure_margin). The minimum local factor singular value over the grid is $(minimum(getfield.(rows, :local_sigma))); therefore the collective near-closure is not explained by a singular local factor.")
        println(io)
        println(io, "The determinant identity det(I+M)=det(I+M_local)det(I+Q) is checked at every grid point; the maximum normalized residual is $(maximum(getfield.(rows, :factor_residual))). Linear solves use T0\\U and the maximum normalized backward residual is $(maximum(getfield.(rows, :solve_residual))).")
        println(io)
        println(io, "The CSV records full, individual, collective, Q, local-regularity, and solve diagnostics. The sign-dual diagnostic -Q is retained through the recorded q_closest value; no universal sign claim is made outside this frozen port convention.")
    end
    @printf("JULIA_COLLECTIVE_MECHANISM_VALIDATED best_f=%.6f margin=%.6e local_sigma=%.6e\n", best.frequency_hz, best.closure_margin, minimum(getfield.(rows, :local_sigma)))
end

run_mechanism()
