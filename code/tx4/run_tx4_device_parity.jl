"""Evaluate the frozen GFL11 device at the same fixed parity point as Python."""

include(joinpath(@__DIR__, "run_tx4_exact_p4_ieee39.jl"))

const OUT_DEVICE = joinpath(CAMPAIGN, "raw", "tx4_exact_p4_julia", "device_parity_julia.csv")
mkpath(dirname(OUT_DEVICE))

const v_device = 0.97 + 0.11im
const x_device = [0.07, 0.003, 0.12, -0.08, 0.01, -0.02, 0.42, -0.11, 0.004, -0.001, -0.08]
const p_device = GFLParams(53.0, 1400.0, 0.03, 0.20, 8.0, 0.20, 8.0, 0.25, 6.0,
    0.15, 0.01, 0.12, -0.08, 0.98, 2.0, 20.0, 0.05)
const derivative_device = gfl_derivatives(p_device, x_device, v_device)
const injection_device = 0.37 * complex(x_device[7], x_device[8]) * cis(x_device[1])

open(OUT_DEVICE, "w") do io
    println(io, "quantity,value")
    for (index, value) in enumerate(derivative_device)
        println(io, "derivative_$(lpad(index, 2, '0')),$(repr(value))")
    end
    println(io, "injection_real,$(repr(real(injection_device)))")
    println(io, "injection_imag,$(repr(imag(injection_device)))")
end

println("JULIA_TX4_GFL11_DEVICE_PARITY_PASS")
