using Test, LinearAlgebra, DataFrames

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "src", "bnd_expQ", "FrequencyMetrics.jl"))
using .FrequencyMetrics
include(joinpath(ROOT, "src", "bnd_expQ", "LinearSecurity.jl"))
using .LinearSecurity

@testset "Frequency definitions" begin
    @test sg_frequency_deviation(1.001) ≈ 0.06 atol=1e-12
    @test pll_frequency_deviation(2pi * 0.1) ≈ 0.1 atol=1e-12
    phase = [3.0, 3.1, -3.083185307179586, -2.983185307179586]
    unwrapped = unwrap_phase(phase)
    @test all(diff(unwrapped) .> 0)
    @test unwrapped[end] ≈ 3.3 atol=1e-12
end

@testset "Linear step transfer and extrema" begin
    model=(;A=reshape([-2.0],1,1),B=[1.0],C=reshape([3.0],1,1),
           output_metadata=DataFrame(bus=[38],kind=["SG"]))
    met=step_metrics(model;disturbance_MW=10.0,horizon_s=5.0,dt_s=0.002)
    @test met.steady_frequency_unit ≈ 1.5 atol=1e-12
    @test met.frequency_peak_unit ≈ 1.5 atol=1e-10
    @test met.rocof_peak_unit ≈ 3.0 atol=1e-10
    @test met.frequency_peak_MW ≈ 15.0 atol=1e-9
    @test met.rocof_peak_MW ≈ 30.0 atol=1e-9
end

@testset "Savitzky-Golay derivative" begin
    dt = 0.01
    t = collect(0.0:dt:1.0)
    y = 2 .+ 3 .* t .- 4 .* t.^2 .+ 0.5 .* t.^3
    exact = 3 .- 8 .* t .+ 1.5 .* t.^2
    got = sg_polynomial_derivative(y, dt; half_window=5, degree=3)
    @test maximum(abs.(got .- exact)) < 1e-9
    f = 0.2
    voltage = cis.(2pi*f .* t)
    fb = bus_frequency_deviation(voltage, dt; half_window=5, degree=3)
    @test maximum(abs.(fb[6:end-5] .- f)) < 1e-9
    cutoff = savgol_cutoff_hz(dt; half_window=5, degree=3)
    @test 0 < cutoff < 0.5/dt
end
