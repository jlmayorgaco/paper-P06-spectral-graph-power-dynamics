using Test, LinearAlgebra
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_design_p","ExpP.jl"))
include(joinpath(ROOT,"src","bnd_expQ","LinearSecurity.jl"))
include(joinpath(@__DIR__,"..","..","src","bnd_expQ2B","FiniteWindow.jl"))
using .FiniteWindow

@testset "causal phase-window estimator" begin
    T=0.5;dt=0.01;times=collect(-1.0:dt:8.0)
    phase=zeros(10,length(times));f0=0.3
    for (j,t) in enumerate(times)
        phase[:,j].=2pi*f0*t
    end
    sig=(;times_s=times,phase_rad=phase,phase_jump_rad=zeros(10),
        F_inf_Hz=fill(f0,10),poles=ComplexF64[],alpha=-0.1,
        A=zeros(1,1),B=zeros(1,1),C_bus=zeros(10,1),D_bus=zeros(10),
        condition_A=1.0,eigenvector_condition=1.0)
    rows=FiniteWindow.window_metrics(sig;windows=(T,),dt_s=dt,horizon_s=8.0,
        disturbance_MW=1.0)
    @test only(rows).F_inf_Hz ≈ f0 atol=1e-12
    @test only(rows).F_peak_Hz ≈ f0 atol=1e-12
    # The phase ramp starts at the event, so the causal window has a finite
    # startup RoCoF; after T seconds it settles to zero. The exact peak is
    # f0/T = 0.6 Hz/s for this fixture.
    @test only(rows).R_peak_Hz_s ≈ f0/T atol=1e-12
end

@testset "finite window bounds a phase jump" begin
    T=0.5;dt=0.01;times=collect(-1.0:dt:4.0)
    phase=zeros(10,length(times));phase[:,findall(times .>= 0)].=0.02
    sig=(;times_s=times,phase_rad=phase,phase_jump_rad=fill(0.02,10),
        F_inf_Hz=zeros(10),poles=ComplexF64[],alpha=-0.1,
        A=zeros(1,1),B=zeros(1,1),C_bus=zeros(10,1),D_bus=zeros(10),
        condition_A=1.0,eigenvector_condition=1.0)
    rows=FiniteWindow.window_metrics(sig;windows=(0.2,0.5,1.0,2.0),
        dt_s=dt,horizon_s=4.0,disturbance_MW=1.0)
    @test length(rows)==4
    @test all(isfinite(r.R_peak_Hz_s) for r in rows)
    @test rows[1].F_peak_Hz > rows[end].F_peak_Hz
end
