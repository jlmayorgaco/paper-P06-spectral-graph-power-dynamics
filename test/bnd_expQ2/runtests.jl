using Test, LinearAlgebra
include(joinpath(@__DIR__,"..","..","src","bnd_expQ","LinearSecurity.jl"))
include(joinpath(@__DIR__,"..","..","src","bnd_expQ","FrequencyMetrics.jl"))
include(joinpath(@__DIR__,"..","..","src","bnd_expQ2","GridFrequency.jl"))
using .LinearSecurity, .FrequencyMetrics, .GridFrequency

@testset "Architecture-invariant grid measurement metrics" begin
    dt=0.01; times=collect(-0.05:dt:0.20)
    phase=zeros(1,length(times)); phase[times.>=0.0].=0.01
    sig=(;times_s=times,phase_rad=phase,f_regular_Hz=zeros(1,length(times)),
        phase_jump_rad=[0.01],F_inf_Hz=[0.0])
    metrics=grid_step_metrics(sig;half_windows=(0,2,5),dt_s=dt,event_time_s=0.0)
    @test length(metrics)==3
    @test metrics[1].estimator=="UNFILTERED_100HZ_FINITE_DIFFERENCE"
    @test metrics[1].phase_jump_peak_rad==0.01
    @test metrics[1].R_peak_Hz_s > metrics[2].R_peak_Hz_s
    @test all(x.includes_event_phase_jump for x in metrics)
end

@testset "Condition-aware continuity comparison" begin
    ref=Dict("F_peak_Hz"=>0.1,"R_peak_Hz_s"=>0.2)
    rows=continuity_rows(30,zeros(10),ref,copy(ref);condition=10.0,
        ref_condition=10.0,rho_value=1e-8,alpha=-0.1)
    @test all(x.pass_within_tolerance for x in rows)
    off=Dict("F_peak_Hz"=>0.2,"R_peak_Hz_s"=>0.2)
    rows=continuity_rows(30,zeros(10),ref,off;condition=10.0,
        ref_condition=10.0,rho_value=1e-8,alpha=-0.1)
    @test !rows[1].pass_within_tolerance
end
