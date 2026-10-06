using Test, LinearAlgebra, CSV, DataFrames, TOML, SHA

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_design_p","ExpP.jl"))
using .ExpP
include(joinpath(ROOT,"src","bnd_design_p","AlgebraicFrontiers.jl"))
using .AlgebraicFrontiers

@testset "ExpP frozen inputs and stage gates" begin
    @test ExpP.verify_expN_freeze(ROOT)=="e2f104608f1eb1705f0ad2c07764e7beab0a5df4deb4d4ce4f359b1e79947f0a"
    cand,csha=ExpP.read_candidate(ROOT)
    @test csha=="3915a8f57da0f552779b11464572f14a86f0b648b532be77124f668eeff9f263"
    p0=CSV.read(joinpath(ROOT,"reports","experiment_P","P0","TABLE_P0_regression.csv"),DataFrame)
    @test all(p0.pass)
    p5=CSV.read(joinpath(ROOT,"reports","experiment_P","P5","TABLE_P5_powerdynamics_validation.csv"),DataFrame)
    @test only(p5.PD_validation)=="PASS"
    tds=CSV.read(joinpath(ROOT,"reports","experiment_P","P5","TABLE_P5_tds_validation.csv"),DataFrame)
    @test nrow(tds)==6 && all(tds.TDS_validation .== "PASS")
    linear=CSV.read(joinpath(ROOT,"reports","experiment_P","P5","TABLE_P5_linear_nonlinear_comparison.csv"),DataFrame)
    @test nrow(linear)==6 && all(linear.match_pass)
    @test maximum(linear.frequency_relative_error)<0.05
    @test maximum(linear.rocof_relative_error)<0.05
    event=CSV.read(joinpath(ROOT,"reports","experiment_P","P5","TABLE_P5_declared_event_validation.csv"),DataFrame)
    @test only(event.event_bus)==16 && only(event.delta_P_MW)==100.0
    @test only(event.pulse_profile)=="sustained Pset step at 1.0 s"
    @test only(event.solver_status)=="Success"
    @test !only(event.overall_pass) && !only(event.rocof_pass) && !only(event.frequency_pass)
    @test only(event.tail_status)=="NOT_SETTLED_WITHIN_60S"
    @test only(event.peak_RoCoF_Hz_s)>only(event.rocof_limit_Hz_s)
    @test only(event.peak_COI_frequency_Hz)>only(event.frequency_limit_Hz)
    result=read(joinpath(ROOT,"reports","experiment_P","P5","P5_RESULTS.json"),String)
    @test occursin("\"status\": \"FAIL_TDS\"",result)
    @test occursin("\"nominal_PD_subgate\": \"PASS\"",result)
    @test occursin("\"small_pulse_scaling_subgate\": \"PASS\"",result)
    guard=CSV.read(joinpath(ROOT,"reports","experiment_P","P2","TABLE_P2_guard_policy.csv"),DataFrame)
    @test only(guard.delta_num_per_s)==1e-6
    @test only(guard.guarded_feasible)
    @test only(guard.all_physical_poles_checked)
    @test only(guard.source_candidate_unchanged) && !only(guard.candidate_was_modified)
    @test 0 < only(guard.retained_MW_increase_vs_frozen) < 1e-5
    p4=read(joinpath(ROOT,"reports","experiment_P","P4","P4_RESULTS.json"),String)
    @test occursin("\"bounded_real_certified\":false",p4)
    @test occursin("\"bounded_real_gamma_upper\":null",p4)
    @test occursin("\"conditional_point_interval_certificate_type\":\"LIPSCHITZ_INTERVAL_FLOAT64_NO_OUTWARD_ROUNDING\"",p4)
    @test occursin("\"joint_robust_KKT_status\":\"NOT_RUN; one-coordinate scalar continuation only\"",p4)
    p4interval=CSV.read(joinpath(ROOT,"reports","experiment_P","P4","TABLE_P4_interval_certificate_conditional.csv"),DataFrame)
    @test only(p4interval.certified) && only(p4interval.nominal_spectrum_pass)
    @test only(p4interval.beta_lower_cert)>only(p4interval.beta_req)
    @test only(p4interval.beta_upper_observed)>=only(p4interval.beta_lower_cert)
    @test only(p4interval.physical_poles)==101 && only(p4interval.nodes)==34333
end

@testset "P1 local PLL rank-one and SG rank-two identities" begin
    ctx=ExpP.PDExactDesignN.design_context(ROOT)
    c,_=ExpP.read_candidate(ROOT)
    rho=Float64.(c["rho"]);kp=fill(ExpP.PDExactDesignN.K0P,10);ki=fill(ExpP.PDExactDesignN.K0I,10)
    pll=AlgebraicFrontiers.pll_pair_identity(ctx,rho,kp,ki)
    @test pll.symbolic.exact
    @test maximum(pll.table.rank2_over_rank1)<1e-10
    @test maximum(pll.table.common_detector_rel_error)<1e-10
    @test maximum(pll.table.matrix_affine_residual)<1e-12
    @test maximum(pll.table.determinant_same_device_cross)<1e-7
    @test maximum(pll.table.determinant_lemma_abs_error)<1e-7
    @test pll.max_between_device_cross.cross>1e-12

    single=AlgebraicFrontiers.single_sg_identity(ctx,kp,ki)
    @test maximum(single.operator_relative_error)<1e-10
    @test maximum(single.determinant_identity_abs_error)<1e-7
    pair=AlgebraicFrontiers.two_sg_identity(ctx,kp,ki)
    @test maximum(pair.operator_relative_error)<1e-10
    @test maximum(pair.determinant_identity_abs_error)<1e-7
    @test maximum(pair.cross_term)>1e-12
end

@testset "P2 exact retention root checks complete physical poles" begin
    ctx=ExpP.PDExactDesignN.design_context(ROOT)
    c,_=ExpP.read_candidate(ROOT)
    allroots=AlgebraicFrontiers.retention_roots(ctx,Float64.(c["Kp"]),Float64.(c["Ki"]),38)
    epsilon0=1-Float64(c["rho"][9])
    rootrow=allroots.roots[argmin(abs.(allroots.roots.epsilon.-epsilon0)),:]
    guardrow=allroots.guarded[argmin(abs.(allroots.guarded.epsilon.-epsilon0)),:]
    @test abs(rootrow.epsilon-epsilon0)/epsilon0<5e-7
    @test Bool(guardrow.guarded_boundary_pass)
    @test maximum(real.(ExpP.design_spectrum(ctx,
        [1,1,1,1,1,1,1,1,1-Float64(guardrow.epsilon),1],
        Float64.(c["Kp"]),Float64.(c["Ki"])).lambda))<=-0.05+2e-7
end

@testset "P3 local KKT and P6 open globality ledger" begin
    kkt=CSV.read(joinpath(ROOT,"reports","experiment_P","P3","TABLE_P3_KKT_audit.csv"),DataFrame)
    @test only(kkt.status)=="LOCAL_KKT_CERTIFIED"
    @test only(kkt.LICQ)
    @test only(kkt.stationarity_residual)<1e-8
    p6=CSV.read(joinpath(ROOT,"reports","experiment_P","P6","TABLE_P6_support_globality_ledger.csv"),DataFrame)
    @test nrow(p6)==1024
    @test count(p6.joint_fixed_support_KKT_evaluated)==1
    @test all(p6.nominal_spectral_objective_lower_bound_MW .== 0.0)
    @test all(.!p6.combined_requirements_incumbent_feasible)
    p6result=read(joinpath(ROOT,"reports","experiment_P","P6","P6_RESULTS.json"),String)
    @test occursin("\"globality_scope\":\"NOMINAL_COMPLETE_SPECTRUM_ONLY\"",p6result)
    @test occursin("\"combined_requirements_feasible_upper_bound_MW\":null",p6result)
    @test occursin("\"declared_transient_pass\":false",p6result)
end
