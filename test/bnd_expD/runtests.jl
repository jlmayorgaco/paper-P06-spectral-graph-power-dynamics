using Test
using LinearAlgebra
using CSV
using DataFrames

include(joinpath(@__DIR__,"..","..","src","bnd_design","BNDDesign.jl"))
using .BNDDesign.AnalyticDeviceModel
using .BNDDesign.AnalyticGFLPLL
using .BNDDesign.LowRankUpdates
using .BNDDesign.WoodburyPLL
using .BNDDesign.ReplacementClosure
using .BNDDesign.AlgebraicBoundary
using .BNDDesign.RhoDesign
using .BNDDesign.AnalyticSG: rating_fractions
import .BNDDesign.AnalyticSG
using .BNDDesign.PortReduction

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_validation","PowerDynamicsValidation.jl"))

@testset "ExpD analytic / validation firewall" begin
    files=filter(f->endswith(f,".jl"),[joinpath(d,f) for (d,_,fs) in walkdir(joinpath(ROOT,"src","bnd_design")) for f in fs])
    imp=r"(?m)^\s*(?:using|import)\s+(?:PowerDynamics|[A-Za-z0-9_\.]*\.PowerDynamics)\b"
    @test !isempty(files)
    @test all(!occursin(imp,read(f,String)) for f in files)
    @test !isfile(joinpath(ROOT,"reports","experiment_D","ANALYTIC_CANDIDATE.json"))
    gate=PowerDynamicsValidation.candidate_gate(ROOT)
    @test !gate.allowed
    validation=PowerDynamicsValidation.run_validation(ROOT)
    @test validation.status=="REFUSED_NO_FROZEN_CANDIDATE"
    @test validation.powerdynamics_imported==false
    @test validation.powerdynamics_calls==0
end

@testset "Independent SimpleGFLDC PLL affine structure" begin
    op=frozen_bus33_operating_point(ROOT)
    x,u,p=op.x,op.u,op.parameters
    J=jacobians(x,u,p)
    gd=gain_derivatives(x,u,p)
    @test length(x)==9
    @test length(rhs(x,u,p))==9
    @test rank_audit(gd.A_kp).rank==1
    @test rank_audit(gd.A_ki).rank==1
    @test rank_audit(hcat(gd.A_kp,gd.A_ki)).rank==2
    @test rank_audit(gd.B_kp).rank==1
    @test rank_audit(gd.B_ki).rank==1
    @test rank_audit(hcat(gd.B_kp,gd.B_ki)).rank==2
    @test norm(gd.C_kp)+norm(gd.C_ki)+norm(gd.D_kp)+norm(gd.D_ki)==0
    for (kp,ki) in ((0.5p.pll_kp,0.7p.pll_ki),(1.7p.pll_kp,1.2p.pll_ki))
        Jp=jacobians(x,u,p;kp,ki)
        @test norm(Jp.A-J.A-(kp-p.pll_kp)*gd.A_kp-(ki-p.pll_ki)*gd.A_ki)/norm(Jp.A)<1e-11
        @test norm(Jp.B-J.B-(kp-p.pll_kp)*gd.B_kp-(ki-p.pll_ki)*gd.B_ki)/norm(Jp.B)<1e-11
        @test Jp.C≈J.C atol=0 rtol=0
        @test Jp.D≈J.D atol=0 rtol=0
    end
end

@testset "Local PLL Woodbury and parameter denominator" begin
    op=frozen_bus33_operating_point(ROOT)
    result=verify_woodbury(op.x,op.u,op.parameters;samples=12)
    @test result.rank_Akp==1
    @test result.rank_Aki==1
    @test result.rank_both_A==2
    @test result.resolvent_median<1e-11
    @test result.resolvent_p95<1e-9
    @test result.port_median<1e-11
    @test result.port_p95<1e-9
    @test result.determinant_factorization_max<1e-9
    F=WoodburyPLL.factorization(op.x,op.u,op.parameters)
    s=0.2+2.5im
    c=denominator_coefficients(F.A0,F.U,F.V,s)
    for (kp,ki) in ((2.0,12.0),(40.0,250.0),(70.0,500.0))
        direct=det(I-Diagonal([kp,ki])*(F.V'*((s*I-F.A0)\F.U)))
        polynomial=c.c00+c.c10*kp+c.c01*ki+c.c11*kp*ki
        @test abs(direct-polynomial)<1e-10
    end
end

@testset "Independent rated SG/GFL local ports at bus 33" begin
    sg=AnalyticSG.frozen_bus_operating_point(ROOT,33)
    sd=AnalyticSG.steady_state_diagnostics(sg.x,sg.u,sg.parameters)
    @test sd.rhs_norm<1e-10
    @test abs(sd.network_power.p-6.32)<1e-9
    @test abs(sd.network_power.q-1.0990598880307398)<1e-8

    gfl_sys=AnalyticGFLPLL.frozen_bus33_operating_point(ROOT)
    gfl=rebase_to_rating(gfl_sys,800.0,100.0)
    @test norm(AnalyticGFLPLL.rhs(gfl.x,gfl.u,gfl.parameters),Inf)<1e-10
    @test norm(rated_output(gfl.x,gfl.u,gfl.parameters,8.0)-AnalyticGFLPLL.output(
        gfl_sys.x,gfl_sys.u,gfl_sys.parameters))<1e-10

    Js=AnalyticSG.jacobians(sg.x,sg.u,sg.parameters)
    Jg_sys=AnalyticGFLPLL.jacobians(gfl_sys.x,gfl_sys.u,gfl_sys.parameters)
    Jg_dev=AnalyticGFLPLL.jacobians(gfl.x,gfl.u,gfl.parameters)
    rel=[]
    for hz in 10.0 .^ range(log10(0.01),log10(10.0),length=41)
        s=0.05+2pi*hz*im
        Ysys=port_transfer(Jg_sys.A,Jg_sys.B,Jg_sys.C,Jg_sys.D,s)
        Yrated=rated_port_transfer(Jg_dev,s,8.0)
        push!(rel,norm(Ysys-Yrated)/max(norm(Ysys),eps(Float64)))
    end
    @test maximum(rel)<1e-12

    Ysg=port_transfer(Js.A,Js.B,Js.C,Js.D,0.1+1.0im)
    Ygfl=rated_port_transfer(Jg_dev,0.1+1.0im,8.0)
    for rho in (0.0,0.1,0.37,0.8,1.0)
        @test affine_rho_residual(Ysg,Ygfl,rho)<1e-14
    end
end

@testset "Local rated-port replication across buses 30/33/35/37" begin
    for bus in (30,33,35,37)
        sg=AnalyticSG.frozen_bus_operating_point(ROOT,bus)
        sd=AnalyticSG.steady_state_diagnostics(sg.x,sg.u,sg.parameters)
        gs=AnalyticGFLPLL.operating_point_for_injection(sg.u,sd.network_power.p,sd.network_power.q)
        ratio=sg.parameters.rating_mva/sg.parameters.system_base_mva
        go=rebase_to_rating(gs,sg.parameters.rating_mva,sg.parameters.system_base_mva)
        @test sd.rhs_norm<1e-10
        @test norm(AnalyticGFLPLL.rhs(go.x,go.u,go.parameters),Inf)<1e-10
        @test norm(sd.port_current-rated_output(go.x,go.u,go.parameters,ratio))<1e-9
        J0=AnalyticGFLPLL.jacobians(gs.x,gs.u,gs.parameters)
        J1=AnalyticGFLPLL.jacobians(go.x,go.u,go.parameters)
        for s in (0.05+0.1im,0.2+3.5im,-0.1+6.2im)
            Y0=port_transfer(J0.A,J0.B,J0.C,J0.D,s)
            @test norm(Y0-rated_port_transfer(J1,s,ratio))/max(norm(Y0),eps(Float64))<1e-11
        end
    end
end

@testset "Generic exact algebra helpers" begin
    A=randn(5,5); U=randn(5,2); V=randn(5,2); Delta=randn(2,2)
    T0=2I+0.3randn(5,5)
    d=determinant_lemma_residual(T0,U,Delta,V)
    @test d.logabsdet_residual<1e-10
    @test d.phase_residual<1e-10
    coeff=rho_determinant_coefficients(Diagonal([-2.0,-3.0]))
    @test coeff≈[1.0,-5.0,6.0]
    roots=rho_candidates(coeff)
    @test roots≈[1/3,1/2] atol=1e-10
    c10=1.0+0im; c01=0.0+2im; c11=1.0-0.5im
    c00=-(2c10+3c01+6c11)
    candidates=bilinear_gain_candidates(c00,c10,c01,c11)
    @test any(c->abs(c.Kp-2)<1e-8 && abs(c.Ki-3)<1e-8,candidates)
    @test rating_fractions(0.0)==(sg=1.0,gfl=0.0)
    @test rating_fractions(1.0)==(sg=0.0,gfl=1.0)
end
